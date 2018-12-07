package org.jcvi.vigor.service;

import com.google.common.base.Predicates;
import com.google.common.collect.Lists;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.springframework.stereotype.Service;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

/**
 * Created by snettem on 5/16/2017.
 */
@Service
public class ModelGenerationService {

    private static final Logger LOGGER = LogManager.getLogger(ModelGenerationService.class);
    private static final int DEFAULT_MAX_ALIGN_MERGE_AA_GAP = 10;
    private static final int DEFAULT_RELAX_MERGE_AA_GAP = 300;
    private static final int DEFAULT_NTOVERLAP_OFFSET = 30;
    private static final int DEFAULT_AAOVERLAP_OFFSET = 10;

    public List<Model> generateModels ( List<Alignment> alignments, VigorConfiguration configuration ) throws ServiceException {

        alignments = mergeIdenticalProteinAlignments(alignments);
        return determineCandidateModels(alignments, configuration);
    }

    /**
     *
     * @param alignments
     * @return
     */
    public List<Alignment> mergeIdenticalProteinAlignments ( List<Alignment> alignments ) {

        List<Alignment> allOutAlignments = new ArrayList<Alignment>();
        Map<String, List<Alignment>> protein2AlignmentsMap = new HashMap<String, List<Alignment>>();

        // sort alignments by protein ID in a map
        Comparator<Alignment> alignmentComparator = ( alignment1, alignment2 )
                -> alignment1.getAlignmentFragments().get(alignment1.getAlignmentFragments().size() - 1).compareTo(alignment2.getAlignmentFragments().get(0));

        //Merge alignments(alignment fragments) belonging to a protein and add up the scores of the alignments
        for (Alignment alignment : alignments) {
            alignment.getAlignmentFragments().sort(AlignmentFragment.Comparators.Ascending);
            String proteinID = alignment.getViralProtein().getProteinID();
            protein2AlignmentsMap.computeIfAbsent(proteinID,(k) -> new ArrayList<>()).add(alignment);
        }

        for (String proteinID : protein2AlignmentsMap.keySet()) {
            Collections.sort(protein2AlignmentsMap.get(proteinID), alignmentComparator);
            List<Alignment> alignmentsTemp = protein2AlignmentsMap.get(proteinID);
            Alignment mergedAlignment = alignmentsTemp.remove(0);

            for (Alignment alignment : alignmentsTemp) {
                mergedAlignment.getAlignmentFragments().addAll(alignment.getAlignmentFragments());
                Map<String, Double> scores = mergedAlignment.getAlignmentScore();
                double score = scores.get(Scores.ALIGNMENT_SCORE);
                score = score + alignment.getAlignmentScore().get(Scores.ALIGNMENT_SCORE);
                scores.put(Scores.ALIGNMENT_SCORE, score);
                mergedAlignment.setAlignmentScore(scores);
            }
            allOutAlignments.add(mergedAlignment);
        }
        return allOutAlignments;
    }

    /**
     * @param alignments
     * @param configuration
     * @return all the possible models of each alignment and after splitting
     * models at the sequence gaps
     */
    public List<Model> determineCandidateModels ( List<Alignment> alignments, VigorConfiguration configuration) throws ServiceException {

        List<Model> initialModels = new ArrayList<Model>();
        List<Model> candidateModels = new ArrayList<Model>();
        for (Alignment alignment: alignments) {
            alignment.getAlignmentFragments().sort(AlignmentFragment.Comparators.Ascending);
            initialModels.addAll(alignmentToModels(alignment, alignment.getViralProtein().getConfiguration()));
        }
        if (configuration.getOrDefault(ConfigurationParameters.Verbose,false))  {
            FormatVigorOutput.printModels(initialModels, "Initial Models");
        }
        List<Range> sequenceGaps = new ArrayList<Range>();
        // get sequence gaps
        if (( !initialModels.isEmpty() ) && initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps() != null) {
            sequenceGaps = initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps();
        }
        int minGapLength = configuration.getOrDefault(ConfigurationParameters.SequenceGapMinimumLength, 0);
        List<Range> validSequenceGaps = sequenceGaps.stream()
                                                    .filter(g -> g.getLength() >= minGapLength)
                                                    .collect(Collectors.toList());
        // split models at sequence gaps
        for (Model model : initialModels) {
            try {
                candidateModels.addAll(splitModelAtSequenceGaps(model, validSequenceGaps));
            } catch (CloneNotSupportedException e) {
                LOGGER.error("problem splitting model {}", model);
                throw new ServiceException(String.format("problem splitting model %s", model), e);
            }
        }
        return candidateModels;
    }

    /**
     * @param alignment
     * @return Models of each alignment.
     */
    // TODO configuration values from viral proteins?
    public List<Model> alignmentToModels ( Alignment alignment, VigorConfiguration defaultConfiguration ) {
        // TODO rely on viral protein having configuration
        VigorConfiguration configuration = alignment.getViralProtein().getConfiguration();
        configuration = configuration != null ? configuration : defaultConfiguration;

        int maxAlignMergeAAGap = configuration.getOrDefault(ConfigurationParameters.MaxAlignMergeAAGap, DEFAULT_MAX_ALIGN_MERGE_AA_GAP);
        int minIntronLength = maxAlignMergeAAGap * 3;
        int relaxMergeAAGap = configuration.getOrDefault(ConfigurationParameters.RelaxAlignMergeAAGap, DEFAULT_RELAX_MERGE_AA_GAP);
        int relaxIntronLength = relaxMergeAAGap * 3;
        //Group alignment fragments based on direction
        Map<Direction, List<AlignmentFragment>> alignmentFragsGroupedList = alignment.getAlignmentFragments().stream()
                .collect(Collectors.groupingBy(w -> w.getDirection()));

        List<Model> models = new ArrayList<Model>();
        for (Direction direction: alignmentFragsGroupedList.keySet()) {
            List<AlignmentFragment> fragments = alignmentFragsGroupedList.get(direction);
            for (List<AlignmentFragment> compatibleFragsList: generateCompatibleFragsChains(fragments, configuration)) {
                int size = 0;
                for (int i = 0; i < 2; i++) {
                    if (i == 0)
                        compatibleFragsList = mergeAlignmentFragments(compatibleFragsList, alignment.getVirusGenome(), minIntronLength, maxAlignMergeAAGap, alignment.getViralProtein());
                    if (i == 1) {
                        compatibleFragsList = mergeAlignmentFragments(compatibleFragsList, alignment.getVirusGenome(), relaxIntronLength, relaxMergeAAGap, alignment.getViralProtein());
                    }
                    if (size != compatibleFragsList.size()) {
                        size = compatibleFragsList.size();
                        Model model = new Model();
                        alignment.setAlignmentFragments(compatibleFragsList);
                        List<Exon> exons = determineVirusGenomeExons(compatibleFragsList);
                        model.getExons().addAll(exons);
                        model.setAlignment(alignment);
                        model.setGeneSymbol(alignment.getViralProtein().getGeneSymbol());
                        model.setDirection(direction);
                        model.getScores().putAll(alignment.getAlignmentScore());
                        model.getStatus().add("Initial Model");
                        List<String> notes = alignment.getViralProtein()
                                                      .getConfiguration()
                                                      .getOrDefault(ConfigurationParameters.Note, Collections.EMPTY_LIST);

                        model.getNotes().addAll(notes.stream().filter(s  -> ! s.isEmpty()).collect(Collectors.toList()));
                        models.add(model);
                    }
                }
            }
        }
        return models;
    }

    /**
     * @return merge two fragments if the intron length is <minIntronLength, if there are no stops in between and if intron length is divisible by 3.
     * Also merge fragments if the missing protein alignment length between two fragments is less than the max align merge aa gap;
     * @param:alignment
     */
    public List<AlignmentFragment> mergeAlignmentFragments ( List<AlignmentFragment> fragments, VirusGenome virusGenome, long minIntronLength, long maxAlignMergeAAGap, ViralProtein viralProtein ) {

        fragments.sort(AlignmentFragment.Comparators.Ascending);
        List<AlignmentFragment> outFragments = new ArrayList<AlignmentFragment>();
        StopTranslationException stopTransExce = viralProtein.getGeneAttributes().getStopTranslationException();
        List<Range> sequenceGaps = virusGenome.getSequenceGaps();
        boolean isPreMerge = false;
        for (int i = 0; i < fragments.size(); i++) {
            if (i == fragments.size() - 1) {
                if (!isPreMerge) {
                    outFragments.add(fragments.get(i));
                }

            } else {
                AlignmentFragment upFragment = fragments.get(i);
                AlignmentFragment downFragment = fragments.get(i + 1);
                Range currentFragment = upFragment.getNucleotideSeqRange();
                Range nextFragment = downFragment.getNucleotideSeqRange();
                Range intronRange;
                if (currentFragment.getEnd() - nextFragment.getBegin() == 0) {
                    intronRange = Range.ofLength(0);
                } else {
                    if (currentFragment.getEnd() > nextFragment.getBegin()) {
                        intronRange = currentFragment.intersection(nextFragment);
                    } else {
                        intronRange = Range.of(currentFragment.getEnd() + 1, nextFragment.getBegin() - 1);
                    }
                }
                Range missingAAalignRange;
                if (upFragment.getProteinSeqRange().getEnd() - downFragment.getProteinSeqRange().getBegin() == 0) {
                    missingAAalignRange = Range.ofLength(0);
                } else {
                    if (downFragment.getProteinSeqRange().getBegin() - 1 > upFragment.getProteinSeqRange().getEnd() + 1) {
                        missingAAalignRange = Range.of(upFragment.getProteinSeqRange().getEnd() + 1, downFragment.getProteinSeqRange().getBegin() - 1);
                    } else {
                        missingAAalignRange = Range.of(downFragment.getProteinSeqRange().getBegin() - 1, upFragment.getProteinSeqRange().getEnd() + 1);
                    }
                }

                if ((intronRange.getLength() <= minIntronLength && missingAAalignRange.getLength() <= maxAlignMergeAAGap && !VigorFunctionalUtils.intheSequenceGap(sequenceGaps, intronRange))) {
                    Map<Frame, List<Long>> intronStops = VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, intronRange);
                    List<Long> upStops = new ArrayList<Long>();
                    List<Long> downStops = new ArrayList<Long>();
                    Frame upSeqFrame = VigorFunctionalUtils.getSequenceFrame(upFragment.getNucleotideSeqRange().getBegin() + upFragment.getFrame().getFrame() - 1);
                    if (intronStops.get(upSeqFrame) != null) {
                        upStops = intronStops.get(upSeqFrame);
                    }
                    Frame downSeqFrame = VigorFunctionalUtils.getSequenceFrame(upFragment.getNucleotideSeqRange().getBegin() + upFragment.getFrame().getFrame() - 1);
                    if (intronStops.get(downSeqFrame) != null) {
                        downStops = intronStops.get(downSeqFrame);
                    }
                    if (stopTransExce.isHasStopTranslationException()) {
                        List<Range> matches = virusGenome.getSequence().findMatches(stopTransExce.getMotif()).distinct().collect(Collectors.toList());
                        int offset = stopTransExce.getOffset();
                        if (offset < 0) {
                            offset = offset + 1;
                        }
                        for (Range match : matches) {
                            long start = match.getEnd() + offset;
                            upStops.remove(start);
                            downStops.remove(start);
                        }
                    }

                    if (upStops.isEmpty() && downStops.isEmpty() && (intronRange.getLength() % 3 == 0)) {
                        if (isPreMerge) {
                            upFragment = outFragments.get(outFragments.size() - 1);
                            outFragments.remove(upFragment);
                            currentFragment = upFragment.getNucleotideSeqRange();
                        }
                        Range adjustedNTrange = Range.of(currentFragment.getBegin(), nextFragment.getEnd());
                        Range adjustedAArange = Range.of(upFragment.getProteinSeqRange().getBegin(), downFragment.getProteinSeqRange().getEnd());
                        AlignmentFragment adjustedFragment = new AlignmentFragment(adjustedAArange, adjustedNTrange, upFragment.getDirection(), upFragment.getFrame());
                        outFragments.add(adjustedFragment);
                        isPreMerge = true;
                    } else {
                        if (!isPreMerge) {
                            outFragments.add(fragments.get(i));
                        }
                        isPreMerge = false;
                    }
                } else {
                    if (!isPreMerge) {
                        outFragments.add(fragments.get(i));
                    }
                    isPreMerge = false;
                }
            }
        }

        if (outFragments.isEmpty()) {
            outFragments.addAll(fragments);
        }
        return outFragments;
    }

    /**
     * @param model
     * @param sequenceGaps
     * @return
     */
    public Model splitExonsAtSequenceGaps ( Model model, List<Range> sequenceGaps ) {

        List<Exon> exons = model.getExons();
        List<Exon> newExons = new ArrayList<Exon>();
        for (Exon exon : exons) {
            for (Range sequenceGap: sequenceGaps ) {
                if (sequenceGap.isSubRangeOf(exon.getRange())) {
                    Exon firstExon = new Exon();
                    firstExon.setRange(Range.of(exon.getRange().getBegin(), sequenceGap.getBegin() - 1));
                    firstExon.setAlignmentFragment(exon.getAlignmentFragment());
                    firstExon.set_3p_adjusted(true);
                    firstExon.setFrame(exon.getFrame());
                    newExons.add(firstExon);
                    int remainder = Math.max(3 - (int) ( sequenceGap.getEnd() - firstExon.getRange().getBegin() + 1 ) % 3,0);
                    exon.setFrame(firstExon.getFrame().shift(remainder));
                    exon.setRange(Range.of(sequenceGap.getEnd() + 1, exon.getRange().getEnd()));
                    exon.set_5p_adjusted(true);
                } else if (sequenceGap.intersects(exon.getRange())) {
                    Range intersection = sequenceGap.intersection(exon.getRange());
                    if (intersection.getBegin() == exon.getRange().getBegin()) {
                        int reminder = Math.max(3 - (int) ( intersection.getLength() ) % 3,0);
                        exon.setFrame(exon.getFrame().shift(reminder));
                        exon.setRange(Range.of(intersection.getEnd() + 1, exon.getRange().getEnd()));
                        exon.set_5p_adjusted(true);
                    } else if (intersection.getEnd() == exon.getRange().getEnd()) {
                        exon.setRange(Range.of(exon.getRange().getBegin(), intersection.getBegin() - 1));
                        exon.set_3p_adjusted(true);
                    }
                }
            }
        }
        exons.addAll(newExons);
        exons.sort(Exon.Comparators.Ascending);
        model.setExons(exons);
        return model;
    }

    /**
     * @param: initModels
     * @param: validSequenceGaps
     * @return Models are split at sequence gaps and new list of models are
     * returned
     */
    public List<Model> splitModelAtSequenceGaps ( Model initModel, List<Range> validSequenceGaps ) throws CloneNotSupportedException {

        initModel = splitExonsAtSequenceGaps(initModel, validSequenceGaps);
        List<Model> newModels = new ArrayList<Model>();
        Model model = initModel.clone();
        if (validSequenceGaps.isEmpty()) {
            newModels.add(initModel);
            return newModels;
        }

        Exon nextExon;
        Exon currentExon;
        Range diffRange;
        List<Exon> firstGroup = new ArrayList<Exon>();
        List<Exon> secondGroup = new ArrayList<Exon>();
        Model firstModel;
        Model secondModel = null;
        List<Exon> modelExons = model.getExons();
        secondGroup.addAll(modelExons);
        boolean startExist = true;
        for (int j = 0; j < modelExons.size(); j++) {
            if (j == modelExons.size() - 1) {
                Exon lastExon = modelExons.get(j);
                long proteinSeqLength = model.getAlignment().getViralProtein().getSequence().getLength();
                Range lastExonAARange = lastExon.getAlignmentFragment().getProteinSeqRange();
                long expectedStart = lastExon.getRange().getEnd();
                if (lastExonAARange.getEnd() < proteinSeqLength - 1) {
                    long difference = ( proteinSeqLength - 1 ) - lastExonAARange.getEnd();
                    expectedStart = lastExon.getRange().getEnd() + difference * 3;
                }
                Range searchRangeTemp = Range.of(lastExon.getRange().getEnd(), expectedStart);
                for (Range gap : validSequenceGaps) {
                    if (gap.intersects(searchRangeTemp)) {
                        model.setPartial3p(true);
                        lastExon.setRange(Range.of(lastExon.getRange().getBegin(), gap.getBegin() - 1));
                        break;
                    }
                }
            } else {
                nextExon = modelExons.get(j + 1);
                currentExon = modelExons.get(j);
                if (nextExon.getRange().getBegin() > currentExon.getRange().getEnd()) {
                    diffRange = Range.of(currentExon.getRange().getEnd(), nextExon.getRange().getBegin());
                    if (diffRange.getLength() >= 20) {
                        firstGroup.add(modelExons.get(j));
                        boolean temp = true;
                        for (int k = 0; k < validSequenceGaps.size(); k++) {
                            if (diffRange.intersects(validSequenceGaps.get(k)) && temp) {
                                secondModel = model.clone();
                                firstModel = model.clone();
                                if (currentExon.getRange().getEnd() < validSequenceGaps.get(k).getBegin() - 1) {
                                    currentExon.setRange(Range.of(currentExon.getRange().getBegin(), validSequenceGaps.get(k).getBegin() - 1));
                                    currentExon.set_3p_adjusted(true);
                                }
                                if (nextExon.getRange().getBegin() > validSequenceGaps.get(k).getEnd() + 1) {
                                    Range seqGapRange = VigorFunctionalUtils.get5pNearestSequenceGap(validSequenceGaps, nextExon.getRange());
                                    int reminder = (int) ( ( nextExon.getRange().getBegin() - 1 - ( seqGapRange.getEnd() + 1 ) ) + 1 ) % 3;
                                    nextExon.setRange(Range.of(seqGapRange.getEnd() + 1, nextExon.getRange().getEnd()));
                                    if (reminder > 0) {
                                        nextExon.setFrame(currentExon.getFrame().shift(reminder));
                                    }
                                    nextExon.set_5p_adjusted(true);
                                }
                                List<Exon> tempFirst = new ArrayList<>();
                                tempFirst.addAll(firstGroup);
                                List<Exon> tempSecond = new ArrayList<>();
                                tempSecond.addAll(secondGroup);
                                firstModel.setExons(tempFirst);
                                firstModel.setPartial3p(true);
                                firstModel.addNote(NoteType.Sequence_Gap);

                                firstModel.getStatus().add("Model split at sequence gaps");
                                if (!startExist) {
                                    firstModel.setPartial5p(true);
                                }
                                secondGroup.removeAll(tempFirst);
                                tempSecond.removeAll(tempFirst);
                                secondModel.setExons(tempSecond);
                                secondModel.setPartial5p(true);
                                secondModel.addNote(NoteType.Sequence_Gap);
                                secondModel.getStatus().add("Model split at sequence gaps");
                                newModels.add(firstModel);
                                startExist = false;
                                firstGroup.clear();
                                temp = false;
                            }
                        }
                    }
                }
            }
        }
        if (secondModel != null && secondModel.getExons().size() > 0) {
            newModels.add(secondModel);
        }
        if (newModels.size() == 0) {
            newModels.add(model);
        }

        return newModels;
    }


    /**
     * @param fragments : alignment fragments that are grouped based on direction
     * @return List<List<AlignmentFragment>>: compatible(check for overlap) list
     * of alignment fragments and their permutations and combinations are
     * grouped
     */
    public List<List<AlignmentFragment>> generateCompatibleFragsChains(List<AlignmentFragment> fragments, VigorConfiguration configuration) {
        int ntOverlap = configuration.getOrDefault(ConfigurationParameters.NTOverlapMaximum, DEFAULT_NTOVERLAP_OFFSET);
        int aaOverlap = configuration.getOrDefault(ConfigurationParameters.AAOverlapMaximum, DEFAULT_AAOVERLAP_OFFSET);

        BiFunction<AlignmentFragment,AlignmentFragment, Boolean> areCompatible = (a, b) ->
                (! a.equals(b)) &&
                (a.getNucleotideSeqRange().getEnd() - ntOverlap) < b.getNucleotideSeqRange().getBegin() &&
                        (a.getProteinSeqRange().getEnd() - aaOverlap) < b.getProteinSeqRange().getBegin();

        List<List<AlignmentFragment>> compatibleFragmentList = new ArrayList<>();
        for (AlignmentFragment fragment: fragments) {
            // starting fragments have no compatible fragments downstream
            if (! fragments.stream()
                           .filter(Predicates.not(fragment::equals))
                           .anyMatch(f -> areCompatible.apply(f,fragment))) {
                compatibleFragmentList.addAll(generateCompatibleFragsList(fragment, fragments, areCompatible));
            }
        }
        return compatibleFragmentList;
    }

    public List<List<AlignmentFragment>> generateCompatibleFragsList(AlignmentFragment currentFragment,
                                                                     List<AlignmentFragment> fragments,
                                                                     BiFunction<AlignmentFragment,AlignmentFragment, Boolean> areCompatible) {


        List<List<AlignmentFragment>> compatibleFragments = new ArrayList<>();

        for (int i=0;i<fragments.size(); i++) {
            AlignmentFragment nextFragment = fragments.get(i);
            if (areCompatible.apply(currentFragment, nextFragment)) {
                // only add consider this fragment by itself if it's not included in some other chain
                if (fragments.stream()
                             .filter(f -> ! (f.equals(currentFragment)|| f.equals(nextFragment)) )
                             .anyMatch(f -> areCompatible.apply(currentFragment, f) && areCompatible.apply(f,nextFragment))) {
                    continue;
                }
                compatibleFragments.addAll(generateCompatibleFragsList(nextFragment,
                                                                       fragments.subList(i+1, fragments.size()),
                                                                       areCompatible));
            }
        }

        if (compatibleFragments.isEmpty()) {
            compatibleFragments.add(new ArrayList<>());
        }

        for (List<AlignmentFragment> compatibleList: compatibleFragments) {
            compatibleList.add(0,currentFragment);
        }
        return compatibleFragments;
    }

    /**
     * @param compatibleFragsList
     * @return List<Exon>: nucleotideSequence range of alignment fragment will
     * be set as range of exon. 5' and 3' edges of exon are not
     * modified.
     */
    public List<Exon> determineVirusGenomeExons ( List<AlignmentFragment> compatibleFragsList ) {

        Iterator<AlignmentFragment> iter = compatibleFragsList.iterator();
        List<Exon> exons = new ArrayList<Exon>();
        while (iter.hasNext()) {
            Exon exon = new Exon();
            AlignmentFragment alignmentFragment = iter.next();
            exon.setRange(alignmentFragment.getNucleotideSeqRange());
            exon.setAlignmentFragment(alignmentFragment);
            exon.setFrame(alignmentFragment.getFrame());
            Frame sequenceFrame = VigorFunctionalUtils.getSequenceFrame(exon.getRange().getBegin() + exon.getFrame().getFrame() - 1);
            exon.setSequenceFrame(sequenceFrame);
            exons.add(exon);
        }
        return exons;
    }
}
