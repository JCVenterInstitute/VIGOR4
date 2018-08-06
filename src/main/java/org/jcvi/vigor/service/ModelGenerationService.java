package org.jcvi.vigor.service;

import org.jcvi.vigor.component.*;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.springframework.stereotype.Service;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by snettem on 5/16/2017.
 */
@Service
public class ModelGenerationService {

    private long minCondensation = 10;
    private long minIntronLength = 30;
    private long relaxIntronLength = 900;
    private long relaxCondensation = 300;
    private static boolean isDebug = false;
    private int AAOverlapOffset = 10;
    private int NTOverlapOffset = 30;
    private static final Logger LOGGER = LogManager.getLogger(ModelGenerationService.class);

    public List<Model> generateModels ( List<Alignment> alignments, VigorForm form ) throws ServiceException {

        VigorConfiguration configuration = form.getConfiguration();
        isDebug = configuration.getOrDefault(ConfigurationParameters.Verbose, false);
        AAOverlapOffset = configuration.getOrDefault(ConfigurationParameters.AAOverlap_offset, AAOverlapOffset);
        NTOverlapOffset = configuration.getOrDefault(ConfigurationParameters.NTOverlap_offset, NTOverlapOffset);
        minCondensation = configuration.getOrDefault(ConfigurationParameters.CondensationMinimum, minCondensation);
        minIntronLength = minCondensation * 3;
        alignments = mergeIdenticalProteinAlignments(alignments);
        return determineCandidateModels(alignments, form);
    }

    public List<Alignment> mergeIdenticalProteinAlignments ( List<Alignment> alignments ) {

        List<Alignment> allOutAlignments = new ArrayList<Alignment>();
        Map<String, List<Alignment>> protein2AlignmentsMap = new HashMap<String, List<Alignment>>();
        for (Alignment alignment : alignments) {
            alignment.getAlignmentFragments().sort(AlignmentFragment.Comparators.Ascending);
            String proteinID = alignment.getViralProtein().getProteinID();
            if (protein2AlignmentsMap.containsKey(proteinID)) {
                List<Alignment> existingAlignments = protein2AlignmentsMap.get(proteinID);
                existingAlignments.add(alignment);
                Collections.sort(existingAlignments, ( alignment1, alignment2 )
                        -> alignment1.getAlignmentFragments().get(alignment1.getAlignmentFragments().size() - 1).compareTo(alignment2.getAlignmentFragments().get(0)));
                protein2AlignmentsMap.put(proteinID, existingAlignments);
            } else {
                protein2AlignmentsMap.put(proteinID, new ArrayList<>(Arrays.asList(alignment)));
            }
        }
        for (String proteinID : protein2AlignmentsMap.keySet()) {
            List<Alignment> alignmentsTemp = protein2AlignmentsMap.get(proteinID);
            Alignment mergedAlignment = alignmentsTemp.get(0);
            alignmentsTemp.remove(mergedAlignment);
            for (Alignment alignment : alignmentsTemp) {
                mergedAlignment.getAlignmentFragments().addAll(alignment.getAlignmentFragments());
                Map<String, Double> scores = mergedAlignment.getAlignmentScore();
                double score = scores.get("alignmentScore");
                score = score + alignment.getAlignmentScore().get("alignmentScore");
                scores.put("alignmentScore", score);
                mergedAlignment.setAlignmentScore(scores);
            }
            allOutAlignments.add(mergedAlignment);
        }
        return allOutAlignments;
    }

    /**
     * @param alignments
     * @param form
     * @return all the possible models of each alignment and after splitting
     * models at the sequence gaps
     */
    public List<Model> determineCandidateModels ( List<Alignment> alignments, VigorForm form ) throws ServiceException {

        List<Model> initialModels = new ArrayList<Model>();
        List<Model> candidateModels = new ArrayList<Model>();
        for (int i = 0; i < alignments.size(); i++) {
            Alignment alignment = alignments.get(i);
            alignment.getAlignmentFragments().sort(AlignmentFragment.Comparators.Ascending);
            initialModels.addAll(alignmentToModels(alignment));
        }
        if (isDebug) {
            FormatVigorOutput.printModels(initialModels, "Initial Models");
        }
        List<Range> sequenceGaps = new ArrayList<Range>();
        // get sequence gaps
        if (( !initialModels.isEmpty() ) && initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps() != null) {
            sequenceGaps = initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps();
        }
        List<Range> validSequenceGaps = new ArrayList<Range>();
        if (sequenceGaps != null && sequenceGaps.size() > 0) {
            long minGapLength = form.getConfiguration().getOrDefault(ConfigurationParameters.SequenceGapMinimumLength, 0);
            for (Range gapRange : sequenceGaps) {
                if (gapRange.getLength() >= minGapLength) {
                    validSequenceGaps.add(gapRange);
                }
            }
        }
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
    public List<Model> alignmentToModels ( Alignment alignment ) {

        Map<Direction, List<AlignmentFragment>> alignmentFragsGroupedList = alignment.getAlignmentFragments().stream()
                .collect(Collectors.groupingBy(w -> w.getDirection()));
        Set<Direction> keyset = alignmentFragsGroupedList.keySet();
        Iterator<Direction> iter = keyset.iterator();
        List<Model> models = new ArrayList<Model>();
        for (Direction direction : keyset) {
            List<List<AlignmentFragment>> ListOfCompatibleFragsList = generateCompatibleFragsChains(
                    alignmentFragsGroupedList.get(iter.next()));
            Iterator<List<AlignmentFragment>> iter1 = ListOfCompatibleFragsList.iterator();
            while (iter1.hasNext()) {
                List<AlignmentFragment> compatibleFragsList = iter1.next();
                int size = 0;
                for (int i = 0; i < 2; i++) {
                    if (i == 0)
                        compatibleFragsList = mergeAlignmentFragments(compatibleFragsList, alignment.getVirusGenome(), minIntronLength, minCondensation, alignment.getViralProtein());
                    if (i == 1) {
                        compatibleFragsList = compatibleFragsList.stream().map(x -> x.clone()).collect(Collectors.toList());
                        compatibleFragsList = mergeAlignmentFragments(compatibleFragsList, alignment.getVirusGenome(), relaxIntronLength, relaxCondensation, alignment.getViralProtein());
                    }
                    if (size != compatibleFragsList.size()) {
                        size = compatibleFragsList.size();
                        alignment.setAlignmentFragments(compatibleFragsList);
                        List<Exon> exons = determineVirusGenomeExons(compatibleFragsList);
                        Model model = new Model();
                        List<NoteType> notes = new ArrayList<>();
                        model.setNotes(notes);
                        Map<String, Double> scores = new HashMap<String, Double>();
                        model.setScores(scores);
                        model.setExons(exons);
                        model.setAlignment(alignment);
                        model.setGeneSymbol(alignment.getViralProtein().getGeneSymbol());
                        model.setDirection(direction);
                        Map<String, Double> alignmentScores = alignment.getAlignmentScore();
                        model.setScores(alignmentScores);
                        List<String> statusList = new ArrayList<String>();
                        statusList.add("Initial Model");
                        model.setStatus(statusList);
                        models.add(model);
                    }
                }
            }
        }
        return models;
    }

    /**
     * @return merge two fragments if the intron length is <minIntronLength and if there are no stops in between and if intron length is divisible by 3.
     * Also merge fragments if the missing protein alignment length between two fragments is less than the minimum_condensation;
     * @param:alignment
     */
    public List<AlignmentFragment> mergeAlignmentFragments ( List<AlignmentFragment> fragments, VirusGenome virusGenome, long thisMinIntronLength, long thisMinCondensation, ViralProtein viralProtein ) {

        fragments.sort(AlignmentFragment.Comparators.Ascending);
        List<AlignmentFragment> outFragments = new ArrayList<AlignmentFragment>();
        StopTranslationException stopTransExce = viralProtein.getGeneAttributes().getStopTranslationException();
        List<Range> sequenceGaps = virusGenome.getSequenceGaps();
        boolean isPreMerge = false;
        if (fragments.size() == 1) {
            outFragments.add(fragments.get(0));
        } else if (fragments.size() > 1) {
            for (int i = 0; i < fragments.size(); i++) {
                if (i != fragments.size() - 1) {
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
                    if (( intronRange.getLength() <= thisMinIntronLength && missingAAalignRange.getLength() <= thisMinCondensation && !VigorFunctionalUtils.intheSequenceGap(sequenceGaps, intronRange) )) {
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
                                if (upStops.contains(start)) upStops.remove(start);
                                if (downStops.contains(start)) downStops.remove(start);
                            }
                        }
                        if (upStops.size() == 0 && downStops.size() == 0 && ( intronRange.getLength() % 3 == 0 )) {
                            if (isPreMerge) {
                                upFragment = outFragments.get(outFragments.size() - 1);
                                outFragments.remove(upFragment);
                                currentFragment = upFragment.getNucleotideSeqRange();
                            }
                            Range adjustedNTrange = Range.of(currentFragment.getBegin(), nextFragment.getEnd());
                            Range adjustedAArange = Range.of(upFragment.getProteinSeqRange().getBegin(), downFragment.getProteinSeqRange().getEnd());
                            upFragment.setNucleotideSeqRange(adjustedNTrange);
                            upFragment.setProteinSeqRange(adjustedAArange);
                            outFragments.add(upFragment);
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
                } else {
                    if (!isPreMerge) {
                        outFragments.add(fragments.get(i));
                    }
                }
            }
        }
        if (outFragments.size() == 0) {
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
            for (int i = 0; i < sequenceGaps.size(); i++) {
                if (sequenceGaps.get(i).isSubRangeOf(exon.getRange())) {
                    Exon firstExon = new Exon();
                    firstExon.setRange(Range.of(exon.getRange().getBegin(), sequenceGaps.get(i).getBegin() - 1));
                    firstExon.setAlignmentFragment(exon.getAlignmentFragment());
                    firstExon.set_3p_adjusted(true);
                    firstExon.setFrame(exon.getFrame());
                    newExons.add(firstExon);
                    int reminder = (int) ( sequenceGaps.get(i).getEnd() - firstExon.getRange().getBegin() + 1 ) % 3;
                    if (reminder > 0) {
                        reminder = 3 - reminder;
                        exon.setFrame(firstExon.getFrame().shift(reminder));
                    } else exon.setFrame(firstExon.getFrame());
                    exon.setRange(Range.of(sequenceGaps.get(i).getEnd() + 1, exon.getRange().getEnd()));
                    exon.set_5p_adjusted(true);
                } else if (sequenceGaps.get(i).intersects(exon.getRange())) {
                    Range intersection = sequenceGaps.get(i).intersection(exon.getRange());
                    if (intersection.getBegin() == exon.getRange().getBegin()) {
                        int reminder = (int) ( intersection.getLength() ) % 3;
                        if (reminder > 0) {
                            reminder = 3 - reminder;
                            exon.setFrame(exon.getFrame().shift(reminder));
                        }
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
     * @return Models are split at sequence gaps and new list of models are
     * returned
     * @param:initModels
     * @param:form
     * @param:genome
     */
    public List<Model> splitModelAtSequenceGaps ( Model initModel, List<Range> validSequenceGaps ) throws CloneNotSupportedException {

        initModel = splitExonsAtSequenceGaps(initModel, validSequenceGaps);
        List<Model> newModels = new ArrayList<Model>();
        Model model = initModel.clone();
        if (validSequenceGaps.size() > 0) {
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
                }
                if (j != modelExons.size() - 1) {
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
                                    List<NoteType> fnotes = firstModel.getNotes();
                                    fnotes.add(NoteType.Sequence_Gap);
                                    firstModel.setNotes(fnotes);
                                    firstModel.getStatus().add("Model split at sequence gaps");
                                    if (!startExist) {
                                        firstModel.setPartial5p(true);
                                    }
                                    secondGroup.removeAll(tempFirst);
                                    tempSecond.removeAll(tempFirst);
                                    secondModel.setExons(tempSecond);
                                    secondModel.setPartial5p(true);
                                    List<NoteType> snotes = secondModel.getNotes();
                                    snotes.add(NoteType.Sequence_Gap);
                                    secondModel.setNotes(snotes);
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
        } else {
            newModels.add(initModel);
            return newModels;
        }
    }

    /**
     * @param alignmentFragments : alignment fragments that are grouped based on direction is
     *                           input to the function
     * @return List<List   <   AlignmentFragment>>: compatible(check for overlap) list
     * of alignment fragments and their permutations and combinations is
     * grouped
     */
    public List<List<AlignmentFragment>> generateCompatibleFragsChains ( List<AlignmentFragment> alignmentFragments ) {

        List<AlignmentFragment> compatibleFragsList = new ArrayList<AlignmentFragment>();
        List<AlignmentFragment> clonedCompatibleFragsList;
        List<List<AlignmentFragment>> ListOfCompatibleFragsList = new ArrayList<List<AlignmentFragment>>();
        compatibleFragsList.add(alignmentFragments.get(0).clone());
        ListOfCompatibleFragsList.add(compatibleFragsList);
        List<List<AlignmentFragment>> tempList;
        if (alignmentFragments.size() > 1) {
            for (int j = 1; j < alignmentFragments.size(); j++) {
                boolean temp = true;
                tempList = new ArrayList<>();
                for (int k = 0; k < ListOfCompatibleFragsList.size(); k++) {
                    List<AlignmentFragment> currentList = ListOfCompatibleFragsList.get(k);
                    AlignmentFragment currentFrag = currentList.get(currentList.size() - 1);
                    long NTEnd = currentFrag.getNucleotideSeqRange().getEnd();
                    long AAEnd = currentFrag.getProteinSeqRange().getEnd();
                    long nextNTStart = alignmentFragments.get(j).getNucleotideSeqRange().getBegin();
                    long nextAAStart = alignmentFragments.get(j).getProteinSeqRange().getBegin();
                    if (nextNTStart >= NTEnd - NTOverlapOffset && nextAAStart >= AAEnd - AAOverlapOffset) {
                        currentList.add(alignmentFragments.get(j).clone());
                    } else if (temp) {
                        clonedCompatibleFragsList = new ArrayList<>();
                        for (int i = 0; i < currentList.size() - 1; i++) {
                            clonedCompatibleFragsList.add(currentList.get(i).clone());
                        }
                        clonedCompatibleFragsList.add(alignmentFragments.get(j).clone());
                        tempList.add(clonedCompatibleFragsList);
                        if (clonedCompatibleFragsList.size() == 1) {
                            temp = false;
                        }
                    }
                }
                if (tempList.size() > 0) {
                    ListOfCompatibleFragsList.addAll(tempList);
                }
            }
        }
        return ListOfCompatibleFragsList;
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
