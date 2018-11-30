package org.jcvi.vigor.service;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.Triplet;
import org.jcvi.vigor.service.exception.ServiceException;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.*;

@Service
public class DetermineStart implements DetermineGeneFeatures {

    private static final Logger LOGGER = LogManager
            .getLogger(DetermineStart.class);

    @Override
    public List<Model> determine (Model model) throws ServiceException {

        VigorConfiguration config = model.getAlignment().getViralProtein().getConfiguration();
        List<Triplet> startCodons = loadStartCodons(model.getAlignment()
                        .getViralProtein().getGeneAttributes()
                        .getStartTranslationException().getAlternateStartCodons(),
                config);
        Integer startCodonWindowParam = config.getOrDefault(ConfigurationParameters.StartCodonSearchWindow, 50);
        LOGGER.trace("determining starts using {}", () -> "TODO");

        try {
            List<Model> models = findStart(startCodons, model, startCodonWindowParam);
            return models;
        } catch (CloneNotSupportedException e) {
            LOGGER.error("for model {} problem finding start using codons {} and search window {}",
                    model,
                    startCodons.stream().map(String:: valueOf).collect(Collectors.joining(",")),
                    startCodonWindowParam);
            throw new ServiceException(e);
        }
    }

    /**
     * @param alternateStartCodons    defined in defline of reference protein
     * @param vigorParameters:Default start codons defined in vigor.ini file
     * @return List of all the possible start codons
     */
    public List<Triplet> loadStartCodons ( List<String> alternateStartCodons,
                                           VigorConfiguration vigorParameters ) {

        alternateStartCodons = alternateStartCodons != null ? alternateStartCodons: Collections.EMPTY_LIST;
        List<String> startCodonStrings = vigorParameters.getOrDefault(ConfigurationParameters.StartCodons,
                                                                      Arrays.asList("ATG"));
        return Stream.concat(startCodonStrings.stream(), alternateStartCodons.stream())
                     .map(s -> Triplet.create(s.charAt(0),s.charAt(1), s.charAt(2)))
                     .collect(Collectors.toList());
    }

    /**
     * @param startCodons
     * @param model
     * @param startCodonWindowParam
     * @return
     * @throws CloneNotSupportedException
     */
    @SuppressWarnings("Duplicates")
    public List<Model> findStart ( List<Triplet> startCodons, Model model,
                                   Integer startCodonWindowParam ) throws CloneNotSupportedException {

        List<Model> newModels = new ArrayList<Model>();
        List<Range> sequenceGaps = model.getAlignment().getVirusGenome().getSequenceGaps();
        long start;
        long end;
        Range startSearchRange;
        Range sequenceGapRange = Range.of(0);
        int windowSize = 50;
        boolean isSequenceMissing = false;
        boolean isSequenceGap = false;
        if (startCodonWindowParam != null) {
                windowSize = startCodonWindowParam;
        }
        Exon firstExon = model.getExons().get(0);
        Frame firstExonFrame = VigorFunctionalUtils.getSequenceFrame(firstExon.getRange().getBegin() + firstExon.getFrame().getFrame() - 1);
        long expectedStart = firstExon.getRange().getBegin()
                - ( ( firstExon.getAlignmentFragment().getProteinSeqRange()
                .getBegin() ) * 3 );
        Map<Range, Double> rangeScoreMap = new HashMap<>();
        if (( expectedStart < 0 && expectedStart > -windowSize ) || expectedStart >= 0) {
            expectedStart = expectedStart >= 0 ? expectedStart : 0;
            start = expectedStart - windowSize;
            //if search window defined above lies outside the sequence, then set isSequenceMissing to true.
            if (start < 0) {
                isSequenceMissing = true;
                start = 0;
            }
            end = expectedStart + windowSize;
            if (end > firstExon.getRange().getEnd()) {
                end = firstExon.getRange().getEnd();
            }
            startSearchRange = Range.of(start, end);
            //find any internal stops
            Map<Frame, List<Long>> internalStops = model.getAlignment().getVirusGenome().getInternalStops();
            List<Long> stopsInFrame = new ArrayList<>();
            if (internalStops != null) {
                for (Frame frame : internalStops.keySet()) {
                    if (frame.equals(firstExonFrame)) {
                        List<Long> tempList = internalStops.get(firstExonFrame);
                        for (Long stop : tempList) {
                            Range tempRange = Range.of(stop);
                            if (tempRange.isSubRangeOf(startSearchRange)) {
                                stopsInFrame.add(stop);
                            }
                        }
                    }
                }
            }
            // Do not allow sequence gaps in the search window
            for (int j = sequenceGaps.size() - 1; j >= 0; j--) {
                Range range = sequenceGaps.get(j);
                if (startSearchRange.intersects(range) || Range.of(startSearchRange.getEnd(), firstExon.getRange().getEnd()).intersects(range)) {
                    long endTemp = ( ( range.getEnd() + 1 + windowSize < firstExon.getRange().getBegin() + 2 ) ? firstExon.getRange().getBegin() + windowSize : range.getEnd() + 1 + windowSize );
                    startSearchRange = Range.of(range.getEnd() + 1, endTemp);
                    sequenceGapRange = range;
                    isSequenceGap = true;
                    break;
                }
            }
            final long tempStart = startSearchRange.getBegin();
            NucleotideSequence NTSequence = model.getAlignment().getVirusGenome()
                    .getSequence().toBuilder(startSearchRange).build();
            // List all starts in frame and assign a score for each match (match closer to expected start scores high)
            for (Triplet triplet : startCodons) {
                Stream<Range> stream = NTSequence.findMatches(triplet.toString());
                List<Range> rangesInFrame = new ArrayList<Range>();
                List<Range> foundRanges = stream.map(
                        range -> {
                            range = Range.of(range.getBegin() + tempStart,
                                    range.getEnd() + tempStart);
                            return range;
                        }).collect(Collectors.toList());
                foundRanges.stream().forEach(x -> {
                    if (VigorFunctionalUtils.getSequenceFrame(x.getBegin()).compareTo(firstExonFrame) == 0) {
                        rangesInFrame.add(x);
                    }
                });
                for (Range range : rangesInFrame) {
                    boolean isValid = true;
                    for (Long stop : stopsInFrame) {
                        if (range.endsBefore(Range.of(stop, stop + 2))) {
                            isValid = false;
                        }
                    }
                    if (isValid) {
                        rangeScoreMap.put(range, VigorFunctionalUtils.generateProximityScore(firstExon.getRange().getBegin(), range.getBegin()));
                    }
                }
            }
            rangeScoreMap
                    .entrySet()
                    .stream()
                    .sorted(Map.Entry.comparingByValue())
                    .collect(
                            Collectors.toMap(Map.Entry:: getKey,
                                    Map.Entry:: getValue, ( e1, e2 ) -> e2,
                                    LinkedHashMap::new));
            for (Range range : rangeScoreMap.keySet()) {
                Model newModel = model.clone();
                Exon fExon = newModel.getExons().get(0);
                fExon.setRange(Range.of(range.getBegin(), fExon.getRange().getEnd()));
                newModel.getScores().put(Scores.START_CODON_SCORE,
                                         rangeScoreMap.get(range));
                newModels.add(newModel);
            }
        } else isSequenceMissing = true;
        //set 5' partial and extend start of the first exon to the beginning of the sequence
        if (( rangeScoreMap.isEmpty() && isSequenceMissing ) || ( rangeScoreMap.isEmpty() && isSequenceGap )) {
            Model newModel = model.clone();
            newModel.setPartial5p(true);
            Exon fExon = newModel.getExons().get(0);
            Range fExonRange = fExon.getRange();
            long bases;
            if (isSequenceMissing) {
                bases = fExonRange.getBegin() - 0;
                fExon.setRange(Range.of(0, fExonRange.getEnd()));
            } else {
                bases = fExonRange.getBegin() - ( sequenceGapRange.getEnd() + 1 );
                fExon.setRange(Range.of(sequenceGapRange.getEnd() + 1, fExonRange.getEnd()));
            }
            int frameshift = (int) bases % 3;
            if (frameshift > 0) fExon.setFrame(fExon.getFrame().shift(frameshift));
            newModels.add(newModel);
        } else if (rangeScoreMap.isEmpty()) {
            Model newModel = model.clone();
            newModel.setPseudogene(true);
            newModels.add(newModel);
        }
        return newModels;
    }
}



