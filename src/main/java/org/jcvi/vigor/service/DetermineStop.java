package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.springframework.stereotype.Service;

@Service
public class DetermineStop implements DetermineGeneFeatures {

    private static Logger LOGGER = LogManager.getLogger(DetermineStop.class);
    @Override
    public List<Model> determine ( Model model ) throws ServiceException {

        VigorConfiguration config = model.getAlignment().getViralProtein().getConfiguration();
        Integer stopCodonWindow = config.getOrDefault(ConfigurationParameters.StopCodonSearchWindow, 50);
        boolean isDebug = config.get(ConfigurationParameters.Verbose).equals("true") ? true : false;
        try {
            return findStop(model, stopCodonWindow, isDebug);
        } catch (CloneNotSupportedException e) {
            throw new ServiceException(String.format("Problem determine stop for model %s", model), e);
        }
    }

    /**
     * @param model
     * @return
     * @throws CloneNotSupportedException
     */
    @SuppressWarnings("Duplicates")
    public List<Model> findStop ( Model model, int stopCodonWindow, boolean isDebug )
            throws CloneNotSupportedException {

        List<Model> newModels = new ArrayList<Model>();
        model.getExons().sort(Exon.Comparators.Ascending);
        List<Exon> exons = model.getExons();
        long start;
        long end;
        Range stopSearchRange;
        NucleotideSequence seq = model.getAlignment().getVirusGenome().getSequence();
        long proteinSeqLength = model.getAlignment().getViralProtein().getSequence().getLength();
        boolean isSequenceMissing = false;
        Exon lastExon = exons.get(exons.size() - 1);
        Range lastExonAARange = lastExon.getAlignmentFragment().getProteinSeqRange();
        long expectedStart;
        // Define search window and search for stops in the search window
        if (lastExonAARange.getEnd() < proteinSeqLength - 1) {
            long difference = proteinSeqLength - 1 - lastExonAARange.getEnd();
            expectedStart = lastExon.getRange().getEnd() + difference * 3;
        } else {
            expectedStart = lastExon.getRange().getEnd();
        }
        Frame lastExonFrame = VigorFunctionalUtils.getSequenceFrame(lastExon.getRange().getBegin() + ( lastExon.getFrame().getFrame() - 1 ));
        start = expectedStart - stopCodonWindow;
        end = expectedStart + stopCodonWindow;
        //if search window defined above lies outside the sequence, then set isSequenceMissing to true.
        if (end > seq.getLength() - 1 || start > seq.getLength() - 1) {
            isSequenceMissing = true;
            end = seq.getLength() - 1;
            start = end - stopCodonWindow;
        }
        if (start < lastExon.getRange().getBegin()) {
            start = lastExon.getRange().getBegin();
        }
        stopSearchRange = Range.of(start, end);
        Map<Range, Double> rangeScoreMap = new HashMap<Range, Double>();
        Map<Frame, List<Long>> stops = model.getAlignment().getVirusGenome().getInternalStops();
        List<Long> stopsInFrame = new ArrayList<Long>();
        if (stops != null) {
            for (Frame frame : stops.keySet()) {
                if (frame.equals(lastExonFrame)) {
                    List<Long> tempList = stops.get(lastExonFrame);
                    for (Long stop : tempList) {
                        Range tempRange = Range.of(stop);
                        if (tempRange.isSubRangeOf(stopSearchRange)) {
                            stopsInFrame.add(stop);
                        }
                    }
                }
            }
        }
        // List all stops in frame and assign a score for each match (match closer to expected start scores high)
        if (stopsInFrame != null && stopsInFrame.size() > 0) {
            for (Long stop : stopsInFrame) {
                rangeScoreMap.put(Range.of(stop, stop + 2), VigorFunctionalUtils.generateProximityScore(lastExon.getRange().getEnd(), stop));
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
                List<Exon> newExons = newModel.getExons();
                Exon lExon = newExons.get(newExons.size() - 1);
                lExon.setRange(Range.of(lExon.getRange().getBegin(), range.getBegin() + 2));
                newModel.getScores().put(Scores.STOP_CODON_SCORE,
                                         rangeScoreMap.get(range));
                newModels.add(newModel);
            }
        }
        //only if stop codon is not found and search window lies outside the sequence (ie: isSequenceMissing=true)
        if (rangeScoreMap.isEmpty() && isSequenceMissing) {
            Model newModel = model.clone();
            newModel.setPartial3p(true);
            newModel.setPseudogene(false);
            //extend the last exon till the end of sequence
            Exon lExon = newModel.getExons().get(newModel.getExons().size() - 1);
            Range lExonRange = lExon.getRange();
            lExon.setRange(Range.of(lExonRange.getBegin(), ( seq.getLength() - 1 )));
            newModels.add(newModel);
        } else if (rangeScoreMap.isEmpty()) {
            Model newModel = model.clone();
            newModel.setPseudogene(true);
            newModels.add(newModel);
        }
        return newModels;
    }
}
