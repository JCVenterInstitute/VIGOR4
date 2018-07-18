package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.stereotype.Service;

@Service
public class DetermineStop implements DetermineGeneFeatures {


    @Override
    public List<Model> determine ( Model model, VigorForm form ) throws ServiceException {

        Integer stopCodonWindow = form.getConfiguration().getOrDefault(ConfigurationParameters.StopCodonSearchWindow, 50);
        boolean isDebug = form.getConfiguration().get(ConfigurationParameters.Verbose).equals("true") ? true : false;
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
        if (lastExonAARange.getEnd() < proteinSeqLength - 1) {
            long difference = proteinSeqLength - 1 - lastExonAARange.getEnd();
            expectedStart = lastExon.getRange().getEnd() + difference * 3;
        } else {
            expectedStart = lastExon.getRange().getEnd();
        }
        Frame lastExonFrame = VigorFunctionalUtils.getSequenceFrame(lastExon.getRange().getBegin() + ( lastExon.getFrame().getFrame() - 1 ));
        start = expectedStart - stopCodonWindow;
        end = expectedStart + stopCodonWindow;
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
            if (!( rangeScoreMap.isEmpty() && rangeScoreMap.size() > 0 )) {
                Set<Range> keys = rangeScoreMap.keySet();
                for (Range range : keys) {
                    Model newModel = model.clone();
                    List<Exon> newExons = newModel.getExons();
                    Exon lExon = newExons.get(newExons.size() - 1);
                    lExon.setRange(Range.of(lExon.getRange().getBegin(), range.getBegin() + 2));
                    newModel.getScores().put("stopCodonScore",
                            rangeScoreMap.get(range));
                    newModels.add(newModel);
                }
            }
        }
        if (rangeScoreMap.isEmpty() && isSequenceMissing) {
            Model newModel = model.clone();
            newModel.setPartial3p(true);
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
