package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.component.Splicing.SpliceSite;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.stereotype.Service;

@Service
public class AdjustUneditedExonBoundaries implements DetermineGeneFeatures {

    private long defaultSearchWindow = 50;
    private long minIntronLength = 20;

    @Override
    public List<Model> determine ( Model model, VigorForm form ) throws ServiceException {

        VigorConfiguration configuration = form.getConfiguration();
        if (VigorUtils.is_Integer(configuration.get(ConfigurationParameters.StopCodonSearchWindow))) {
            defaultSearchWindow = Integer.parseInt(configuration.get(ConfigurationParameters.StopCodonSearchWindow));
        }
        if (VigorUtils.is_Integer(configuration.get(ConfigurationParameters.IntronMinimumSize))) {
            minIntronLength = Integer.parseInt(configuration.get(ConfigurationParameters.IntronMinimumSize));
        }
        List<Model> models = null;
        try {
            models = adjustSpliceSites(model);
        } catch (CloneNotSupportedException e) {
            throw new ServiceException(String.format("Problem adjusting exon boundaries for model %s", model), e);
        }
        return models;
    }

    /**
     * @param model
     * @return : models with permutations and combinations of different splice sites found for each splice region. Exon boundaries are adjusted to the splice region.
     * @throws CloneNotSupportedException
     */
    public List<Model> adjustSpliceSites ( Model model ) throws CloneNotSupportedException {

        List<Model> models = new ArrayList<Model>();
        if (model.getAlignment().getViralProtein().getIntrons().size() >= 1) {
            List<Model> tempModels = new ArrayList<Model>();
            VirusGenome virusGenome = model.getAlignment().getVirusGenome();
            List<SpliceSite> splicePairs = new ArrayList<SpliceSite>();
            if (model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getNonCanonical_spliceSites() != null) {
                splicePairs.addAll(model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getNonCanonical_spliceSites());
            }
            for (int i = 0; i < model.getExons().size() - 1; i++) {
                if (i != model.getExons().size() - 1) {
                    Exon upExon = model.getExons().get(i);
                    Exon downExon = model.getExons().get(i + 1);
                    Range currentExon = upExon.getRange();
                    Range nextExon = downExon.getRange();
                    if (nextExon.getBegin() - currentExon.getEnd() <= minIntronLength) {
                        upExon.set_3p_adjusted(true);
                        downExon.set_5p_adjusted(true);
                    }
                    if (upExon.is_3p_adjusted() != true && downExon.is_5p_adjusted() != true) {
                        //Check if donor and acceptor are found at start and end of intron respectively
                        boolean foundSplicePair = false;
                        String expectedDonor = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(currentExon.getEnd() + 1, currentExon.getEnd() + 2)).build().toString();
                        String expectedAcceptor = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(nextExon.getBegin() - 2, nextExon.getBegin() - 1)).build().toString();
                        for (SpliceSite var : splicePairs) {
                            if (var.donor.equals(expectedDonor) && var.acceptor.equals(expectedAcceptor)) {
                                boolean isCompatible = checkSplicePairCompatibility(currentExon, nextExon, currentExon, nextExon, upExon.getFrame(), downExon.getFrame());
                                if (isCompatible) {
                                    foundSplicePair = true;
                                    tempModels.add(model);
                                }
                            }
                        }
                        boolean isBoundaryAdjusted = false;
                        boolean isPesudogene = false;
                        if (!foundSplicePair && !isBoundaryAdjusted && !isPesudogene) {
                            //determine Donor search window
                            long donorStart = currentExon.getEnd() - defaultSearchWindow;
                            if (donorStart < 0) {
                                donorStart = 0;
                            }
                            long donorEnd = currentExon.getEnd() + defaultSearchWindow;
                            if (donorEnd > model.getAlignment().getVirusGenome().getSequence().getLength()) {
                                donorEnd = model.getAlignment().getVirusGenome().getSequence().getLength();
                            }
                            if (donorStart < currentExon.getBegin()) {
                                donorStart = currentExon.getBegin();
                            }
                            Range donorSearchWindow = Range.of(donorStart, donorEnd);
                            Map<Frame, List<Long>> donorInternalStopsMap = VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, donorSearchWindow);
                            if (donorInternalStopsMap.get(upExon.getSequenceFrame()) != null) {
                                List<Long> donorInternalStops = donorInternalStopsMap.get(upExon.getSequenceFrame());
                                if (donorInternalStops.size() > 0) {
                                    Collections.sort(donorInternalStops);
                                    Long stopCoordinate = donorInternalStops.get(0);
                                    if (stopCoordinate >= currentExon.getBegin()) {
                                        donorSearchWindow = Range.of(donorStart, stopCoordinate);
                                    }
                                }
                            }
                            NucleotideSequence donorSearchSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(donorSearchWindow)
                                    .build();
                            //determine acceptor search window
                            long acceptorStart = nextExon.getBegin() - defaultSearchWindow;
                            if (acceptorStart < 0) {
                                acceptorStart = 0;
                            }
                            long acceptorEnd = nextExon.getBegin() + defaultSearchWindow;
                            if (acceptorEnd > model.getAlignment().getVirusGenome().getSequence().getLength()) {
                                acceptorEnd = model.getAlignment().getVirusGenome().getSequence().getLength();
                            }
                            if (acceptorEnd > nextExon.getEnd()) {
                                acceptorEnd = nextExon.getEnd();
                            }
                            Range acceptorSearchWindow = Range.of(acceptorStart, acceptorEnd);
                            Map<Frame, List<Long>> acceptorInternalStopsMap = VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, acceptorSearchWindow);
                            if (acceptorInternalStopsMap.get(downExon.getSequenceFrame()) != null) {
                                List<Long> acceptorInternalStops = acceptorInternalStopsMap.get(downExon.getSequenceFrame());
                                if (acceptorInternalStops.size() > 0) {
                                    Collections.sort(acceptorInternalStops, Collections.reverseOrder());
                                    Long stopCoordinate = acceptorInternalStops.get(0);
                                    if (stopCoordinate < nextExon.getEnd()) {
                                        acceptorSearchWindow = Range.of(stopCoordinate, acceptorEnd);
                                    }
                                }
                            }
                            NucleotideSequence acceptorSearchSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(acceptorSearchWindow)
                                    .build();
                            for (SpliceSite splicePair : splicePairs) {
                                String donor = splicePair.donor;
                                String acceptor = splicePair.acceptor;
                                List<Range> donorRanges = donorSearchSeq.findMatches(donor, true).distinct().collect(Collectors.toList());
                                List<Range> acceptorRanges = acceptorSearchSeq.findMatches(acceptor, true).distinct().collect(Collectors.toList());
                                final long acceptorStartTemp1 = acceptorSearchWindow.getBegin();
                                final long donorStartTemp1 = donorSearchWindow.getBegin();
                                donorRanges = donorRanges.stream().map(range -> Range.of(range.getBegin() + donorStartTemp1, range.getEnd() + donorStartTemp1)).collect(Collectors.toList());
                                acceptorRanges = acceptorRanges.stream().map(range -> Range.of(range.getBegin() + acceptorStartTemp1, range.getEnd() + acceptorStartTemp1)).collect(Collectors.toList());
                                boolean isNewSpliceSite = true;
                                tempModels.add(model);
                                for (Range donorRange : donorRanges) {
                                    for (Range acceptorRange : acceptorRanges) {
                                        Range foundUpExonRange = Range.of(currentExon.getBegin(), donorRange.getBegin() - 1);
                                        Range foundDownExonRange = Range.of(acceptorRange.getEnd() + 1, nextExon.getEnd());
                                        boolean isCompatible = checkSplicePairCompatibility(currentExon, nextExon, foundUpExonRange, foundDownExonRange, upExon.getFrame(), downExon.getFrame());
                                        if (isCompatible) {
                                            if (isNewSpliceSite && models.size() > 0) {
                                                tempModels.clear();
                                                tempModels.addAll(models);
                                                models.clear();
                                            }
                                            for (Model newModelprev : tempModels) {
                                                Model newModel = newModelprev.clone();
                                                newModel.getExons().get(i).setRange(foundUpExonRange);
                                                newModel.getExons().get(i + 1).setRange(foundDownExonRange);
                                                double donorScore = VigorFunctionalUtils.generateScore(currentExon.getEnd(), donorRange.getBegin());
                                                double acceptorScore = VigorFunctionalUtils.generateScore(nextExon.getBegin(), acceptorRange.getBegin());
                                                double spliceScore = donorScore + acceptorScore;
                                                Map<String, Double> scores = newModel.getScores();
                                                if (scores.containsKey("spliceScore")) {
                                                    double existingScore = scores.get("spliceScore");
                                                    double spliceScoreSum = existingScore + spliceScore;
                                                    scores.replace("spliceScore", spliceScoreSum);
                                                } else {
                                                    scores.put("spliceScore", spliceScore);
                                                }
                                                newModel.setScores(scores);
                                                models.add(newModel);
                                            }
                                            isNewSpliceSite = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (models.size() == 0) {
            models.add(model);
        }
        return models;
    }

    public boolean checkSplicePairCompatibility ( Range upExonRange, Range downExonRange, Range foundUpExonRange, Range foundDownExonRange, Frame upExonFrame, Frame downExonFrame ) {

        int upNucleotides = (int) ( upExonRange.getLength() - ( upExonFrame.getFrame() - 1 ) ) % 3;
        int downNucleotides = downExonFrame.getFrame() - 1;
        int extendedUpNucleotides = Math.abs((int) ( upExonRange.getEnd() - foundUpExonRange.getEnd() ));
        int extendedDownNucleotides = Math.abs((int) ( downExonRange.getBegin() - foundDownExonRange.getBegin() ));
        return ( upNucleotides + downNucleotides + extendedDownNucleotides + extendedUpNucleotides ) % 3 == 0;
    }
}
