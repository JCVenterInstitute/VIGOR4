package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.component.SpliceSite;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.springframework.stereotype.Service;

/**
 *
 */
@Service
public class AdjustUneditedExonBoundaries implements DetermineGeneFeatures {

    private static Logger LOGGER = LogManager.getLogger(AdjustUneditedExonBoundaries.class);
    private static final int DEFAULT_STOP_CODON_SEARCH_WINDOW = 50;
    private static final int DEFAULT_MIN_INTRON_SIZE = 20;

    @Override
    public List<Model> determine ( Model model ) throws ServiceException {

        VigorConfiguration configuration = model.getAlignment().getViralProtein().getConfiguration();
        int defaultSearchWindow = configuration.getOrDefault(ConfigurationParameters.StopCodonSearchWindow,
                                                             DEFAULT_STOP_CODON_SEARCH_WINDOW);
        int minIntronLength = configuration.getOrDefault(ConfigurationParameters.IntronMinimumSize,
                                                         DEFAULT_MIN_INTRON_SIZE);

        LOGGER.trace("adjusting unedited exon boundaries using", () -> {return "TODO print values and sources"; });

        try {
            return adjustSpliceSites(model, defaultSearchWindow, minIntronLength);
        } catch (CloneNotSupportedException e) {
            throw new ServiceException(String.format("Problem adjusting exon boundaries for model %s", model), e);
        }
    }

    /**
     * @param model
     * @return : models with permutations and combinations of different splice sites found for each splice region. Exon boundaries are adjusted to the splice region.
     * @throws CloneNotSupportedException
     */
    public List<Model> adjustSpliceSites ( Model model, int defaultSearchWindow, int minIntronLength ) throws CloneNotSupportedException {

        // no introns means no place to splice?
        if (model.getAlignment().getViralProtein().getIntrons().isEmpty()) {
            return Lists.newArrayList(model);
        }

        List<Model> models = new ArrayList<>();
        List<Model> tempModels = new ArrayList<>();
        VirusGenome virusGenome = model.getAlignment().getVirusGenome();
        NucleotideSequence virusGenomeSequence = virusGenome.getSequence();

        List<SpliceSite> splicePairs = model.getAlignment().getViralProtein().getGeneAttributes().getSpliceSites();
        // -2 to skip the last exon
        for (int i = 0; i < model.getExons().size() - 2; i++) {
            Exon upExon = model.getExons().get(i);
            Exon downExon = model.getExons().get(i + 1);
            Range currentExon = upExon.getRange();
            Range nextExon = downExon.getRange();
            if (nextExon.getBegin() - currentExon.getEnd() <= minIntronLength) {
                upExon.set_3p_adjusted(true);
                downExon.set_5p_adjusted(true);
            }
            if ( ! (upExon.is_3p_adjusted() || downExon.is_5p_adjusted()) ) {
                //Check if donor and acceptor are found at start and end of intron respectively
                boolean foundSplicePair = false;
                // expected donor is two nucleotides after end of upstream exon. expected acceptor is two nucleotides before start of the next exon
                String expectedDonor = virusGenomeSequence.toBuilder(Range.of(currentExon.getEnd() + 1, currentExon.getEnd() + 2)).build().toString();
                String expectedAcceptor = virusGenomeSequence.toBuilder(Range.of(nextExon.getBegin() - 2, nextExon.getBegin() - 1)).build().toString();
                for (SpliceSite spliceSite : splicePairs) {
                    if (spliceSite.donor.equals(expectedDonor) && spliceSite.acceptor.equals(expectedAcceptor)) {
                        boolean isCompatible = checkSplicePairCompatibility(currentExon, nextExon, currentExon, nextExon, upExon.getFrame(), downExon.getFrame());
                        if (isCompatible) {
                            foundSplicePair = true;
                            //This has to be added to tempModels, as next exons may enter below loop to find splicesites.
                            tempModels.add(model);
                        }
                    }
                }

                //If we find compatible splice sites above (in the expected location as above), we dont look for any other compatible splice pairs as below
                if (!foundSplicePair ) {
                    //determine Donor search window
                    long donorStart = Math.max(currentExon.getEnd() - defaultSearchWindow, currentExon.getBegin());
                    long donorEnd = Math.min(currentExon.getEnd() + defaultSearchWindow,
                                             model.getAlignment().getVirusGenome().getSequence().getLength());
                    Range donorSearchWindow = Range.of(donorStart, donorEnd);
                    Map<Frame, List<Long>> donorInternalStopsMap = VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, donorSearchWindow);
                    List<Long> donorInternalStops = donorInternalStopsMap.getOrDefault(upExon.getSequenceFrame(), Collections.EMPTY_LIST);
                    if (donorInternalStops.size() > 0) {
                        Collections.sort(donorInternalStops);
                        Long stopCoordinate = donorInternalStops.get(0);
                        if (stopCoordinate >= currentExon.getBegin()) {
                            donorSearchWindow = Range.of(donorStart, stopCoordinate);
                        }
                    }
                    NucleotideSequence donorSearchSeq = virusGenomeSequence.toBuilder(donorSearchWindow)
                                                                           .build();
                    //determine acceptor search window
                    long acceptorStart = Math.max(nextExon.getBegin() - defaultSearchWindow,0);
                    long acceptorEnd = Math.min(nextExon.getBegin() + defaultSearchWindow, nextExon.getEnd());
                    Range acceptorSearchWindow = Range.of(acceptorStart, acceptorEnd);

                    Map<Frame, List<Long>> acceptorInternalStopsMap = VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, acceptorSearchWindow);
                    List<Long> acceptorInternalStops = acceptorInternalStopsMap.getOrDefault(downExon.getSequenceFrame(), Collections.EMPTY_LIST);
                    if (acceptorInternalStops.size() > 0) {
                        Collections.sort(acceptorInternalStops, Collections.reverseOrder());
                        Long stopCoordinate = acceptorInternalStops.get(0);
                        if (stopCoordinate < nextExon.getEnd()) {
                            acceptorSearchWindow = Range.of(stopCoordinate, acceptorEnd);
                        }
                    }
                    NucleotideSequence acceptorSearchSeq = virusGenomeSequence.toBuilder(acceptorSearchWindow)
                                                                              .build();
                    for (SpliceSite splicePair : splicePairs) {
                        String donor = splicePair.donor;
                        String acceptor = splicePair.acceptor;
                        List<Range> donorRanges = donorSearchSeq.findMatches(donor, true).distinct().collect(Collectors.toList());
                        List<Range> acceptorRanges = acceptorSearchSeq.findMatches(acceptor, true).distinct().collect(Collectors.toList());
                        final long acceptorStartTemp1 = acceptorSearchWindow.getBegin();
                        final long donorStartTemp1 = donorSearchWindow.getBegin();
                        donorRanges = donorRanges.stream()
                                                 .map(range -> range.toBuilder().shift(donorStartTemp1).build())
                                                 .collect(Collectors.toList());
                        acceptorRanges = acceptorRanges.stream()
                                                       .map(range -> range.toBuilder().shift(acceptorStartTemp1).build())
                                                       .collect(Collectors.toList());
                        boolean isNewSpliceSite = true;
                        tempModels.add(model);
                        //Upstream exon should be trimmed/extended till the donor splice site and downstream exon should be trimmed/extended till the acceptor splice site.
                        // If there are multiple splice pairs , For each Splice pair determine compatible donor and acceptor pairs as below
                        for (Range donorRange : donorRanges) {
                            for (Range acceptorRange : acceptorRanges) {
                                Range foundUpExonRange = Range.of(currentExon.getBegin(), donorRange.getBegin() - 1);
                                Range foundDownExonRange = Range.of(acceptorRange.getEnd() + 1, nextExon.getEnd());
                                if (! checkSplicePairCompatibility(currentExon, nextExon, foundUpExonRange, foundDownExonRange, upExon.getFrame(), downExon.getFrame())) {
                                    continue;
                                }
                                if (isNewSpliceSite && models.size() > 0) {
                                    tempModels.clear();
                                    tempModels.addAll(models);
                                    models.clear();
                                }

                                //Here in the end we have cloned models where each Model has spliceScore from different spliceSites & for all splice pairs exon coordinates are adjusted in a model.
                                //for each compatible donor and acceptor pair clone a new model
                                for (Model newModelprev : tempModels) {
                                    Model newModel = newModelprev.clone();
                                    newModel.getExons().get(i).setRange(foundUpExonRange);
                                    newModel.getExons().get(i + 1).setRange(foundDownExonRange);
                                    //Score each pair based on the location where donor and acceptor are found.
                                    // Donor close to the end of the first exon and acceptor close to the start of the second exon will score high
                                    //If spliceScore already exists it is the model cloned from the previous splice pair, so sum up to scores.
                                    double spliceScore = VigorFunctionalUtils.generateProximityScore(currentExon.getEnd(), donorRange.getBegin()) +
                                            VigorFunctionalUtils.generateProximityScore(nextExon.getBegin(), acceptorRange.getBegin());
                                    Map<String, Double> scores = newModel.getScores();
                                    scores.put("spliceScore",scores.getOrDefault("spliceScore", 0d) + spliceScore);
                                    newModel.setScores(scores);
                                    //In the end (closure of main for loop) we have a all the splice sites adjusted in a model
                                    models.add(newModel);
                                }
                                isNewSpliceSite = false;
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

    //Leftover number of nucleotides at the end of the upstream exon + leftover number of nucleotides at start of the downstream exon = 3 or 0 then splice pair is compatible
    public boolean checkSplicePairCompatibility ( Range upExonRange, Range downExonRange,
                                                  Range foundUpExonRange, Range foundDownExonRange,
                                                  Frame upExonFrame, Frame downExonFrame ) {

        int upNucleotides = (int) ( upExonRange.getLength() - ( upExonFrame.getFrame() - 1 ) ) % 3;
        int downNucleotides = downExonFrame.getFrame() - 1;
        int extendedUpNucleotides = Math.abs((int) ( upExonRange.getEnd() - foundUpExonRange.getEnd() ));
        int extendedDownNucleotides = Math.abs((int) ( downExonRange.getBegin() - foundDownExonRange.getBegin() ));
        return ( upNucleotides + downNucleotides + extendedDownNucleotides + extendedUpNucleotides ) % 3 == 0;
    }
}
