package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.RNA_Editing;
import org.jcvi.vigor.component.Ribosomal_Slippage;
import org.jcvi.vigor.component.StopTranslationException;
import org.springframework.stereotype.Service;

@Service
public class AdjustViralTricks implements DetermineGeneFeatures {

    private static final Logger LOGGER = LogManager.getLogger(AdjustViralTricks.class);
    private static double DEFAULT_LEAKYSTOP_NOTFOUND_SCORE = .80d;
    @Override
    public List<Model> determine ( Model model) throws ServiceException {

        List<Model> outputModels = new ArrayList<>();
        List<Model> rnaEditedModels = new ArrayList<>();
        double leakyStopNotFoundScore = model.getAlignment().getViralProtein().getConfiguration().getOrDefault(ConfigurationParameters.ScoreFactorLeakyStopNotFound, DEFAULT_LEAKYSTOP_NOTFOUND_SCORE);
        try {
            List<Model> riboAdjustedmodels = adjustRibosomalSlippage(model);
            for (Model riboAdjustedModel : riboAdjustedmodels) {
                rnaEditedModels.addAll(adjustRNAEditing(riboAdjustedModel));
            }
            for (Model rnaEditeddModel : rnaEditedModels) {
                LOGGER.trace("checking for leaky stop using {}", () -> {
                    VigorConfiguration.ValueWithSource val = model.getAlignment().getViralProtein().getConfiguration().getWithSource(ConfigurationParameters.ScoreFactorLeakyStopNotFound).orElse(VigorConfiguration.ValueWithSource.of(String.valueOf(DEFAULT_LEAKYSTOP_NOTFOUND_SCORE),"service default"));
                    return String.format("%s=%s from %s", ConfigurationParameters.ScoreFactorLeakyStopNotFound.configKey, val.source, val.value);
                });
                outputModels.addAll(checkForLeakyStop(rnaEditeddModel, leakyStopNotFoundScore));
            }
        } catch (CloneNotSupportedException e) {
            throw new ServiceException(String.format("Problem adjusting model %s for viral tricks", model), e);
        }
        return outputModels;
    }

    /**
     * @param model
     * @return list of models.Cloned model for each match found.
     * @throws CloneNotSupportedException
     */
    public List<Model> adjustRibosomalSlippage ( Model model ) throws CloneNotSupportedException {

        Ribosomal_Slippage riboSlippage = model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
        List<Model> models = new ArrayList<Model>();
        if (riboSlippage.isHas_ribosomal_slippage()) {
            List<Range> sequenceGaps = model.getAlignment().getVirusGenome().getSequenceGaps();
            long CDSStart = model.getExons().get(0).getRange().getBegin();
            long CDSEnd = model.getExons().get(model.getExons().size() - 1).getRange().getEnd();
            NucleotideSequence cds = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(CDSStart, CDSEnd))
                    .build();
            int offset = riboSlippage.getSlippage_offset();
            //+1 is added if offset is negative, this is to start count from the point where match is found.
            if (offset < 0) {
                offset = offset + 1;
            }
            //Once the matches are found, get the coordinates relative to complete sequence
            List<Range> matches = cds.findMatches(riboSlippage.getSlippage_motif()).map(x -> x.toBuilder().shift(CDSStart).build())
                    .sequential()
                    .collect(Collectors.toList());
            //If no match found and model is not partial then model is marked is pseudogene
            if (matches.isEmpty() && !( model.isPartial3p() || model.isPartial5p() )) {
                model.setPseudogene(true);
            }
            for (Range match : matches) {
                if (!VigorFunctionalUtils.intheSequenceGap(sequenceGaps, match)) {
                    Model newModel = model.clone();
                    Range slippagePoint = Range.of(match.getEnd() + offset);
                    newModel.setRibosomalSlippageRange(slippagePoint);
                    for (int i = 0; i < newModel.getExons().size(); i++) {
                        Range exonRange = newModel.getExons().get(i).getRange();
                        if (exonRange.intersects(slippagePoint)) {
                            PointOfOccurrence pointOfOccurance = determineLocation(exonRange, slippagePoint, true);
                            //if match found in the middle of an exon , adjust 3' end of upstream part exon and  5' end of downstream part of exon (results in two exons)
                            if (pointOfOccurance == PointOfOccurrence.MIDDLE) {
                                Exon exon = newModel.getExons().get(i).clone();
                                exon.set_5p_adjusted(true);
                                exon.setRange(Range.of(slippagePoint.getBegin() + riboSlippage.getSlippage_frameshift(), exonRange.getEnd()));
                                exon.setFrame(Frame.ONE);
                                newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(), slippagePoint.getBegin() - 1));
                                newModel.getExons().get(i).set_3p_adjusted(true);
                                newModel.getExons().add(exon);
                            //if match found at the start of the exon, adjust 5' end of the exon and 3' end of previous exon, if any.
                            } else if (pointOfOccurance == PointOfOccurrence.START) {
                                newModel.getExons().get(i).setRange(Range.of(slippagePoint.getBegin() + riboSlippage.getSlippage_frameshift(), exonRange.getEnd()));
                                newModel.getExons().get(i).set_5p_adjusted(true);
                                if (i != 0) {
                                    Range prevExonRange = newModel.getExons().get(i - 1).getRange();
                                    newModel.getExons().get(i - 1).setRange(Range.of(prevExonRange.getBegin(), slippagePoint.getBegin() - 1));
                                    newModel.getExons().get(i - 1).set_3p_adjusted(true);
                                }
                            //if match is found at the end of the exon, adjust 3' end of the exon and the 5' end of next exon, if any
                            } else if (pointOfOccurance == PointOfOccurrence.END) {
                                newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(), slippagePoint.getBegin() - 1));
                                newModel.getExons().get(i).set_3p_adjusted(true);
                                if (i != newModel.getExons().size() - 1) {
                                    Range nextExonRange = newModel.getExons().get(i + 1).getRange();
                                    newModel.getExons().get(i + 1).setRange(Range.of(slippagePoint.getBegin() + riboSlippage.getSlippage_frameshift(), nextExonRange.getEnd()));
                                    newModel.getExons().get(i + 1).set_5p_adjusted(true);
                                }
                            }
                        //Below scenario is when slippage point lies in the intron, extend upstream and downstream exons till the slippage point
                        } else if (i != newModel.getExons().size() - 1) {
                            Range nextExonRange = newModel.getExons().get(i + 1).getRange();
                            if (slippagePoint.intersects(Range.of(exonRange.getEnd() + 1, nextExonRange.getBegin() - 1))) {
                                newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(), slippagePoint.getBegin() - 1));
                                newModel.getExons().get(i).set_3p_adjusted(true);
                                newModel.getExons().get(i + 1).setRange(Range.of(slippagePoint.getBegin() + riboSlippage.getSlippage_frameshift(), nextExonRange.getEnd()));
                                newModel.getExons().get(i + 1).set_5p_adjusted(true);
                            }
                        }
                    }
                    models.add(newModel);
                }
            }
        }
        if (models.size() == 0) {
            models.add(model);
        }
        return models;
    }

    /**
     * @param model
     * @return List of models. Cloned model for each match found.
     * @throws CloneNotSupportedException
     */
    public List<Model> adjustRNAEditing ( Model model ) throws CloneNotSupportedException {

        List<Model> models = new ArrayList<Model>();
        if (model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing().isHas_RNA_editing()) {
            RNA_Editing rna_editing = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing();
            Range modelRange = model.getRange();
            long CDSStart = modelRange.getBegin();
            NucleotideSequence cds = model.getAlignment().getVirusGenome().getSequence().toBuilder(modelRange)
                    .build();
            List<Range> matches = cds.findMatches(rna_editing.getRegExp())
                    .distinct()
                    .map(x -> x = x.toBuilder().shift(CDSStart).build())
                    .sequential()
                    .collect(Collectors.toList());
            List<Range> sequenceGaps = model.getAlignment().getVirusGenome().getSequenceGaps();
            int offset = rna_editing.getOffset();
            //+1 is added if offset is negative, this is to start count from the point where match is found.
            if (rna_editing.getOffset() < 0) {
                offset = offset + 1;
            }
            for (Range match : matches) {
                if (!VigorFunctionalUtils.intheSequenceGap(sequenceGaps, match)) {
                    Model newModel = model.clone();
                    Range pointOfInsertion = Range.of(match.getEnd() + offset, match.getEnd() + rna_editing.getInsertionString().length() - 1);
                    newModel.setInsertRNAEditingRange(pointOfInsertion);
                    for (int i = 0; i < newModel.getExons().size(); i++) {
                        Range exonRange = newModel.getExons().get(i).getRange();
                        if (exonRange.intersects(pointOfInsertion)) {
                            PointOfOccurrence pointOfOccurance = determineLocation(exonRange, pointOfInsertion);
                            // If the position where insertionString has to be inserted lies at the start, adjust previous exon and current exon ranges till the point of insertion
                            if (pointOfOccurance == PointOfOccurrence.START) {
                                newModel.getExons().get(i).setRange(Range.of(pointOfInsertion.getBegin(), exonRange.getEnd()));
                                if (i != 0) {
                                    Range prevExonRange = newModel.getExons().get(i - 1).getRange();
                                    newModel.getExons().get(i - 1).setRange(Range.of(prevExonRange.getBegin(), pointOfInsertion.getBegin() - 1));
                                }
                            // If the position where insertionString has to be inserted lies in the end of the exon , adjust current exon and next exon ranges till the point of insertion
                            } else if (pointOfOccurance == PointOfOccurrence.END) {
                                newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(), pointOfInsertion.getBegin() - 1));
                                if (i != newModel.getExons().size() - 1) {
                                    Range nextExonRange = newModel.getExons().get(i + 1).getRange();
                                    newModel.getExons().get(i + 1).setRange(Range.of(pointOfInsertion.getBegin(), nextExonRange.getEnd()));
                                } else {
                                    Exon exon = new Exon(Range.of(pointOfInsertion.getEnd() + 1, exonRange.getEnd()), Frame.ONE);
                                    Frame sequenceFrame = VigorFunctionalUtils.getSequenceFrame(exon.getRange().getBegin() + exon.getFrame().getFrame() - 1);
                                    exon.setSequenceFrame(sequenceFrame);
                                    exon.setAlignmentFragment(newModel.getExons().get(i).getAlignmentFragment());
                                    newModel.getExons().add(exon);
                                }
                            }
                        } else if (i != newModel.getExons().size() - 1) {
                            Range nextExonRange = newModel.getExons().get(i + 1).getRange();
                            if (pointOfInsertion.intersects(Range.of(exonRange.getEnd() + 1, nextExonRange.getBegin() - 1))) {
                                newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(), pointOfInsertion.getBegin() - 1));
                                newModel.getExons().get(i + 1).setRange(Range.of(pointOfInsertion.getBegin(), nextExonRange.getEnd()));
                            }
                        }
                    }
                    newModel.addNote(NoteType.RNA_Editing);
                    models.add(newModel);
                }
            }
        }
        if (models.size() == 0) {
            models.add(model);
        }
        return models;
    }

    /**
     * @param model
     * @return List of models. Cloned model for each match found.
     * @throws CloneNotSupportedException
     */
    public List<Model> checkForLeakyStop (Model model, double leakyStopNotFoundScore ) throws CloneNotSupportedException {

        List<Model> newModels = new ArrayList<Model>();
        List<Range> sequenceGaps = model.getAlignment().getVirusGenome().getSequenceGaps();
        StopTranslationException stopTransExce = model.getAlignment().getViralProtein().getGeneAttributes().getStopTranslationException();
        if (stopTransExce.isHasStopTranslationException()) {
            NucleotideSequence cds = VigorFunctionalUtils.getCDS(model);
            List<Range> matches = cds.findMatches(stopTransExce.getMotif()).distinct().collect(Collectors.toList());
            int offset = stopTransExce.getOffset();
            //+1 is added if offset is negative, this is to start count from the point where match is found.
            if (offset < 0) {
                offset = offset + 1;
            }
            if (matches != null) {
                for (Range match : matches) {
                    long leakyStopStart = VigorFunctionalUtils.getNTRange(model.getExons(), match.getBegin());
                    Range leakyStopRange = Range.of(leakyStopStart, leakyStopStart + match.getLength() - 1);
                    long start = leakyStopRange.getEnd() + offset;
                    //match in frame has to considered. For every match found, clone a model and save replacement stop codon range
                    if (VigorFunctionalUtils.isInFrameWithExon(model.getExons(), start) && !VigorFunctionalUtils.intheSequenceGap(sequenceGaps, leakyStopRange)) {
                        Model newModel;
                        newModel = model.clone();
                        LOGGER.trace("Stop codon read through " + start);
                        LOGGER.trace("Sequence {}", () -> model.getAlignment().getVirusGenome().getSequence().toBuilder().trim(Range.of(start, start + 2)).build());
                        Map<String, Double> scores = newModel.getScores();
                        scores.put(Scores.LEAKY_STOP_SCORE, 100.00);
                        newModel.setScores(scores);

                        newModel.setReplaceStopCodonRange(Range.of(start, start + 2));
                        newModel.addNote(NoteType.StopCodonReadThrough);
                        newModels.add(newModel);
                    }
                }
            }
            if (newModels.isEmpty()) {
                Map<String, Double> scores = model.getScores();
                scores.put(Scores.LEAKY_STOP_SCORE, leakyStopNotFoundScore);
                model.setScores(scores);
            }
        }

        if (newModels.isEmpty()) {
            newModels.add(model);
        }
        return newModels;
    }

    private enum PointOfOccurrence {START, MIDDLE, END, NONE}

    ;

    private PointOfOccurrence determineLocation ( Range searchRange, Range inputRange ) {

        return determineLocation(searchRange, inputRange, false);
    }

    private PointOfOccurrence determineLocation ( Range searchRange, Range inputRange, boolean checkMiddle ) {

        PointOfOccurrence returnVal = PointOfOccurrence.NONE;
        long length = searchRange.getLength();
        if (inputRange.isSubRangeOf(searchRange)) {
            Range start;
            Range middle = null;
            if (checkMiddle) {
                start = Range.of(searchRange.getBegin(), searchRange.getBegin() + length / 4);
                middle = Range.of(start.getEnd() + 1, start.getEnd() + 1 + length / 2);
            } else {
                start = Range.of(searchRange.getBegin(), searchRange.getBegin() + length / 2);
            }
            if (inputRange.isSubRangeOf(start)) {
                returnVal = PointOfOccurrence.START;
            } else if (middle != null && inputRange.isSubRangeOf(middle)) {
                returnVal = PointOfOccurrence.MIDDLE;
            } else {
                returnVal = PointOfOccurrence.END;
            }
        }
        //LOGGER.debug("Returning {} for searchRange {}, inputRange {}", returnVal, searchRange, inputRange);
        return returnVal;
    }
}
