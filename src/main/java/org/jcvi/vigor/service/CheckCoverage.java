package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.*;
import org.springframework.stereotype.Service;

@Service
public class CheckCoverage implements EvaluateModel {

    @Override
    public Model evaluate ( Model model, VigorConfiguration configuration ) {

        NucleotideSequence cds = determineCDS(model);
        List<Range> internalStops = getInternalStops(model);
        if (internalStops.size() > 0) {
            model.setPseudogene(true);
            List<NoteType> notes = model.getNotes();
            if (internalStops.size() > 1) notes.add(NoteType.StopCodonsInterruption);
            else notes.add(NoteType.StopCodonInterruption);
            model.setNotes(notes);
        }
        model = determineHomology(model, cds);
        return model;
    }

    /**
     * @param model
     * @param cds
     * @return
     */
    public Model determineHomology ( Model model, NucleotideSequence cds ) {

        long replacementOffset = 0;
        if (model.getReplaceStopCodonRange() != null) {
            replacementOffset = model.getReplaceStopCodonRange().getBegin();
        }
        //replace stopcodon with an amino acid as per viral protein specifications
        if (replacementOffset != 0) {
            replacementOffset = getTranslatedProteinCooridnate(model.getExons(), replacementOffset, model.getInsertRNAEditingRange());
        }
        Frame fFrame = model.getExons().get(0).getFrame();
        AminoAcid replacementAA = model.getAlignment().getViralProtein().getGeneAttributes().getStopTranslationException().getReplacementAA();
        ProteinSequence translatedSeq = IupacTranslationTables.STANDARD.translate(cds, fFrame);
        ProteinSequenceBuilder proteinSeqBuilder = new ProteinSequenceBuilder(translatedSeq);
        if (replacementOffset != 0 && replacementAA != null) {
            proteinSeqBuilder.replace((int) replacementOffset, replacementAA);
        }
        ProteinSequence querySeq = proteinSeqBuilder.build();
        model.setTranslatedSeq(querySeq);
        ProteinSequence subSeq = model.getAlignment().getViralProtein().getSequence();
        AminoAcidSubstitutionMatrix blosom62 = BlosumMatrices.blosum62();
        ProteinPairwiseSequenceAlignment actual = PairwiseAlignmentBuilder
                .createProtienAlignmentBuilder(querySeq,
                        subSeq, blosom62).gapPenalty(-8, -8)
                .build();
        Map<String, Double> scores = new HashMap<String, Double>();
        if (model.getScores() != null) {
            scores.putAll(model.getScores());
        }
        double percentIdentity = actual.getPercentIdentity() * 100;
        long maxSeqLength = Long.max(querySeq.getLength(), subSeq.getLength());
        double percentSimilarity = SequenceUtils.computePercentSimilarity(actual.getGappedQueryAlignment(), actual.getGappedSubjectAlignment(), maxSeqLength, blosom62);
        double maxAlignmentLength = Long.max(actual.getQueryRange().getLength(), actual.getSubjectRange().getLength());
        double percentCoverage = ( maxAlignmentLength / maxSeqLength ) * 100;
        if (percentCoverage > 100) percentCoverage = 100;
        scores.put("%identity", percentIdentity);
        scores.put("%similarity", percentSimilarity);
        scores.put("%coverage", percentCoverage);
        double modelScore = percentIdentity + percentSimilarity + percentCoverage;
        scores.put("modelScore", modelScore);
        model.setScores(scores);
        return model;
    }

    /**
     * @param model
     * @return coding sequence of model
     */
    public NucleotideSequence determineCDS ( Model model ) {

        model.getExons().sort(Exon.Comparators.Ascending);
        List<Exon> exons = model.getExons();
        NucleotideSequence virusGenomeSeq = model.getAlignment().getVirusGenome().getSequence();
        boolean inserted = false;
        NucleotideSequenceBuilder NTSeqBuilder = new NucleotideSequenceBuilder("");
        NucleotideSequence NTSeq;
        for (int i = 0; i < exons.size(); i++) {
            Range exonRange = exons.get(i).getRange();
            // trim off stop codon unless the model is partial
            if (i == exons.size() - 1 && !model.isPartial3p()) {
                exonRange = Range.of(exonRange.getBegin(), exonRange.getEnd() - 3);
            }
            if (!inserted && model.getInsertRNAEditingRange() != null && model.getInsertRNAEditingRange().getBegin() == ( exonRange.getEnd() + 1 )) {
                NTSeqBuilder.append(virusGenomeSeq.toBuilder(exonRange));
                String insertionString = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing().getInsertionString();
                NTSeqBuilder.append(insertionString);
                inserted = true;
            } else {
                NTSeqBuilder.append(virusGenomeSeq.toBuilder(exonRange));
            }
        }
        NTSeq = NTSeqBuilder.build();
        return NTSeq;
    }

    /**
     * @param exons
     * @param NTOffset
     * @param insertionRange
     * @return
     */
    public long getTranslatedProteinCooridnate ( List<Exon> exons, long NTOffset, Range insertionRange ) {

        long translatedProteinLength = 0;
        long proteinCoordinate = 0;
        Frame adjustedFrame = null;
        for (Exon exon : exons) {
            long insertionLength = 0;
            Range exonRange = exon.getRange();
            long length = exonRange.getLength();
            if (insertionRange != null) {
                if (insertionRange.getBegin() == exon.getRange().getEnd() + 1) {
                    insertionLength = insertionRange.getLength();
                    long tempLength = length - exon.getFrame().getFrame() + 1 + insertionLength;
                    long leftOvers = tempLength % 3;
                    if (leftOvers == 1) {
                        adjustedFrame = Frame.THREE;
                    } else if (leftOvers == 2) {
                        adjustedFrame = Frame.TWO;
                    } else if (leftOvers == 0) {
                        adjustedFrame = Frame.ONE;
                    }
                }
            }
            if (exonRange.intersects(Range.of(NTOffset))) {
                Frame exonFrame = exon.getFrame();
                if (adjustedFrame != null) {
                    exonFrame = adjustedFrame;
                }
                long startCoordinate = exonRange.getBegin() + exonFrame.getFrame() - 1;
                long difference = NTOffset - startCoordinate;
                if (difference < 0) {
                    proteinCoordinate = proteinCoordinate + translatedProteinLength;
                    break;
                } else {
                    float temp = ( (float) difference ) / 3;
                    proteinCoordinate = (long) Math.ceil(temp);
                    proteinCoordinate = proteinCoordinate + translatedProteinLength;
                    break;
                }
            }
            length = length - exon.getFrame().getFrame() + 1 + insertionLength;
            float temp = ( (float) length ) / 3;
            translatedProteinLength = (long) ( translatedProteinLength + Math.ceil(temp) );
        }
        return proteinCoordinate;
    }

    /**
     * @param model
     * @return list of stop codon ranges
     */
    public List<Range> getInternalStops ( Model model ) {

        List<Range> internalStops = new ArrayList<Range>();
        NucleotideSequence cds = VigorFunctionalUtils.getCDS(model);
        Map<Frame, List<Long>> stops = IupacTranslationTables.STANDARD.findStops(cds);
        Frame fFrame = model.getExons().get(0).getFrame();
        for (Map.Entry<Frame, List<Long>> pair : stops.entrySet()) {
            if (pair.getKey().equals(fFrame)) {
                List<Long> cdsStops = pair.getValue();
                for (Long stop : cdsStops) {
                    long NTStop = VigorFunctionalUtils.getNTRange(model.getExons(), stop);
                    Range NTStopRange = Range.of(NTStop, NTStop + 2);
                    boolean internalStop = true;
                    if (model.getReplaceStopCodonRange() != null && NTStopRange.equals(model.getReplaceStopCodonRange())) {
                        internalStop = false;
                    }
                    if (Range.of(stop).equals(Range.of(cds.getLength() - 3))) {
                        internalStop = false;
                    }
                    if (internalStop) internalStops.add(NTStopRange);
                }
            }
        }
        return internalStops;
    }
}
