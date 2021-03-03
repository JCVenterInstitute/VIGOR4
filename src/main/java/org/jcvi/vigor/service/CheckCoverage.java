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
import org.jcvi.vigor.component.MappedNucleotideSequence;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.*;
import org.springframework.stereotype.Service;

@Service
public class CheckCoverage implements EvaluateModel {

    @Override
    public Model evaluate ( Model model, VigorConfiguration configuration ) {

        List<Range> internalStops = getInternalStops(model);
        if (internalStops.size() > 0) {
            model.setPseudogene(true);
            model.addNote(internalStops.size() > 1 ? NoteType.StopCodonsInterruption : NoteType.StopCodonInterruption);
        }
        model = determineHomology(model);
        return model;
    }

    /**
     * @param model
     * @return
     */
    public Model determineHomology ( Model model) {

        MappedNucleotideSequence mappedNucleotideSequence = model.getCDS();
        Frame fFrame = model.getExons().get(0).getFrame();
        ProteinSequence translatedSeq = IupacTranslationTables.STANDARD.translate(mappedNucleotideSequence.getSequence(), fFrame);

        //replace stopcodon with an amino acid as per viral protein specifications
        if (model.getReplaceStopCodonRange() != null) {
            AminoAcid replacementAA = model.getAlignment().getViralProtein().getGeneAttributes().getStopTranslationException().getReplacementAA();
            if (replacementAA != null) {
                long replacementOffset = model.getReplaceStopCodonRange().getBegin();
                if (replacementOffset != 0) {
                    replacementOffset = getTranslatedProteinCoordinate(model.getExons(), replacementOffset, model.getInsertRNAEditingRange());
                }
                if (replacementOffset != 0) {
                    ProteinSequenceBuilder proteinSeqBuilder = new ProteinSequenceBuilder(translatedSeq);
                    proteinSeqBuilder.replace((int) replacementOffset, replacementAA);
                    translatedSeq = proteinSeqBuilder.build();
                }
            }
        }
        model.setTranslatedSeq(translatedSeq);
        ProteinSequence subSeq = model.getAlignment().getViralProtein().getSequence();
        AminoAcidSubstitutionMatrix blosom62 = BlosumMatrices.blosum62();
        ProteinPairwiseSequenceAlignment actual = PairwiseAlignmentBuilder
                .createProteinAlignmentBuilder(translatedSeq,
                        subSeq, blosom62).gapPenalty(-8, -8)
                .build();
        Map<String, Double> scores = new HashMap<String, Double>();
        if (model.getScores() != null) {
            scores.putAll(model.getScores());
        }
        double percentIdentity = actual.getPercentIdentity() * 100;
        long maxSeqLength = Long.max(translatedSeq.getLength(), subSeq.getLength());
        double percentSimilarity = SequenceUtils.computePercentSimilarity(actual.getGappedQueryAlignment(), actual.getGappedSubjectAlignment(), maxSeqLength, blosom62);
        double maxAlignmentLength = Long.max(actual.getQueryRange().getLength(), actual.getSubjectRange().getLength());
        double percentCoverage = Double.min( (maxAlignmentLength / maxSeqLength ) * 100, 100d);

        scores.put(Scores.IDENTITY_SCORE, percentIdentity);
        scores.put(Scores.SIMILARITY_SCORE, percentSimilarity);
        scores.put(Scores.COVERAGE_SCORE, percentCoverage);
        double modelScore = percentIdentity + percentSimilarity + percentCoverage;
        scores.put(Scores.MODEL_SCORE, modelScore);
        model.setScores(scores);
        return model;
    }


    /**
     * @param exons
     * @param NTOffset
     * @param insertionRange
     * @return
     */
    public long getTranslatedProteinCoordinate(List<Exon> exons, long NTOffset, Range insertionRange ) {

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
        MappedNucleotideSequence cds = model.getCDS(false);
        Map<Frame, List<Long>> stops = IupacTranslationTables.STANDARD.findStops(cds.getSequence());
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
                    if (Range.of(stop).equals(Range.of(cds.getSequence().getLength() - 3))) {
                        internalStop = false;
                    }
                    if (internalStop) internalStops.add(NTStopRange);
                }
            }
        }
        return internalStops;
    }
}
