package org.jcvi.vigor.service;

import java.util.*;

import com.google.common.collect.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.*;
import org.jcvi.jillion.core.Direction;

@Service
public class DetermineMissingExons implements DetermineGeneFeatures {

    private static Logger LOGGER = LogManager.getLogger(DetermineMissingExons.class);

    @Override
    public List<Model> determine ( Model model ) {

        List<Model> outModels = new ArrayList<Model>();
        VigorConfiguration config = model.getAlignment().getViralProtein().getConfiguration();
        int maxIntronSize = config.getOrDefault(ConfigurationParameters.IntronMaximumSize, 2500);
        int min_missing_AA_size = config.getOrDefault(ConfigurationParameters.MinimumMissingAASize, 10);
        String proteinID = model.getProteinID();
        LOGGER.trace(() ->
                     {
                         VigorConfiguration.ValueWithSource serviceDefault = VigorConfiguration.ValueWithSource.of("",
                                                                                                                   "service default");
                         return String.format("Determining missing exons for %s using %s=%s (%s) %s=%s (%s)",
                                       proteinID,
                                       ConfigurationParameters.IntronMaximumSize.configKey,
                                       maxIntronSize,
                                       config.getWithSource(ConfigurationParameters.IntronMaximumSize)
                                             .orElse(serviceDefault).source,
                                       ConfigurationParameters.MinimumMissingAASize.configKey,
                                       min_missing_AA_size,
                                       config.getWithSource(ConfigurationParameters.MinimumMissingAASize)
                         .orElse(serviceDefault).source);

                     });
        model.getExons().sort(Exon.Comparators.Ascending);
        model.getExons().addAll(findMissingExons(model, maxIntronSize, min_missing_AA_size));
        model.getExons().sort(Exon.Comparators.Ascending);
        outModels.add(model);
        return outModels;
    }

    /**
     * @param NTRange
     * @param AARange
     * @param NTSequence
     * @param AASequence
     * @param modelDirection
     * @return
     */
    public Optional<Exon> performJillionPairWiseAlignment (Range NTRange, Range AARange,
                                                           NucleotideSequence NTSequence, ProteinSequence AASequence, Direction modelDirection ) {

        Exon exon = null;

        NucleotideSequence NTSubSequence = NTSequence.toBuilder(NTRange)
                                                     .build();
        ProteinPairwiseSequenceAlignment actual = null;
        ProteinSequence subjectAASequence = AASequence.toBuilder(AARange)
                                                      .build();
        Map<Frame, ProteinPairwiseSequenceAlignment> alignments = new HashMap<Frame, ProteinPairwiseSequenceAlignment>();
        ProteinPairwiseSequenceAlignment bestAlignment = null;
        for (Frame frame : Frame.forwardFrames()) {
            ProteinSequence queryAASequence = IupacTranslationTables.STANDARD.translate(NTSubSequence, frame);
            AminoAcidSubstitutionMatrix blosom62 = BlosumMatrices.blosum62();
            actual = PairwiseAlignmentBuilder
                    .createProtienAlignmentBuilder(queryAASequence,
                                                   subjectAASequence, blosom62).gapPenalty(-8, -8)
                    .build();
            alignments.put(frame, actual);
        }

        for (Frame myFrame : alignments.keySet()) {
            if (bestAlignment == null) {
                bestAlignment = alignments.get(myFrame);
            }
            if (bestAlignment.getScore() < alignments.get(myFrame).getScore()) {
                bestAlignment = alignments.get(myFrame);
                if (myFrame == Frame.TWO) {
                    NTRange = Range.of(NTRange.getBegin() + 1, NTRange.getEnd());
                } else if (myFrame == Frame.THREE) {
                    NTRange = Range.of(NTRange.getBegin() + 2, NTRange.getEnd());
                }
            }
        }

       // A model should have all the alignment fragments in the same direction. Hence the bestAlignment has to be in the same direction as model.
        if (bestAlignment != null && bestAlignment.getQueryRange().getDirection().equals(modelDirection)) {
            Range modelExonAARange = Range.of(bestAlignment.getSubjectRange()
                                                           .getRange().getBegin()
                                                      + AARange.getBegin(), bestAlignment.getSubjectRange()
                                                                                         .getRange().getEnd()
                                                      + AARange.getBegin());
            exon = new Exon();
            Range range = bestAlignment.getQueryRange().getRange();
            Range modelExonNTRange = Range.of(
                    (range.getBegin() * 3) + NTRange.getBegin(),
                    (((range.getEnd() + 1) * 3) - 1) + NTRange.getBegin());
            exon.setRange(modelExonNTRange);
            AlignmentFragment alignmentFragment = new AlignmentFragment();
            alignmentFragment.setDirection(bestAlignment.getQueryRange()
                                                        .getDirection());
            alignmentFragment.setNucleotideSeqRange(modelExonNTRange);
            alignmentFragment.setProteinSeqRange(modelExonAARange);
            exon.setAlignmentFragment(alignmentFragment);
            exon.setFrame(Frame.ONE);

        }

        return Optional.ofNullable(exon);
    }

    /**
     * @param model
     * @return
     */
    public List<Exon> findMissingExons ( Model model, int maxIntronSize, int min_missing_AA_size ) {

        List<Exon> exons = model.getExons();
        List<Exon> missingExons = new ArrayList<>();
        if (exons.isEmpty()) {
            return missingExons;
        }
        long proteinLength = model.getAlignment().getViralProtein().getSequence().getLength();
        NucleotideSequence NTSeq = model.getAlignment().getVirusGenome().getSequence();
        ProteinSequence AASeq = model.getAlignment().getViralProtein().getSequence();
        List<Range> sequenceGaps = model.getAlignment().getVirusGenome().getSequenceGaps();
        long NTSeqLength = NTSeq.getLength();
        Range preNTRange = exons.get(0).getRange();
        Range preAARange = exons.get(0).getAlignmentFragment().getProteinSeqRange();
        Exon prevExon = null;

        for (int i = 0; i <= exons.size(); i++) {

            if ( (i == 0 && model.isPartial5p()) || (i == exons.size() && model.isPartial3p()) ) {
                continue;
            }

            int j = i;
            if (i == exons.size()) {
                j = exons.size() - 1;
            }
            Exon exon = exons.get(j);
            Range aaRange = exon.getAlignmentFragment().getProteinSeqRange();
            Range ntRange = exon.getRange();
            Range missingAARange = Range.ofLength(0);
            Range missingNTRange = Range.ofLength(0);
            if (i == 0) {
                missingAARange = Range.of(0, aaRange.getBegin());
                missingNTRange = Range.of(0, ntRange.getBegin());
            } else if (i == exons.size()) {
                missingAARange = Range.of(aaRange.getEnd() + 1, proteinLength - 1);
                missingNTRange = Range.of(ntRange.getEnd() + 1, NTSeqLength - 1);
            } else if (preAARange.getEnd() - aaRange.getBegin() > min_missing_AA_size &&
                    preNTRange.getEnd() - ntRange.getBegin() > min_missing_AA_size * 3) {
                missingAARange = Range.of(preAARange.getEnd() + 1, aaRange.getBegin() - 1);
                missingNTRange = Range.of(preNTRange.getEnd() + 1, ntRange.getBegin() - 1);
            }

            if (missingAARange.getLength() > min_missing_AA_size) {

                long temp = maxIntronSize + (missingAARange.getLength() * 3);
                if (missingNTRange.getLength() > temp) {
                    missingNTRange = Range.of(missingNTRange.getBegin(), missingNTRange.getBegin() + temp);
                }
                boolean sequenceGap = false;
                Optional<Range> gapRange = Lists.reverse(sequenceGaps).stream().filter(missingNTRange::intersects).findFirst();
                if (gapRange.isPresent()) {
                    Range intersection = gapRange.get().intersection(missingNTRange);
                    Range leftOver;
                    if (i == 0) {
                        leftOver = Range.of(intersection.getEnd() + 1, missingNTRange.getEnd());
                    } else {
                        leftOver = Range.of(missingNTRange.getBegin(), intersection.getBegin() - 1);
                    }
                    if (leftOver != null && leftOver.getLength() >= 20) {
                        missingNTRange = leftOver;
                    } else {
                        sequenceGap = true;
                    }
                }

                //do not check for missing exons in the sequence gap
                if (!sequenceGap && missingNTRange.getLength() > min_missing_AA_size * 3) {
                    Optional<Exon> determinedExon = performJillionPairWiseAlignment(missingNTRange,
                                                                                    missingAARange, NTSeq, AASeq, model.getDirection());
                    if (determinedExon.isPresent()) {
                        missingExons.add(determinedExon.get());

                        // Unset 3' /5' adjusted flags if there is any missing protein alignment up/down stream of the exon
                        if (prevExon != null) {
                            prevExon.set_3p_adjusted(false);
                        }
                        if (i == exons.size()) {
                            exon.set_3p_adjusted(false);
                        } else {
                            exon.set_5p_adjusted(false);
                        }
                    }
                }
            }
            preNTRange = ntRange;
            preAARange = aaRange;
            prevExon = exon;

        }
        return missingExons;
    }
}
