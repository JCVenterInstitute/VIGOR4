package org.jcvi.vigor.component;

import lombok.Data;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

@Component
@Scope("prototype")
@Data
@SuppressWarnings("serial")
public class MaturePeptideMatch {
    /**
     * Protein from model
     */
    private ProteinSequence protein;

    /**
     * Reference peptide
     */
    private ViralProtein reference;

    /**
     * position of alignment of reference to protein starting from beginning of protein.
     */
    private Range proteinRange;

    /**
     * position of alignment of reference to protein starting from beginning of reference.
     */
    private Range referenceRange;

    /**
     * gaps, truncations at beginning of alignment or misalignment to previous peptide
     */
    private boolean fuzzyBegin= false;
    /**
     * gaps, truncations at end of alignment or misalignment to next peptide
     */
    private boolean fuzzyEnd = false;

    private double identity;
    private double coverage;
    private double similarity;

    public static MaturePeptideMatch of(ProteinSequence protein, ViralProtein reference, Range proteinRange, Range referenceRange) {
        return of(protein, reference, proteinRange, referenceRange, false, false, 0, 0, 0);
    }

    public static MaturePeptideMatch of(ProteinSequence protein, ViralProtein reference, Range proteinRange, Range referenceRange, boolean fuzzyBegin, boolean fuzzyEnd, double identity, double similarity, double coverage) {
        MaturePeptideMatch mpMatch = new MaturePeptideMatch();
        mpMatch.setProtein(protein);
        mpMatch.setReference(reference);
        mpMatch.setProteinRange(proteinRange);
        mpMatch.setReferenceRange(referenceRange);
        mpMatch.setFuzzyBegin(fuzzyBegin);
        mpMatch.setFuzzyEnd(fuzzyEnd);
        mpMatch.setIdentity(identity);
        mpMatch.setSimilarity(similarity);
        mpMatch.setCoverage(coverage);
        return mpMatch;
    }
}
