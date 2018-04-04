package org.jcvi.vigor.component;

import lombok.Data;
import org.jcvi.jillion.core.Range;
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
    private ViralProtein protein;

    /**
     * Reference peptide
     */
    private ViralProtein reference;

    /**
     * position of alignment of reference to protein starting from beginning of protein (1 based coordinates).
     */
    private Range proteinRange;

    /**
     * position of alignment of reference to protein starting from beginning of reference (1 based coordinates).
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


    public static MaturePeptideMatch of(ViralProtein protein, ViralProtein reference, Range proteinRange, Range referenceRange) {
        return of(protein, reference, proteinRange, referenceRange, false, false);
    }

    public static MaturePeptideMatch of(ViralProtein protein, ViralProtein reference, Range proteinRange, Range referenceRange, boolean fuzzyBegin, boolean fuzzyEnd) {
        MaturePeptideMatch mpMatch = new MaturePeptideMatch();
        mpMatch.setProtein(protein);
        mpMatch.setReference(reference);
        mpMatch.setProteinRange(proteinRange);
        mpMatch.setReferenceRange(referenceRange);
        mpMatch.setFuzzyBegin(fuzzyBegin);
        mpMatch.setFuzzyEnd(fuzzyEnd);
        return mpMatch;
    }
}
