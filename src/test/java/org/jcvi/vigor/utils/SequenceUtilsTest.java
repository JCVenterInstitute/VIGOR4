package org.jcvi.vigor.utils;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.junit.Test;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.MatcherAssert.assertThat;

public class SequenceUtilsTest {

    @Test
    public void testcomputePercentSimilarity () {

        String proteinString = "MSLLTEVETYTLSIIPSGPL";
        int proteinStringLength = proteinString.length();
        ProteinSequence first = new ProteinSequenceBuilder(proteinString).build();
        for (AminoAcidSubstitutionMatrix matrix : new AminoAcidSubstitutionMatrix[] {
                BlosumMatrices.blosum30(),
                BlosumMatrices.blosum40(),
                BlosumMatrices.blosum50(),
                BlosumMatrices.blosum62(),
                BlosumMatrices.blosum90() }) {
            ProteinSequence second = first.toBuilder().build();
            double similarity = SequenceUtils.computePercentSimilarity(first, second, proteinStringLength, matrix);
            assertThat("identical sequences should be 100% similar", similarity, equalTo(100.0d));
            AminoAcid firstAcid = first.get(0);
            for (AminoAcid acid : AminoAcid.values()) {
                second = second.toBuilder().replace(0, acid).build();
                similarity = SequenceUtils.computePercentSimilarity(first, second, proteinStringLength, matrix);
                if (matrix.getValue(firstAcid, acid) > 0) {
                    assertThat(String.format("replacing %s with an related amino acid %s should be 100%% similar", firstAcid, acid), similarity, equalTo(100.0d));
                } else {
                    assertThat(String.format("replacing %s with an unrelated amino acid %s should be 95%% similar", firstAcid, acid), similarity, equalTo(95.0d));
                }
            }
            second = first.toBuilder().trim(Range.of(0, 9)).build();
            similarity = SequenceUtils.computePercentSimilarity(first, second, proteinStringLength, matrix);
            assertThat("An identical sequence 1/2 as long should be 50% similar", similarity, equalTo(50.0d));
            // the order shouldn't matter
            similarity = SequenceUtils.computePercentSimilarity(second, first, proteinStringLength, matrix);
            assertThat("A sequence compared to an identical sequence 1/2 as long should be 50% similar", similarity, equalTo(50.0d));
            second = first.toBuilder().replace(0, AminoAcid.Gap).build();
            similarity = SequenceUtils.computePercentSimilarity(first, second, proteinStringLength, matrix);
            assertThat("replace one amino acid with a gap should be 95% similar", similarity, equalTo(95.0d));
        }
    }
}
