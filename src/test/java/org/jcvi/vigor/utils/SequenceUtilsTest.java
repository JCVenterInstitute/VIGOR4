package org.jcvi.vigor.utils;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.junit.Test;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.MatcherAssert.assertThat;

public class SequenceUtilsTest {

    @Test
    public void testComputeSimilarity() {
        AminoAcidSubstitutionMatrix blosum62 = BlosumMatrices.blosum62();

        ProteinSequence first = new ProteinSequenceBuilder("MSLLTEVETYTLSIIPSGPL").build();

        ProteinSequence second = first.toBuilder().build();
        double similarity = SequenceUtils.computePercentSimilarity(first, second,20, blosum62);
        assertThat("identical sequences should be 100% similary", similarity, equalTo(100d));

        second = second.toBuilder().replace(0,AminoAcid.Arginine).build();
        similarity = SequenceUtils.computePercentSimilarity(first, second,20, blosum62);
        assertThat("replacing one amino acid with an unrelated amino acid should be 95% similary", similarity, equalTo(95d));

        second = first.toBuilder().replace(0, AminoAcid.Isoleucine).build();
        similarity = SequenceUtils.computePercentSimilarity(first, second,20, blosum62);
        assertThat("replacing one amino acid with an related amino acid should be 100% similary", similarity, equalTo(100d));

    }
}
