package org.jcvi.vigor.service;

import static junit.framework.TestCase.assertTrue;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.notNullValue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.testutils.assembly.cas.FastaRecordWriter;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.testing.category.Fast;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.testing.category.ReferenceDatabase;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.component.Alignment;
import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

@Category({Fast.class, ReferenceDatabase.class})
@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class AlignmentGenerationServiceTest {

    @Autowired
    private VigorInitializationService initializationService;

    @Test
    public void generateAlignmentsTest () throws VigorException, IOException {

        File virusGenomeSeqFile = new File(this.getClass().getResource("/vigorUnitTestInput/sequence_flua.fasta").getFile());
        File alignmentOutput = new File(this.getClass().getResource("/vigorUnitTestInput/sequence_flua.txt").getFile());Files.createTempFile("test-alignments",".text").toFile();
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String refereceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(refereceDBPath);
        assertThat("reference database path is required", refereceDBPath, is(notNullValue()));
        String referenceDB = Paths.get(refereceDBPath, "flua_db").toString();
        List<Alignment> alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile, referenceDB, alignmentOutput, config);
        assertEquals(11, alignments.size());
        assertTrue("all alignment fragments are on the positive strand",
                   alignments.stream()
                             .flatMap(a -> a.getAlignmentFragments().stream())
                   .allMatch(af -> af.getDirection() == Direction.FORWARD)
        );
    }

    @Test
    public void testReverseComplement() throws VigorException, IOException {
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        File virusGenomeSeqFile = new File(this.getClass().getResource("/vigorUnitTestInput/sequence_flua.fasta").getFile());
        File alignmentOutput = new File(this.getClass().getResource("/vigorUnitTestInput/sequence_flua.txt").getFile());Files.createTempFile("test-alignments",".text").toFile();
        File reverseAlignmentOutput = new File(this.getClass().getResource("/vigorUnitTestInput/sequence_flua-reverse.txt").getFile());Files.createTempFile("test-alignments",".text").toFile();
        String refereceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(refereceDBPath);
        assertThat("reference database path is required", refereceDBPath, is(notNullValue()));
        String referenceDB = Paths.get(refereceDBPath, "flua_db").toString();
        Path reversedFile = Files.createTempFile("sequence_flua-rev", ".fasta");
        //Files.deleteIfExists(reversedFile);
        try (NucleotideFastaDataStore input = new NucleotideFastaFileDataStoreBuilder(virusGenomeSeqFile).build();
             FastaRecordWriter output = new FastaRecordWriter(reversedFile.toFile().getParentFile(),
                                                              reversedFile.toFile().getName(),
                                                              Integer.MAX_VALUE)) {
            input.records().throwingForEach(r -> output.write(r.getId(),
                                                              r.getSequence().toBuilder().reverseComplement().build()));
       }
        List<Alignment> alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile, referenceDB, alignmentOutput, config);
        List<Alignment> alignmentsReverse = VigorTestUtils.getAlignments(reversedFile.toFile(), referenceDB, reverseAlignmentOutput, config);
        assertEquals("expected 11 alignments", 11, alignments.size());
        assertEquals("complemented sequences should find the same number of alignments as the original",
                     alignments.size(),
                     alignmentsReverse.size());

        assertTrue("all alignment fragments are on the reverse strand",
                   alignmentsReverse.stream()
                                    .flatMap(a -> a.getAlignmentFragments().stream())
                                    .allMatch(af -> af.getDirection() == Direction.REVERSE)
        );

    }
}
