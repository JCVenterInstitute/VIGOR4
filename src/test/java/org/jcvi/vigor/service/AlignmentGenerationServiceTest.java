package org.jcvi.vigor.service;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.notNullValue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;

import java.io.File;
import java.nio.file.Paths;
import java.util.List;

import org.jcvi.vigor.Application;
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
    public void generateAlignmentsTest () throws VigorException {

        File resources = new File("src/test/resources");
        File virusGenomeSeqFile = new File(resources.getAbsolutePath() + File.separator + "vigorUnitTestInput/sequence_flua.fasta");
        File alignmentOutput = new File(resources.getAbsolutePath() + File.separator + "vigorUnitTestInput/sequence_flua.txt");
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String refereceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(refereceDBPath);
        assertThat("reference database path is required", refereceDBPath, is(notNullValue()));
        String referenceDB = Paths.get(refereceDBPath, "flua_db").toString();
        List<Alignment> alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile, referenceDB, alignmentOutput, config);
        assertEquals(11, alignments.size());
    }
}
