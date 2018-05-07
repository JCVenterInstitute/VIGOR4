package org.jcvi.vigor.service;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.notNullValue;
import static org.jcvi.vigor.service.CommandLineParameters.referenceDB;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;

import java.io.File;
import java.nio.file.Paths;
import java.util.List;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.component.Alignment;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;


@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class AlignmentGenerationServiceTest {

	@Autowired
	private VigorInitializationService initializationService;

	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader();



	@Test
	public void generateAlignmentsTest() throws VigorException {
		File virusGenomeSeqFile = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua.fasta"). getFile());
		File alignmentOutput = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua.txt"). getFile());
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String refereceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        assertThat("reference database path is required", refereceDBPath, is(notNullValue()));
        String referenceDB = Paths.get(refereceDBPath, "flua_db").toString();
        List<Alignment> alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile,referenceDB,alignmentOutput, config);
		assertEquals(11,alignments.size());
	}

}
