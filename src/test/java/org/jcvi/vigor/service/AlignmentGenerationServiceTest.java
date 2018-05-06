package org.jcvi.vigor.service;

import static org.junit.Assert.assertEquals;
import java.io.File;
import java.util.List;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.component.Alignment;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;


@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class AlignmentGenerationServiceTest {

	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader();



	@Test
	public void generateAlignmentsTest() throws VigorException {
		File virusGenomeSeqFile = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua.fasta"). getFile());
		File alignmentOutput = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua.txt"). getFile());
		String referenceDB = classLoader.getResource("vigorResources/data3/flua_db").getFile().toString();
        List<Alignment> alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile,referenceDB,alignmentOutput);
		assertEquals(11,alignments.size());
	}

}
