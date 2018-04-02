package org.jcvi.vigor.service;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.AppConfig;
import org.jcvi.vigor.service.AlignmentGenerationService;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.VirusGenome;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;


@RunWith(SpringRunner.class)
@ContextConfiguration(classes = AppConfig.class)
public class AlignmentGenerationServiceTest {

	@Autowired
	private AlignmentGenerationService alignmentGenerationService;
	private VirusGenome virusGenome = new VirusGenome();
	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader();
	private File file = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua.fasta"). getFile());

	@Before
	public void getVirusGenome() throws IOException{

		NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(
				file).hint(
				DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
				.build();
		Optional<NucleotideFastaRecord> record = dataStore.records().findFirst();
		NucleotideFastaRecord seqRecord = record.get();
		virusGenome.setDefline(seqRecord.getComment());
		virusGenome.setId(seqRecord.getId());
		virusGenome.setSequence(seqRecord.getSequence());

	}

	@Test
	public void generateAlignmentsTest() {
		List<Alignment> alignments = alignmentGenerationService.generateExonerateAlignment(virusGenome,new AlignmentEvidence("flua_db"));
		assertEquals(4,alignments.size());
	}

}
