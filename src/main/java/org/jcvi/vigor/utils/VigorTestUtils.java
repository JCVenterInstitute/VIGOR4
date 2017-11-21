package org.jcvi.vigor.utils;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.ExonerateService;
import org.jcvi.vigor.service.ViralProteinService;
import org.jcvi.vigor.service.VirusGenomeService;

public class VigorTestUtils {

	private static final Logger LOGGER = LogManager.getLogger(VigorTestUtils.class);
	
	
	
	public static List<Alignment> getAlignments(String inputFilePath, String refDB,String workspace,String proteinID) {
       
		NucleotideFastaDataStore dataStore;
		List<Alignment> alignments = new ArrayList<Alignment>();
		ExonerateService exonerateService = new ExonerateService();
		ViralProteinService viralProteinService = new ViralProteinService();
		try {
			File file = new File(inputFilePath);
			dataStore = new NucleotideFastaFileDataStoreBuilder(file)
					.hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED).build();
			Stream<NucleotideFastaRecord> records = dataStore.records();
			Iterator<NucleotideFastaRecord> i = records.iterator();
			NucleotideFastaRecord record = i.next();
			VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(), false,
					false);
			List<Range> sequenceGaps = VirusGenomeService.findSequenceGapRanges("20",
					virusGenome.getSequence());
			virusGenome.setSequenceGaps(sequenceGaps);
			// create alignment evidence
			AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
			alignmentEvidence.setReference_db(refDB);
			
			String fileName = GenerateExonerateOutput.queryExonerate(
					virusGenome, refDB, workspace,proteinID);
			File outputFile = new File(fileName);
			alignments = exonerateService
					.parseExonerateOutput(outputFile,
							alignmentEvidence, virusGenome);
			alignments = alignments
					.stream()
					.map(alignment -> viralProteinService
							.setViralProteinAttributes(alignment,new VigorForm()))
					.collect(Collectors.toList());
						
			//System.out.println("Number of alignments are :" + alignments.size());

		} catch (Exception e) {
			LOGGER.error(e.getMessage(), e);
		}
		return alignments;
	}
	
		
	
	

}
