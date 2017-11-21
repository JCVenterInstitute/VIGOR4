package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.utils.GenerateExonerateOutput;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.ExonerateService;
import org.jcvi.vigor.service.ModelGenerationService;
import org.jcvi.vigor.service.ViralProteinService;
import org.jcvi.vigor.utils.FormatVigorOutput;

public class GenerateVigor4Models {

	public Map<String, List<Model>> generateModels(String workspace,
			String inputFilePath, String refDB) 
			 {
		Map<String, List<Model>> vigor4Models = new HashMap<String, List<Model>>();
		try{
		File file = new File(workspace);
		if (!file.isDirectory())
			file = file.getParentFile();
		if (file.exists()) {
			File inputFile = new File(inputFilePath);
			if (inputFile.exists()) {
				ExonerateService exonerateService = new ExonerateService();
				ViralProteinService viralProteinService = new ViralProteinService();
				ModelGenerationService modelGenerationService = new ModelGenerationService();
				NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(
						inputFile).hint(
						DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
						.build();
				Stream<NucleotideFastaRecord> records = dataStore.records();
				Iterator<NucleotideFastaRecord> i = records.iterator();
				while (i.hasNext()) {
					NucleotideFastaRecord record = i.next();
					VirusGenome virusGenome = new VirusGenome(
							record.getSequence(), record.getComment(),
							record.getId(), false, false);

					String fileName = GenerateExonerateOutput.queryExonerate(
							virusGenome, refDB, file.getAbsolutePath(),null);
					File outputFile = new File(fileName);
					AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
					alignmentEvidence.setReference_db(refDB);
					List<Alignment> alignments = exonerateService
							.parseExonerateOutput(outputFile,
									alignmentEvidence, virusGenome);
					alignments = alignments
							.stream()
							.map(alignment -> viralProteinService
									.setViralProteinAttributes(alignment, new VigorForm()))
							.collect(Collectors.toList());
					List<Model> candidateModels = modelGenerationService
							.determineCandidateModels(alignments,
									new VigorForm());
					vigor4Models.put(virusGenome.getId(), candidateModels);
				//	FormatVigorOutput.printModels(candidateModels);
					/*System.out.println("outputFile" + outputFile);
					System.out.println("actual one in the end:"
							+ virusGenome.getId());
					System.out.println("counter" + count);*/
				
				}

				/*
				 * Set<String> keys = vigor4Models.keySet(); for(String key :
				 * keys ){ System.out.println("Reference ID" + key
				 * +" "+vigor4Models.get(key).size()); }
				 */
			} else {
				System.out.println("InputFile does not exist");
			}
		}

		else {
			System.out.println("Workspace folder does not exit");
		}}
		catch(IOException e){
			System.out.println(e.getMessage());		
		}
		catch(Exception e){
			System.out.println(e.getMessage());	
		}
		
		return vigor4Models;

	}

}
