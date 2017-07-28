package com.vigor.RegressionTest;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;

import com.vigor.component.Alignment;
import com.vigor.component.AlignmentEvidence;
import com.vigor.component.Model;
import com.vigor.component.VirusGenome;
import com.vigor.forms.VigorForm;
import com.vigor.service.ExonerateService;
import com.vigor.service.ModelGenerationService;
import com.vigor.service.ViralProteinService;
import com.vigor.utils.VigorTestUtils;

public class generateVigor4Models {
		
	public static void main(String[] args){
		try{
		generateModels("C:/git/VIGOR4/src/test/resources/VigorRegressionTestOutput","C:/git/VIGOR4/src/test/resources/VigorRegressionTestInput/flua.fasta");
		}
		catch(Exception e){
			System.out.println(e.getMessage());
		}
	}
	
	public static void generateModels(String workspace,String inputFilePath) throws IOException, InterruptedException{
		File file = new File(workspace);
		if (!file.isDirectory())
		   file = file.getParentFile();
		if (file.exists()){
		 File inputFile = new File(inputFilePath) ;
		 if(inputFile.exists()){
		int count =0;
	HashMap<String,List<Model>> vigor4Models = new HashMap<String,List<Model>>();
    ExonerateService exonerateService = new ExonerateService();
    ViralProteinService viralProteinService = new ViralProteinService();
    ModelGenerationService modelGenerationService = new ModelGenerationService();
	GenerateExonerateOutputTest generateExonerateOutputTest = new GenerateExonerateOutputTest();
	NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(
			inputFile).hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
					.build();
	Stream<NucleotideFastaRecord> records = dataStore.records();
	Iterator<NucleotideFastaRecord> i = records.iterator();
	while (i.hasNext()) {
		NucleotideFastaRecord record = i.next();
		VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
				false, false);
		String refDB = "flua_db";
		String fileName = generateExonerateOutputTest.queryExonerate(virusGenome, refDB,file.getAbsolutePath());
		Thread.sleep(10000);
		File outputFile = new File(fileName);
		AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
		alignmentEvidence.setReference_db("flua_db");
		List<Alignment> alignments = exonerateService.parseExonerateOutput(outputFile, alignmentEvidence,virusGenome);
		alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment))
				.collect(Collectors.toList());
		List<Model> candidateModels = modelGenerationService.determineCandidateModels(alignments, new VigorForm());
		vigor4Models.put(virusGenome.getId(), candidateModels);
		System.out.println("counter"+count);
		      
	}
		
	Set<String> keys = vigor4Models.keySet();
	for(String key : keys ){
	System.out.println("Reference ID" + key +" "+vigor4Models.get(key).size());
	}}
		 else{
			 System.out.println("InputFile does not exist");
		 }
	}
	
	else{
		System.out.println("Workspace folder does not exit");
	}
	}
	
}
