package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.utils.TBLParser.TBLFileParser;
import org.jcvi.vigor.utils.TBLParser.TBLModel;

public class GenerateVigor3Models {

	public Map<String,List<Model>> generateModels(String Vigor3OutputFolder,String referenceDBFolder){
		Map<String, List<Model>> vigor3Models = new HashMap<String, List<Model>>();
		String TBLFilePath = Vigor3OutputFolder+File.separator+referenceDBFolder+File.separator+"PRODUCTION"+File.separator+referenceDBFolder+".tbl";
		try{
		NucleotideFastaDataStore datastore = new NucleotideFastaFileDataStoreBuilder(new File("/home/snettem/git/Vigor4/src/test/resources/VigorRegressionTestInput/regr_testing/flu.fasta")).build();
		System.out.println("Number of records in the fasta file are : "+ datastore.getNumberOfRecords());
		}
		catch(Exception e)
		{
			System.out.println(e.getMessage());
		}
		String PEPFilePath = Vigor3OutputFolder+File.separator+referenceDBFolder+File.separator+"PRODUCTION"+File.separator+referenceDBFolder+".pep";
		TBLFileParser TBLParser = new TBLFileParser();
		List<TBLModel> TBLModels = TBLParser.getModels(TBLFilePath, PEPFilePath);
		System.out.println("Total Number of models are :"+TBLModels.size() );
		for(TBLModel tblModel : TBLModels){
			Model model = new Model();
			List<Model> models = new ArrayList<Model>();
			Alignment alignment = new Alignment();
			ViralProtein viralProtein = new ViralProtein();
			VirusGenome virusGenome = new VirusGenome();
			model.setGeneSymbol(tblModel.getGene());
			viralProtein.setProteinID(tblModel.getViralProteinID());
			alignment.setViralProtein(viralProtein);
			virusGenome.setId(tblModel.getVirusGenomeID());
			alignment.setVirusGenome(virusGenome);
			model.setAlignment(alignment);
			model.setExons(tblModel.getExons());
			String virusGenomeID = model.getAlignment().getVirusGenome().getId();
			if(vigor3Models.containsKey(virusGenomeID)){
			models = vigor3Models.get(virusGenomeID);
			}
		    models.add(model);
		    vigor3Models.put(virusGenomeID, models);
		    
		}
		System.out.println(vigor3Models.entrySet().size());
		vigor3Models.entrySet().forEach(entry -> {System.out.println("key:"+ entry.getKey()+"value:"+entry.getValue().size());
			
		});
		return vigor3Models;
	}
	
	
}
