package org.jcvi.vigor.RegressionTest;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.jcvi.vigor.RegressionTest.*;
import org.jcvi.vigor.component.Model;
import org.junit.Before;
import org.junit.Test;

public class ValidateVigor4Models {
	
	private Map<String,List<Model>> allVigor4Models = new HashMap<String,List<Model>>();
	private Map<String,List<Model>> allVigor3Models = new HashMap<String,List<Model>>();
	
	@Before
	public void getModels() throws IOException, InterruptedException{
		GenerateVigor3Models generateVigor3Models = new GenerateVigor3Models();
		allVigor3Models = generateVigor3Models.generateModels("/home/snettem/git/Vigor4/src/test/resources/vigor3Output","veev");
				
		GenerateVigor4Models generateVigor4Models = new GenerateVigor4Models();
		allVigor4Models = generateVigor4Models.generateModels("/home/snettem/workspace/vigor4RegressionOutput",
					"/home/snettem/git/Vigor4/src/test/resources/VigorRegressionTestInput/veev.fasta",
					"veev_db");
		
		
	}
	
	
	@Test
	public void validate(){
	  allVigor4Models.entrySet().forEach(entry -> { if(allVigor3Models.containsKey(entry.getKey())){
			  List<Model> vigor4Models = entry.getValue();
			  List<Model> vigor3Models = allVigor3Models.get(entry.getKey());
			  for(Model vigor3Model : vigor3Models){
				  String vigor3ProteinID = vigor3Model.getAlignment().getViralProtein().getProteinID();
				  String vigor3GenomeID = vigor3Model.getAlignment().getVirusGenome().getId();
				  boolean flag = false;
				 for(Model vigor4Model : vigor4Models ){
				  String vigor4ProteinID = vigor4Model.getAlignment().getViralProtein().getProteinID();
				  if(vigor3ProteinID.equals(vigor4ProteinID)){
					  flag= true;
				  }
				 }
				if(!flag){
					 System.out.println(vigor3ProteinID+" Model missing for VirusGenome ID"+vigor3GenomeID);
				 }
				 flag = false;
			  }
		      
	          }
		  		  
	  });
	  				
	}

}