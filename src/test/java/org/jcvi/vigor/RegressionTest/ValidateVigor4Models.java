package org.jcvi.vigor.RegressionTest;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jcvi.vigor.AppConfig;
import org.jcvi.vigor.component.Model;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

@RunWith(SpringRunner.class)
@ContextConfiguration(classes = AppConfig.class)
public class ValidateVigor4Models {

	@Autowired
	GenerateVigor4GeneModels generateVigor4GeneModels;
	private Map<String,List<Model>> allVigor4Models = new HashMap<String,List<Model>>();
	private Map<String,List<Model>> allVigor3Models = new HashMap<String,List<Model>>();
	
	@Test
	public void getModels() throws IOException, InterruptedException{
	/*	GenerateVigor3Models generateVigor3Models = new GenerateVigor3Models();
		allVigor3Models = generateVigor3Models.generateModels("/Users/snettem/git/Vigor4/src/test/resources/vigor3Output","veev");
				*/
		allVigor4Models = generateVigor4GeneModels.generateModels("/home/snettem/workspace/vigor4_regression",
					"/home/snettem/git/VIGOR4/src/test/resources/VigorRegressionTestInput/flu.fasta",
					"flua_db");
		
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
