/*
package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.utils.LoadDefaultParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class ValidateVigor4Models {

	@Autowired
	GenerateVigor4GeneModels generateVigor4GeneModels;
	private Map<String,List<Model>> allVigor4Models = new HashMap<String,List<Model>>();
	private Map<String,List<Model>> allVigor3Models = new HashMap<String,List<Model>>();
    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4Models.class);
	@Test
	public void getModels() throws IOException,VigorException{
		//Todo; inputs from user/developer to run regression test
        String vigor3TblFile = "/home/snettem/git/VIGOR4/src/test/resources/vigor3Output/flu/PRODUCTION/flu.tbl";
        File vigor3OutputTbl = new File(vigor3TblFile);
        if(!vigor3OutputTbl.exists() || !vigor3OutputTbl.isFile()){
           LOGGER.error("{} file does not exist. Please input valid vigor3 output tbl file", vigor3TblFile);
           System.exit(1);
        }
        String vigor3PepFile = "/home/snettem/git/VIGOR4/src/test/resources/vigor3Output/flu/PRODUCTION/flu.pep";
        File vigor3OutputPep =  new File(vigor3PepFile);
        if(!vigor3OutputPep.exists() || !vigor3OutputPep.isFile()){
            LOGGER.error("{} file does not exist. Please input valid vigor3 output pep file");
            System.exit(1);
        }
        String vigor4InputFile = "/home/snettem/git/VIGOR4/src/test/resources/VigorRegressionTestInput/flu.fasta";
        String workSpace = "/home/snettem/workspace/vigor4_regression";
        File workSpacDir = new File(workSpace);
        if(!workSpacDir.exists()|| !workSpacDir.isDirectory()){
            LOGGER.error("{} workspace directory does not exists",workSpacDir);
            System.exit(1);
        }
        File inputFASTAFile = new File(vigor4InputFile);
        if(!inputFASTAFile.exists()||!inputFASTAFile.isFile()){
            LOGGER.error("{} input fasta file does not exist",vigor4InputFile);
            System.exit(1);
        }
        String referenceDB = "flua_db";
        String referenceDB_Path = "";
        File refDBPath = new File(referenceDB_Path);
        if(!refDBPath.isDirectory() || referenceDB.isEmpty() ){
            ClassLoader classLoader = VigorTestUtils.class.getClassLoader();
            referenceDB = classLoader.getResource("vigorResources/data3/"+referenceDB).getFile().toString();
        }else{
            referenceDB=refDBPath+File.separator+referenceDB;
        }
        File referenceDBFile = new File(referenceDB);
        if(!referenceDBFile.exists()){
            LOGGER.error("{} reference db does not exist",referenceDBFile);
        }

		GenerateVigor3Models generateVigor3Models = new GenerateVigor3Models();
		allVigor3Models = generateVigor3Models.generateModels(vigor3OutputTbl.getAbsolutePath(),vigor3OutputPep.getAbsolutePath(),inputFASTAFile.getAbsolutePath());
		allVigor4Models = generateVigor4GeneModels.generateModels(workSpacDir.getAbsolutePath(),
					inputFASTAFile.getAbsolutePath(),
					referenceDBFile.getAbsolutePath());
		
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
*/
