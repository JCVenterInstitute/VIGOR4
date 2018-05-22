package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.DetermineMissingExonsTest;
import org.jcvi.vigor.service.VigorInitializationService;
import org.jcvi.vigor.utils.*;
import org.junit.ClassRule;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ErrorCollector;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.test.context.junit4.rules.SpringClassRule;
import org.springframework.test.context.junit4.rules.SpringMethodRule;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertFalse;
import static org.hamcrest.CoreMatchers.equalTo;

@RunWith(Parameterized.class)
@ContextConfiguration(classes = Application.class)
public class ValidateVigor4Models {

	@Autowired
	GenerateVigor4GeneModels generateVigor4GeneModels;

	private Map<String,List<Model>> allVigor3Models = new HashMap<String,List<Model>>();
    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4Models.class);
    @Autowired
    private VigorInitializationService initializationService;
    private String vigor3OutputTBL;
    private String vigor3OutputPep;
    private String inputFasta;
    private String workspace;
    private String refDB;
    @ClassRule
    public static final SpringClassRule springClassRule = new SpringClassRule();

    @Rule
    public final SpringMethodRule springMethodRule = new SpringMethodRule();

    @Rule
    public ErrorCollector collector = new ErrorCollector();

    public ValidateVigor4Models(String vigor3OutputTBL, String vigor3OutputPep, String inputFasta, String workspace,String refDB) {
        this.vigor3OutputTBL = vigor3OutputTBL;
        this.vigor3OutputPep = vigor3OutputPep;
        this.inputFasta = inputFasta;
        this.workspace = workspace;
        this.refDB = refDB;
    }

    private void isAppPackagedCorrectly(String res) throws VigorException{

        if(ValidateVigor4Models.class.getClassLoader().getResource(res)==null){
            throw new VigorException("Missing regression test resource  "+res);
        }
    }

    @Parameterized.Parameters(name="ValidateVigor4Models[#{index} {0}]")
    public static Collection<Object[]> getTestData() {
        List<Object[]> testData = new ArrayList<>();
        testData.add(
                new Object[]{
                        "vigor3Output/flu/PRODUCTION/flu.tbl",
                        "vigor3Output/flu/PRODUCTION/flu.pep",
                        "vigorRegressionTestInput/flua.fasta",
                        "vigor4Output","flua_db"
                });

     return testData;

    }
    private void setAbsolutePathOfTestResource() throws VigorException{
        File resources = new File("src/test/resources");
        this.vigor3OutputTBL = resources.getAbsolutePath()+File.separator+vigor3OutputTBL;
        this.vigor3OutputPep = resources.getAbsolutePath()+File.separator+vigor3OutputPep;
        this.inputFasta = resources.getAbsolutePath()+File.separator+inputFasta;
        this.workspace = resources.getAbsolutePath()+File.separator+workspace;
    }

    public Map<String,List<Model>> getVigor4Models() throws VigorException{
         VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
         String referenceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
         String referenceDB = Paths.get(referenceDBPath, refDB).toString();
         if(referenceDB==null){
             LOGGER.error("{} reference db does not exist",referenceDB);
         }
         refDB = referenceDB;
        Map<String,List<Model>> vigor4Models = new HashMap<String,List<Model>>();
        isAppPackagedCorrectly(inputFasta);
        isAppPackagedCorrectly(workspace);
        setAbsolutePathOfTestResource();
        config.put(ConfigurationParameters.OutputDirectory,workspace);
        config.put(ConfigurationParameters.Verbose,"false");
        vigor4Models = generateVigor4GeneModels.generateModels(inputFasta,refDB,config);
        return vigor4Models;
    }

   public Map<String,List<Model>> getVigor3Models() throws IOException{
        GenerateVigor3Models generateVigor3Models = new GenerateVigor3Models();
		allVigor3Models = generateVigor3Models.generateModels(vigor3OutputTBL,vigor3OutputPep,inputFasta);
        return allVigor3Models;
	}

	@Test
	public void validate() throws IOException,VigorException{
       Map<String,List<Model>> allVigor4Models = getVigor4Models();
       Map<String,List<Model>> allVigor3Models = getVigor3Models();
	  allVigor4Models.entrySet().forEach(entry -> { if(allVigor3Models.containsKey(entry.getKey())){
			  List<Model> vigor4Models = entry.getValue();
			  List<Model> vigor3Models = allVigor3Models.get(entry.getKey());
			  for(Model vigor3Model : vigor3Models){
				  String vigor3ProteinID = vigor3Model.getAlignment().getViralProtein().getProteinID();
				  String vigor3GenomeID = vigor3Model.getAlignment().getVirusGenome().getId();
				  List<Exon> vigor3Exons = vigor3Model.getExons();
				  boolean flag = false;
				 for(Model vigor4Model : vigor4Models ){
				  String vigor4ProteinID = vigor4Model.getAlignment().getViralProtein().getProteinID();
				  List<Exon> vigor4Exons = vigor4Model.getExons();
				  if(vigor3ProteinID.equals(vigor4ProteinID)){
					  flag= true;
					  if(vigor4Exons.size()!=vigor3Exons.size()){
					      System.out.println("Break");
                      }
					  collector.checkThat(String.format("Exon count differs for Vigor3 & Vigor4 model with gene ID %s of VirusGenome Sequence %s",vigor3ProteinID,vigor3GenomeID),vigor4Exons.size(),equalTo(vigor3Exons.size()));
					  for(int i=0;i<vigor3Exons.size();i++){
					      Range temp = vigor4Exons.get(i).getRange();
					      Range vigor4ExonRange = Range.of(temp.getBegin()+1,temp.getEnd()+1);
					      collector.checkThat(String.format("Exon range differs for Vigor3 & Vigor4 model with gene ID %s of VirusGenome Sequence %s",vigor3ProteinID,vigor3GenomeID),vigor4ExonRange,equalTo(vigor3Exons.get(i).getRange()));
                      }
				  }
				 }

				collector.checkThat(String.format("Vigor3 & Vigor4 gene models do not match for VirusGenome Sequence %s.Expected gene ID %s ",vigor3GenomeID,vigor3ProteinID),flag,equalTo(true) );

			  }
		      
	          }
		  		  
	  });
	  				
	}

}
