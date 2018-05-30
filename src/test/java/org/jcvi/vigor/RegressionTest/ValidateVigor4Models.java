package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.assertj.core.api.SoftAssertions;
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
                        "vigor4Output","flua_db",
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
        File errorReport = new File(workspace+File.separator+"ErrorReport.txt");
        final FileWriter writer=new FileWriter(errorReport);
        writer.write("***********************Error Report**************************\n");
        allVigor3Models.entrySet().forEach(entry -> { if(allVigor4Models.containsKey(entry.getKey())) {
             try{
            List<Model> vigor3Models = entry.getValue();
            List<Model> vigor4Models = allVigor4Models.get(entry.getKey());
            boolean failed=false;
            for (Model vigor3Model : vigor3Models) {
                SoftAssertions softly = new SoftAssertions();
                String vigor3GeneID = vigor3Model.getGeneID();
                String vigor3GenomeID = vigor3Model.getAlignment().getVirusGenome().getId();
                List<Exon> vigor3Exons = vigor3Model.getExons();
                boolean flag = false;
                for (Model vigor4Model : vigor4Models) {
                    List<Exon> vigor4Exons = vigor4Model.getExons();
                    if (vigor3Model.getGeneID().equals(vigor4Model.getGeneID())) {
                        flag = true;
                        if (vigor4Model.isPartial3p() != vigor3Model.isPartial3p()) {
                            writer.write(String.format("Partial 3' gene feature mismatch found for gene symbol %s of VirusGenome Sequence %s. partial 3' %s expected \n", vigor3GeneID, vigor3GenomeID, vigor3Model.isPartial3p()));
                            failed = true;
                        }
                        /*if(vigor3Model.getAlignment().getViralProtein().getProduct() != vigor4Model.getAlignment().getViralProtein().getProduct()) {
                            writer.write(String.format("product id mismatch found for gene symbol %s of VirusGenome Sequence %s. partial 5' %s expected \n", vigor3GeneID, vigor3GenomeID, vigor3Model.isPartial5p()));
                            failed = true;
                        }*/

                        if(vigor3Model.isPartial5p() != vigor4Model.isPartial5p()) {
                            writer.write(String.format("Partial 5' gene feature mismatch found for gene symbol %s of VirusGenome Sequence %s. partial 5' %s expected \n", vigor3GeneID, vigor3GenomeID, vigor3Model.isPartial5p()));
                            failed = true;
                        }
                        if((vigor3Model.getRibosomalSlippageRange()!=null && vigor4Model.getRibosomalSlippageRange()==null)||(vigor3Model.getRibosomalSlippageRange()==null && vigor4Model.getRibosomalSlippageRange()!=null)){
                            writer.write(String.format("Ribosomal Slippage mismatch found for gene %s of VirusGenome Sequence \n", vigor3GeneID, vigor3GenomeID));
                            failed=true;
                        }
                        if((vigor3Model.getRibosomalSlippageRange()!=null && vigor4Model.getRibosomalSlippageRange()==null)||(vigor3Model.getRibosomalSlippageRange()==null && vigor4Model.getRibosomalSlippageRange()!=null)){
                            writer.write(String.format("Ribosomal Slippage mismatch found for gene %s of VirusGenome Sequence \n", vigor3GeneID, vigor3GenomeID));
                            failed=true;
                        }
                        if((vigor3Model.getReplaceStopCodonRange()!= null && vigor4Model.getReplaceStopCodonRange()!=null) &&(vigor3Model.getReplaceStopCodonRange()!=vigor4Model.getReplaceStopCodonRange())){
                            writer.write(String.format("Stop Codon Readthrough feature mismatch for gene %s of VirusGenome Sequence \n", vigor3GeneID, vigor3GenomeID));
                            failed=true;
                        }
                        if(vigor3Model.isPseudogene()!=vigor4Model.isPseudogene()) {
                            writer.write(String.format("Pseudogene feature mismatch found for gene %s of VirusGenome Sequence %s. pseudogene %s expected \n", vigor3GeneID, vigor3GenomeID, vigor3Model.isPseudogene()));
                            failed = true;
                        }
                        if(vigor4Exons.size()!=vigor3Exons.size()) {
                            writer.write(String.format("Exon count differs for Vigor3 & Vigor4 model with gene symbol %s of VirusGenome Sequence %s \n", vigor3GeneID, vigor3GenomeID));
                            failed = true;
                        }
                        if(vigor4Exons.size()==vigor3Exons.size()){
                        for (int i = 0; i < vigor3Exons.size(); i++) {
                            Range temp = vigor4Exons.get(i).getRange();
                            Range vigor4ExonRange = Range.of(temp.getBegin() + 1, temp.getEnd() + 1);
                            if(vigor4ExonRange!=vigor3Exons.get(i).getRange()) {
                                writer.write(String.format("Exon range differs for Vigor3 & Vigor4 model with gene symbol %s of VirusGenome Sequence %s \n", vigor3GeneID, vigor3GenomeID));
                                failed = true;
                            }
                        }}
                        break;
                    }
                }
                if(!flag) {
                    writer.write(String.format("Vigor3 & Vigor4 gene models do not match for VirusGenome Sequence %s.Expected gene symbol %s \n", vigor3GenomeID, vigor3GeneID));
                    failed=true;
                }

            } if(failed) {
                writer.write(">Genome "+entry.getKey());
                writer.write("\n***************************************************************************************************************************************************" + "\n");
                 }
            }catch(IOException e){e.printStackTrace();}
        }

	  });
        writer.close();
	  				
	}

}
