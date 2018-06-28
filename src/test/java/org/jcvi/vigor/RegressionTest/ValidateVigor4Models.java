package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.io.FileWriter;
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
import org.jcvi.vigor.service.VigorInitializationService;
import org.jcvi.vigor.utils.*;
import org.junit.ClassRule;
import org.junit.Rule;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.rules.SpringClassRule;
import org.springframework.test.context.junit4.rules.SpringMethodRule;

@RunWith(Parameterized.class)
@ContextConfiguration(classes = Application.class)
public class ValidateVigor4Models {

    private static String outputDirectory;
    private String virusSpecifcDirectory;
    private String vigor3OutputTBL;
    private String vigor3OutputPep;
    private String inputFasta;
    private String refDB;
    @Autowired
    GenerateVigor4GeneModels generateVigor4GeneModels;
    private Map<String, List<Model>> allVigor3Models = new HashMap<String, List<Model>>();
    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4Models.class);
    @Autowired
    private VigorInitializationService initializationService;
    @ClassRule
    public static final SpringClassRule springClassRule = new SpringClassRule();
    @Rule
    public final SpringMethodRule springMethodRule = new SpringMethodRule();

    public ValidateVigor4Models ( String vigor3OutputTBL, String vigor3OutputPep, String inputFasta, String refDB ) {

        this.vigor3OutputTBL = vigor3OutputTBL;
        this.vigor3OutputPep = vigor3OutputPep;
        this.inputFasta = inputFasta;
        this.refDB = refDB;
    }

    public static void prepare ( String outputDir ) {

        ValidateVigor4Models.outputDirectory = outputDir;
    }

    private void isAppPackagedCorrectly ( String res ) throws VigorException {

        if (ValidateVigor4Models.class.getClassLoader().getResource(res) == null) {
            throw new VigorException("Missing regression test resource  " + res);
        }
    }

    @Parameterized.Parameters(name = "ValidateVigor4Models[#{index} {0}]")
    public static Collection<Object[]> getTestData () {

        List<Object[]> testData = new ArrayList<>();
        testData.add(
                new Object[] {
                        "vigor3Output/flu/PRODUCTION/flu.tbl",
                        "vigor3Output/flu/PRODUCTION/flu.pep",
                        "vigorRegressionTestInput/flua.fasta",
                        "flua_db",
                });
        testData.add(
                new Object[] {
                        "vigor3Output/veev/PRODUCTION/veev.tbl",
                        "vigor3Output/veev/PRODUCTION/veev.pep",
                        "vigorRegressionTestInput/veev.fasta",
                        "veev_db",
                });
        testData.add(
                new Object[] {
                        "vigor3Output/veev/PRODUCTION/rsv.tbl",
                        "vigor3Output/veev/PRODUCTION/rsv.pep",
                        "vigorRegressionTestInput/rsv.fasta",
                        "rsv_db",
                });
        return testData;
    }

    private void setAbsolutePathOfTestResource () {

        File resources = new File("src/test/resources");
        this.vigor3OutputTBL = resources.getAbsolutePath() + File.separator + vigor3OutputTBL;
        this.vigor3OutputPep = resources.getAbsolutePath() + File.separator + vigor3OutputPep;
        this.inputFasta = resources.getAbsolutePath() + File.separator + inputFasta;
    }

    public Map<String, List<Model>> getVigor4Models () throws VigorException {

        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String referenceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        String referenceDB = Paths.get(referenceDBPath, refDB).toString();
        if (referenceDB == null) {
            LOGGER.error("{} reference db does not exist", referenceDB);
        }
        refDB = referenceDB;
        Map<String, List<Model>> vigor4Models;
        isAppPackagedCorrectly(inputFasta);
        isAppPackagedCorrectly(vigor3OutputPep);
        isAppPackagedCorrectly(vigor3OutputTBL);
        setAbsolutePathOfTestResource();
        File refDBFile = new File(refDB);
        String outputPrefix = refDBFile.getName().replace("_db", "");
        config.put(ConfigurationParameters.OutputPrefix, outputPrefix);
        String virusSpecificPath = outputDirectory + File.separator + config.get(ConfigurationParameters.OutputPrefix);
        File virusSpecificDir = new File(virusSpecificPath);
        if (!virusSpecificDir.exists()) virusSpecificDir.mkdir();
        virusSpecifcDirectory = virusSpecificDir.getAbsolutePath();
        config.put(ConfigurationParameters.OutputDirectory, virusSpecifcDirectory);
        config.put(ConfigurationParameters.Verbose, "false");
        vigor4Models = generateVigor4GeneModels.generateModels(inputFasta, refDB, config);
        return vigor4Models;
    }

    public Map<String, List<Model>> getVigor3Models () throws IOException {

        GenerateVigor3Models generateVigor3Models = new GenerateVigor3Models();
        allVigor3Models = generateVigor3Models.generateModels(vigor3OutputTBL, vigor3OutputPep, inputFasta);
        return allVigor3Models;
    }

    @Test
    public void validate () {

        try {
            Map<String, List<Model>> allVigor4Models = getVigor4Models();
            Map<String, List<Model>> allVigor3Models = getVigor3Models();
            File errorReport = new File(virusSpecifcDirectory, "differencesReport.txt");
            final FileWriter writer = new FileWriter(errorReport, false);
            writer.write("***********************Differences Report**************************\n");
            allVigor3Models.entrySet().forEach(entry -> {
                if (allVigor4Models.containsKey(entry.getKey())) {
                    try {
                        List<Model> vigor3Models = entry.getValue();
                        List<Model> vigor4Models = allVigor4Models.get(entry.getKey());
                        boolean failed = false;
                        for (Model vigor3Model : vigor3Models) {
                            String vigor3GeneID = vigor3Model.getGeneSymbol();
                            String vigor3GenomeID = vigor3Model.getAlignment().getVirusGenome().getId();
                            List<Exon> vigor3Exons = vigor3Model.getExons();
                            boolean flag = false;
                            for (Model vigor4Model : vigor4Models) {
                                List<Exon> vigor4Exons = vigor4Model.getExons();
                                if (vigor3Model.getGeneSymbol().equals(vigor4Model.getGeneSymbol())) {
                                    flag = true;
                                    if (vigor4Model.isPartial3p() != vigor3Model.isPartial3p()) {
                                        writer.write(String.format("Partial 3' gene feature mismatch found for gene symbol %s of VirusGenome Sequence %s. partial 3' %s expected \n", vigor3GeneID, vigor3GenomeID, vigor3Model.isPartial3p()));
                                        failed = true;
                                    }
                                    if (vigor3Model.isPartial5p() != vigor4Model.isPartial5p()) {
                                        writer.write(String.format("Partial 5' gene feature mismatch found for gene symbol %s of VirusGenome Sequence %s. partial 5' %s expected \n", vigor3GeneID, vigor3GenomeID, vigor3Model.isPartial5p()));
                                        failed = true;
                                    }
                                    if (( vigor3Model.getRibosomalSlippageRange() != null && vigor4Model.getRibosomalSlippageRange() == null ) || ( vigor3Model.getRibosomalSlippageRange() == null && vigor4Model.getRibosomalSlippageRange() != null )) {
                                        writer.write(String.format("Ribosomal Slippage mismatch found for gene %s of VirusGenome Sequence \n", vigor3GeneID, vigor3GenomeID));
                                        failed = true;
                                    }
                                    if (( vigor3Model.getRibosomalSlippageRange() != null && vigor4Model.getRibosomalSlippageRange() == null ) || ( vigor3Model.getRibosomalSlippageRange() == null && vigor4Model.getRibosomalSlippageRange() != null )) {
                                        writer.write(String.format("Ribosomal Slippage mismatch found for gene %s of VirusGenome Sequence \n", vigor3GeneID, vigor3GenomeID));
                                        failed = true;
                                    }
                                    if (( vigor3Model.getReplaceStopCodonRange() != null && vigor4Model.getReplaceStopCodonRange() != null ) && ( vigor3Model.getReplaceStopCodonRange() != vigor4Model.getReplaceStopCodonRange() )) {
                                        writer.write(String.format("Stop Codon Readthrough feature mismatch for gene %s of VirusGenome Sequence \n", vigor3GeneID, vigor3GenomeID));
                                        failed = true;
                                    }
                                    if (vigor3Model.isPseudogene() != vigor4Model.isPseudogene()) {
                                        writer.write(String.format("Pseudogene feature mismatch found for gene %s of VirusGenome Sequence %s. pseudogene %s expected \n", vigor3GeneID, vigor3GenomeID, vigor3Model.isPseudogene()));
                                        failed = true;
                                    }
                                    if (vigor4Exons.size() != vigor3Exons.size()) {
                                        writer.write(String.format("Exon count differs for Vigor3 & Vigor4 model with gene symbol %s of VirusGenome Sequence %s \n", vigor3GeneID, vigor3GenomeID));
                                        failed = true;
                                    }
                                    if (vigor4Exons.size() == vigor3Exons.size()) {
                                        for (int i = 0; i < vigor3Exons.size(); i++) {
                                            Range temp = vigor4Exons.get(i).getRange();
                                            Range vigor4ExonRange = Range.of(temp.getBegin() + 1, temp.getEnd() + 1);
                                            if (vigor4ExonRange != vigor3Exons.get(i).getRange()) {
                                                writer.write(String.format("Exon range differs for Vigor3 & Vigor4 model with gene symbol %s of VirusGenome Sequence %s \n", vigor3GeneID, vigor3GenomeID));
                                                failed = true;
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                            if (!flag) {
                                writer.write(String.format("Vigor3 & Vigor4 gene models do not match for VirusGenome Sequence %s.Expected gene symbol %s \n", vigor3GenomeID, vigor3GeneID));
                                failed = true;
                            }
                        }
                        if (failed) {
                            writer.write(">Genome " + entry.getKey());
                            writer.write("\n***************************************************************************************************************************************************" + "\n");
                        }
                    } catch (IOException e) {
                        LOGGER.error(e);
                        System.exit(1);
                    }
                }
            });
            writer.close();
        } catch (VigorException e) {
            LOGGER.error(e);
            System.exit(1);
        } catch (Exception e) {
            LOGGER.error(e);
            System.exit(1);
        }
    }
}
