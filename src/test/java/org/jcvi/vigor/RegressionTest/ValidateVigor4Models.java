package org.jcvi.vigor.RegressionTest;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
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
                        "vigor3Output/flua/flua.tbl",
                        "vigor3Output/flua/flua.pep",
                        "vigor3Output/flua/flua.fasta",
                        "flua_db",
                });
        testData.add(
                new Object[] {
                        "vigor3Output/veev/veev.tbl",
                        "vigor3Output/veev/veev.pep",
                        "vigor3Output/veev/veev.fasta",
                        "veev_db",
                });
        testData.add(
                new Object[] {
                        "vigor3Output/rsv/rsv.tbl",
                        "vigor3Output/rsv/rsv.pep",
                        "vigor3Output/rsv/rsv.fasta",
                        "rsv_db",
                });
        return testData;
    }

    private void setAbsolutePathOfTestResource () {
        this.vigor3OutputTBL = ValidateVigor4Models.class.getClassLoader().getResource(vigor3OutputTBL).getFile();
        this.vigor3OutputPep = ValidateVigor4Models.class.getClassLoader().getResource(vigor3OutputPep).getFile();
        this.inputFasta = ValidateVigor4Models.class.getClassLoader().getResource(inputFasta).getFile();
    }

    public Map<String, List<Model>> getVigor4Models () throws VigorException, IOException {

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
        if (config.get(ConfigurationParameters.OverwriteOutputFiles) == null) {
            config.put(ConfigurationParameters.OverwriteOutputFiles, "true");
        }
        // TODO allow user to set this
        if (config.get(ConfigurationParameters.TemporaryDirectory) == null) {
            final Path tempDir = Files.createTempDirectory(Paths.get(outputDirectory), "vigor-tmp");
            // delete on shutdown
            Runtime.getRuntime().addShutdownHook(new Thread(() -> VigorUtils.deleteDirectory(tempDir)));
            config.put(ConfigurationParameters.TemporaryDirectory, tempDir.toString());
        }
        vigor4Models = generateVigor4GeneModels.generateModels(inputFasta, refDB, config);
        return vigor4Models;
    }

    public Map<String, List<Model>> getVigor3Models () throws IOException {

        GenerateVigor3Models generateVigor3Models = new GenerateVigor3Models();
        allVigor3Models = generateVigor3Models.generateModels(vigor3OutputTBL, vigor3OutputPep, inputFasta);
        return allVigor3Models;
    }

    @Test
    public void validate() throws IOException, VigorException {
        validate(getVigor4Models(), getVigor3Models());
    }

    public void validate (Map<String, List<Model>> allVigor4Models, Map<String, List<Model>> allVigor3Models) {

        File errorReport = new File(virusSpecifcDirectory, "differencesReport.txt");
        try (FileWriter writer = new FileWriter(errorReport, false)) {
            writer.write("***********************Differences Report**************************\n");
            for (String genome: allVigor3Models.keySet()) {
                List<Model> vigor4Models = allVigor4Models.getOrDefault(genome, Collections.EMPTY_LIST);
                if (vigor4Models.isEmpty()) {
                    writer.write(String.format("\nNo vigor4 models found for genome %s\n", genome));
                    continue;
                }
                List<Model> vigor3Models = allVigor3Models.get(genome);

                for (Model vigor3Model : vigor3Models) {
                    List<String> errors = new ArrayList<>();
                    String vigor3GeneID = vigor3Model.getGeneSymbol();
                    String vigor3GenomeID = vigor3Model.getAlignment().getVirusGenome().getId();
                    Optional<Model> vigor4Model = vigor4Models.stream().filter(m -> vigor3GeneID.equals(m.getGeneSymbol())).findFirst();
                    if (!vigor4Model.isPresent()) {
                        errors.add(String.format("Vigor3 & Vigor4 gene models do not match for VirusGenome Sequence %s.Expected gene symbol %s \n", vigor3GenomeID, vigor3GeneID));
                    } else {
                        errors.addAll(compareModels(vigor3Model, vigor4Model.get()));
                    }

                    if (!errors.isEmpty()) {
                        writer.write("\n***************************************************************************************************************************************************" + "\n");
                        writer.write(String.format(">Genome %s\n", genome));
                        for (String error : errors) {
                            writer.write(error);
                            writer.write("\n");
                        }
                        writer.write("\n***************************************************************************************************************************************************" + "\n");
                    }
                }
            }
        } catch (Exception e) {
            LOGGER.error(e);
            System.exit(1);
        }
    }
    
    List<String> compareModels(Model vigor3Model, Model vigor4Model) {
        List<String> errors = new ArrayList<>();

        List<Exon> vigor4Exons = vigor4Model.getExons();
        List<Exon> vigor3Exons = vigor3Model.getExons();
        String vigor3GeneID = vigor3Model.getGeneSymbol();
        String vigor3GenomeID = vigor3Model.getAlignment().getVirusGenome().getId();
        
        if (vigor4Model.isPartial3p() != vigor3Model.isPartial3p()) {
            errors.add(String.format("Partial 3' gene feature mismatch found for gene symbol %s of VirusGenome Sequence %s. partial 3' %s expected", vigor3GeneID, vigor3GenomeID, vigor3Model.isPartial3p()));
        }
        if (vigor3Model.isPartial5p() != vigor4Model.isPartial5p()) {
            errors.add(String.format("Partial 5' gene feature mismatch found for gene symbol %s of VirusGenome Sequence %s. partial 5' %s expected", vigor3GeneID, vigor3GenomeID, vigor3Model.isPartial5p()));
        }
        if (vigor3Model.getRibosomalSlippageRange() == null ^ vigor4Model.getRibosomalSlippageRange() == null ) {
            errors.add(String.format("Ribosomal Slippage mismatch found for gene %s of VirusGenome Sequence %s", vigor3GeneID, vigor3GenomeID));
        }
        if (vigor3Model.getReplaceStopCodonRange() == null ^ vigor4Model.getReplaceStopCodonRange() == null ) {
            errors.add(String.format("Stop Codon Readthrough feature mismatch for gene %s of VirusGenome Sequence %s", vigor3GeneID, vigor3GenomeID));
        }
        if (vigor3Model.isPseudogene() != vigor4Model.isPseudogene()) {
            errors.add(String.format("Pseudogene feature mismatch found for gene %s of VirusGenome Sequence %s. pseudogene %s expected", vigor3GeneID, vigor3GenomeID, vigor3Model.isPseudogene()));
        }
        if (vigor4Exons.size() != vigor3Exons.size()) {
            errors.add(String.format("Exon count differs for Vigor3 (%s) & Vigor4 (s) models with gene symbol %s of VirusGenome Sequence %s",
                                     vigor3Exons.size(), vigor4Exons.size(), vigor3GeneID, vigor3GenomeID));
        } else {
            for (int i = 0; i < vigor3Exons.size(); i++) {
                Range vigor3ExonRange = vigor3Exons.get(i).getRange(); 
                Range vigor4ExonRange = vigor4Exons.get(i).getRange();
                if (! vigor4ExonRange.equals(vigor3ExonRange)) {
                    errors.add(String.format("Exon range differs for Vigor3 (%s) & Vigor4 (%s) model with gene symbol %s of VirusGenome Sequence %s",
                                             vigor3ExonRange.toString(Range.CoordinateSystem.RESIDUE_BASED),
                                             vigor4ExonRange.toString(Range.CoordinateSystem.RESIDUE_BASED),
                                             vigor3GeneID, vigor3GenomeID));
                }
            }
        }
        return errors;
    }
}
