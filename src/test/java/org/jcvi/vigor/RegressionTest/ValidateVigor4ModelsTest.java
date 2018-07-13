package org.jcvi.vigor.RegressionTest;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.testing.category.Regression;
import org.jcvi.vigor.testing.category.Slow;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.VigorInitializationService;
import org.jcvi.vigor.utils.*;
import org.junit.ClassRule;
import org.junit.Rule;
import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.rules.SpringClassRule;
import org.springframework.test.context.junit4.rules.SpringMethodRule;

import static org.junit.Assert.fail;

@Category({Slow.class, Regression.class})
@RunWith(Parameterized.class)
@ContextConfiguration(classes = Application.class)
public class ValidateVigor4ModelsTest {
    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4ModelsTest.class);

    private String vigor3OutputTBL;
    private String vigor3OutputPep;
    private String inputFasta;
    private String referenceDatabaseName;

    @Autowired
    GenerateVigor4GeneModels generateVigor4GeneModels;

    @Autowired
    private VigorInitializationService initializationService;

    @ClassRule
    public static final SpringClassRule springClassRule = new SpringClassRule();
    @Rule
    public final SpringMethodRule springMethodRule = new SpringMethodRule();

    public ValidateVigor4ModelsTest(String vigor3OutputTBL, String vigor3OutputPep, String inputFasta, String referenceDatabaseName ) {

        this.vigor3OutputTBL = getResource(vigor3OutputTBL).orElseThrow(
                () -> new IllegalArgumentException(String.format("%s not found", vigor3OutputTBL)));
        this.vigor3OutputPep = getResource(vigor3OutputPep).orElseThrow(
                () -> new IllegalArgumentException(String.format("%s not found", vigor3OutputPep)));
        this.inputFasta = getResource(inputFasta).orElseThrow(
                () -> new IllegalArgumentException(String.format("%s not found", inputFasta)));
        this.referenceDatabaseName = referenceDatabaseName;
    }

    private Optional<String> getResource (String resource) {

        URL url = ValidateVigor4ModelsTest.class.getClassLoader().getResource(resource);
        if ( url == null) {
            return Optional.empty();
        }
        return Optional.of(Paths.get(url.getFile()).toAbsolutePath().toString());
    }

    @Parameterized.Parameters(name = "ValidateVigor4ModelsTest[#{index} {0}]")
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

    public Map<String, List<Model>> getVigor4Models (VigorConfiguration config, String inputFasta) throws VigorException, IOException {
        checkConfig(config);
        String referenceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        String referenceDB = Paths.get(referenceDBPath, referenceDatabaseName).toString();
        VigorUtils.checkFilePath("reference database", referenceDB,
                                 VigorUtils.FileCheck.EXISTS,
                                 VigorUtils.FileCheck.READ,
                                 VigorUtils.FileCheck.FILE);
        return generateVigor4GeneModels.generateModels(inputFasta, referenceDB, config);
    }

    public Map<String, List<Model>> getVigor3Models () throws IOException {
        return new GenerateVigor3Models().generateModels(vigor3OutputTBL, vigor3OutputPep, inputFasta);
    }

    @Test
    public void validate() throws IOException, VigorException {
        VigorConfiguration config = getConfiguration();
        Map<String, List<String>> errors = validate(getVigor4Models(config, inputFasta), getVigor3Models());
        boolean hasErrors = errors.entrySet()
                                  .stream()
                                  .filter( e -> ! e.getValue().isEmpty())
                                  .findAny().isPresent();
        if (hasErrors) {
            StringBuilder sb = new StringBuilder();
            for (String genome: errors.keySet()) {
                if (errors.get(genome).isEmpty()) {
                    continue;
                }
                sb.append("\n*************\n");
                sb.append("For Genome: ");
                sb.append(genome);
                sb.append("\n\n");
                sb.append(String.join("\n", errors.get(genome)));
                sb.append("\n\n*************\n\n");
            }

            if ("true".equals(System.getProperty("vigor.regression_test.write_report")) ||
                    "true".equals(System.getenv("VIGOR_REGRESSION_TEST_WRITE_REPORT"))) {
                Path reportFile  = Paths.get(config.get(ConfigurationParameters.OutputDirectory), "differencesReport.txt");
                boolean overwrite = "true".equals(config.get(ConfigurationParameters.OverwriteOutputFiles));

                List<OpenOption> openOptionsList = new ArrayList<>();
                if (overwrite) {
                    openOptionsList.add(StandardOpenOption.CREATE);
                    openOptionsList.add(StandardOpenOption.TRUNCATE_EXISTING);
                } else {
                    openOptionsList.add(StandardOpenOption.CREATE_NEW);
                }
                try (BufferedWriter writer = Files.newBufferedWriter(reportFile,
                                                                     Charset.forName("UTF-8"),
                                                                     openOptionsList.toArray(new OpenOption[] {}) )) {
                    writer.write(sb.toString());
                    writer.flush();
                }
            }

            fail(sb.toString());
        }
    }

    private VigorConfiguration getConfiguration() throws IOException, VigorException {
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String outputPrefix = new File(referenceDatabaseName).getName().replace("_db", "");
        config.put(ConfigurationParameters.OutputPrefix, outputPrefix);
        String outDir = config.get(ConfigurationParameters.OutputDirectory);
        String tmpDir = config.get(ConfigurationParameters.TemporaryDirectory);
        if (outDir == null) {
            Path outDirPath;
            if ( tmpDir != null) {
                outDirPath = Files.createTempDirectory(Paths.get(tmpDir), "vigor4_test");
            } else {
                outDirPath = Files.createTempDirectory("vigor4_test");
            }
            outDir = outDirPath.toString();
            Runtime.getRuntime().addShutdownHook(new Thread(() -> VigorUtils.deleteDirectory(outDirPath)));
        }
        String virusSpecificPath = Paths.get(outDir,config.get(ConfigurationParameters.OutputPrefix)).toString();
        File virusSpecificDir = new File(virusSpecificPath);
        if (! (virusSpecificDir.exists() || virusSpecificDir.mkdir()))  {
            throw new VigorException(String.format("virus specific output directory %s doesn't exist and could not be created", virusSpecificPath));
        }
        config.put(ConfigurationParameters.OutputDirectory, virusSpecificDir.getAbsolutePath());
        config.put(ConfigurationParameters.Verbose, "false");
        if (config.get(ConfigurationParameters.OverwriteOutputFiles) == null) {
            config.put(ConfigurationParameters.OverwriteOutputFiles, "true");
        }
        // TODO allow user to set this
        if (tmpDir == null) {
            final Path tempDir = Files.createTempDirectory(Paths.get(outDir), "vigor-tmp");
            // delete on shutdown
            Runtime.getRuntime().addShutdownHook(new Thread(() -> VigorUtils.deleteDirectory(tempDir)));
            config.put(ConfigurationParameters.TemporaryDirectory, tempDir.toString());
        }

        checkConfig(config);
        return config;
    }

    private void checkConfig(VigorConfiguration config) throws VigorException {
        VigorUtils.checkFilePath("reference database path", config.get(ConfigurationParameters.ReferenceDatabasePath),
                                 VigorUtils.FileCheck.EXISTS, VigorUtils.FileCheck.READ);
        VigorUtils.checkFilePath("temporary directory path", config.get(ConfigurationParameters.TemporaryDirectory),
                                 VigorUtils.FileCheck.EXISTS, VigorUtils.FileCheck.READ, VigorUtils.FileCheck.WRITE);
        // must be set, (but will be created if it doesn't exist?)
        VigorUtils.checkFilePath("output directory", config.get(ConfigurationParameters.OutputDirectory));
    }

    public Map<String, List<String>> validate (Map<String, List<Model>> allVigor4Models, Map<String, List<Model>> allVigor3Models) {

        Map<String, List<String>> allErrors = new HashMap<>();

        for (String genome: allVigor3Models.keySet()) {
            List<Model> vigor4Models = allVigor4Models.getOrDefault(genome, Collections.EMPTY_LIST);
            List<String> errors = allErrors.computeIfAbsent(genome, k-> new ArrayList<>());
            if (vigor4Models.isEmpty()) {
                errors.add(String.format("No vigor4 models found for genome %s", genome));
                continue;
            }
            List<Model> vigor3Models = allVigor3Models.get(genome);

            for (Model vigor3Model : vigor3Models) {
                String vigor3GeneID = vigor3Model.getGeneSymbol();
                String vigor3GenomeID = vigor3Model.getAlignment().getVirusGenome().getId();
                Optional<Model> vigor4Model = vigor4Models.stream().filter(m -> vigor3GeneID.equals(m.getGeneSymbol())).findFirst();
                if (!vigor4Model.isPresent()) {
                    errors.add(String.format("Vigor3 & Vigor4 gene models do not match for VirusGenome Sequence %s. Expected gene symbol %s", vigor3GenomeID, vigor3GeneID));
                } else {
                    errors.addAll(compareModels(vigor3Model, vigor4Model.get()));
                }

            }
            LOGGER.debug("{} errors for genome {}", errors.size(), genome);
        }
        return allErrors;
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
            errors.add(String.format("Exon count differs for Vigor3 (%s) & Vigor4 (%s) models with gene symbol %s of VirusGenome Sequence %s",
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
