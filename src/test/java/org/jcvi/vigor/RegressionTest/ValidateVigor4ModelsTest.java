package org.jcvi.vigor.RegressionTest;

import java.io.*;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;
import com.opencsv.CSVReader;
import com.opencsv.CSVReaderBuilder;
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

@Category( { Slow.class, Regression.class })
@RunWith(com.googlecode.junittoolbox.ParallelParameterized.class)
@ContextConfiguration(classes = Application.class)
public class ValidateVigor4ModelsTest {

    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4ModelsTest.class);
    private String referenceOutputTBL;
    private String inputFasta;
    private String referenceDatabaseName;
    private String referenceType;
    @Autowired
    GenerateVigor4GeneModels generateVigor4GeneModels;
    @Autowired
    private VigorInitializationService initializationService;
    @ClassRule
    public static final SpringClassRule springClassRule = new SpringClassRule();
    @Rule
    public final SpringMethodRule springMethodRule = new SpringMethodRule();

    public ValidateVigor4ModelsTest ( String referenceOutputTBL, String inputFasta, String referenceDatabaseName, String referenceType ) throws Exception{

        this.referenceOutputTBL = getResource(referenceOutputTBL).orElseThrow(
                () -> new IllegalArgumentException(String.format("%s not found", referenceOutputTBL)));
        this.inputFasta = getResource(inputFasta).orElseThrow(
                () -> new IllegalArgumentException(String.format("%s not found", inputFasta)));
        this.referenceDatabaseName = referenceDatabaseName;
        this.referenceType = referenceType;
    }

    private Optional<String> getResource ( String resource ) {

        URL url = ValidateVigor4ModelsTest.class.getClassLoader().getResource(resource);
        if (url == null) {
            return Optional.empty();
        }
        return Optional.of(Paths.get(url.getFile()).toAbsolutePath().toString());
    }

    @Parameterized.Parameters(name = "ValidateVigor4ModelsTest[#{index} {0}]")
    public static Collection<Object[]> getTestData () throws IOException {

        List<Object[]> testData = new ArrayList<>();
        boolean hasHeader = true;
        Reader dataReader;
        String configData = System.getProperty("vigor.regression_test.config_csv");
        if (! NullUtil.isNullOrEmpty(configData)) {
            dataReader = new StringReader(configData);
            hasHeader = false;
        } else {
            String configFile = System.getProperty("vigor.regression_test.config_file",
                                                   ValidateVigor4ModelsTest.class.getClassLoader().getResource("config/RegressionTestConfig.csv").getFile());
            dataReader = new FileReader(new File(configFile));
        }

        final CSVReader reader = new CSVReaderBuilder(dataReader).withSkipLines(hasHeader ? 1 : 0).build();
        String[] nextLine;
        // format is: expected TBL file, input fasta, reference database, vigor version
        while((nextLine= reader.readNext()) !=null){
             testData.add(Arrays.copyOfRange(nextLine,0,4));
         }
        LOGGER.info("Will run DBs {}", testData.stream().map(t -> (String) t[2]).collect(Collectors.joining(",")));
        return testData;
    }

    public Map<String, List<Model>> getVigor4Models ( VigorConfiguration config, String inputFasta ) throws VigorException, IOException {

        checkConfig(config);
        String referenceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        String referenceDB = Paths.get(referenceDBPath, referenceDatabaseName).toString();
        VigorUtils.checkFilePath("reference database", referenceDB,
                VigorUtils.FileCheck.EXISTS,
                VigorUtils.FileCheck.READ,
                VigorUtils.FileCheck.FILE);
        return generateVigor4GeneModels.generateModels(inputFasta, referenceDB, config);
    }

    public Map<String, List<Model>> getReferenceModels () throws IOException {

        return new GenerateReferenceModels().generateModels(referenceOutputTBL, inputFasta);
    }

    @Test
    public void validate () throws IOException, VigorException {

        LOGGER.info("Validating DB {} using input {} and TBL file {}", this.referenceDatabaseName, this.inputFasta, this.referenceOutputTBL);
        VigorConfiguration config = getConfiguration();
        Map<String, List<String>> errors = compareWithReferenceModels(getVigor4Models(config, inputFasta), getReferenceModels());
        String errorReport = String.format(config.get(ConfigurationParameters.OutputPrefix)+"_differencesReport_%sRef.txt", referenceType);
        boolean hasErrors = errors.entrySet()
                .stream()
                .filter(e -> !( e.getValue() == null || e.getValue().isEmpty() ))
                .findAny().isPresent();
        if (hasErrors) {
            StringBuilder sb = new StringBuilder();
            for (String genome : errors.keySet()) {
                if (errors.get(genome).isEmpty()) {
                    continue;
                }
                sb.append("\n**********************************************************************\n");
                sb.append("For Genome: ");
                sb.append(genome);
                sb.append("\n\n");
                sb.append(String.join("\n", errors.get(genome)));
                sb.append("\n\n*********************************************************************\n\n");
            }
            if ("true".equals(System.getProperty("vigor.regression_test.write_report")) ||
                    "true".equals(System.getenv("VIGOR_REGRESSION_TEST_WRITE_REPORT"))) {
                Path reportFile = Paths.get(config.get(ConfigurationParameters.OutputDirectory), errorReport);
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
                        openOptionsList.toArray(new OpenOption[] {}))) {
                    writer.write(sb.toString());
                    writer.flush();
                }
            }
            fail(sb.toString());
        }
    }

    private VigorConfiguration getConfiguration () throws IOException, VigorException {

        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());

        String outputPrefix = new File(referenceDatabaseName).getName().replace("_db", "");
        config.putString(ConfigurationParameters.OutputPrefix, outputPrefix);
        String outDir = config.get(ConfigurationParameters.OutputDirectory);
        String tmpDir = config.get(ConfigurationParameters.TemporaryDirectory);
        if (outDir == null) {
            Path outDirPath;
            if (tmpDir != null) {
                outDirPath = Files.createTempDirectory(Paths.get(tmpDir), "vigor4_test");
            } else {
                outDirPath = Files.createTempDirectory("vigor4_test");
            }
            outDir = outDirPath.toString();
            Runtime.getRuntime().addShutdownHook(new Thread(() -> VigorUtils.deleteDirectory(outDirPath)));
        }
        String virusSpecificPath = Paths.get(outDir, outputPrefix).toString();
        File virusSpecificDir = new File(virusSpecificPath);
        if (!( virusSpecificDir.exists() || virusSpecificDir.mkdir() )) {
            throw new VigorException(String.format("virus specific output directory %s doesn't exist and could not be created", virusSpecificPath));
        }
        config.putString(ConfigurationParameters.OutputDirectory, virusSpecificDir.getAbsolutePath());
        config.putString(ConfigurationParameters.Verbose, "false");
        if (config.get(ConfigurationParameters.OverwriteOutputFiles) == null) {
            config.putString(ConfigurationParameters.OverwriteOutputFiles, "true");
        }
        // TODO allow user to set this
        if (tmpDir == null) {
            final Path tempDir = Files.createTempDirectory(Paths.get(outDir), "vigor-tmp");
            // delete on shutdown
            Runtime.getRuntime().addShutdownHook(new Thread(() -> VigorUtils.deleteDirectory(tempDir)));
            config.putString(ConfigurationParameters.TemporaryDirectory, tempDir.toString());
        }
        checkConfig(config);

        String referenceDatabasePath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        File virusSpecificConfig = new File(referenceDatabasePath, referenceDatabaseName + ".ini").getAbsoluteFile();
        if (virusSpecificConfig.exists()) {
            VigorConfiguration virusConfig = LoadDefaultParameters.loadVigorConfiguration(virusSpecificConfig.toString(),
                                                                                          virusSpecificConfig,
                                                                                          section -> EnumSet.of(ConfigurationParameters.Flags.GENE_SET,
                                                                                                                ConfigurationParameters.Flags.VIRUS_SET));
            virusConfig.setDefaults(config);
            config = virusConfig;
            checkConfig(config);
        }
        return config;
    }

    private void checkConfig ( VigorConfiguration config ) throws VigorException {

        VigorUtils.checkFilePath("reference database path", config.get(ConfigurationParameters.ReferenceDatabasePath),
                VigorUtils.FileCheck.EXISTS, VigorUtils.FileCheck.READ);
        VigorUtils.checkFilePath("temporary directory path", config.get(ConfigurationParameters.TemporaryDirectory),
                VigorUtils.FileCheck.EXISTS, VigorUtils.FileCheck.READ, VigorUtils.FileCheck.WRITE);
        // must be set, but will be created if it doesn't exist.
        VigorUtils.checkFilePath("output directory", config.get(ConfigurationParameters.OutputDirectory), VigorUtils.FileCheck.SET);
    }

    public Map<String, List<String>> compareWithReferenceModels ( Map<String, List<Model>> allVigor4Models, Map<String, List<Model>> allReferenceModels ) {

        Map<String, List<String>> allErrors = new HashMap<>();
        for (String genome : allReferenceModels.keySet()) {
            List<Model> vigor4Models = allVigor4Models.getOrDefault(genome, Collections.EMPTY_LIST);
            List<String> errors = allErrors.computeIfAbsent(genome, k -> new ArrayList<>());
            if (vigor4Models.isEmpty()) {
                errors.add(String.format("No vigor4 models found for genome %s", genome));
                continue;
            }
            List<Model> refModels = allReferenceModels.get(genome);
            List<String> actualGeneSymbs = vigor4Models.stream().map(Model::getGeneSymbol).collect(Collectors.toCollection(ArrayList::new));
            List<String> refGeneSymbs = refModels.stream().map(Model::getGeneSymbol).collect(Collectors.toCollection(ArrayList::new));
            Set<String> vigor4Difference = actualGeneSymbs.stream().filter(t->!refGeneSymbs.contains(t)).collect(Collectors.toSet());
            boolean reportAllModels=false;
            if(vigor4Difference.size()>0) {
                errors.add(String.format("Vigor4 reported additional/different geneModels for VirusGenome Sequence %s ", genome));
                vigor4Difference.stream().forEach(g->{
                    Model diffModel = vigor4Models.stream().filter(m -> g.equals(m.getGeneSymbol())).findFirst().get();
                    errors.add("\n"+diffModel+"\n");
                });
                reportAllModels=true;

            }
            for (Model refModel : refModels) {
                boolean errorFound = false;
                String refGeneID = refModel.getGeneSymbol();
                String refGenomeID = refModel.getAlignment().getVirusGenome().getId();

                Model vigor4Model=null;
                for(Model tempVigo4Model : vigor4Models){
                    if(tempVigo4Model.getGeneSymbol().equals(refModel.getGeneSymbol())) {
                        vigor4Model = tempVigo4Model;
                        break;
                    }
                }
                if (vigor4Model==null) {
                    errors.add(String.format(referenceType + " reference models & Vigor4 models do not match for VirusGenome Sequence %s. Expected gene symbol %s", refGenomeID, refGeneID));
                    reportAllModels=true;
                } else {
                    List<String> outErrors = compareModels(refModel, vigor4Model);
                    if (outErrors.size() > 0) {
                        errorFound = true;
                    }
                    errors.addAll(outErrors);
                }
                if (errorFound) {
                    errors.add(String.format("\nReferenceModel :%s \n\nVigor4Model :%s", refModel.toString(), vigor4Model.toString()));
                }
            }
            if(reportAllModels){
                errors.add(String.format("\nReferenceModels:\n %s \n\nVigor4Models:\n %s", refModels.stream().map(Object::toString)
                        .collect(Collectors.joining("\n")),vigor4Models.stream().map(Object::toString)
                        .collect(Collectors.joining("\n"))));
            }
            LOGGER.debug("{} errors for genome {}", errors.size(), genome);
        }
        return allErrors;
    }

    List<String> compareModels ( Model refModel, Model model ) {

        List<String> errors = new ArrayList<>();
        List<Exon> foundExons = model.getExons();
        List<Exon> refExons = refModel.getExons();
        String refGenomeID = refModel.getAlignment().getVirusGenome().getId();
        String expectedGeneSymbol = refModel.getGeneSymbol();
        if (model.isPartial3p() != refModel.isPartial3p()) {
            errors.add(String.format("\nPartial 3' gene feature mismatch found for gene %s of VirusGenome Sequence %s", expectedGeneSymbol, refGenomeID));
        }
        if (refModel.isPartial5p() != model.isPartial5p()) {
            errors.add(String.format("\nPartial 5' gene feature mismatch found for gene %s of VirusGenome Sequence %s", expectedGeneSymbol, refGenomeID));
        }
        if (refModel.getRibosomalSlippageRange() == null ^ model.getRibosomalSlippageRange() == null) {
            errors.add(String.format("\nRibosomal Slippage mismatch found for gene %s of VirusGenome Sequence %s", expectedGeneSymbol, refGenomeID));
        }
        if (refModel.getReplaceStopCodonRange() == null ^ model.getReplaceStopCodonRange() == null) {
            errors.add(String.format("\nStop Codon Readthrough feature mismatch found for gene %s of VirusGenome Sequence %s", expectedGeneSymbol, refGenomeID));
        }
        if (refModel.isPseudogene() != model.isPseudogene()) {
            errors.add(String.format("\nPseudogene feature mismatch found for gene %s of VirusGenome Sequence %s", expectedGeneSymbol, refGenomeID));
        }
        if (foundExons.size() != refExons.size()) {
            errors.add(String.format("\nExon count differs for models with gene symbol %s of VirusGenome Sequence %s", expectedGeneSymbol, refGenomeID));
        } else {
            for (int i = 0; i < refExons.size(); i++) {
                Range vigor3ExonRange = refExons.get(i).getRange();
                Range vigor4ExonRange = foundExons.get(i).getRange();
                if (!vigor4ExonRange.equals(vigor3ExonRange)) {
                    errors.add(String.format("\nExon range differs for models with gene symbol %s of VirusGenome Sequence %s. Expected %s , found %s ",
                            expectedGeneSymbol, refGenomeID, vigor3ExonRange,
                            vigor4ExonRange
                    ));
                }
            }
        }
        return errors;
    }
}
