package org.jcvi.vigor;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.Filter;
import org.apache.logging.log4j.core.Layout;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.filter.LevelRangeFilter;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.*;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.service.exception.UserFacingException;
import org.jcvi.vigor.utils.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;

@Component
public class Vigor {

    private static Logger LOGGER = LogManager.getLogger(Vigor.class);
    @Autowired
    private VigorInitializationService initializationService;
    @Autowired
    private AlignmentGenerationService alignmentGenerationService;
    @Autowired
    ModelGenerationService modelGenerationService;
    @Autowired
    private VigorInputValidationService inputValidationService;
    @Autowired
    private GeneModelGenerationService geneModelGenerationService;
    @Autowired
    private PeptideMatchingService peptideMatchingService;


    public void run ( String... args ) {

        Namespace parsedArgs = parseArgs(args);
        int verbosity = parsedArgs.getInt(CommandLineParameters.verbose);
        setVerboseLogging(verbosity);
        try {
            VigorConfiguration vigorConfiguration = getVigorConfiguration(parsedArgs);
            if (parsedArgs.getBoolean(CommandLineParameters.listDatabases)) {
                String referenceDatabasePath = vigorConfiguration.get(ConfigurationParameters.ReferenceDatabasePath);
                VigorUtils.checkFilePath("reference database path", referenceDatabasePath,
                                         VigorUtils.FileCheck.EXISTS,
                                         VigorUtils.FileCheck.DIRECTORY,
                                         VigorUtils.FileCheck.READ);
                List<VigorInitializationService.DatabaseInfo> databases = initializationService.getDatabaseInfo(referenceDatabasePath);
                printDatabaseInfo(referenceDatabasePath, databases);
                System.exit(0);
            }
            checkConfig(vigorConfiguration);
            File outputDirectory = new File((String) vigorConfiguration.get(ConfigurationParameters.OutputDirectory));

            String outputPrefix = vigorConfiguration.get(ConfigurationParameters.OutputPrefix);
            initiateReportFile(outputDirectory.getAbsolutePath(), outputPrefix, verbosity);

            String referenceDB = vigorConfiguration.get(ConfigurationParameters.ReferenceDatabaseFile);
            LOGGER.info("Command line arguments: {}", String.join(" ", args));
            LOGGER.info("Current working directory: {}", Paths.get("").toAbsolutePath().normalize().toString());
            String inputFileName = parsedArgs.getString("input_fasta");
            generateAnnotations(inputFileName, referenceDB, vigorConfiguration);
        } catch (UserFacingException e) {
            System.err.println(e.getMessage());
            System.exit(1);
        } catch (Exception e) {
            LOGGER.error("Exception encountered: exiting 1",e);
            System.exit(1);
        }
    }

    private void printDatabaseInfo(String referenceDatabasePath, List<VigorInitializationService.DatabaseInfo> databases) {
        LOGGER.info("Databases found under {}", referenceDatabasePath);
        List<VigorInitializationService.DatabaseInfo> sortedDatabases = databases.stream()
                                                                                 .sorted(Comparator.comparing(d -> d.databaseFile.getName(), String.CASE_INSENSITIVE_ORDER))
                                                                                 .collect(Collectors.toList());
        for (VigorInitializationService.DatabaseInfo db: sortedDatabases) {
            List<String> dbInfo = new ArrayList<>(8);
            dbInfo.add(String.format("\nDatabase file: %s", db.databaseFile.getName()));
            dbInfo.add(String.format("  Config file: %s", db.configFile.orElse(new File("no config found")).getName()));
            if (db.configFile.isPresent()) {
                try {
                    LOGGER.trace("loading file {}", db.configFile.get().getAbsolutePath());
                    VigorConfiguration vigorConfiguration = initializationService.mergeConfigurations(initializationService.loadVirusConfiguration(db.configFile.get().getAbsoluteFile()));
                    if (vigorConfiguration.hasSection(VigorConfiguration.METADATA_SECTION)) {
                        String virusName = vigorConfiguration.get(VigorConfiguration.METADATA_SECTION, ConfigurationParameters.VirusName);
                        if (!NullUtil.isNullOrEmpty(virusName)) {
                            dbInfo.add(String.format("  Virus: %s", virusName));
                        }
                        List<String> aliases = vigorConfiguration.get(VigorConfiguration.METADATA_SECTION, ConfigurationParameters.Alias);
                        if (! (aliases == null || aliases.isEmpty()) ) {
                            dbInfo.add(String.format("  Alias(es): %s", String.join(", ", aliases)));
                        }
                        String description = vigorConfiguration.get(VigorConfiguration.METADATA_SECTION, ConfigurationParameters.Description);
                        if (!NullUtil.isNullOrEmpty(description)) {
                            dbInfo.add(String.format("  Description: %s", description));
                        }
                        List<String> notes = vigorConfiguration.get(VigorConfiguration.METADATA_SECTION, ConfigurationParameters.Note);
                        if (!(notes == null || notes.isEmpty())) {
                            dbInfo.add(String.format("  Note(s): %s", String.join(", ", notes)));
                        }
                    }
                } catch (VigorException e) {
                    LOGGER.error(String.format("Unable to load config file %s", db.configFile.get().getName()), e);
                }
            }
            LOGGER.info(String.join("\n", dbInfo));
        }
    }



    /**
     * Check that required parameters are set
     * @param config
     * @throws VigorException
     */
    private void checkConfig(VigorConfiguration config) throws VigorException {
        List<String> errors =  new ArrayList<>();

        String outputDirectoryPath = config.get(ConfigurationParameters.OutputDirectory);
        String outputPrefix = config.get(ConfigurationParameters.OutputPrefix);
        if (! NullUtil.isNullOrEmpty(outputDirectoryPath)) {
            File outputDirectory = new File(outputDirectoryPath);
            if (!(outputDirectory.exists() || outputDirectory.mkdirs())) {
                errors.add(String.format("unable to create directory %s", outputDirectory));
            }
            if (!(outputDirectory.exists() && outputDirectory.isDirectory())) {
                errors.add(String.format("Invalid output prefix %s/%s. Please provide a directory followed by a file prefix", outputDirectory, outputPrefix));
            }
        } else {
            errors.add("output directory is not set");
        }

        String reference_db = config.get(ConfigurationParameters.ReferenceDatabaseFile);

        try {
            VigorUtils.checkFilePath("Reference database file", reference_db,
                                     VigorUtils.FileCheck.EXISTS,
                                     VigorUtils.FileCheck.FILE,
                                     VigorUtils.FileCheck.READ);
        } catch (VigorException e) {
            errors.add(e.getMessage());
        }

        String temporaryDirectory = config.get(ConfigurationParameters.TemporaryDirectory);
        if ( NullUtil.isNullOrEmpty(temporaryDirectory)) {
            errors.add("temporary directory not set");
        }

        File tempDir = Paths.get(temporaryDirectory).toFile();
        if (tempDir.exists()) {
            try {
                VigorUtils.checkFilePath("temporary directory", temporaryDirectory,
                                         VigorUtils.FileCheck.DIRECTORY,
                                         VigorUtils.FileCheck.READ,
                                         VigorUtils.FileCheck.WRITE);
            } catch (VigorException e) {
                errors.add(e.getMessage());
            }
        } else {
            if (! tempDir.mkdirs()) {
                errors.add(String.format("unable to create temporary directory %s", tempDir));
            }
        }

        if (! errors.isEmpty()) {
            throw new UserFacingException(String.join("\n", errors));
        }

    }

    private void printConfiguration(VigorConfiguration vigorParameters){
        VigorConfiguration.ValueWithSource unset = VigorConfiguration.ValueWithSource.of("NOTSET", "unknown");
        LOGGER.info(() -> vigorParameters.entrySet()
                                         .stream()
                                         .sorted(Comparator.comparing(es -> es.getKey().configKey, String.CASE_INSENSITIVE_ORDER))
                                         .map(e -> String.format("%-50s%s (%s)",
                                                                 e.getKey().configKey,
                                                                 e.getValue(),
                                                                 vigorParameters.getWithSource(e.getKey()).orElse(unset).source))
                                         .collect(Collectors.joining("\n")));
    }

    public void generateAnnotations(String inputFileName, String referenceDB, VigorConfiguration vigorParameters) throws VigorException {
        try {
            VigorUtils.checkFilePath("input file", inputFileName, VigorUtils.FileCheck.EXISTS, VigorUtils.FileCheck.READ);
        } catch (VigorException e) {
            throw new UserFacingException(e.getMessage());
        }
        printConfiguration(vigorParameters);
        String outputDir = vigorParameters.get(ConfigurationParameters.OutputDirectory);
        String outputPrefix = vigorParameters.get(ConfigurationParameters.OutputPrefix);
        VigorUtils.checkFilePath("output directory", outputDir,
                                 VigorUtils.FileCheck.EXISTS,
                                 VigorUtils.FileCheck.WRITE,
                                 VigorUtils.FileCheck.DIRECTORY);
        List<IOutputWriter> writers = getWriters(vigorParameters);
        try (NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(new File(inputFileName))
                .hint(DataStoreProviderHint.ITERATION_ONLY)
                .build();
             Outfiles outfiles = getOutfiles(vigorParameters);
        ) {
            // TODO move all this file handling to method
            // TODO checkout output earlier.
            writeEffectiveConfig(outputDir, outputPrefix, vigorParameters);
            // initialize the writers. This will fail if the files exist and we're not overwriting
            for (IOutputWriter writer: writers) {
                writer.getWriter(outfiles, new OutputContext());
            }
            Iterator<NucleotideFastaRecord> recordIterator = dataStore.records().iterator();
            while (recordIterator.hasNext()) {
                NucleotideFastaRecord record = recordIterator.next();
                LOGGER.debug("processing {}", record.getId());
                List<Model> geneModels = modelsFromNucleotideRecord(record, referenceDB, vigorParameters);
                if (geneModels.isEmpty()) {
                    LOGGER.warn("No gene models generated for sequence {}", record.getId());
                    continue;
                }
                outputModels(writers, outfiles, geneModels);
            }
        } catch (DataStoreException e) {
            throw new VigorException(String.format("problem reading input file %s", inputFileName), e);
        } catch (FileAlreadyExistsException e) {
            throw new UserFacingException(String.format("File already exists %s", e.getMessage()));
        } catch (IOException e) {
            throw new VigorException(String.format("File issue. Got %s: %s", e.getClass().getSimpleName(), e.getMessage()), e);
        }
    }

    private List<IOutputWriter> getWriters(VigorConfiguration config) throws VigorException {
        // TODO get writer preferences from config
        List<IOutputWriter> writers = new ArrayList<>();
        Set<String> selectedWriters = config.getOrDefault(ConfigurationParameters.OutputFormats, Collections.EMPTY_SET);
        for (String selectedWriter: selectedWriters) {
            Supplier<IOutputWriter> writerSupplier = OutputWriters.Writers.get(selectedWriter);
            if (writerSupplier == null) {
                throw new VigorException("Unknown format " + selectedWriter);
            }
            writers.add(writerSupplier.get());
        }

        if (writers.isEmpty()) {
            throw new UserFacingException("No writers configured");
        }

        for (IOutputWriter writer: writers) {
            if (IConfigurable.class.isAssignableFrom(writer.getClass())) {
                ((IConfigurable) writer).configure(config);
            }
        }
        return writers;
    }

    public void outputModels(List<IOutputWriter> writers, Outfiles outfiles, List<Model> geneModels) throws IOException, VigorException {
        for (IOutputWriter writer: writers) {
            writer.writeModels(outfiles, geneModels);
            outfiles.flush();
        }
        FormatVigorOutput.printSequenceFeatures(geneModels, "GeneModels");
        outfiles.flush();
    }

    public List<Model> modelsFromNucleotideRecord(NucleotideFastaRecord record, String referenceDB, VigorConfiguration vigorParameters) throws VigorException {
        LOGGER.info("Getting alignments for {}", record.getId());
        VirusGenome virusGenome = VirusGenomeService.fastaRecordToVirusGenome(record, vigorParameters);
        List<Alignment> alignments = generateAlignments(virusGenome, referenceDB, vigorParameters);
        alignments = handleReverseAlignments(vigorParameters, alignments);
        LOGGER.info("{} alignment(s) found for sequence {}", alignments.size(), record.getId());
        List<Model> candidateModels = generateModels(alignments, vigorParameters);
        LOGGER.info("{} candidate model(s) found for sequence {}", candidateModels.size(), record.getId());
        List<Model> geneModels = generateGeneModels(candidateModels, vigorParameters);
        LOGGER.info("{} gene model(s) found for sequence {}", geneModels.size(), record.getId());
        geneModels = findPeptides(vigorParameters, geneModels);
        LOGGER.debug("Found {} peptides for {} models for sequence {}",
                     geneModels.stream()
                               .map( m -> m.getMaturePeptides().size() )
                               .reduce(Integer::sum).orElse(0),
                     geneModels.size(),
                     record.getId());

        // sort by begin,end
        return geneModels.stream()
                         .sorted(Comparator.comparing(m -> VigorFunctionalUtils.getDirectionBasedRange(m.getRange(),
                                                                                                       m.getAlignment().getVirusGenome().getSequence().getLength(),
                                                                                                       m.getDirection()),
                                                      Range.Comparators.ARRIVAL))
                         .collect(Collectors.toList());

    }

    /**
     * Handle alignments on the opposite strand by reverse complementing the sequence and recalculating stops and gaps
     * @param vigorParameters
     * @param alignments
     * @return
     */
    private List<Alignment> handleReverseAlignments(VigorConfiguration vigorParameters, List<Alignment> alignments) {
        for (Alignment alignment: alignments) {
            if (alignment.getDirection() == Direction.REVERSE) {
                VirusGenome complement = VirusGenomeService.reverseComplementVirusGenome(alignment.getVirusGenome(),vigorParameters);
                alignment.setVirusGenome(complement);
            }
        }
        return alignments;
    }

    private void setVerboseLogging ( int verbosity ) {

        if (verbosity == 0) {
            return;
        }
        Level verboseLevel = verbosity == 1 ? Level.DEBUG : Level.TRACE;
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        lc.getConfiguration().getLoggerConfig("org.jcvi.vigor").setLevel(verboseLevel);
        lc.updateLoggers();
        LOGGER.debug("verbose logging enabled at level {}", verboseLevel);
    }

    private PeptideMatchingService.Scores getPeptideScores ( VigorConfiguration config ) {

        double minIdentity = config.get(ConfigurationParameters.MaturePeptideMinimumIdentity);
        minIdentity = minIdentity / 100.0d;
        double minCoverage = config.get(ConfigurationParameters.MaturePeptideMinimumCoverage);
        minCoverage = minCoverage / 100.0d;
        double minSimilarity = config.get(ConfigurationParameters.MaturePeptideMinimumSimilarity);
        minSimilarity = minSimilarity / 100.0d;
        return PeptideMatchingService.Scores.of(minIdentity, minCoverage, minSimilarity);
    }

    private List<Model> findPeptides ( VigorConfiguration config, List<Model> geneModels ) throws VigorException {

        PeptideMatchingService.Scores scores = getPeptideScores(config);
        for (Model model : geneModels) {
            String maturePeptideDB = model.getAlignment().getAlignmentEvidence().getMatpep_db();
            // TODO check peptides for psuedogenes?
            if (!( maturePeptideDB == null || maturePeptideDB.isEmpty() )) {
                LOGGER.debug("finding mature peptides for {} using db {}", model.getGeneID(), maturePeptideDB);
                model.setMaturePeptides(peptideMatchingService.findPeptides(model,
                        new File(maturePeptideDB),
                        scores));
                LOGGER.debug("for {} found {} peptides.", model.getGeneID(), model.getMaturePeptides().size());
            }
        }
        return geneModels;
    }

    public Namespace parseArgs ( String[] args ) {

        return inputValidationService.processInput(args);
    }

    public VigorConfiguration getVigorConfiguration ( Namespace args ) throws VigorException {

        return initializationService.initializeVigor(args);
    }

    public List<Alignment> generateAlignments ( VirusGenome genome, String referenceDB, VigorConfiguration configuration ) throws VigorException {

        return alignmentGenerationService.generateAlignment(genome, referenceDB, configuration);
    }

    public List<Model> generateModels ( List<Alignment> alignments, VigorConfiguration configuration ) throws ServiceException {

        return modelGenerationService.generateModels(alignments, configuration);
    }

    public List<Model> generateGeneModels ( List<Model> models, VigorConfiguration configuration ) throws ServiceException {

        return geneModelGenerationService.generateGeneModel(models, configuration);
    }

    private Outfiles getOutfiles (VigorConfiguration config) throws IOException, VigorException {
        String outputDir = config.get(ConfigurationParameters.OutputDirectory);
        VigorUtils.checkFilePath("output directory", outputDir,
                                 VigorUtils.FileCheck.EXISTS,
                                 VigorUtils.FileCheck.WRITE,
                                 VigorUtils.FileCheck.DIRECTORY);
        boolean overwrite = config.getOrDefault(ConfigurationParameters.OverwriteOutputFiles, false);
        String fileBase = config.get(ConfigurationParameters.OutputPrefix);

        Outfiles outfiles = new Outfiles(Paths.get(outputDir), fileBase, overwrite);
        return outfiles;
    }

    private void writeEffectiveConfig(String outputDir, String outputPrefix, VigorConfiguration configuration) throws IOException {

        String dateString = new SimpleDateFormat("yyyyMMdd-HHmmss").format(new Date());
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(outputDir,
                                                                       String.format("%s-%s.ini", outputPrefix, dateString)),
                                                             Charset.forName("UTF-8"),
                                                             StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)) {
            writer.write(String.format("; Effective configuration %s\n\n", dateString));
            VigorConfiguration.ValueWithSource val;
            List<ConfigurationParameters> parameters = configuration.keySet()
                                                                    .stream()
                                                                    .sorted(Comparator.comparing(p -> p.configKey,
                                                                                                 String.CASE_INSENSITIVE_ORDER))
                                                                    .collect(Collectors.toList());
            for (ConfigurationParameters param: parameters) {
                val = configuration.getWithSource(param).get();
                writer.write("; source: ");
                writer.write(val.source);
                writer.newLine();
                writer.write(String.format("%s = \"%s\"", param.configKey, param.valueToString(val.value)));
                writer.newLine();
                writer.newLine();
            }
        }

    }

    public void initiateReportFile(String outputDir, String outputPrefix, int verbose){
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        Configuration config = lc.getConfiguration();

        FileAppender fa = FileAppender.newBuilder()
                                      .withName("mylogger")
                                      .withAppend(false)
                                      .withFileName(new File(outputDir, outputPrefix+".rpt").toString())
                                      .build();
        fa.start();
        config.addAppender(fa);

        Layout warningLayout = PatternLayout.newBuilder()
                                            .withConfiguration(config)
                                            .withPattern("%level %msg %exception{full}\n")
                                            .build();

        Filter warningFilter = LevelRangeFilter.createFilter(Level.FATAL, Level.WARN, Filter.Result.ACCEPT, Filter.Result.DENY);

        FileAppender warnings = FileAppender.newBuilder()
                                            .withName("__warnings")
                                            .withAppend(false)
                                            .withLayout(warningLayout)
                                            .withFileName(new File(outputDir, outputPrefix + ".warnings").toString())
                                            .build();
        warningFilter.start();
        warnings.addFilter(warningFilter);
        warnings.start();
        config.addAppender(warnings);

        lc.getLogger("org.jcvi.vigor").addAppender(warnings);
        lc.getLogger("org.jcvi.vigor").addAppender(fa);

        if (verbose > 0) {
            lc.getConfiguration().getLoggerConfig("org.jcvi.vigor").setLevel(verbose == 1 ? Level.DEBUG: Level.TRACE);
        }
        lc.updateLoggers();
    }


}

