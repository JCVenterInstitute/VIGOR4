package org.jcvi.vigor;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
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
    private GenerateVigorOutput generateVigorOutput;
    @Autowired
    private GenerateGFF3Output generateGFF3Output;
    @Autowired
    private PeptideMatchingService peptideMatchingService;
    @Autowired
    private GenerateAlignmentOuput generateAlignmentOuput;

    public void run ( String... args ) {

        Namespace parsedArgs = parseArgs(args);
        int verbosity = parsedArgs.getInt(CommandLineParameters.verbose);
        if (verbosity > 0) {
            setVerboseLogging(verbosity);
            LOGGER.debug("verbose logging enabled");
        }
        String inputFileName = parsedArgs.getString("input_fasta");
        try {
            VigorConfiguration vigorConfiguration = getVigorConfiguration(parsedArgs);
            checkConfig(vigorConfiguration);
            String referenceDB = vigorConfiguration.get(ConfigurationParameters.ReferenceDatabaseFile);
            generateAnnotations(inputFileName, referenceDB, vigorConfiguration);
        } catch (Exception e) {
            LOGGER.error("Exception encountered: exiting 1",e);
            System.exit(1);
        }
    }

    /**
     * Check that required parameters are set
     * @param config
     * @throws VigorException
     */
    private void checkConfig(VigorConfiguration config) throws VigorException {
        List<String> errors =  new ArrayList<>();
        Object val;
        for (ConfigurationParameters parameter: Arrays.stream(ConfigurationParameters.values())
                                                      .filter(p -> p.hasFlag(ConfigurationParameters.Flags.REQUIRED))
                                                      .collect(Collectors.toList())) {
            val = config.get(parameter);
            if (val == null || (val instanceof String && ((String) val).isEmpty())) {
                errors.add(String.format("parameter %s is not set but is required", parameter.configKey));
            }
        }

        if (! errors.isEmpty()) {
            throw new VigorException(String.join("\n", errors));
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
        VigorUtils.checkFilePath("input file", inputFileName, VigorUtils.FileCheck.EXISTS, VigorUtils.FileCheck.READ);
        printConfiguration(vigorParameters);
        String outputDir = vigorParameters.get(ConfigurationParameters.OutputDirectory);
        String outputPrefix = vigorParameters.get(ConfigurationParameters.OutputPrefix);
        VigorUtils.checkFilePath("output directory", outputDir,
                                 VigorUtils.FileCheck.EXISTS,
                                 VigorUtils.FileCheck.WRITE,
                                 VigorUtils.FileCheck.DIRECTORY);
        try (NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(new File(inputFileName))
                .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
                .build();
             GenerateVigorOutput.Outfiles outfiles = getOutfiles(outputDir,
                                                                 outputPrefix,
                                                                 vigorParameters.getOrDefault(ConfigurationParameters.OverwriteOutputFiles, false))
        ) {
            // TODO move all this file handling to method
            // TODO checkout output earlier.
            writeEffectiveConfig(outputDir, vigorParameters);
            outfiles.get(GenerateVigorOutput.Outfile.GFF3).write("##gff-version 3\n");
            Iterator<NucleotideFastaRecord> recordIterator = dataStore.records().iterator();
            while (recordIterator.hasNext()) {
                NucleotideFastaRecord record = recordIterator.next();
                LOGGER.debug("processing {}", record.getId());
                List<Model> geneModels = modelsFromNucleotideRecord(record, referenceDB, vigorParameters);
                if (geneModels.isEmpty()) {
                    LOGGER.warn("No gene models generated for sequence {}", record.getId());
                    continue;
                }
                outputModels(vigorParameters, outfiles, geneModels);
            }
        } catch (DataStoreException e) {
            throw new VigorException(String.format("problem reading input file %s", inputFileName), e);
        } catch (IOException e) {
            throw new VigorException(String.format("File issue. Got %s: %s", e.getClass().getSimpleName(), e.getMessage()), e);
        }
    }

    public void outputModels(VigorConfiguration vigorParameters, GenerateVigorOutput.Outfiles outfiles, List<Model> geneModels) throws IOException {
        generateAlignmentOutput(geneModels, outfiles);
        generateOutput(vigorParameters, geneModels, outfiles);
        generateGFF3Output(geneModels, outfiles);
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

        Level verboseLevel = verbosity == 1 ? Level.DEBUG : Level.TRACE;
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        lc.getConfiguration().getLoggerConfig("org.jcvi.vigor").setLevel(verboseLevel);
        lc.updateLoggers();
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

    public void generateOutput ( VigorConfiguration config, List<Model> models, GenerateVigorOutput.Outfiles outfiles ) throws IOException {

        generateVigorOutput.generateOutputFiles(config, outfiles, models);
    }

    public void generateGFF3Output ( List<Model> models, GenerateVigorOutput.Outfiles outfiles ) throws IOException {

        generateGFF3Output.generateOutputFile(outfiles, models);
    }

    public void generateAlignmentOutput ( List<Model> models, GenerateVigorOutput.Outfiles outfiles) {

        generateAlignmentOuput.generateOutputFile(outfiles, models);
    }

    private GenerateVigorOutput.Outfiles getOutfiles ( String outputDir, String outputPrefix, boolean overwrite ) throws IOException {

        GenerateVigorOutput.Outfiles outfiles = new GenerateVigorOutput.Outfiles();
        List<OpenOption> openOptionsList = new ArrayList<>();
        if (overwrite) {
            openOptionsList.add(StandardOpenOption.CREATE);
            openOptionsList.add(StandardOpenOption.TRUNCATE_EXISTING);
        } else {
            openOptionsList.add(StandardOpenOption.CREATE_NEW);
        }
        OpenOption[] openOptions = openOptionsList.toArray(new OpenOption[] {});
        for (GenerateVigorOutput.Outfile outfile : GenerateVigorOutput.Outfile.values()) {
            outfiles.put(outfile, Files.newBufferedWriter(Paths.get(outputDir, outputPrefix + "." + outfile.extension),
                    Charset.forName("UTF-8"), openOptions));
        }
        return outfiles;
    }

    private void writeEffectiveConfig(String outputDir, VigorConfiguration configuration) throws IOException {

        String dateString = new SimpleDateFormat("yyyyMMdd-HHmmss").format(new Date());
        try (BufferedWriter writer = Files.newBufferedWriter(Paths.get(outputDir, String.format("vigor-%s.ini", dateString)),
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
                writer.write(String.format("%s = \"%s\"", param.configKey, val.value));
                writer.newLine();
                writer.newLine();
            }
        }

    }


}

