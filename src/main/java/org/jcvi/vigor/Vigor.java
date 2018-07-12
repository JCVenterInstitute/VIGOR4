package org.jcvi.vigor;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.PartialProteinSequence;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
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
            VigorForm vigorForm = getVigorForm(parsedArgs);
            generateAnnotations(inputFileName, vigorForm);
        } catch (VigorException e) {
            LOGGER.error(e);
            System.exit(1);
        }
    }

    public void generateAnnotations(String inputFileName, VigorForm vigorForm) throws VigorException {
        File inputFile = new File(inputFileName);
        if (!inputFile.exists()) {
            LOGGER.error("input file {} doesn't exists.", inputFileName);
            System.exit(1);
        } else if (!inputFile.canRead()) {
            LOGGER.error("input file {} isn't readable.", inputFileName);
            System.exit(1);
        }
        try {
            VigorConfiguration vigorParameters = vigorForm.getConfiguration();
            VigorConfiguration.ValueWithSource unset = VigorConfiguration.ValueWithSource.of("NOTSET", "unknown");
            LOGGER.info(() -> vigorParameters.entrySet()
                                             .stream()
                                             .sorted(Comparator.comparing(es -> es.getKey().configKey, String.CASE_INSENSITIVE_ORDER))
                                             .map(e -> String.format("%-50s%s (%s)",
                                                                     e.getKey().configKey,
                                                                     e.getValue(),
                                                                     vigorParameters.getWithSource(e.getKey()).orElse(unset).source))
                                             .collect(Collectors.joining("\n")));
            // TODO check file exists and is readable
            // TODO check output directory and permissions
            NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(inputFile)
                    .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
                    .build();
            // TODO move all this file handling to method
            // TODO checkout output earlier.
            String outputDir = vigorParameters.get(ConfigurationParameters.OutputDirectory);
            String outputPrefix =vigorParameters.get(ConfigurationParameters.OutputPrefix);
            writeEffectiveConfig(outputDir, vigorParameters);
            try (GenerateVigorOutput.Outfiles outfiles = getOutfiles(outputDir,
                                                                     outputPrefix,
                                                                     vigorParameters.get(ConfigurationParameters.OverwriteOutputFiles) == "true")) {
                outfiles.get(GenerateVigorOutput.Outfile.GFF3).write("##gff-version 3\n");
                Iterator<NucleotideFastaRecord> i = dataStore.records().iterator();
                while (i.hasNext()) {
                    NucleotideFastaRecord record = i.next();
                    List<Model> geneModels = modelsFromNucleotideRecord(record, vigorForm, vigorParameters);
                    if (geneModels.isEmpty()) {
                        LOGGER.warn("No gene models generated for sequence {}", record.getId());
                        continue;
                    }
                    generateAlignmentOutput(outfiles, vigorForm);
                    generateOutput(vigorParameters, geneModels, outfiles);
                    generateGFF3Output(geneModels, outfiles);
                    FormatVigorOutput.printSequenceFeatures(geneModels, "GeneModels");
                    outfiles.flush();
                }
            }
        } catch (DataStoreException e) {
            throw new VigorException(String.format("problem reading input file %s", inputFileName), e);
        } catch (IOException e) {
            throw new VigorException("File issue", e);
        }
    }

    public List<Model> modelsFromNucleotideRecord(NucleotideFastaRecord record, VigorForm vigorForm, VigorConfiguration vigorParameters) throws VigorException {
        VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
                                                  "1".equals(vigorParameters.get(ConfigurationParameters.CompleteGene)),
                                                  "1".equals(vigorParameters.get(ConfigurationParameters.CircularGene)));
        List<Alignment> alignments = generateAlignments(virusGenome, vigorForm);
        LOGGER.info("{} alignment(s) found for sequence {}", alignments.size(), record.getId());
        List<Model> candidateModels = generateModels(alignments, vigorForm);
        LOGGER.info("{} candidate model(s) found for sequence {}", candidateModels.size(), record.getId());
        List<Model> geneModels = generateGeneModels(candidateModels, vigorForm);
        LOGGER.info("{} gene model(s) found for sequence {}", geneModels.size(), record.getId());
        geneModels = findPeptides(vigorParameters, geneModels);
        // sort by begin,end
        return geneModels.stream()
                         .sorted(Comparator.comparing(g -> g.getRange(), Range.Comparators.ARRIVAL))
                         .collect(Collectors.toList());

    }

    private void setVerboseLogging ( int verbosity ) {

        Level verboseLevel = verbosity == 1 ? Level.DEBUG : Level.TRACE;
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        lc.getConfiguration().getLoggerConfig("org.jcvi.vigor").setLevel(verboseLevel);
        lc.updateLoggers();
    }

    private PeptideMatchingService.Scores getPeptideScores ( VigorConfiguration config ) {

        String minIdentityString = config.get(ConfigurationParameters.MaturePeptideMinimumIdentity);
        double minIdentity = Double.parseDouble(minIdentityString) / 100.0d;
        String minCoverageString = config.get(ConfigurationParameters.MaturePeptideMinimumCoverage);
        double minCoverage = Double.parseDouble(minCoverageString) / 100.0d;
        String minSimilarityString = config.get(ConfigurationParameters.MaturePeptideMinimumSimilarity);
        double minSimilarity = Double.parseDouble(minSimilarityString) / 100.0d;
        return PeptideMatchingService.Scores.of(minIdentity, minCoverage, minSimilarity);
    }

    private List<Model> findPeptides ( VigorConfiguration config, List<Model> geneModels ) throws VigorException {

        PeptideMatchingService.Scores scores = getPeptideScores(config);
        for (Model model : geneModels) {
            String maturePeptideDB = model.getAlignment().getAlignmentEvidence().getMatpep_db();
            // TODO check peptides for psuedogenes?
            if (!( maturePeptideDB == null || maturePeptideDB.isEmpty() )) {
                LOGGER.debug("finding mature peptides for {} using db {}", model.getGeneID(), maturePeptideDB);
                model.setMaturePeptides(peptideMatchingService.findPeptides(
                        PartialProteinSequence.of(model.getTranslatedSeq(), model.isPartial3p(), model.isPartial5p()),
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

    public VigorForm getVigorForm ( Namespace args ) throws VigorException {

        return initializationService.initializeVigor(args);
    }

    public List<Alignment> generateAlignments ( VirusGenome genome, VigorForm form ) throws VigorException {

        return alignmentGenerationService.generateAlignment(genome, form);
    }

    public List<Model> generateModels ( List<Alignment> alignments, VigorForm form ) throws ServiceException {

        return modelGenerationService.generateModels(alignments, form);
    }

    public List<Model> generateGeneModels ( List<Model> models, VigorForm form ) throws ServiceException {

        return geneModelGenerationService.generateGeneModel(models, form);
    }

    public void generateOutput ( VigorConfiguration config, List<Model> models, GenerateVigorOutput.Outfiles outfiles ) throws IOException {

        generateVigorOutput.generateOutputFiles(config, outfiles, models);
    }

    public void generateGFF3Output ( List<Model> models, GenerateVigorOutput.Outfiles outfiles ) throws IOException {

        generateGFF3Output.generateOutputFile(outfiles, models);
    }

    public void generateAlignmentOutput ( GenerateVigorOutput.Outfiles outfiles, VigorForm form ) {

        generateAlignmentOuput.generateOutputFile(outfiles, form);
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
            writer.write(String.format("# Effective configuration %s\n\n", dateString));
            VigorConfiguration.ValueWithSource val;
            for (ConfigurationParameters param: configuration.keySet()) {
                val = configuration.getWithSource(param).get();
                writer.write("# source: ");
                writer.write(val.source);
                writer.newLine();
                writer.write(String.format("%s = \"%s\"", param.configKey, val.value));
                writer.newLine();
                writer.newLine();
            }
        }

    }
}

