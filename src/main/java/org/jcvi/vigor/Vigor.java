package org.jcvi.vigor;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.*;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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


	public void run(String ... args) {

        Namespace parsedArgs = parseArgs(args);

        String inputFileName = parsedArgs.getString("input_fasta");
        File inputFile = new File(inputFileName);
        if (! inputFile.exists()) {
            LOGGER.error("input file {} doesn't exists.", inputFileName);
            System.exit(1);
        } else if (! inputFile.canRead()) {
            LOGGER.error("input file {} isn't readable.", inputFileName);
            System.exit(1);
        }
        try{
            VigorForm vigorForm = getVigorForm(parsedArgs);
            VigorConfiguration vigorParameters = vigorForm.getConfiguration();
            VigorConfiguration.ValueWithSource unset = VigorConfiguration.ValueWithSource.of("NOTSET", "unknown");
            LOGGER.info( () ->  vigorParameters.entrySet()
                                               .stream()
                                               .sorted(Comparator.comparing(es -> es.getKey().configKey, String.CASE_INSENSITIVE_ORDER))
                                               .map( e -> String.format("%-50s%s (%s)",
                                                       e.getKey().configKey,
                                                       e.getValue(),
                                                       vigorParameters.getWithSource(e.getKey()).orElse(unset).source))
                                               .collect(Collectors.joining("\n")) );
            // TODO check file exists and is readable
            // TODO check output directory and permissions
            NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(inputFile)
                    .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
                    .build();

            // TODO move all this file handling to method
            // TODO checkout output earlier.
            String outputDir = vigorForm.getConfiguration().get(ConfigurationParameters.OutputDirectory);
            String outputPrefix = vigorForm.getConfiguration().get(ConfigurationParameters.OutputPrefix);
            try (GenerateVigorOutput.Outfiles outfiles = getOutfiles(outputDir,
                    outputPrefix,
                    vigorForm.getConfiguration().get(ConfigurationParameters.OverwriteOutputFiles) == "true")) {
                outfiles.get(GenerateVigorOutput.Outfile.GFF3).write("##gff-version 3\n");
                PeptideMatchingService.Scores peptideScores = getPeptideScores(vigorParameters);
                Iterator<NucleotideFastaRecord> i = dataStore.records().iterator();
                while (i.hasNext()) {
                    NucleotideFastaRecord record = i.next();
                    VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
                            "1".equals(vigorParameters.get(ConfigurationParameters.CompleteGene)),
                            "1".equals(vigorParameters.get(ConfigurationParameters.CircularGene)));

                    // Call referenceDBGenerationService methods to generate alignmentEvidence.
                    List<Alignment> alignments = generateAlignments(virusGenome, vigorForm);
                    LOGGER.info("{} alignment(s) found for sequence {}", alignments.size(), record.getId());
                    List<Model> candidateModels = generateModels(alignments, vigorForm);
                    LOGGER.info("{} candidate model(s) found for sequence {}", candidateModels.size(), record.getId());
                    List<Model> geneModels = generateGeneModels(candidateModels, vigorForm);
                    LOGGER.info("{} gene model(s) found for sequence {}", geneModels.size(), record.getId());
                    geneModels = findPeptides(vigorParameters, geneModels, vigorForm);
                    if (geneModels.isEmpty()) {
                        LOGGER.warn("No gene models generated for sequence {}", record.getId());
                        continue;
                    }
                    // sort by begin,end
                    geneModels = geneModels.stream()
                                           .sorted(Comparator.comparing(g -> g.getRange(), Range.Comparators.ARRIVAL))
                                           .collect(Collectors.toList());
                    generateOutput(vigorParameters, geneModels, outfiles);
                    generateGFF3Output(vigorParameters, geneModels, outfiles);
                    FormatVigorOutput.printSequenceFeatures(geneModels,"GeneModels");
                    outfiles.flush();
                }
            }


        } catch (DataStoreException e) {
            LOGGER.error(String.format("problem reading input file %s", inputFileName), e);
            System.exit(1);
        } catch (IOException e) {
            LOGGER.error("file problem", e);
            System.exit(1);
        } catch (VigorException e) {
            LOGGER.error(e);
            System.exit(1);
        }

    }

    private PeptideMatchingService.Scores getPeptideScores(VigorConfiguration config) {
	    String minIdentityString = config.get(ConfigurationParameters.MaturePeptideMinimumIdentity);
	    double minIdentity = Double.parseDouble(minIdentityString)/ 100.0d;

        String minCoverageString = config.get(ConfigurationParameters.MaturePeptideMinimumCoverage);
        double minCoverage = Double.parseDouble(minCoverageString)/ 100.0d;

        String minSimilarityString = config.get(ConfigurationParameters.MaturePeptideMinimumSimilarity);
        double minSimilarity = Double.parseDouble(minSimilarityString) /100.0d;

        return PeptideMatchingService.Scores.of(minIdentity, minCoverage, minSimilarity);
    }

    private List<Model> findPeptides(VigorConfiguration config, List<Model> geneModels, VigorForm vigorForm) throws VigorException {

	    PeptideMatchingService.Scores scores = getPeptideScores(config);

	    for (Model model: geneModels) {
            String maturePeptideDB = model.getAlignment().getAlignmentEvidence().getMatpep_db();
            // TODO check peptides for psuedogenes?
            if (! (maturePeptideDB == null || maturePeptideDB.isEmpty()) ) {
                model.setMaturePeptides(peptideMatchingService.findPeptides(model.getTanslatedSeq(),
                        new File(maturePeptideDB), scores));
            }
        }
        return geneModels;
    }



    public Namespace parseArgs(String[] args) {
        return inputValidationService.processInput(args);
    }

    public VigorForm getVigorForm(Namespace args) throws VigorException {
        return initializationService.initializeVigor(args);
    }

    public List<Alignment> generateAlignments(VirusGenome genome, VigorForm form) throws VigorException{
        return alignmentGenerationService.generateAlignment(genome, form);
    }

    public List<Model> generateModels(List<Alignment> alignments, VigorForm form) throws ServiceException {
        return modelGenerationService.generateModels(alignments, form);
    }

    public List<Model> generateGeneModels(List<Model> models, VigorForm form) throws ServiceException {
        return geneModelGenerationService.generateGeneModel(models, form);
    }

    public void generateOutput(VigorConfiguration config, List<Model> models, GenerateVigorOutput.Outfiles outfiles) throws ServiceException, IOException{
        generateVigorOutput.generateOutputFiles(config, outfiles, models);
    }
    public void generateGFF3Output(VigorConfiguration config, List<Model> models, GenerateVigorOutput.Outfiles outfiles) throws IOException {
        generateGFF3Output.generateOutputFile(config, outfiles, models);
    }


    private GenerateVigorOutput.Outfiles getOutfiles(String outputDir, String outputPrefix, boolean overwrite) throws IOException {

        GenerateVigorOutput.Outfiles outfiles = new GenerateVigorOutput.Outfiles();
        List<OpenOption> openOptionsList = new ArrayList<>();
        if (overwrite) {
            openOptionsList.add(StandardOpenOption.CREATE);
            openOptionsList.add(StandardOpenOption.TRUNCATE_EXISTING);
        } else {
            openOptionsList.add(StandardOpenOption.CREATE_NEW);
        }

        OpenOption[] openOptions = openOptionsList.toArray(new OpenOption[]{});
            for (GenerateVigorOutput.Outfile outfile : GenerateVigorOutput.Outfile.values()) {
                outfiles.put(outfile, Files.newBufferedWriter(Paths.get(outputDir, outputPrefix + "." + outfile.extension),
                        Charset.forName("UTF-8"), openOptions));
            }
        return outfiles;
    }

}

