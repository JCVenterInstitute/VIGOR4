package org.jcvi.vigor.RegressionTest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.Vigor;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.CommandLineParameters;
import org.jcvi.vigor.service.VigorInitializationService;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.*;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.stream.Collectors;

@Service
@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class GenerateVigor4GeneModels {

    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4Models.class);
    @Autowired
    Vigor vigor;
    @Autowired
    VigorInitializationService vigorInitializationService;

    public Map<String, List<Model>> generateModels ( String inputFASTA, String refDB, VigorConfiguration config ) {

        Map<String, List<Model>> vigor4Models = new HashMap<String, List<Model>>();
        String outputPrefix = config.get(ConfigurationParameters.OutputPrefix);
        VigorForm vigorForm = new VigorForm(config);
        String outputDir = config.get(ConfigurationParameters.OutputDirectory);
        try {
            GenerateVigorOutput.Outfiles outfiles = new GenerateVigorOutput.Outfiles();
            List<OpenOption> openOptionsList = new ArrayList<>();
            openOptionsList.add(StandardOpenOption.CREATE_NEW);
            OpenOption[] openOptions = openOptionsList.toArray(new OpenOption[] {});
            for (GenerateVigorOutput.Outfile outfile : GenerateVigorOutput.Outfile.values()) {
                outfiles.put(outfile, Files.newBufferedWriter(Paths.get(outputDir, outputPrefix + "." + outfile.extension),
                        Charset.forName("UTF-8"), openOptions));
            }
            outfiles.get(GenerateVigorOutput.Outfile.GFF3).write("##gff-version 3\n");
            VigorConfiguration.ValueWithSource unset = VigorConfiguration.ValueWithSource.of("NOTSET", "unknown");
            LOGGER.info(() -> config.entrySet()
                    .stream()
                    .sorted(Comparator.comparing(es -> es.getKey().configKey, String.CASE_INSENSITIVE_ORDER))
                    .map(e -> String.format("%-50s%s (%s)",
                            e.getKey().configKey,
                            e.getValue(),
                            config.getWithSource(e.getKey()).orElse(unset).source))
                    .collect(Collectors.joining("\n")));
            NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(new File(inputFASTA))
                    .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
                    .build();
            Iterator<NucleotideFastaRecord> i = dataStore.records().iterator();
            AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
            alignmentEvidence.setReference_db(refDB);
            vigorForm.setAlignmentEvidence(alignmentEvidence);
            while (i.hasNext()) {
                NucleotideFastaRecord record = i.next();
                vigorInitializationService.initiateReportFile(outputDir, outputPrefix, 0);
                VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
                        "1".equals(config.get(ConfigurationParameters.CompleteGene)),
                        "1".equals(config.get(ConfigurationParameters.CircularGene)));
                List<Alignment> alignments = vigor.generateAlignments(virusGenome, vigorForm);
                LOGGER.info("{} alignment(s) found for sequence {}", alignments.size(), record.getId());
                List<Model> candidateModels = vigor.generateModels(alignments, vigorForm);
                LOGGER.info("{} candidate model(s) found for sequence {}", candidateModels.size(), record.getId());
                List<Model> geneModels = vigor.generateGeneModels(candidateModels, vigorForm);
                LOGGER.info("{} gene model(s) found for sequence {}", geneModels.size(), record.getId());
                VigorUtils.deleteTempFiles(vigorForm.getTempDirectoryPath());
                if (geneModels.isEmpty()) {
                    LOGGER.warn("No gene models generated for sequence {}", record.getId());
                    continue;
                }
                geneModels = geneModels.stream()
                        .sorted(Comparator.comparing(g -> g.getRange(), Range.Comparators.ARRIVAL))
                        .collect(Collectors.toList());
                vigor.generateAlignmentOutput(outfiles, vigorForm);
                vigor.generateOutput(vigorForm.getConfiguration(), geneModels, outfiles);
                vigor.generateGFF3Output(geneModels, outfiles);
                FormatVigorOutput.printSequenceFeatures(geneModels, "GeneModels");
                outfiles.flush();
                vigor4Models.put(virusGenome.getId(), geneModels);
            }
        } catch (DataStoreException e) {
            LOGGER.error(String.format("problem reading input file %s", inputFASTA)
                    , e);
            System.exit(1);
        } catch (IOException e) {
            LOGGER.error("file problem", e);
            System.exit(1);
        } catch (ServiceException e) {
            LOGGER.error(e);
            System.exit(1);
        } catch (VigorException e) {
            LOGGER.error(e);
            System.exit(1);
        }
        return vigor4Models;
    }
}
