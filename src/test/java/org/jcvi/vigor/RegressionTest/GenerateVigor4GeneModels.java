/*
package org.jcvi.vigor.RegressionTest;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
import org.jcvi.vigor.service.*;
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
    public Map<String, List<Model>> generateModels(String workspace,
                                                   String inputFilePath, String refDB)
    {
        String[] args = new String[] {};
        args.
        Namespace parsedArgs = vigor.parseArgs(args);

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
            VigorForm vigorForm = vigor.getVigorForm(parsedArgs);
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
            GenerateVigorOutput.Outfiles outfiles = new GenerateVigorOutput.Outfiles();
            List<OpenOption> openOptionsList = new ArrayList<>();
            if (vigorForm.getConfiguration().get(ConfigurationParameters.OverwriteOutputFiles) == "true") {
                openOptionsList.add(StandardOpenOption.CREATE);
                openOptionsList.add(StandardOpenOption.TRUNCATE_EXISTING);
            } else {
                openOptionsList.add(StandardOpenOption.CREATE_NEW);
            }

            OpenOption[] openOptions =  openOptionsList.toArray(new OpenOption[] {});
            for (GenerateVigorOutput.Outfile outfile: GenerateVigorOutput.Outfile.values()) {
                outfiles.put(outfile, Files.newBufferedWriter(Paths.get(outputDir, outputPrefix + "." + outfile.extension),
                        Charset.forName("UTF-8"), openOptions));
            }

            Iterator<NucleotideFastaRecord> i = dataStore.records().iterator();
            while (i.hasNext()) {
                NucleotideFastaRecord record = i.next();
                VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
                        "1".equals(vigorParameters.get(ConfigurationParameters.CompleteGene)),
                        "1".equals(vigorParameters.get(ConfigurationParameters.CircularGene)));

                // Call referenceDBGenerationService methods to generate alignmentEvidence.
                List<Alignment> alignments = vigor.generateAlignments(virusGenome, vigorForm);
                LOGGER.info("{} alignment(s) found for sequence {}", alignments.size(), record.getId());
                List<Model> candidateModels = generateModels(alignments, vigorForm);
                LOGGER.info("{} candidate model(s) found for sequence {}", candidateModels.size(), record.getId());
                List<Model> geneModels = vigor.generateGeneModels(candidateModels, vigorForm);
                LOGGER.info("{} gene model(s) found for sequence {}", geneModels.size(), record.getId());
                geneModels = vigor.findPeptides(geneModels, vigorForm);
                if (geneModels.isEmpty()) {
                    LOGGER.warn("No gene models generated for sequence {}", record.getId());
                    continue;
                }

            }

        } catch (DataStoreException e) {
            LOGGER.error(String.format("problem reading input file %s", inputFileName)
                    , e);
            System.exit(1);
        } catch (IOException e) {
            LOGGER.error("file problem", e);
            System.exit(1);
        } catch (ServiceException e) {
            LOGGER.error(e);
            System.exit(1);
        }
        catch(VigorException e){
            LOGGER.error(e);
            System.exit(1);
        }


    }


    public VigorForm setVigorParameters(VigorForm form, String workspace){
        VigorConfiguration configuration = new VigorConfiguration("test");
        configuration.put(ConfigurationParameters.OutputDirectory,workspace+"/TEST");
        form.setConfiguration(configuration);
        return form;
    }
}
*/
