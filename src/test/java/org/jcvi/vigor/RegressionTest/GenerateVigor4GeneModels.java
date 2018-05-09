package org.jcvi.vigor.RegressionTest;

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
import java.util.*;



@Service
@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class GenerateVigor4GeneModels {
    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4Models.class);

    @Autowired
    Vigor vigor;
    public Map<String, List<Model>> generateModels(String inputFASTA, String refDB,VigorConfiguration config)
    {
        Map<String, List<Model>> vigor4Models = new HashMap<String, List<Model>>();
       try{
           NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(new File(inputFASTA))
                   .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
                   .build();
           Iterator<NucleotideFastaRecord> i = dataStore.records().iterator();
           VigorForm vigorForm = new VigorForm(config);
           AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
           alignmentEvidence.setReference_db(refDB);
           vigorForm.setAlignmentEvidence(alignmentEvidence);
           while (i.hasNext()) {
               NucleotideFastaRecord record = i.next();
               VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
                       "1".equals(config.get(ConfigurationParameters.CompleteGene)),
                       "1".equals(config.get(ConfigurationParameters.CircularGene)));

               // Call referenceDBGenerationService methods to generate alignmentEvidence.
               List<Alignment> alignments = vigor.generateAlignments(virusGenome, vigorForm);
               LOGGER.info("{} alignment(s) found for sequence {}", alignments.size(), record.getId());
               List<Model> candidateModels = vigor.generateModels(alignments, vigorForm);
               LOGGER.info("{} candidate model(s) found for sequence {}", candidateModels.size(), record.getId());
               List<Model> geneModels = vigor.generateGeneModels(candidateModels, vigorForm);
               LOGGER.info("{} gene model(s) found for sequence {}", geneModels.size(), record.getId());
               if (geneModels.isEmpty()) {
                   LOGGER.warn("No gene models generated for sequence {}", record.getId());
                   continue;
               }
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
        }
        catch(VigorException e){
            LOGGER.error(e);
            System.exit(1);
        }
    return vigor4Models;

    }


    public VigorForm setVigorParameters(VigorForm form, String workspace){
        VigorConfiguration configuration = new VigorConfiguration("test");
        configuration.put(ConfigurationParameters.OutputDirectory,workspace+"/TEST");
        form.setConfiguration(configuration);
        return form;
    }
}
