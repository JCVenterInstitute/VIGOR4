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

    public Map<String, List<Model>> generateModels ( String inputFASTA, String refDB, VigorConfiguration config) throws VigorException {

        try {
            VigorForm vigorForm = new VigorForm(config);
            AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
            alignmentEvidence.setReference_db(refDB);
            vigorForm.setAlignmentEvidence(alignmentEvidence);
            vigor.generateAnnotations(inputFASTA, vigorForm);
            String outputDir = config.get(ConfigurationParameters.OutputDirectory);
            String outputPrefix = config.get(ConfigurationParameters.OutputPrefix);
            return modelsFromResults(outputDir, outputPrefix);
        } catch (IOException e) {
            throw new VigorException("Error generating vigor4 models", e);
        }
    }

    public  Map<String, List<Model>> modelsFromResults(String outputDirectory, String outputPrefix) throws IOException {
        return modelsFromResults(outputDirectory, outputPrefix, null);
    }
    public  Map<String, List<Model>> modelsFromResults(String outputDirectory, String outputPrefix, String inputFasta) throws IOException {
        String tblFile = Paths.get(outputDirectory, outputPrefix + ".tbl").toString();
        String pepFile = Paths.get(outputDirectory, outputPrefix + ".pep").toString();
        return new GenerateVigor3Models().generateModels(tblFile, pepFile, inputFasta);
    }
}
