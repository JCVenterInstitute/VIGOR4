package org.jcvi.vigor.RegressionTest;

import javafx.util.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.Vigor;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.*;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

@Service
@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class GenerateVigor4GeneModels {

    private final static Logger LOGGER = LogManager.getLogger(ValidateVigor4ModelsTest.class);
    @Autowired
    Vigor vigor;

    public Map<String, List<Model>> generateModels ( String inputFASTA, String refDB, VigorConfiguration config ) throws VigorException {

        try {
            VigorForm vigorForm = new VigorForm(config);
            AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
            alignmentEvidence.setReference_db(refDB);
            vigorForm.setAlignmentEvidence(alignmentEvidence);
            String outputDir = config.get(ConfigurationParameters.OutputDirectory);
            String outputPrefix = config.get(ConfigurationParameters.OutputPrefix);
            Pair<String, Boolean> outputFile = getVigor4OutputFiles(outputDir, outputPrefix);
            if (!( outputFile.getValue() )) {
                vigor.generateAnnotations(inputFASTA, vigorForm);
            }
            return modelsFromResults(outputFile.getKey());
        } catch (IOException e) {
            throw new VigorException("Error generating vigor4 models", e);
        }
    }

    public Map<String, List<Model>> modelsFromResults ( String tblFile ) throws IOException {

        return modelsFromResults(tblFile, null);
    }

    public Map<String, List<Model>> modelsFromResults ( String tblFile, String inputFasta ) throws IOException {

        return new GenerateReferenceModels().generateModels(tblFile, inputFasta);
    }

    private Pair<String, Boolean> getVigor4OutputFiles ( String outputDirectory, String outputPrefix ) {

        String tblFilePath = Paths.get(outputDirectory, outputPrefix + ".tbl").toString();
        Path tblPath = Paths.get(tblFilePath);
        return new Pair<>(tblFilePath, Files.exists(tblPath));
    }
}
