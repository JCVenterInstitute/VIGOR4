package org.jcvi.vigor.RegressionTest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.Vigor;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.utils.*;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

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

    private class Tuple<K,V> {
        public final K first;
        public final V second;

        Tuple(K first, V second) {
            this.first = first;
            this.second = second;
        }

    }
    @Autowired
    Vigor vigor;

    public Map<String, List<Model>> generateModels ( String inputFASTA, String refDB, VigorConfiguration config ) throws VigorException {

        try {

            config.put(ConfigurationParameters.ReferenceDatabaseFile, refDB);
            String outputDir = config.get(ConfigurationParameters.OutputDirectory);
            String outputPrefix = config.get(ConfigurationParameters.OutputPrefix);
            Tuple<String, Boolean> outputFile = getVigor4OutputFiles(outputDir, outputPrefix);
            boolean overwrite = config.getOrDefault(ConfigurationParameters.OverwriteOutputFiles, false);
            if (! outputFile.second  || overwrite) {
                vigor.generateAnnotations(inputFASTA, refDB, config);
            }
            return modelsFromResults(outputFile.first);
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

    private Tuple<String, Boolean> getVigor4OutputFiles ( String outputDirectory, String outputPrefix ) {

        String tblFilePath = Paths.get(outputDirectory, outputPrefix + ".tbl").toString();
        Path tblPath = Paths.get(tblFilePath);
        return new Tuple<>(tblFilePath, Files.exists(tblPath));
    }
}
