package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Vigor;
import org.jcvi.vigor.component.Model;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.util.List;
import java.util.stream.Collectors;

@Service
public class GenerateAlignmentOuput {

    private static Logger LOGGER = LogManager.getLogger(Vigor.class);

    public void generateOutputFile ( GenerateVigorOutput.Outfiles outfiles, List<Model> models ) {
        List<File> raw_files = models.stream()
                                    .map(m -> m.getAlignment().getAlignmentEvidence().getRaw_alignment())
                                    .distinct().
                collect(Collectors.toList());
        for (File raw_alignment: raw_files) {
            printAlignment(outfiles.get(GenerateVigorOutput.Outfile.ALN), raw_alignment);
        }
        List<File> temp_directories = models.stream()
                .map(m -> m.getAlignment().getAlignmentEvidence().getResults_directory())
                .distinct()
                .collect(Collectors.toList());
        for (File temp_directory: temp_directories) {
            VigorUtils.deleteTempFiles(temp_directory.getAbsolutePath().toString());
        }
    }

    public void printAlignment ( BufferedWriter bw, File inputFile ) {

        try {
            FileInputStream fRead = new FileInputStream(inputFile);
            int c;
            while (( c = fRead.read() ) != -1) {
                bw.write((char) c);
            }
            fRead.close();
        } catch (Exception e) {
            LOGGER.warn("Error reading temporory alignment file", inputFile.getAbsolutePath());
        }
    }
}

