package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Vigor;
import org.jcvi.vigor.forms.VigorForm;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;

@Service
public class GenerateAlignmentOuput {
    private static Logger LOGGER = LogManager.getLogger(Vigor.class);
    public void generateOutputFile(GenerateVigorOutput.Outfiles outfiles,VigorForm form) {
            printAlignment(outfiles.get(GenerateVigorOutput.Outfile.ALN), new File(form.getAlignmentOutputTempFile()));
            deleteTempFiles(form.getTempDirectoryPath());
    }

    public void printAlignment(BufferedWriter bw, File inputFile){
        try {
            FileInputStream fRead = new FileInputStream(inputFile);
            int c;
            while ((c = fRead.read()) != -1) {
                bw.write((char) c);
            }
            fRead.close();
        }
        catch(Exception e ){
            LOGGER.warn("Error reading temparory alignment file", inputFile.getAbsolutePath());
        }
    }
    public static void deleteTempFiles(String path) {
        Path workspace = new File(path).toPath();
        try {
            // clean up
            if (workspace != null) {
                Files.walk(workspace)
                        .map(Path::toFile)
                        .sorted(Comparator.reverseOrder())
                        .forEach(File::delete);
            }
        } catch (IOException e) {
            LOGGER.warn("Error deleting temporary working directory {}", workspace);
        }
    }


}

