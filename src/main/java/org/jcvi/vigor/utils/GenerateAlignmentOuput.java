package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Vigor;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;
import org.springframework.stereotype.Service;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

@Service
public class GenerateAlignmentOuput {

    private static Logger LOGGER = LogManager.getLogger(Vigor.class);

    public void generateOutputFile ( Outfiles outfiles, List<Model> models ) throws IOException, VigorException {
        if (models.isEmpty()) {
            LOGGER.warn("No models to print to ALN file");
        }
        String genomeID = models.get(0).getAlignment().getVirusGenome().getId();
        Path genomeAlignmentPath = Paths.get(GenerateVigorOutput.getSequenceFilePath(genomeID) + ".aln");
        BufferedWriter alignmentWriter = outfiles.getWriter( GenerateVigorOutput.Outfile.ALN);
        try (BufferedWriter genomeWriter = outfiles.getWriter( genomeAlignmentPath)) {
            List<File> raw_files = models.stream()
                                         .map(m -> m.getAlignment().getAlignmentEvidence().getRaw_alignment())
                                         .distinct().
                                                 collect(Collectors.toList());

            GenerateVigorOutput.WriterBundle bw = new GenerateVigorOutput.WriterBundle(alignmentWriter, genomeWriter);
            for (File raw_alignment : raw_files) {
                printAlignment(bw, raw_alignment);
            }
        } catch (IOException e ) {
            String message = String.format("Issue opening/closing genome file %s for alignment", genomeAlignmentPath);
            LOGGER.error(message, e);
            throw e;
        }
        List<File> temp_directories = models.stream()
                .map(m -> m.getAlignment().getAlignmentEvidence().getResults_directory())
                .distinct()
                .collect(Collectors.toList());
        for (File temp_directory: temp_directories) {
            VigorUtils.deleteTempFiles(temp_directory.getAbsolutePath().toString());
        }
    }

    // TODO don't write alignments one character at a time.
    public void printAlignment (GenerateVigorOutput.WriterBundle bw, File inputFile ) {

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

