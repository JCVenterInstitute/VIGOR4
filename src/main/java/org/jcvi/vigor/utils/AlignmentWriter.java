package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

public class AlignmentWriter extends BaseOutputWriter {
    private static Logger LOGGER = LogManager.getLogger(AlignmentWriter.class);

    @Override
    public void writeModels(Outfiles outfiles, List<Model> models) throws VigorException, IOException {
        if (models.isEmpty()) {
            LOGGER.warn("No models to print to ALN file");
        }

        String genomeID = models.get(0).getAlignment().getVirusGenome().getId();
        OutputContext context = new OutputContext();
        context.addContext(OutputContext.Key.GENOME, genomeID);

        try (WriterBundle bw  = getWriter(outfiles, context, OutputContext.Key.GENOME)) {
            List<File> raw_files = models.stream()
                                         .map(m -> m.getAlignment().getAlignmentEvidence().getRaw_alignment())
                                         .distinct()
                                         .collect(Collectors.toList());

            for (File raw_alignment : raw_files) {
                printAlignment(bw, raw_alignment);
            }
        }
        List<File> temp_directories = models.stream()
                .map(m -> m.getAlignment().getAlignmentEvidence().getResults_directory())
                .distinct()
                .collect(Collectors.toList());
        for (File temp_directory: temp_directories) {
            VigorUtils.deleteTempFiles(temp_directory.getAbsolutePath());
        }
    }

    @Override
    public String getExtension() {
        return "aln";
    }

        // TODO don't write alignments one character at a time.
    public void printAlignment (WriterBundle bw, File inputFile ) {

        try {
            FileInputStream fRead = new FileInputStream(inputFile);
            int c;
            while (( c = fRead.read() ) != -1) {
                bw.write((char) c);
            }
            fRead.close();
        } catch (Exception e) {
            LOGGER.warn("Error reading temporory alignment file {}", inputFile.getAbsolutePath());
        }
    }

}
