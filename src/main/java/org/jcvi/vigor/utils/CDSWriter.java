package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.jcvi.vigor.utils.OutputWriterUtils.getDefline;

public class CDSWriter extends BaseOutputWriter {

    private static Logger LOGGER = LogManager.getLogger(CDSWriter.class);

    @Override
    public String getExtension() {
        return "cds";
    }

    @Override
    public void writeModels(Outfiles outfiles, List<Model> geneModels) throws IOException, VigorException {
        if (geneModels.isEmpty()) {
            LOGGER.warn("No gene models to write CDS reports");
            return;
        }
        OutputContext context = new OutputContext();
        context.addContext(OutputContext.Key.GENOME, geneModels.get(0).getGeneID());
        try (WriterBundle unused = getWriter(outfiles, context, OutputContext.Key.GENOME)) {
            for (Model model : geneModels) {
                context.addContext(OutputContext.Key.GENE, model.getGeneID());
                try (WriterBundle bw = getWriter(outfiles, context, OutputContext.Key.GENE)) {
                    bw.write(getDefline(model));
                    bw.newLine();
                    NucleotideSequenceBuilder builder = new NucleotideSequenceBuilder();
                    NucleotideSequence virusGenome = model.getAlignment().getVirusGenome().getSequence();
                    long seqLength = virusGenome.getLength();
                    List<Range> translatedRanges = model.getExons()
                                                        .stream()
                                                        .map(e -> VigorFunctionalUtils.getDirectionBasedRange(e.getRange(),
                                                                                                              seqLength,
                                                                                                              model.getDirection()))
                                                        .collect(Collectors.toList());
                    Collections.sort(translatedRanges, Range.Comparators.ARRIVAL);
                    for (Range exonRange : translatedRanges) {
                        builder.append(virusGenome.toBuilder(exonRange).build());
                    }

                    OutputWriterUtils.writeSequence(bw, builder.build());
                }
            }
        }
    }
}
