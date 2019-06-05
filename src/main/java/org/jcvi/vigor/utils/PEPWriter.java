package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.MaturePeptideMatch;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;

import static org.jcvi.vigor.utils.OutputWriterUtils.formatMaturePeptideRange;

public class PEPWriter extends BaseOutputWriter {

    private static Logger LOGGER = LogManager.getLogger(PEPWriter.class);

    @Override
    public void writeModels(Outfiles outfiles, List<Model> models) throws VigorException, IOException {
        StringBuilder defline;
        long seqLength = models.get(0).getAlignment().getVirusGenome().getSequence().getLength();
        OutputContext context = new OutputContext();
        context.addContext(OutputContext.Key.GENOME, models.get(0).getGeneID());
        // this ensures that the genome writer is closed
        try (WriterBundle unused = getWriter(outfiles, context, OutputContext.Key.GENOME)) {
            for (Model model : models) {
                context.addContext(OutputContext.Key.GENE, model.getGeneID());
                try (WriterBundle sequenceWriter = getWriter(outfiles, context, OutputContext.Key.GENE)) {
                    sequenceWriter.write(OutputWriterUtils.getDefline(model));
                    sequenceWriter.newLine();
                    OutputWriterUtils.writeSequence(sequenceWriter, model.getTranslatedSeq());
                }
                //
                context.removeContext(OutputContext.Key.GENE);
                IDGenerator idGenerator = IDGenerator.of(model.getGeneID());
                for (MaturePeptideMatch match : model.getMaturePeptides()) {
                    String pepID = idGenerator.next();
                    context.addContext(OutputContext.Key.PEP, pepID);
                    try (WriterBundle bw = getWriter(outfiles, context, OutputContext.Key.PEP)) {
                        defline = new StringBuilder();
                        defline.append(">" + pepID);
                        if (model.isPseudogene()) {
                            defline.append(" pseudogene");
                        }
                        defline.append(" mat_peptide");
                        List<Range> cdsRanges = VigorFunctionalUtils.proteinRangeToCDSRanges(model, match.getProteinRange());
                        // TODO handle truncation etc
                        Exon initialExon = model.getExons().get(0);
                        long startCoordinate = VigorFunctionalUtils.getDirectionBasedCoordinate(initialExon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());
                        long endCoordinate = VigorFunctionalUtils.getDirectionBasedCoordinate(model.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());
                        defline.append(String.format(" location=%s", formatMaturePeptideRange(model,
                                                                                              match,
                                                                                              cdsRanges,
                                                                                              Range.CoordinateSystem.RESIDUE_BASED,
                                                                                              "..",
                                                                                              // start_codon adjustment
                                                                                              startCoordinate + initialExon.getFrame().getFrame() - 1,
                                                                                              endCoordinate, true)));
                        defline.append(String.format(" gene=\"%s\"", model.getGeneSymbol()));
                        String product = match.getReference().getProduct();
                        if (!NullUtil.isNullOrEmpty(product)) {
                            defline.append(String.format(" product=\"%s\"", VigorUtils.putativeName(product, match.isFuzzyEnd(), match.isFuzzyBegin())));
                        } else {
                            LOGGER.warn("Missing product for {}", idGenerator.getCurrent());
                        }
                        String refDB = model.getAlignment().getAlignmentEvidence().getMatpep_db();
                        if (!NullUtil.isNullOrEmpty(refDB)) {
                            defline.append(String.format(" ref_db=\"%s\"", Paths.get(refDB).getFileName().toString()));
                        }
                        defline.append(String.format(" ref_id=\"%s\"", match.getReference().getProteinID()));
                        bw.write(defline.toString());
                        bw.newLine();
                        OutputWriterUtils.writeSequence(bw, match.getProtein().toBuilder().trim(match.getProteinRange()).build());
                    }
                }
            }
        }
    }

    @Override
    public String getExtension() {
        return "pep";
    }
}
