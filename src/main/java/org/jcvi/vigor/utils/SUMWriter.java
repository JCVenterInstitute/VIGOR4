package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.Scores;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public class SUMWriter extends BaseOutputWriter {

    private static Logger LOGGER = LogManager.getLogger(SUMWriter.class);

    @Override
    public void writeModels(Outfiles outfiles, List<Model> models) throws VigorException, IOException {
        if (models.isEmpty()) {
            LOGGER.warn("no gene models to write to file");
            return;
        }

        OutputContext context = new OutputContext();
        context.addContext(OutputContext.Key.GENOME, models.get(0).getAlignment().getVirusGenome().getId());

        try (WriterBundle bw = getWriter(outfiles, context, OutputContext.Key.GENOME)) {

            double identityAvg = 0;
            double similarityAvg = 0;
            double coverageAvg = 0;
            long totalCDSBases = 0;
            long totalPepBases = 0;
            VirusGenome virusGenome = models.get(0).getAlignment().getVirusGenome();
            long seqLength = virusGenome.getSequence().getLength();
            bw.write("gene_id\t%identity\t%similarity\t%coverage\tstart..stop\tpep_size\tref_size\tref_id\tgene\tgene_product\n");
            for (Model model : models) {
                ViralProtein viralProtein = model.getAlignment().getViralProtein();
                Map<String, Double> scores = model.getScores();
                long cdsBases = 0;
                bw.write(model.getGeneID());
                bw.write("\t" + String.format("%.02f", scores.get(Scores.IDENTITY_SCORE)));
                bw.write("\t" + String.format("%.02f", scores.get(Scores.SIMILARITY_SCORE)));
                bw.write("\t" + String.format("%.02f", scores.get(Scores.COVERAGE_SCORE)) + "\t");
                for (int i = 0; i < model.getExons().size(); i++) {
                    Exon exon = model.getExons().get(i);
                    String start = Long.toString(exon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED));
                    String end = Long.toString(exon.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED));
                    if (i == 0 && model.isPartial5p()) {
                        start = "<" + start;
                    }
                    if (i == model.getExons().size() - 1 && model.isPartial3p()) {
                        end = ">" + end;
                    }
                    if (model.getExons().size() > 1 && i != 0) {
                        bw.write(", ");
                    }
                    bw.write(start + ".." + end);
                    cdsBases = cdsBases + exon.getRange().getLength();
                }
                if (!model.isPartial3p()) {
                    cdsBases = cdsBases - 3;
                }
                bw.write("\t" + (cdsBases) / 3);
                bw.write("\t" + viralProtein.getSequence().getLength());
                bw.write("\t" + viralProtein.getProteinID());
                bw.write("\t" + model.getGeneSymbol());
                bw.write("\t" + viralProtein.getProduct());

                bw.write(System.lineSeparator());
                totalCDSBases = totalCDSBases + cdsBases;
                totalPepBases = cdsBases + totalPepBases;
                identityAvg = identityAvg + scores.get(Scores.IDENTITY_SCORE);
                similarityAvg = similarityAvg + scores.get(Scores.SIMILARITY_SCORE);
                coverageAvg = coverageAvg + scores.get(Scores.COVERAGE_SCORE);

                IDGenerator idGenerator = IDGenerator.of(model.getGeneID());
                for (MaturePeptideMatch match : model.getMaturePeptides()) {
                    bw.write(idGenerator.next());
                    bw.write("\t" + String.format("%.02f", match.getIdentity() * 100));
                    bw.write("\t" + String.format("%.02f", match.getSimilarity() * 100));
                    bw.write("\t" + String.format("%.02f", match.getCoverage() * 100) + "\t");
                    List<Range> cdRanges = VigorFunctionalUtils.proteinRangeToCDSRanges(model, match.getProteinRange());
                    long start = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());
                    long end = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqLength, model.getDirection());

                    bw.write(OutputWriterUtils.formatMaturePeptideRange(model, match, cdRanges, Range.CoordinateSystem.RESIDUE_BASED,
                                                                        "..", start + model.getExons().get(0).getFrame().getFrame() - 1, end, true));
                    bw.write("\t" + match.getProteinRange().getLength());
                    bw.write("\t" + match.getReference().getSequence().getLength());
                    bw.write("\t" + match.getReference().getProteinID());
                    bw.write("\t" + model.getGeneSymbol());
                    bw.write("\t"
                                           + VigorUtils.putativeName(match.getReference().getProduct(), match.isFuzzyEnd(), match.isFuzzyBegin()));
                    bw.write(System.lineSeparator());
                }
            }
        }

    }

    @Override
    public String getExtension() {
        return "sum";
    }
}
