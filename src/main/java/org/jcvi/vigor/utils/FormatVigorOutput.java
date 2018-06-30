package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.component.*;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;

/**
 * Created by snettem on 5/24/2017.
 */
public class FormatVigorOutput {

    private static final Logger LOGGER = LogManager.getLogger(FormatVigorOutput.class);
    private static Range.CoordinateSystem oneBased = Range.CoordinateSystem.RESIDUE_BASED;

    public static void printModels ( List<Model> models, String message ) {

        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
                "*********************************" + message + "*****************************************************");
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-20s%-20s%-20s%-20s%-20s%-30s", "Protein_ID", "geneSymbol", "Direction", "spliceform", "NTSeqRange", "AASeqRange", "Scores"));
        content.append(System.lineSeparator());
        for (Model model : models) {
            Map<String, Double> scores = model.getScores();
            List<Exon> exons = model.getExons();
            content.append(String.format("%-32s", model.getAlignment().getViralProtein().getProteinID()));
            content.append(String.format("%-20s", model.getGeneSymbol()));
            content.append(String.format("%-20s", model.getDirection()));
            content.append(String.format("%-20s", model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform()));
            List<Range> NTranges = exons.stream().map(e -> e.getRange()).collect(Collectors.toList());
            List<Range> AAranges = exons.stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
                    .collect(Collectors.toList());
            for (int i = 0; i < NTranges.size(); i++) {
                String start = Long.toString(NTranges.get(i).getBegin(oneBased));
                String end = Long.toString(NTranges.get(i).getEnd(oneBased));
                if (i == 0 && model.isPartial5p()) start = "<" + start;
                if (i == NTranges.size() - 1 && model.isPartial3p()) end = ">" + end;
                content.append(String.format("%-20s", start + "-" + end));
                content.append(String.format("%-20s", AAranges.get(i).getBegin(oneBased) + "-" + AAranges.get(i).getEnd(oneBased)));
                content.append(System.lineSeparator());
                content.append(String.format("%-92s", ""));
            }
            content.append(String.format("%-40s", ""));
            for (String key : scores.keySet()) {
                content.append(String.format("%-30s", key + "=" + scores.get(key)));
                content.append(System.lineSeparator());
                content.append(String.format("%-132s", ""));
            }
            content.append(System.lineSeparator());
        }
        LOGGER.info(content);
    }

    public static void printAllGeneModelsWithScores ( List<Model> geneModels, String msg ) {

        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
                "*********************************" + msg + "*****************************************************");
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-10s%-10s%-10s%-10s%-10s%-10s%-20s%-10s%-10s%-20s", "Reference_ID", "%id", "%sim", "%cov", "%t5", "%gap", "%t3", "start..stop", "pep_size", "ref_size", "definition"));
        content.append(System.lineSeparator());
        for (Model model : geneModels) {
            ViralProtein viralProtein = model.getAlignment().getViralProtein();
            Map<String, Double> scores = model.getScores();
            long cdsBases = 0;
            content.append(String.format("%-32s", viralProtein.getProteinID()));
            content.append(String.format("%-10s", String.format("%.02f", scores.get("%identity"))));
            content.append(String.format("%-10s", String.format("%.02f", scores.get("%similarity"))));
            content.append(String.format("%-10s", String.format("%.02f", scores.get("%coverage"))));
            content.append(String.format("%-10s", "0.0"));
            content.append(String.format("%-10s", "0.0"));
            content.append(String.format("%-10s", "0.0"));
            for (int i = 0; i < model.getExons().size(); i++) {
                Exon exon = model.getExons().get(i);
                String start = Long.toString(exon.getRange().getBegin(oneBased));
                String end = Long.toString(exon.getRange().getEnd(oneBased));
                if (i == 0 && model.isPartial5p()) start = "<" + start;
                if (i == model.getExons().size() - 1 && model.isPartial3p()) end = ">" + end;
                if (i != 0) {
                    content.append(System.lineSeparator());
                    content.append(String.format("%-92s", ""));
                }
                content.append(String.format("%-20s", start + ".." + end));
                cdsBases = cdsBases + exon.getRange().getLength();
            }
            if (model.isPartial3p()) {
                content.append(String.format("%-10s", ( cdsBases ) / 3));
            } else {
                content.append(String.format("%-10s", ( cdsBases - 3 ) / 3));
            }
            content.append(String.format("%-10s", viralProtein.getSequence().getLength()));
            content.append(String.format("%-20s", model.getGeneSymbol() + " | " + viralProtein.getProduct()));
            content.append(System.lineSeparator());
        }
        LOGGER.info(content);
    }

    public static void printAlignments ( List<Alignment> alignments ) {

        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
                "*********************************Initial list of Alignments*****************************************************");
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-20s%-20s%-20s%-20s", "Protein_ID", "Alignment_Algorithm", "NTSeqRange", "AASeqRange", "Frame"));
        content.append(System.lineSeparator());
        for (Alignment alignment : alignments) {
            content.append(String.format("%-32s", alignment.getViralProtein().getProteinID()));
            content.append(String.format("%-20s", alignment.getAlignmentTool().getToolName()));
            List<Range> NTranges = alignment.getAlignmentFragments().stream().map(e -> e.getNucleotideSeqRange())
                    .collect(Collectors.toList());
            List<Frame> frames = alignment.getAlignmentFragments().stream().map(e -> e.getFrame()).collect(Collectors.toList());
            List<Range> AAranges = alignment.getAlignmentFragments().stream().map(e -> e.getProteinSeqRange())
                    .collect(Collectors.toList());
            for (int i = 0; i < NTranges.size(); i++) {
                if (i != 0) {
                    content.append(System.lineSeparator());
                    content.append(String.format("%-52s", ""));
                }
                content.append(String.format("%-20s", NTranges.get(i).getBegin(oneBased) + "-" + NTranges.get(i).getEnd(oneBased)));
                content.append(String.format("%-20s", AAranges.get(i).getBegin(oneBased) + "-" + AAranges.get(i).getEnd(oneBased)));
                content.append(String.format("%-20s", frames.get(i).getFrame()));
            }
            content.append(System.lineSeparator());
        }
        LOGGER.info(content);
    }

    public static void printSequenceFeatures ( List<Model> geneModels, String msg ) {

        double identityAvg = 0;
        double similarityAvg = 0;
        double coverageAvg = 0;
        long totalCDSBases = 0;
        long totalPepBases = 0;
        VirusGenome virusGenome = geneModels.get(0).getAlignment().getVirusGenome();
        String refDb = geneModels.get(0).getAlignment().getAlignmentEvidence().getReference_db();
        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(System.lineSeparator());
        content.append(String.format("%-20s%-10s%-10s%-10s%-10s%-10s%-10s%-20s%-10s%-10s%-20s%-20s", "gene_id", "%id", "%sim", "%cov", "%t5", "%gap", "%t3", "start..stop", "pep_size", "ref_size", "ref_id", "definition"));
        content.append(System.lineSeparator());
        for (Model model : geneModels) {
            ViralProtein viralProtein = model.getAlignment().getViralProtein();
            Map<String, Double> scores = model.getScores();
            long cdsBases = 0;
            content.append(String.format("%-20s", model.getGeneID()));
            content.append(String.format("%-10s", String.format("%.02f", scores.get("%identity"))));
            content.append(String.format("%-10s", String.format("%.02f", scores.get("%similarity"))));
            content.append(String.format("%-10s", String.format("%.02f", scores.get("%coverage"))));
            content.append(String.format("%-10s", "0.0"));
            content.append(String.format("%-10s", "0.0"));
            content.append(String.format("%-10s", "0.0"));
            for (int i = 0; i < model.getExons().size(); i++) {
                Exon exon = model.getExons().get(i);
                String start = Long.toString(exon.getRange().getBegin(oneBased));
                String end = Long.toString(exon.getRange().getEnd(oneBased));
                if (i == 0 && model.isPartial5p()) start = "<" + start;
                if (i == model.getExons().size() - 1 && model.isPartial3p()) end = ">" + end;
                if (i != 0) {
                    content.append(System.lineSeparator());
                    content.append(String.format("%-80s", ""));
                }
                content.append(String.format("%-20s", start + ".." + end));
                cdsBases = cdsBases + exon.getRange().getLength();
            }
            if (!model.isPartial3p()) cdsBases = cdsBases - 3;
            content.append(String.format("%-10s", ( cdsBases ) / 3));
            content.append(String.format("%-10s", viralProtein.getSequence().getLength()));
            content.append(String.format("%-20s", viralProtein.getProteinID()));
            content.append(String.format("%-20s", model.getGeneSymbol() + " | " + viralProtein.getProduct()));
            content.append(System.lineSeparator());
            totalCDSBases = totalCDSBases + cdsBases;
            totalPepBases = cdsBases + totalPepBases;
            identityAvg = identityAvg + scores.get("%identity");
            similarityAvg = similarityAvg + scores.get("%similarity");
            coverageAvg = coverageAvg + scores.get("%coverage");

            IDGenerator idGenerator = IDGenerator.of(model.getGeneID());
            for (MaturePeptideMatch match : model.getMaturePeptides()) {
                content.append(String.format("%-20s", idGenerator.next()));
                content.append(String.format("%-10s",String.format("%.02f",match.getIdentity() * 100)));
                content.append(String.format("%-10s",String.format("%.02f",match.getSimilarity() * 100)));
                content.append(String.format("%-10s",String.format("%.02f",match.getCoverage() * 100)));
                content.append(String.format("%-10s","0.0"));
                content.append(String.format("%-10s","0.0"));
                content.append(String.format("%-10s","0.0"));
                List<Range> cdRanges = VigorFunctionalUtils.proteinRangeToCDSRanges(model, match.getProteinRange());

                content.append(String.format("%-20s",GenerateVigorOutput.formatMaturePeptideRange(model,
                                                                                                  match,
                                                                                                  cdRanges,
                                                                                                  Range.CoordinateSystem.RESIDUE_BASED,
                                                                                                  "..",
                                                                                                  model.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED) + model.getExons().get(0).getFrame().getFrame() - 1 ,
                                                                                                  model.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED)
                )));
                content.append(String.format("%-10s",match.getProteinRange().getLength()));
                content.append(String.format("%-10s",match.getReference().getSequence().getLength()));
                content.append(String.format("%-20s",match.getReference().getProteinID()));
                content.append(String.format("%-20s",String.join("|", model.getGeneSymbol(), VigorUtils.putativeName(match.getReference().getProduct(), match.isFuzzyEnd(), match.isFuzzyBegin()))));
                content.append(System.lineSeparator());
            }
        }
        totalPepBases = totalPepBases / 3;
        identityAvg = identityAvg / geneModels.size();
        similarityAvg = similarityAvg / geneModels.size();
        coverageAvg = coverageAvg / geneModels.size();
        StringBuffer contentSummary = new StringBuffer("");
        contentSummary.append(System.lineSeparator());
        contentSummary.append(String.format("%-32s%-10s%-10s%-15s%-15s%-15s%-15s%-17s%-15s%-20s", "Sequence", "Length", "Genes", "Pseudogenes", "CDS_Bases", "Peptide_Bases", "%Ref_Identity", "%Ref_Similarity", "%Ref_Coverage", "Ref_DB"));
        contentSummary.append(System.lineSeparator());
        contentSummary.append(String.format("%-32s", virusGenome.getId()));
        contentSummary.append(String.format("%-10s", virusGenome.getSequence().getLength()));
        contentSummary.append(String.format("%-10s", geneModels.size()));
        contentSummary.append(String.format("%-15s", "0"));
        contentSummary.append(String.format("%-15s", totalCDSBases));
        contentSummary.append(String.format("%-15s", totalPepBases));
        contentSummary.append(String.format("%-15s", String.format("%.02f", identityAvg)));
        contentSummary.append(String.format("%-17s", String.format("%.02f", similarityAvg)));
        contentSummary.append(String.format("%-15s", String.format("%.02f", coverageAvg)));
        contentSummary.append(String.format("%-20s", refDb));
        contentSummary.append(System.lineSeparator());
        LOGGER.info(contentSummary);
        LOGGER.info(content);
    }
}
