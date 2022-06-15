package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;

import java.io.IOException;
import java.util.*;



public class TBLWriter extends BaseOutputWriter {
    private static Logger LOGGER = LogManager.getLogger(TBLWriter.class);
    private String locusPrefix = "";

    @Override
    public void configure(VigorConfiguration config) {
        super.configure(config);
        locusPrefix = config.getOrDefault(ConfigurationParameters.Locustag,"");
    }

    @Override
    public String getExtension() {
        return "tbl";
    }

    @Override
    public void writeModels(Outfiles outfiles, List<Model> models) throws VigorException, IOException {
        if (models.isEmpty()) {
            LOGGER.warn("no gene models to write to file");
            return;
        }

        VirusGenome virusGenome = models.get(0).getAlignment().getVirusGenome();
        String genomeID = virusGenome.getId();
        OutputContext context = new OutputContext().addContext(OutputContext.Key.GENOME, models.get(0).getGeneID());

        boolean writeLocus = ! NullUtil.isNullOrEmpty(locusPrefix);
        long seqlength = virusGenome.getSequence().getLength();

        try (WriterBundle bw = getWriter(outfiles, context, OutputContext.Key.GENOME)) {

            bw.write(">Features " + genomeID + "\n");
            String proteinID = "";
            List<Map> previousExonDirections = new ArrayList<>(0);
            for (int i = 0; i < models.size(); i++) {
                Model model = models.get(i);
                List<String> modelNotes = model.getNotes();
                Ribosomal_Slippage riboSlippage = model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
                RNA_Editing rna_editing = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing();
                List<SpliceSite> spliceSites = model.getAlignment().getViralProtein().getGeneAttributes().getSpliceSites();
                StringBuilder notes = new StringBuilder();
                List<Exon> exons = model.getExons();
                Collections.sort(exons, Comparator.comparing(e -> VigorFunctionalUtils.getDirectionBasedRange(e.getRange(), seqlength, model.getDirection()), Range.Comparators.ARRIVAL));
                Exon firstExon = exons.get(0);
                int codon_start = firstExon.getFrame().getFrame();
                if (!model.getAlignment().getViralProtein().getProteinID().equals(proteinID) || !isExonDirectionsIdentical(previousExonDirections, exons, seqlength, model.getDirection())) {
                    bw.write(OutputWriterUtils.getGeneCoordinatesString(model, models));
                    bw.write("\tgene\n");
                    if (writeLocus) {
                        bw.write("\t\t\tlocus_tag\t" + VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
                    }
                    bw.write("\t\t\tgene\t" + model.getGeneSymbol() + "\n");
                }
                proteinID = model.getAlignment().getViralProtein().getProteinID();
                previousExonDirections.clear();
                String geneSynonym = model.getAlignment().getViralProtein().getGeneSynonym();
                if (!NullUtil.isNullOrEmpty(geneSynonym)) {
                    bw.write("\t\t\tgene_syn\t" + geneSynonym + "\n");
                }
                for (int j = 0; j < exons.size(); j++) {
                    Exon exon = exons.get(j);
                    String Cstart = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqlength, model.getDirection()));
                    String Cend = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqlength, model.getDirection()));
                    Map<String, String> exonDirections = new HashMap<>(2);
                    exonDirections.put("start", Cstart);
                    exonDirections.put("end", Cend);
                    previousExonDirections.add(exonDirections);
                    if (j == 0 && model.isPartial5p()) {
                        Cstart = "<" + Cstart;
                    }
                    if (j == exons.size() - 1 && model.isPartial3p()) {
                        Cend = ">" + Cend;
                    }
                    if (j == 0) {
                        bw.write(String.join("\t",
                                             Cstart,
                                             Cend,
                                             model.isPseudogene() ? "misc_feature" : "CDS"));
                        bw.newLine();
                    } else {
                        bw.write(Cstart + "\t" + Cend + "\n");
                    }
                }
                bw.write("\t\t\tcodon_start\t" + codon_start + "\n");
                if (model.getReplaceStopCodonRange() != null) {
                    long replaceStopBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getReplaceStopCodonRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED),
                            seqlength,
                            model.getDirection());
                    long replaceStopEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getReplaceStopCodonRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED),
                            seqlength,
                            model.getDirection());
                    bw.write("\t\t\ttransl_except\t" + String.format("(pos:%s..%s,aa:R)", replaceStopBegin, replaceStopEnd) + "\n");
                }
                bw.write("\t\t\tprotein_id\t" + model.getGeneID() + "\n");
                if (writeLocus) {
                    bw.write("\t\t\tlocus_tag\t" + VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
                }
                bw.write("\t\t\tgene\t" + model.getGeneSymbol() + "\n");
                String product = model.getAlignment().getViralProtein().getProduct();
                if (!NullUtil.isNullOrEmpty(product)) {
                    bw.write("\t\t\tproduct\t" + VigorUtils.putativeName(product, model.isPartial3p(), model.isPartial5p()) + "\n");
                } else {
                    LOGGER.warn("Missing product for {}", genomeID);
                }
                if (riboSlippage.isHas_ribosomal_slippage()) {
                    bw.write("\t\t\tribosomal_slippage\n");
                }
                if (rna_editing.isHas_RNA_editing()) {
                    bw.write("\t\t\texception\tRNA editing\n");
                    notes = notes.append(rna_editing.getNote() + ";");
                }
                if (spliceSites != SpliceSite.DEFAULT_SPLICE_SITES) {
                    notes.append("non-canonical splicing");
                }
                if (model.getInsertRNAEditingRange() != null) {
                    bw.write("\t\t\tnote\t" + notes + "\n");
                    long insertBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getInsertRNAEditingRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED),
                            seqlength,
                            model.getDirection());
                    long insertEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(
                            model.getInsertRNAEditingRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED),
                            seqlength,
                            model.getDirection());
                    // TODO coordinate system?
                    bw.write(insertBegin + "\t" + insertEnd + "\t" + "misc_feature\n");
                    NucleotideSequence subSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(model.getInsertRNAEditingRange()).build();
                    bw.write("\t\t\tnote\tlocation of RNA editing (" + subSeq + "," + rna_editing.getInsertionString() + ") in " + model.getAlignment().getViralProtein().getProduct() + "\n");
                }
                if (modelNotes.size() > 0) {
                    String notesText = String.join(",", modelNotes);
                    bw.write("\t\t\tnote\t" + notesText + "\n");
                }
            }

            for (Model model : models) {
                if (model.getMaturePeptides() != null && !model.getMaturePeptides().isEmpty()) {
                    bw.write(">Features " + model.getGeneID());
                    long proteinLength = model.getAlignment().getViralProtein().getSequence().getLength();
                    bw.newLine();
                    String product;
                    for (MaturePeptideMatch match : model.getMaturePeptides()) {
                        long start = VigorFunctionalUtils.getDirectionBasedCoordinate(1, proteinLength, model.getDirection());
                        long end = VigorFunctionalUtils.getDirectionBasedCoordinate(proteinLength, proteinLength, model.getDirection());
                        bw.write(OutputWriterUtils.formatMaturePeptideRange(model,
                                                                            match,
                                                                            Arrays.asList(match.getProteinRange()),
                                                                            Range.CoordinateSystem.RESIDUE_BASED,
                                                                            "\t",
                                                                            start,
                                                                            end, false));
                        bw.write("\t");
                        product = match.getReference().getProduct();
                        if (!NullUtil.isNullOrEmpty(product)) {
                            if (product.contains("signal")) {
                                // TODO pre-classify type
                                bw.write("sig_peptide");
                            } else {
                                bw.write("mat_peptide");
                                bw.newLine();
                                bw.write("\t\t\tproduct\t");
                                // TODO check that there aren't other factors here.
                                bw.write(VigorUtils.putativeName(product, match.isFuzzyEnd(), match.isFuzzyBegin()));
                            }
                        } else {
                            LOGGER.warn("Missing product for mature peptide {}", match.getReference().getProteinID());
                        }
                        bw.newLine();
                        String geneSymbol = model.getGeneSymbol();
                        if (!NullUtil.isNullOrEmpty(geneSymbol)) {
                            if (writeLocus) {
                                bw.write("\t\t\tlocus_tag\t");
                                bw.write(VigorUtils.nameToLocus(geneSymbol, locusPrefix, model.isPseudogene()));
                                bw.newLine();
                            }
                            bw.write("\t\t\tgene\t");
                            bw.write(geneSymbol);
                            bw.newLine();
                        }
                    }
                }
            }

        } catch (IOException e) {
            throw new VigorException(e);
        }
    }

    private boolean isExonDirectionsIdentical(List<Map> previousExonDirections, List<Exon> exons, long seqlength, Direction direction) {
        if (previousExonDirections.size() != exons.size()) {
            return false;
        }

        for (int i = 0; i < exons.size(); i++) {
            Exon exon = exons.get(i);
            Map previousExonDirection = previousExonDirections.get(i);

            String start = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqlength, direction));
            String end = Long.toString(VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqlength, direction));

            if (!(start.equals(previousExonDirection.get("start")) && start.equals(previousExonDirection.get("end")))) {
                return false;
            }
        }

        return true;
    }
}
