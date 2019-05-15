package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.exception.VigorRuntimeException;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.function.Consumer;

public class GFF3Writer extends BaseOutputWriter {
    
    private static Logger LOGGER = LogManager.getLogger(GFF3Writer.class);
    private static Consumer<BufferedWriter> onOpen = b -> {
            try {
                b.write("##gff-version 3");
                b.newLine();
            } catch (IOException e) {
                throw new VigorRuntimeException("Problem writing GFF3 header");
            }
        };


    @Override
    public Consumer<BufferedWriter> getOnOpen() {
        return onOpen;
    }

    @Override
    public void writeModels(Outfiles outfiles, List<Model> models) throws VigorException, IOException {
            if (models.isEmpty()) {
            LOGGER.warn("no gene models for GFF3 output");
            return;
        }
        String genomeID = models.get(0).getAlignment().getVirusGenome().getId();

        OutputContext context = new OutputContext();
        context.addContext(OutputContext.Key.GENOME, genomeID);
        try (WriterBundle bw = getWriter(outfiles, context, OutputContext.Key.GENOME)) {

            for (Model geneModel : models) {
                List<String> notes = geneModel.getNotes();
                int i = 1;
                VirusGenome virusGenome = geneModel.getAlignment().getVirusGenome();
                long seqlength = virusGenome.getSequence().getLength();
                String geneomeSeqID = virusGenome.getId();
                List<Exon> exons = geneModel.getExons();
                String geneName = geneModel.getAlignment().getViralProtein().getGeneSymbol();
                long CDSStart = VigorFunctionalUtils.getDirectionBasedCoordinate
                        (exons.get(0).getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                long CDSEnd = VigorFunctionalUtils.getDirectionBasedCoordinate
                        (exons.get(exons.size() - 1).getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                String mRnaID = geneModel.getGeneID() + "." + i;
                bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
                //gene
                if (geneModel.isPseudogene()) {
                    bw.write("pseudogene" + "\t");
                } else {
                    bw.write("gene" + "\t");
                }
                bw.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
                bw.write("." + "\t");
                if (geneModel.getDirection().equals(Direction.FORWARD)) {
                    bw.write("+" + "\t");
                } else bw.write("-" + "\t");
                bw.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
                bw.write("ID=" + geneModel.getGeneID() + ";" + "Name=" + geneName + ";");
                if (geneModel.isPartial3p() || geneModel.isPartial5p()) {
                    bw.write("Partial" + ",");
                }
                if (notes.contains(NoteType.Gene.toString())) {
                    bw.write(String.format("Note=%s;", NoteType.Gene));
                }
                bw.write("\n");
                //mRNA
                bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
                bw.write("mRNA" + "\t");
                bw.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
                bw.write("." + "\t");
                if (geneModel.getDirection().equals(Direction.FORWARD)) {
                    bw.write("+" + "\t");
                } else bw.write("-" + "\t");
                bw.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
                bw.write("ID=" + mRnaID + ";" + "Parent=" + geneModel.getGeneID() + ";");
                if (geneModel.isPartial3p() || geneModel.isPartial5p()) {
                    bw.write("Partial;");
                }
                bw.write("\n");
                //exon
                IDGenerator idGenerator = new IDGenerator(mRnaID);
                for (int j = 0; j < exons.size(); j++) {
                    Exon exon = exons.get(j);
                    bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
                    bw.write("exon" + "\t");
                    long exonBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    long exonEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(exon.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    bw.write(Math.min(exonBegin, exonEnd) + "\t" + Math.max(exonBegin, exonEnd) + "\t");
                    bw.write("." + "\t");
                    bw.write((geneModel.getDirection() == Direction.FORWARD ? "+" : '-') + "\t");
                    bw.write(exon.getFrame().getFrame() - 1 + "\t");
                    bw.write(String.format("ID=%s;Parent=%s;", idGenerator.next(), mRnaID));
                    if ((j == 0 && geneModel.isPartial5p()) || (j == exons.size() - 1 && geneModel.isPartial3p())) {
                        bw.write("Partial" + ";");
                    }
                    bw.write("\n");
                }
                //CDS
                bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
                bw.write("CDS" + "\t");
                bw.write(Math.min(CDSStart, CDSEnd) + "\t" + Math.max(CDSStart, CDSEnd) + "\t");
                bw.write("." + "\t");
                if (geneModel.getDirection().equals(Direction.FORWARD)) {
                    bw.write("+" + "\t");
                } else bw.write("-" + "\t");
                bw.write(exons.get(0).getFrame().getFrame() - 1 + "\t");
                bw.write(String.format("ID=%s;Parent=%s;", idGenerator.next(), geneModel.getGeneID()));
                if (geneModel.isPartial3p() || geneModel.isPartial5p()) {
                    bw.write("Partial" + ";");
                }
                bw.write("\n");
                //insertion
                if (geneModel.getInsertRNAEditingRange() != null) {
                    bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
                    bw.write("insertion" + "\t");
                    long insertBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(geneModel.getInsertRNAEditingRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    long insertEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(geneModel.getInsertRNAEditingRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    bw.write(Math.min(insertBegin, insertEnd) + "\t" + Math.min(insertBegin, insertEnd) + "\t");
                    bw.write("." + "\t");
                    if (geneModel.getDirection().equals(Direction.FORWARD)) {
                        bw.write("+" + "\t");
                    } else bw.write("-" + "\t");
                    bw.write("." + "\t");
                    bw.write(String.format("ID=%s;Parent=%s;", idGenerator.next(), mRnaID));
                    if (notes.contains(NoteType.RNA_Editing.toString())) {
                        bw.write(String.format("Note=%s;", NoteType.RNA_Editing));
                    }
                    bw.write("\n");
                }
                //Stop_codon_read_through
                if (geneModel.getReplaceStopCodonRange() != null) {
                    bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
                    bw.write("stop_codon_read_through" + "\t");
                    long replaceStopBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(geneModel.getReplaceStopCodonRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    long replaceStopEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(geneModel.getReplaceStopCodonRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    bw.write(Math.min(replaceStopBegin, replaceStopEnd) + "\t" + Math.max(replaceStopBegin, replaceStopEnd) + "\t");
                    bw.write("." + "\t");
                    if (geneModel.getDirection().equals(Direction.FORWARD)) {
                        bw.write("+" + "\t");
                    } else bw.write("-" + "\t");
                    bw.write("." + "\t");
                    bw.write(String.format("ID=%s;Parent=%s;", idGenerator.next(), mRnaID));
                    if (notes.contains(NoteType.StopCodonReadThrough.toString()))
                        bw.write(String.format("Note=%s;", (NoteType.StopCodonReadThrough)));
                    bw.write("\n");
                }
                //plus/minus_1_translationally_frameshifted
                if (geneModel.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage().isHas_ribosomal_slippage() && geneModel.getRibosomalSlippageRange() != null) {
                    int frameshift = geneModel.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage().getSlippage_frameshift();
                    long slippageBegin = VigorFunctionalUtils.getDirectionBasedCoordinate(geneModel.getRibosomalSlippageRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    long slippageEnd = VigorFunctionalUtils.getDirectionBasedCoordinate(geneModel.getRibosomalSlippageRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqlength, geneModel.getDirection());
                    if (frameshift == -1 || frameshift == 1) {
                        bw.write(geneomeSeqID + "\t" + "vigor" + "\t");
                        bw.write(frameshift == -1 ? "mRNA_with_minus_1_frameshift" : "mRNA_with_plus_1_frameshift" + "\t");
                        bw.write(Math.min(slippageBegin, slippageEnd) + "\t" + Math.max(slippageBegin, slippageEnd) + "\t");
                        bw.write("." + "\t");
                        if (geneModel.getDirection().equals(Direction.FORWARD)) {
                            bw.write("+" + "\t");
                        } else bw.write("-" + "\t");
                        bw.write("." + "\t");
                        bw.write(String.format("ID=%s;Parent=%s;", idGenerator.next(), mRnaID));
                        bw.write("\n");
                    }
                }
            }
        }
    }

    @Override
    public String getExtension() {
        return "gff3";
    }
}
