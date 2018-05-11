package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.*;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

@Service
public class GenerateVigorOutput {

    public enum Outfile {
        TBL("tbl"),
        CDS("cds"),
        GFF3("gff3"),
        PEP("pep");

        final public String extension;
        Outfile(String extension) {
            this.extension = extension;
        }
    }

    public static class Outfiles extends EnumMap<Outfile, BufferedWriter> implements AutoCloseable {

        public Outfiles() {
            super(Outfile.class);
        }

        public void close() throws IOException {
            List<IOException> exceptions = new ArrayList<>();
            for (BufferedWriter writer: values()) {
                try {
                    writer.close();
                } catch (IOException e) {
                    exceptions.add(e);
                }

            }

            if (! exceptions.isEmpty()) {
                // TODO join all exceptions somehow meaningfully
                throw exceptions.get(0);
            }
        }

        public void flush() throws IOException{
            for (BufferedWriter writer: values()) {
                writer.flush();
            }
        }
    }

    private static final Logger LOGGER = LogManager.getLogger(GenerateVigorOutput.class);

    public void generateOutputFiles(VigorConfiguration config, Outfiles outfiles ,List<Model> geneModels) throws IOException {
        generateTBLReport(config, outfiles.get(Outfile.TBL),geneModels);
        generateCDSReport(config, outfiles.get(Outfile.CDS),geneModels);
        generatePEPReport(config, outfiles.get(Outfile.PEP),geneModels);
    }


    public void generateTBLReport(VigorConfiguration config, BufferedWriter bw, List<Model> geneModels) throws IOException {
        if (geneModels.isEmpty()) {
            LOGGER.warn("no gene models to write to file");
            return;
        }

        String locusPrefix = config.get(ConfigurationParameters.Locustag);
        boolean writeLocus = !(locusPrefix == null || locusPrefix.isEmpty());

        String genomeID = geneModels.get(0).getAlignment().getVirusGenome().getId();
        String[] genomeIDParts = genomeID.split(Pattern.quote("|"));
        String proteinIDOfGenome;
        if (genomeIDParts.length >= 2) {
            proteinIDOfGenome = genomeIDParts[0] + "_" + genomeIDParts[1];
        } else {
            proteinIDOfGenome = genomeIDParts[0];
        }

        IDGenerator idGenerator = IDGenerator.of(proteinIDOfGenome);
        bw.write(">Features " + genomeID + "\n");
        for (int i = 0; i < geneModels.size(); i++) {
            Model model = geneModels.get(i);
            Ribosomal_Slippage riboSlippage = model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
            RNA_Editing rna_editing = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing();
            Splicing splicing = model.getAlignment().getViralProtein().getGeneAttributes().getSplicing();
            StringBuilder notes = new StringBuilder("");
            List<Exon> exons = model.getExons();
            long start = exons.get(0).getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED);
            long end = exons.get(exons.size() - 1).getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED);
            bw.write(Long.toString(start));
            bw.write("\t");
            bw.write(Long.toString(end));
            bw.write("\tgene\n");
            if (writeLocus) {
                bw.write("\t\t\tlocus_tag\t" + VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
            }
            bw.write("\t\t\tgene\t" + model.getAlignment().getViralProtein().getGeneSymbol() + "\n");
            for (int j = 0; j < exons.size(); j++) {
                Exon exon = exons.get(j);
                if (j == 0) {
                    bw.write(String.join("\t",
                            Long.toString(exon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED)),
                            Long.toString(exon.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED)),"CDS"));
                    bw.newLine();
                } else {
                    bw.write(Long.toString(exon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED)) + "\t" + Long.toString(exon.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED)) + "\n");
                }
            }
            bw.write("\t\t\tcodon_start\t" + start + "\n");
            bw.write("\t\t\tprotein_id\t" + model.getGeneID() + "\n");
            if (writeLocus) {
                bw.write("\t\t\tlocus_tag\t" +  VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
            }
            bw.write("\t\t\tgene\t" + model.getAlignment().getViralProtein().getGeneSymbol() + "\n");
            bw.write("\t\t\tproduct\t" + VigorUtils.putativeName(model.getAlignment().getViralProtein().getProduct(), model.isPartial3p(), model.isPartial5p()) + "\n");
            if (riboSlippage.isHas_ribosomal_slippage()) {
                bw.write("\t\t\tribosomal_slippage\n");
            }
            if (rna_editing.isHas_RNA_editing()) {
                bw.write("\t\t\texception\tRNA editing\n");
                notes = notes.append(rna_editing.getNote() + ";");
            }
            if (splicing.getNonCanonical_spliceSites() != null && splicing.getNonCanonical_spliceSites().size() > 1) {
                notes.append("non-canonical splicing");
            }

            if (model.getInsertRNAEditingRange() != null) {
                bw.write("\t\t\tnote\t" + notes + "\n");
                // TODO coordinate system?
                bw.write(model.getInsertRNAEditingRange().getBegin() + "\t" + model.getInsertRNAEditingRange().getEnd() + "\t" + "misc_feature\n");
                NucleotideSequence subSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(model.getInsertRNAEditingRange()).build();
                bw.write("\t\t\tnote\tlocation of RNA editing (" + subSeq + "," + rna_editing.getInsertionString() + ") in " + model.getAlignment().getViralProtein().getProduct() + "\n");
            }

            if (model.getMaturePeptides()!=null && !model.getMaturePeptides().isEmpty()) {
                bw.write(">Features " + idGenerator.next());
                bw.newLine();
                String product;
                for (MaturePeptideMatch match : model.getMaturePeptides()) {

                    bw.write(String.format("%s\t%s\t", formatMaturePeptideRange(match)));
                    product = match.getReference().getProduct();
                    if (product.contains("signal")) {
                        // TODO pre-classify type
                        bw.write("sig_peptide");
                    } else {
                        bw.write("mat_peptide");
                        bw.newLine();
                        bw.write("\t\t\tproduct\t");
                        // TODO check that there aren't other factors here.
                        bw.write(VigorUtils.putativeName(product, model.isPartial3p(), model.isPartial5p()));
                    }
                    bw.newLine();
                    String geneSymbol = model.getGeneSymbol();
                    if (!(geneSymbol == null || geneSymbol.isEmpty())) {
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
    }

    private void writeDefline(BufferedWriter bw, Model model) throws IOException {
        String reference_db = model.getAlignment().getAlignmentEvidence().getReference_db();
        ViralProtein refProtein = model.getAlignment().getViralProtein();
        List<Exon> exons = model.getExons();
        Range cdsRange = Range.of(exons.get(0).getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED),
                exons.get(exons.size() - 1).getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED));
        StringBuilder defline = new StringBuilder();
        defline.append(">" + model.getGeneID());

        if (model.isPseudogene()) {
            defline.append(" pseudogene");
        }
        defline.append(" location=" + cdsRange);
        defline.append(" codon_start=" + cdsRange.getBegin());
        defline.append(" gene=" + refProtein.getGeneSymbol());
        defline.append(" product=" + refProtein.getProduct());
        defline.append(" ref_db=\"" + reference_db+"\"");
        defline.append(" ref_id=\"" + refProtein.getProteinID()+"\"");
        bw.write(defline.toString());
        bw.newLine();

    }

    private void writeSequence(BufferedWriter bw, Sequence seq) throws IOException {
        Iterator<String> sequenceLineIter = SequenceUtils.steamOf(seq, 70).iterator();
        while(sequenceLineIter.hasNext()) {
            bw.write(sequenceLineIter.next());
            bw.newLine();
        }
        bw.newLine();

    }

    public void generateCDSReport(VigorConfiguration config, BufferedWriter bw ,List<Model> geneModels) throws IOException {
        for (Model model: geneModels) {
            writeDefline(bw, model);
            writeSequence(bw, model.getTanslatedSeq());
        }
    }



    public void generatePEPReport(VigorConfiguration config, BufferedWriter bw, List<Model> geneModels) throws IOException {

        StringBuilder defline;
        for (Model model : geneModels) {

            writeDefline(bw, model);
            writeSequence(bw, model.getTanslatedSeq());

            IDGenerator idGenerator = IDGenerator.of(model.getGeneID());

            if (model.getMaturePeptides() != null && !model.getMaturePeptides().isEmpty()) {
                for (MaturePeptideMatch match : model.getMaturePeptides()) {
                    defline = new StringBuilder();
                    defline.append(">" + idGenerator.next());
                    if (model.isPseudogene()) {
                        defline.append(" pseudogene");
                    }
                    defline.append(" mat_peptide");
                    defline.append(String.format(" location=%s..%s", formatMaturePeptideRange(match)));
                    // TODO make sure this is correct
                    defline.append(String.format(" gene=\"%s\"", model.getGeneSymbol()));
                    // TODO this needs some formatting
                    defline.append(String.format(" product=\"%s\"", match.getReference().getProduct()));
                    defline.append(String.format(" ref_db=\"%s\"", model.getAlignment().getAlignmentEvidence().getMatpep_db()));
                    defline.append(String.format(" ref_id=\"%s\"", match.getReference().getProteinID()));
                    bw.write(defline.toString());
                    bw.newLine();
                    writeSequence(bw, match.getProtein().toBuilder().trim(match.getProteinRange()).build());
                }

            }
        }
    }

    private static Object[] formatMaturePeptideRange(MaturePeptideMatch match) {
        String start = String.valueOf(match.getProteinRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED));
        if (match.isFuzzyBegin()) {
            start = "<" + start;
        }
        String end = String.valueOf(match.getProteinRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED));
        if (match.isFuzzyEnd()) {
            end = ">" + end;
        }
        return new Object[] {start, end};
    }
}
