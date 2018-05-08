package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.*;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.IOException;
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

    public static class Outfiles extends EnumMap<Outfile, BufferedWriter> {

        public Outfiles() {
            super(Outfile.class);
        }
    }

    private static final Logger LOGGER = LogManager
            .getLogger(GenerateVigorOutput.class);

    public void generateOutputFiles(Outfiles outfiles ,List<Model> geneModels) throws IOException {
        generateTBLReport(outfiles.get(Outfile.TBL),geneModels);
        generateCDSReport(outfiles.get(Outfile.CDS),geneModels);
        generatePEPReport(outfiles.get(Outfile.PEP),geneModels);
    }


    public void generateTBLReport(BufferedWriter bw, List<Model> geneModels) throws IOException {
        if (geneModels.isEmpty()) {
            LOGGER.warn("no gene models to write to file");
            return;
        }
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
            long start = exons.get(0).getRange().getBegin() + 1;
            long end = exons.get(exons.size() - 1).getRange().getEnd() + 1;
            bw.write(String.format("%-13s", Long.toString(start) + "\t" + Long.toString(end)) + String.format("%-5s", "gene") + "\n");
            bw.write("\t\t\t\t\t" + String.format("%-13s", "locus_tag") + "\tvigor_" + model.getGeneSymbol() + "\n");
            bw.write("\t\t\t\t\t" + String.format("%-13s", "gene") + "\t" + model.getAlignment().getViralProtein().getGeneSymbol() + "\n");
            for (int j = 0; j < exons.size(); j++) {
                Exon exon = exons.get(j);
                if (j == 0) {
                    bw.write(String.format("%-13s", Long.toString(exon.getRange().getBegin() + 1) + "\t" + Long.toString(exon.getRange().getEnd() + 1)) + String.format("%-5s", "CDS") + "\n");
                } else {
                    bw.write(Long.toString(exon.getRange().getBegin() + 1) + "\t" + Long.toString(exon.getRange().getEnd() + 1) + "\n");
                }
            }
            bw.write("\t\t\t\t\t" + String.format("%-13s", "codon_start") + "\t" + start + "\n");
            bw.write("\t\t\t\t\t" + String.format("%-13s", "protein_id") + "\t" + model.getGeneID() + "\n");
            bw.write("\t\t\t\t\t" + String.format("%-13s", "locus_tag") + "\tvigor_" + model.getGeneSymbol() + "\n");
            bw.write("\t\t\t\t\t" + String.format("%-13s", "gene") + "\t" + model.getAlignment().getViralProtein().getGeneSymbol() + "\n");
            bw.write("\t\t\t\t\t" + String.format("%-13s", "product") + "\t" + model.getAlignment().getViralProtein().getProduct() + "\n");
            if (riboSlippage.isHas_ribosomal_slippage()) {
                bw.write("\t\t\t\t\t" + String.format("%-13s", "ribosomal_slippage") + "\n");
            }
            if (rna_editing.isHas_RNA_editing()) {
                bw.write("\t\t\t\t\t" + String.format("%-13s", "exception") + "\tRNA editing\n");
                notes = notes.append(rna_editing.getNote() + ";");
            }
            if (splicing.getNonCanonical_spliceSites() != null && splicing.getNonCanonical_spliceSites().size() > 1) {
                notes.append("non-canonical splicing");
            }

            if (model.getInsertRNAEditingRange() != null) {
                bw.write("\t\t\t\t\t" + String.format("%-13s", "note") + "\t" + notes + "\n");
                bw.write(String.format("%-13s", model.getInsertRNAEditingRange().getBegin()) + "\t" + model.getInsertRNAEditingRange().getEnd() + "\t" + "misc_feature\n");


                if (model.getInsertRNAEditingRange() != null) {
                    NucleotideSequence subSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(model.getInsertRNAEditingRange()).build();
                    bw.write("\t\t\t\t\t" + String.format("%-13s", "note") + "\tlocation of RNA editing (" + subSeq + "," + rna_editing.getInsertionString() + ") in " + model.getAlignment().getViralProtein().getProduct() + "\n");
                }

                if (!model.getMaturePeptides().isEmpty()) {
                    bw.write(">Features " + idGenerator.next());
                    bw.newLine();
                    String product;
                    boolean isSignalPeptide;
                    for (MaturePeptideMatch match : model.getMaturePeptides()) {

                        bw.write(String.format("%s\t%s\t", formatMaturePeptideRange(match)));
                        product = match.getReference().getProduct();
                        isSignalPeptide = product.contains("signal");
                        if (isSignalPeptide) {
                            // TODO pre-classify type
                            bw.write("sig_peptide");
                        } else {
                            bw.write("mat_peptide");
                            bw.newLine();
                            bw.write("\t\t\tproduct\t");
                            bw.write(product);
                        }
                        bw.newLine();
                        String geneSymbol = model.getGeneSymbol();
                        if (!(geneSymbol == null || geneSymbol.isEmpty())) {
                            bw.write("\t\t\tgene\t");
                            bw.write(geneSymbol);
                            bw.newLine();
                        }
                        // TODO other attributes (locus_tag, note)
                    }
                }
            }

        }
    }

    public void generateCDSReport(BufferedWriter bw ,List<Model> geneModels) throws IOException {
        for (Model model: geneModels) {
            String reference_db = model.getAlignment().getAlignmentEvidence().getReference_db();
            ViralProtein refProtein = model.getAlignment().getViralProtein();
            List<Exon> exons = model.getExons();
            Range cdsRange = Range.of(exons.get(0).getRange().getBegin()+1, exons.get(exons.size() - 1).getRange().getEnd()+1);
            StringBuilder defline = new StringBuilder("");
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
            //  bw.newLine();
            bw.write(defline.toString());
            bw.newLine();

            Iterator<String> sequenceIter = SequenceUtils.steamOf(model.getCds(),70).iterator();
            while(sequenceIter.hasNext() ){
                bw.write(sequenceIter.next());
                bw.newLine();
            }
            bw.newLine();
        }
    }



    public void generatePEPReport(BufferedWriter bw, List<Model> geneModels) throws IOException {

        for (Model model: geneModels) {
            String reference_db = model.getAlignment().getAlignmentEvidence().getReference_db();
            ViralProtein refProtein = model.getAlignment().getViralProtein();
            List<Exon> exons = model.getExons();
            Range cdsRange = Range.of(exons.get(0).getRange().getBegin() + 1, exons.get(exons.size() - 1).getRange().getEnd() + 1);
            StringBuilder defline = new StringBuilder();
            defline.append(">" + model.getGeneID());

            if (model.isPseudogene()) {
                defline.append(" pseudogene");
            }
            defline.append(" location=" + cdsRange);
            defline.append(" codon_start=" + cdsRange.getBegin());
            defline.append(" gene=" + refProtein.getGeneSymbol());
            defline.append(" product=" + refProtein.getProduct());
            defline.append(" ref_db=\"" + reference_db + "\"");
            defline.append(" ref_id=\"" + refProtein.getProteinID() + "\"");
            bw.write(defline.toString());
            bw.newLine();
            Iterator<String> sequenceLineIter = SequenceUtils.steamOf(model.getTanslatedSeq(), 70).iterator();
            while (sequenceLineIter.hasNext()) {
                bw.write(sequenceLineIter.next());
                bw.newLine();
            }
            bw.newLine();

            IDGenerator idGenerator = IDGenerator.of(model.getGeneID());
            ProteinSequence matchSequence;

            if (model.getMaturePeptides()!=null && !model.getMaturePeptides().isEmpty()){
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


                    matchSequence = match.getProtein().toBuilder().trim(match.getProteinRange()).build();
                    Iterator<String> lineIter = SequenceUtils.steamOf(matchSequence, 70).iterator();
                    while (lineIter.hasNext()) {
                        bw.write(lineIter.next());
                        bw.newLine();
                    }
                }
        }
        }



    }

    private static Object[] formatMaturePeptideRange(MaturePeptideMatch match) {
        String start = String.valueOf(match.getProteinRange().getBegin());
        if (match.isFuzzyBegin()) {
            start = "<" + start;
        }
        String end = String.valueOf(match.getProteinRange().getEnd());
        if (match.isFuzzyEnd()) {
            end = ">" + end;
        }
        return new Object[] {start, end};
    }
}
