package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.*;
import org.springframework.stereotype.Service;

import java.io.*;
import java.nio.file.Paths;
import java.util.*;

@Service
public class GenerateVigorOutput {

    public enum Outfile {
        TBL("tbl"),
        CDS("cds"),
        GFF3("gff3"),
        PEP("pep"),
        ALN("aln");

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
    private static Range.CoordinateSystem oneBased = Range.CoordinateSystem.RESIDUE_BASED;

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
        bw.write(">Features " + genomeID + "\n");
        String proteinID ="";
        for (int i = 0; i < geneModels.size(); i++) {
            Model model = geneModels.get(i);
            List<NoteType> modelNotes = model.getNotes();
            Ribosomal_Slippage riboSlippage = model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
            RNA_Editing rna_editing = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing();
            Splicing splicing = model.getAlignment().getViralProtein().getGeneAttributes().getSplicing();
            StringBuilder notes = new StringBuilder("");
            List<Exon> exons = model.getExons();
            Exon firstExon = exons.get(0);
            String start = Long.toString(firstExon.getRange().getBegin(oneBased));
            int codon_start = firstExon.getFrame().getFrame();
            String end = Long.toString(exons.get(exons.size() - 1).getRange().getEnd(oneBased));
            if (model.isPartial5p()) start = "<" + start;
            if (model.isPartial3p()) end = ">" + end;
            if (!model.getAlignment().getViralProtein().getProteinID().equals(proteinID)) {
                /* bw.write(start);
                bw.write("\t");
                bw.write(end);*/
               bw.write(getGeneCoordinates(model,geneModels));
                bw.write("\tgene\n");
                if (writeLocus) {
                    bw.write("\t\t\tlocus_tag\t" + VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
                }
                bw.write("\t\t\tgene\t" + model.getGeneSymbol() + "\n");

            }
            proteinID = model.getAlignment().getViralProtein().getProteinID();
            String geneSynonym = model.getAlignment().getViralProtein().getGeneSynonym();
            if (geneSynonym != null && geneSynonym != "") {
                bw.write("\t\t\tgene_syn\t" + geneSynonym + "\n");
            }
            for (int j = 0; j < exons.size(); j++) {
                Exon exon = exons.get(j);
                String Cstart = Long.toString(exon.getRange().getBegin(oneBased));
                String Cend = Long.toString(exon.getRange().getEnd(oneBased));
                if (j == 0 && model.isPartial5p()) Cstart = "<" + Cstart;
                if (j == exons.size() - 1 && model.isPartial3p()) Cend = ">" + Cend;
                if (j == 0) {
                    bw.write(String.join("\t",
                            Cstart,
                            Cend, "CDS"));
                    bw.newLine();
                } else {
                    bw.write(Cstart + "\t" + Cend + "\n");
                }
            }
            bw.write("\t\t\tcodon_start\t" + codon_start + "\n");
            if (model.getReplaceStopCodonRange() != null) {
                bw.write("\t\t\ttransl_except\t" + String.format("(pos:%s..%s,aa:R)", model.getReplaceStopCodonRange().getBegin(oneBased), model.getReplaceStopCodonRange().getEnd(oneBased)) + "\n");
            }
            bw.write("\t\t\tprotein_id\t" + model.getGeneID() + "\n");
            if (writeLocus) {
                bw.write("\t\t\tlocus_tag\t" + VigorUtils.nameToLocus(model.getGeneSymbol(), locusPrefix, model.isPseudogene()) + "\n");
            }
            bw.write("\t\t\tgene\t" + model.getGeneSymbol() + "\n");
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
                bw.write(model.getInsertRNAEditingRange().getBegin(oneBased) + "\t" + model.getInsertRNAEditingRange().getEnd(oneBased) + "\t" + "misc_feature\n");
                NucleotideSequence subSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(model.getInsertRNAEditingRange()).build();
                bw.write("\t\t\tnote\tlocation of RNA editing (" + subSeq + "," + rna_editing.getInsertionString() + ") in " + model.getAlignment().getViralProtein().getProduct() + "\n");
            }
            if (modelNotes.contains(NoteType.Sequence_Gap)) {
                bw.write("\t\t\tnote\t" + NoteType.Sequence_Gap + "\n");
            }
        }
       for(Model model : geneModels) {
              if (model.getMaturePeptides() != null && !model.getMaturePeptides().isEmpty()) {
                    bw.write(">Features " + model.getGeneID());
                    bw.newLine();
                    String product;
                    for (MaturePeptideMatch match : model.getMaturePeptides()) {

                        bw.write(formatMaturePeptideRange(model,
                                match,
                                Arrays.asList(match.getProteinRange()),
                                Range.CoordinateSystem.RESIDUE_BASED,
                                "\t",
                                1,
                                match.getProtein().getLength()));
                        bw.write("\t");
                        product = match.getReference().getProduct();
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

    private String getGeneCoordinates(Model model, List<Model> geneModels){
        long start=model.getExons().get(0).getRange().getBegin()+1;
        long end =0;
        String startString="";
        String endString ="";
        Model endGeneModel=model;
        for(Model geneModel : geneModels){
            if(model.getAlignment().getViralProtein().getProteinID().equals(geneModel.getAlignment().getViralProtein().getProteinID())){
                endGeneModel=geneModel;
            }
        }
        end = endGeneModel.getExons().get(endGeneModel.getExons().size()-1).getRange().getEnd()+1;
        if(model.isPartial5p()){
            startString="<"+start;
        }else startString=""+start;
        if(endGeneModel.isPartial3p()){
            endString=">"+end;
        }else endString=""+end;

        return startString+"\t"+endString;
    }

    private void writeDefline(BufferedWriter bw, Model model) throws IOException {

        ViralProtein refProtein = model.getAlignment().getViralProtein();
        StringBuilder defline = new StringBuilder();
        defline.append(">" + model.getGeneID());
        if (model.isPseudogene()) defline.append(" pseudogene");
        int codon_start = model.getExons().get(0).getFrame().getFrame();
        defline.append(" location=");
        for(int i=0;i<model.getExons().size();i++) {
            Exon exon = model.getExons().get(i);
            String start = Long.toString(exon.getRange().getBegin(oneBased));
            String end = Long.toString(exon.getRange().getEnd(oneBased));
            if (model.isPartial5p() && i==0) {
                start = "<" + start;
            }
            if (model.isPartial3p() && i==model.getExons().size()-1) {
                end = ">" + end;
            }
            if(i!=0) defline.append(",");
            defline.append(String.format("%s..%s",start,end));
        }
        defline.append(" codon_start=" + codon_start);
        defline.append(String.format(" gene=\"%s\"" , refProtein.getGeneSymbol()));
        defline.append(String.format(" product=\"%s\"" , VigorUtils.putativeName(refProtein.getProduct(), model.isPartial3p(), model.isPartial5p())));
        String reference_db = model.getAlignment().getAlignmentEvidence().getReference_db();
        if (! (reference_db == null || reference_db.isEmpty())) {
            defline.append(String.format(" ref_db=\"%s\"", Paths.get(reference_db).getFileName().toString()));
        }
        defline.append(String.format(" ref_id=\"%s\"" , refProtein.getProteinID()));
        bw.write(defline.toString());
        bw.newLine();

    }

    private void writeSequence(BufferedWriter bw, Sequence seq) throws IOException {
        Iterator<String> sequenceLineIter = SequenceUtils.steamOf(seq, 60).iterator();
        while (sequenceLineIter.hasNext()) {
            bw.write(sequenceLineIter.next());
            bw.newLine();
        }
    }

    public void generateCDSReport(VigorConfiguration config, BufferedWriter bw ,List<Model> geneModels) throws IOException {
        for (Model model: geneModels) {
            writeDefline(bw, model);
            writeSequence(bw, model.getTranslatedSeq());
        }
    }



    public void generatePEPReport(VigorConfiguration config, BufferedWriter bw, List<Model> geneModels) throws IOException {

        StringBuilder defline;
        for (Model model : geneModels) {

            writeDefline(bw, model);
            writeSequence(bw, model.getTranslatedSeq());

            IDGenerator idGenerator = IDGenerator.of(model.getGeneID());

            for (MaturePeptideMatch match : model.getMaturePeptides()) {
                defline = new StringBuilder();
                defline.append(">" + idGenerator.next());
                if (model.isPseudogene()) {
                    defline.append(" pseudogene");
                }
                defline.append(" mat_peptide");
                List<Range> cdsRanges = VigorFunctionalUtils.proteinRangeToCDSRanges(model, match.getProteinRange());
                // TODO handle truncation etc
                Exon initialExon = model.getExons().get(0);
                long endCoordinate = model.getRange().getEnd(Range.CoordinateSystem.RESIDUE_BASED);
                defline.append(String.format(" location=%s", formatMaturePeptideRange(model,
                        match,
                        cdsRanges,
                        Range.CoordinateSystem.RESIDUE_BASED,
                        "..",
                        // start_codon adjustment
                        initialExon.getRange().getBegin(Range.CoordinateSystem.RESIDUE_BASED) + initialExon.getFrame().getFrame() - 1,
                        endCoordinate)));
                defline.append(String.format(" gene=\"%s\"", model.getGeneSymbol()));
                defline.append(String.format(" product=\"%s\"", VigorUtils.putativeName(match.getReference().getProduct(), match.isFuzzyEnd(), match.isFuzzyBegin())));
                String refDB = model.getAlignment().getAlignmentEvidence().getMatpep_db();
                if (! (refDB == null || refDB.isEmpty() )) {
                    defline.append(String.format(" ref_db=\"%s\"", Paths.get(refDB).getFileName().toString()));
                }
                defline.append(String.format(" ref_id=\"%s\"", match.getReference().getProteinID()));
                bw.write(defline.toString());
                bw.newLine();
                writeSequence(bw, match.getProtein().toBuilder().trim(match.getProteinRange()).build());
            }

        }

    }

    private static String formatMaturePeptideRange(Model model,
                                                   MaturePeptideMatch match,
                                                   List<Range> ranges,
                                                   Range.CoordinateSystem coordinateSystem,
                                                   String rangeDelimiter,
                                                   long startCoordinate,
                                                   long endCoordinate) {
        List<String> rangeStrings = new ArrayList<>(ranges.size());
        Exon initialExon = model.getExons().get(0);
        Exon lastExon = model.getExons().get(model.getExons().size() - 1);

        for (Range range: ranges) {
            LOGGER.trace("range {}-{} start {} sframe {} end {} eframe {} 5p {} 3p {}",
                         range.getBegin(coordinateSystem),
                         range.getEnd(coordinateSystem),
                         startCoordinate,
                         initialExon.getFrame().getFrame(),
                         endCoordinate,
                         lastExon.getFrame().getFrame(),
                         model.isPartial5p(),
                         model.isPartial3p());
            long start = range.getBegin(coordinateSystem);
            String startStr = String.valueOf(start);
            if (match.isFuzzyBegin() || model.isPartial5p() && start == startCoordinate) {
                startStr = "<" + startStr;
            }
            long end = range.getEnd(coordinateSystem);
            String endStr = String.valueOf(end);
            if (match.isFuzzyEnd() || (model.isPartial3p() && end == endCoordinate) ) {
                endStr = ">" + endStr;
            }
            rangeStrings.add(String.format("%s%s%s", startStr, rangeDelimiter, endStr));
        }
        return String.join(",", rangeStrings);
    }

}
