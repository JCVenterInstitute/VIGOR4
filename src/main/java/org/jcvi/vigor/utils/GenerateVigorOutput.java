package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.service.DetermineStart;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

@Service
public class GenerateVigorOutput {
    private static final Logger LOGGER = LogManager
            .getLogger(GenerateVigorOutput.class);
    public void generateOutputFiles(String outputFile,List<Model> geneModels){
        File file = new File(outputFile+".tbl");
        generateTBLReport(file,geneModels);
        File CDSFile = new File(outputFile+".cds");
        generateCDSReport(CDSFile,geneModels);
        File PEPFile = new File(outputFile+".pep");
        generatePEPReport(PEPFile,geneModels);
    }

    public void generateTBLReport(File TBLFile,List<Model> geneModels){
        BufferedWriter bw = null;
        FileWriter fw = null;
        String genomeID = geneModels.get(0).getAlignment().getVirusGenome().getId();
        String regex = Pattern.quote("|");
        String[] genomeIDParts = genomeID.split(regex);
        String proteinIDOfGenome=genomeID;
        if(genomeIDParts.length>=2) {
            proteinIDOfGenome = genomeIDParts[0] + "_" + genomeIDParts[1];
        }else{
            proteinIDOfGenome=genomeIDParts[0];
        }
        try{
            fw = new FileWriter(TBLFile,true);
            bw = new BufferedWriter(fw);
            bw.write(">Features "+genomeID+"\n");
            for(int i=0;i<geneModels.size();i++){
                int temp = i+1;
                proteinIDOfGenome = proteinIDOfGenome+"."+temp;
                Model model = geneModels.get(i);
                model.setProteinID(proteinIDOfGenome);
                Ribosomal_Slippage riboSlippage= model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
                RNA_Editing rna_editing = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing();
                Splicing splicing = model.getAlignment().getViralProtein().getGeneAttributes().getSplicing();
                StringBuilder notes = new StringBuilder("");
                List<Exon> exons = model.getExons();
                long start = exons.get(0).getRange().getBegin()+1;
                long end = exons.get(exons.size()-1).getRange().getEnd()+1;
                bw.write(Long.toString(start)+"\t"+Long.toString(end)+"\t"+"gene"+"\n");
                bw.write("\t\t\tlocus_tag\t\n");
                bw.write("\t\t\tgene\t"+model.getAlignment().getViralProtein().getGeneSymbol()+"\n");
                for(int j=0;j<exons.size();j++){
                    Exon exon = exons.get(j);
                    if(j==0){
                        bw.write(Long.toString(exon.getRange().getBegin()+1)+"\t"+Long.toString(exon.getRange().getEnd()+1)+"\t"+"CDS\n");
                    }else{
                        bw.write(Long.toString(exon.getRange().getBegin()+1)+"\t"+Long.toString(exon.getRange().getEnd()+1)+"\n");
                    }
                }
                bw.write("\t\t\tcodon_start\t"+start+"\n");
                bw.write("\t\t\tprotein_id\t"+proteinIDOfGenome+"\n");
                bw.write("\t\t\tlocus_tag\t\n");
                bw.write("\t\t\tgene\t"+model.getAlignment().getViralProtein().getGeneSymbol()+"\n");
                bw.write("\t\t\tproduct\t"+model.getAlignment().getViralProtein().getProduct()+"\n");
                                if(riboSlippage.isHas_ribosomal_slippage()){
                    bw.write("\t\t\tribosomal_slippage\n");
                }
                if(rna_editing.isHas_RNA_editing()){
                    bw.write("\t\t\texception\tRNA editing\n");
                    notes = notes.append(rna_editing.getNote()+";");
                }
                if(splicing.getNonCanonical_spliceSites()!=null && splicing.getNonCanonical_spliceSites().size()>1){
                    notes.append("non-canonical splicing");
                }
                bw.write("\t\t\tnote\t"+notes+"\n");
                if(model.getInsertRNAEditingRange()!=null){
                    bw.write(model.getInsertRNAEditingRange().getBegin()+"\t"+model.getInsertRNAEditingRange().getEnd()+"\t"+"misc_feature\n");

                }
                if(model.getInsertRNAEditingRange()!=null) {
                    NucleotideSequence subSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(model.getInsertRNAEditingRange()).build();
                    bw.write("\t\t\tnote\tlocation of RNA editing (" + subSeq + "," + rna_editing.getInsertionString() + ") in " + model.getAlignment().getViralProtein().getProduct() + "\n");
                }
            }

        }catch(IOException e ){
            LOGGER.error(e.getMessage(),e);
        }
        finally{
            try{
                if(bw!=null){
                    bw.close();
                }
                if(fw!=null){
                    fw.close();
                }
            }catch(IOException e )
            {
                LOGGER.error(e.getMessage(),e);
            }
        }



    }


    public void generateCDSReport(File CDSFile,List<Model> geneModels) {
       try {
            FileWriter fw = new FileWriter(CDSFile, true);
            BufferedWriter bw = new BufferedWriter(fw);

            for (int i = 0; i < geneModels.size(); i++) {
                Model model = geneModels.get(i);
                String reference_db = model.getAlignment().getAlignmentEvidence().getReference_db();
                ViralProtein refProtein = model.getAlignment().getViralProtein();
                List<Exon> exons = model.getExons();
                Range cdsRange = Range.of(exons.get(0).getRange().getBegin()+1, exons.get(exons.size() - 1).getRange().getEnd()+1);
                StringBuilder defline = new StringBuilder("");
                defline.append(">" + model.getProteinID());
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
                List<String> sequenceLines = java.util.Arrays.asList(model.getCds().toString().split("(?<=\\G.{70})"));
                sequenceLines.stream().forEach(line -> {
                    try {
                        bw.write(line);
                        bw.newLine();
                    } catch (IOException e) {
                        LOGGER.error(e.getMessage(),e);
                    }
                });
                bw.newLine();
                bw.close();
                fw.close();
            }
        }catch(IOException e){
            LOGGER.error(e.getMessage(),e);
        }


    }

    public void generatePEPReport(File PEPFile,List<Model> geneModels){
      try {
            FileWriter fw = new FileWriter(PEPFile, true);
            BufferedWriter bw = new BufferedWriter(fw);

            for (int i = 0; i < geneModels.size(); i++) {
                Model model = geneModels.get(i);
                String reference_db = model.getAlignment().getAlignmentEvidence().getReference_db();
                ViralProtein refProtein = model.getAlignment().getViralProtein();
                List<Exon> exons = model.getExons();
                Range cdsRange = Range.of(exons.get(0).getRange().getBegin()+1, exons.get(exons.size() - 1).getRange().getEnd()+1);
                StringBuilder defline = new StringBuilder();
                defline.append(">" + model.getProteinID());
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
                List<String> sequenceLines = java.util.Arrays.asList(model.getTanslatedSeq().toString().split("(?<=\\G.{70})"));
                sequenceLines.stream().forEach(line -> {
                    try {
                        bw.write(line);
                        bw.newLine();
                    } catch (IOException e) {
                        LOGGER.error(e.getMessage(),e);
                    }
                });
                bw.newLine();
                bw.close();
                fw.close();
            }
        }catch(IOException e){
            LOGGER.error(e.getMessage(),e);
        }


    }
}
