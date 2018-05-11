package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.component.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.vigor.exception.VigorException;

/**
 * Created by snettem on 5/24/2017.
 */
public class FormatVigorOutput {

    private static final Logger LOGGER = LogManager.getLogger(FormatVigorOutput.class);

	public static void printModels(List<Model> models, String message) {
        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
				"*********************************"+message+"*****************************************************");
        content.append(System.lineSeparator());
		content.append(String.format("%-32s%-20s%-20s%-20s%-20s", "Gene_Symbol", "Direction","spliceform", "NTSeqRange", "AASeqRange"));
        content.append(System.lineSeparator());
		for (Model model : models) {
			List<Exon> exons = model.getExons();
            content.append(String.format("%-32s", model.getAlignment().getViralProtein().getProteinID()));
            content.append(String.format("%-20s", model.getDirection()));
            content.append(String.format("%-20s",model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform()));
			List<Range> NTranges =exons.stream().map(e -> e.getRange()).collect(Collectors.toList());
			List<Range> AAranges = exons.stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
                content.append(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
                content.append(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
                content.append(System.lineSeparator());
                content.append(String.format("%-72s", ""));
			}
            content.append(System.lineSeparator());
		}
        LOGGER.debug(content);
	}

	public static void printModelsWithAllFeatures(List<Model> models) throws FileNotFoundException, VigorException {
	    String genomeID = models.get(0).getAlignment().getVirusGenome().getId();
	    PrintStream o = new PrintStream(new File(VigorUtils.getVigorWorkSpace()+"/Unit_Test_Output/"+File.separator+genomeID));
	    PrintStream console = System.out;
	    System.setOut(o);
        System.out.println("Genomic Sequence: "+genomeID+"  Genomic Sequence Length: "+ models.get(0).getAlignment().getVirusGenome().getSequence().getLength());
	    System.out
				.println(String.format("%-32s%-20s%-20s%-20s%-20s%-20s%-32s%-20s%-20s%-20s", "Gene_Symbol","ProteinLength", "Direction","spliceform", "NTSeqRange", "AASeqRange", "Scores","partial5'","partial3'","isPseudogene"));
		for (Model model : models) {
			List<Exon> exons = model.getExons();
			System.out.print(String.format("%-32s", model.getAlignment().getViralProtein().getProteinID()));
			System.out.print(String.format("%-20s",model.getAlignment().getViralProtein().getSequence().getLength()));
			System.out.print(String.format("%-20s", model.getDirection()));
			System.out.print(String.format("%-20s",model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform()));
			List<Range> NTranges =exons.stream().map(e -> e.getRange()).collect(Collectors.toList());
			List<Range> AAranges = exons.stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
				System.out.print(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
				System.out.println(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
				System.out.print(String.format("%-92s", ""));
			}
			Map<String,Double> scores = model.getScores();
			Set<String> keys = scores.keySet();
            System.out.print(String.format("%-40s",""));
			for(String key : keys){
			     System.out.println(String.format("%-32s",key+":"+scores.get(key)));
			     System.out.print(String.format("%-132s",""));
			}
			System.out.print(String.format("%-32s",""));
            System.out.print(String.format("%-20s",model.isPartial5p()));
            System.out.print(String.format("%-20s",model.isPartial3p()));
            System.out.print(String.format("%-20s",model.isPseudogene()));
			System.out.println("");

		}
	}
    public static void printAllGeneModelsWithScores(List<Model> geneModels,String msg){

        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
                ""+msg);
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-10s%-10s%-10s%-10s%-10s%-10s%-30s%-10s%-10s%-20s", "Reference_ID", "%id","%sim","%cov", "%t5","%gap","%t3","start..stop","pep_sz","ref_sz","definition"));
        content.append(System.lineSeparator());
        for(Model model : geneModels){
            ViralProtein viralProtein = model.getAlignment().getViralProtein();
            Map<String,Double> scores  =  model.getScores();
            long cdsBases =0;
            content.append(String.format("%-32s", viralProtein.getProteinID()));
            content.append(String.format("%-10s",String.format("%.02f",scores.get("%identity"))));
            content.append(String.format("%-10s",String.format("%.02f",scores.get("%similarity"))));
            content.append(String.format("%-10s",String.format("%.02f",scores.get("%coverage"))));
            content.append(String.format("%-10s","0.0"));
            content.append(String.format("%-10s","0.0"));
            content.append(String.format("%-10s","0.0"));
            for (int i=0;i<model.getExons().size();i++) {
                Exon exon = model.getExons().get(i);
                if(i!=0){
                    content.append(System.lineSeparator());
                    content.append(String.format("%-92s", ""));
                }
                content.append(String.format("%-30s",exon.getRange().getBegin()+".."+exon.getRange().getEnd()));
                cdsBases=cdsBases+exon.getRange().getLength();
            }
            content.append(String.format("%-10s",cdsBases/3));
            content.append(String.format("%-10s",viralProtein.getSequence().getLength()));
            content.append(String.format("%-20s",model.getGeneSymbol() +" | " +viralProtein.getProduct()));
            content.append(System.lineSeparator());
        }
        LOGGER.info(content);

    }

	public static void printModels2(List<Model> models,String message) {
        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
				"*********************************"+message+"*****************************************************");
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-20s%-20s%-20s", "Protein_ID","Direction","NTSeqRange", "AASeqRange"));
        content.append(System.lineSeparator());
		for (Model model : models) {
            content.append(String.format("%-32s", model.getAlignment().getViralProtein().getProteinID()));
            content.append(String.format("%-20s", model.getDirection()));
			List<Range> NTranges = model.getExons().stream().map(e -> e.getRange()).collect(Collectors.toList());
						List<Range> AAranges = model.getExons().stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
			    if(i!=0){
                    content.append(System.lineSeparator());
                    content.append(String.format("%-52s", ""));
                }
                content.append(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
                content.append(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
                }
            content.append(System.lineSeparator());
		}
        LOGGER.debug(content);
	}

	/*public static void printAlignments(List<Alignment> alignments) {
		System.out.println(
				"*********************************Initial list of Alignments*****************************************************");
		System.out.println(String.format("%-32s%-20s%-20s%-20s%-20s", "Protein_ID","Alignment_Tool","NTSeqRange", "AASeqRange","Frame"));

		for (Alignment alignment : alignments) {
			System.out.print(String.format("%-32s", alignment.getViralProtein().getProteinID()));
			System.out.print(String.format("%-20s", alignment.getAlignmentTool_name()));
			List<Range> NTranges = alignment.getAlignmentFragments().stream().map(e -> e.getNucleotideSeqRange())
					.collect(Collectors.toList());
			List<Frame> frames = alignment.getAlignmentFragments().stream().map(e->e.getFrame()).collect(Collectors.toList());
			List<Range> AAranges = alignment.getAlignmentFragments().stream().map(e -> e.getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
				System.out.print(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
				System.out.print(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
				System.out.println(String.format("%-20s", frames.get(i).getFrame()));
				System.out.print(String.format("%-52s", ""));
			}
			System.out.println("");
		}
	 }*/
    public static void printAlignments(List<Alignment> alignments) {
        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
                "*********************************Initial list of Alignments*****************************************************");
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-20s%-20s%-20s%-20s", "Protein_ID","Alignment_Tool","NTSeqRange", "AASeqRange","Frame"));
        content.append(System.lineSeparator());
        for (Alignment alignment : alignments) {
            content.append(String.format("%-32s", alignment.getViralProtein().getProteinID()));
            content.append(String.format("%-20s", alignment.getAlignmentTool_name()));
            List<Range> NTranges = alignment.getAlignmentFragments().stream().map(e -> e.getNucleotideSeqRange())
                    .collect(Collectors.toList());
            List<Frame> frames = alignment.getAlignmentFragments().stream().map(e->e.getFrame()).collect(Collectors.toList());
            List<Range> AAranges = alignment.getAlignmentFragments().stream().map(e -> e.getProteinSeqRange())
                    .collect(Collectors.toList());
            for (int i = 0; i < NTranges.size(); i++) {
                if(i!=0){
                    content.append(System.lineSeparator());
                    content.append(String.format("%-52s", ""));
                }
                content.append(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
                content.append(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
                content.append(String.format("%-20s", frames.get(i).getFrame()));

            }
           content.append(System.lineSeparator());
        }
        LOGGER.debug(content);
    }
	
	public static void printModelsWithStart(List<Model> models, String message) {
        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
				"*********************************"+message+"*****************************************************");
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-20s%-20s%-20s%-20s", "Gene_Symbol", "Direction","spliceform","NTSeqRange", "AASeqRange"));
        content.append(System.lineSeparator());
		for (Model model : models) {
            content.append(String.format("%-32s", model.getAlignment().getViralProtein().getProteinID()));
            content.append(String.format("%-20s", model.getDirection()));
            content.append(String.format("%-20s",model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform()));
			List<Range> NTranges = model.getExons().stream().map(e -> e.getRange()).collect(Collectors.toList());
			List<Range> AAranges = model.getExons().stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
			    if(i!=0){
                    content.append(System.lineSeparator());
                    content.append(String.format("%-72s", ""));
                }
                content.append(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
                content.append(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
			}
            content.append(System.lineSeparator());
		}
        LOGGER.debug(content);
	}

	public static void printSequenceFeatures(List<Model> geneModels,String msg){
        double identityAvg=0;
        double similarityAvg=0;
        double coverageAvg=0;
        long totalCDSBases=0;
        long totalPepBases=0;
        VirusGenome virusGenome = geneModels.get(0).getAlignment().getVirusGenome();
        String refDb = geneModels.get(0).getAlignment().getAlignmentEvidence().getReference_db();
        StringBuffer content = new StringBuffer("");
        content.append(System.lineSeparator());
        content.append(
                ""+msg);
        content.append(System.lineSeparator());
        content.append(String.format("%-32s%-10s%-10s%-10s%-10s%-10s%-10s%-30s%-10s%-10s%-20s%-20s", "gene_id", "%id","%sim","%cov", "%t5","%gap","%t3","start..stop","pep_sz","ref_sz","ref_id","definition"));
        content.append(System.lineSeparator());
        for(Model model : geneModels){
           ViralProtein viralProtein = model.getAlignment().getViralProtein();
           Map<String,Double> scores  =  model.getScores();
           long cdsBases =0;
            content.append(String.format("%-32s", model.getGeneID()));
            content.append(String.format("%-10s",String.format("%.02f",scores.get("%identity"))));
            content.append(String.format("%-10s",String.format("%.02f",scores.get("%similarity"))));
            content.append(String.format("%-10s",String.format("%.02f",scores.get("%coverage"))));
            content.append(String.format("%-10s","0.0"));
            content.append(String.format("%-10s","0.0"));
            content.append(String.format("%-10s","0.0"));
            for (int i=0;i<model.getExons().size();i++) {
                Exon exon = model.getExons().get(i);
                if(i!=0){
                    content.append(System.lineSeparator());
                    content.append(String.format("%-92s", ""));
                }
                content.append(String.format("%-30s",exon.getRange().getBegin()+".."+exon.getRange().getEnd()));
                cdsBases=cdsBases+exon.getRange().getLength();
            }
            content.append(String.format("%-10s",cdsBases/3));
            content.append(String.format("%-10s",viralProtein.getSequence().getLength()));
            content.append(String.format("%-20s",viralProtein.getProteinID()));
            content.append(String.format("%-20s",model.getGeneSymbol() +" | " +viralProtein.getProduct()));
            content.append(System.lineSeparator());
            totalCDSBases = totalCDSBases+cdsBases;
            totalPepBases = totalPepBases + cdsBases/3;
            identityAvg = identityAvg+scores.get("%identity");
            similarityAvg = similarityAvg +scores.get("%similarity");
            coverageAvg = coverageAvg+scores.get("%coverage");
        }
        identityAvg = identityAvg/geneModels.size();
        similarityAvg = similarityAvg/geneModels.size();
        coverageAvg = coverageAvg/geneModels.size();
        StringBuffer contentSummary = new StringBuffer("");
        contentSummary.append(System.lineSeparator());
        contentSummary.append(String.format("%-32s%-10s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-20s", "Sequence", "Length","Genes","Pseudogenes", "CDS_Bases","Peptide_Bases","%Ref_Identity","%Ref_Similarity","%Ref_Coverage","Ref_DB"));
        contentSummary.append(System.lineSeparator());
        contentSummary.append(String.format("%-32s", virusGenome.getId()));
        contentSummary.append(String.format("%-10s", virusGenome.getSequence().getLength()));
        contentSummary.append(String.format("%-15s",geneModels.size()));
        contentSummary.append(String.format("%-15s","0"));
        contentSummary.append(String.format("%-15s",totalCDSBases));
        contentSummary.append(String.format("%-15s",totalPepBases));
        contentSummary.append(String.format("%-15s",String.format("%.02f",identityAvg)));
        contentSummary.append(String.format("%-15s",String.format("%.02f",similarityAvg)));
        contentSummary.append(String.format("%-15s",String.format("%.02f",coverageAvg)));
        contentSummary.append(String.format("%-20s",refDb));
        contentSummary.append(System.lineSeparator());
        LOGGER.info(contentSummary);
        LOGGER.info(content);

    }


}
