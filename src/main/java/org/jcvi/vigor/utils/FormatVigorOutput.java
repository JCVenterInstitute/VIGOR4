package org.jcvi.vigor.utils;

import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;

/**
 * Created by snettem on 5/24/2017.
 */
public class FormatVigorOutput {

	public static void printModels(List<Model> models, String message) {
		System.out.println(
				"********************************"+message+"**************************************");
		System.out
				.println(String.format("%-32s%-20s%-20s%-20s%-20s", "Gene_Symbol", "Direction","spliceform", "NTSeqRange", "AASeqRange"));

		for (Model model : models) {
			List<Exon> exons = model.getExons();
		    System.out.print(String.format("%-32s", model.getGeneSymbol()));
			System.out.print(String.format("%-20s", model.getDirection()));
			System.out.print(String.format("%-20s",model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform()));
			List<Range> NTranges =exons.stream().map(e -> e.getRange()).collect(Collectors.toList());
			List<Range> AAranges = exons.stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
				System.out.print(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
				System.out.println(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
				System.out.print(String.format("%-72s", ""));
			}
			System.out.println("");
		}
	}

	public static void printModelsWithAllFeatures(List<Model> models) throws FileNotFoundException{
	    String genomeID = models.get(0).getAlignment().getVirusGenome().getId();
	    PrintStream o = new PrintStream(new File(VigorUtils.getVigorWorkSpace()+"/Unit_Test_Output/"+File.separator+genomeID));
	    PrintStream console = System.out;
	    System.setOut(o);
        System.out.println("Genomic Sequence: "+genomeID+"  Genomic Sequence Length: "+ models.get(0).getAlignment().getVirusGenome().getSequence().getLength());
	    System.out
				.println(String.format("%-32s%-20s%-20s%-20s%-20s%-20s%-32s%-20s%-20s%-20s", "Gene_Symbol","ProteinLength", "Direction","spliceform", "NTSeqRange", "AASeqRange", "Scores","partial5'","partial3'","isPseudogene"));
		for (Model model : models) {
			List<Exon> exons = model.getExons();
			System.out.print(String.format("%-32s", model.getGeneSymbol()));
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

	public static void printModels2(List<Model> models) {
		for (Model model : models) {

			System.out.print(String.format("%-32s", model.getGeneSymbol()));
			System.out.print(String.format("%-20s", model.getDirection()));
			List<Range> NTranges = model.getExons().stream().map(e -> e.getRange()).collect(Collectors.toList());
						List<Range> AAranges = model.getExons().stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
				System.out.print(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
				System.out.println(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
				System.out.print(String.format("%-52s", ""));
			}
			System.out.println("");
		}
	}

	public static void printAlignments(List<Alignment> alignments) {
		System.out.println(
				"*********************************Intial list of Alignments*****************************************************");
		System.out.println(String.format("%-32s%-20s%-20s%-20s%-20s", "Protein_ID","Alignment_Tool","Score","NTSeqRange", "AASeqRange","Frame"));

		for (Alignment alignment : alignments) {
			System.out.print(String.format("%-32s", alignment.getViralProtein().getProteinID()));
			System.out.print(String.format("%-20s", alignment.getAlignmentTool_name()));
			System.out.print(String.format("%-20s", alignment.getAlignmentScore().get("ExonerateScore")));
			List<Range> NTranges = alignment.getAlignmentFragments().stream().map(e -> e.getNucleotideSeqRange())
					.collect(Collectors.toList());
			List<Frame> frames = alignment.getAlignmentFragments().stream().map(e->e.getFrame()).collect(Collectors.toList());
			List<Range> AAranges = alignment.getAlignmentFragments().stream().map(e -> e.getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
				System.out.print(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
				System.out.print(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
				System.out.println(String.format("%-20s", frames.get(i).getFrame()));
				System.out.print(String.format("%-72s", ""));
			}
			System.out.println("");
		}
	 }
	
	public static void printModelsWithStart(List<Model> models, String message) {
		System.out.println(
				"********************************"+message+"**************************************");
		System.out
				.println(String.format("%-32s%-20s%-20s%-20s%-20s", "Gene_Symbol", "Direction","spliceform","NTSeqRange", "AASeqRange"));

		for (Model model : models) {
			System.out.print(String.format("%-32s", model.getGeneSymbol()));
			System.out.print(String.format("%-20s", model.getDirection()));
			System.out.print(String.format("%-20s",model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform()));
			List<Range> NTranges = model.getExons().stream().map(e -> e.getRange()).collect(Collectors.toList());
			List<Range> AAranges = model.getExons().stream().map(e -> e.getAlignmentFragment().getProteinSeqRange())
					.collect(Collectors.toList());
			for (int i = 0; i < NTranges.size(); i++) {
				System.out.print(String.format("%-20s", NTranges.get(i).getBegin() + "-" + NTranges.get(i).getEnd()));
				System.out.println(String.format("%-20s", AAranges.get(i).getBegin() + "-" + AAranges.get(i).getEnd()));
				System.out.print(String.format("%-92s", ""));
			}
			System.out.println("");
		}
	}

}
