package com.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.stereotype.Service;

import com.vigor.component.Exon;
import com.vigor.component.Model;
import com.vigor.component.ViralProtein;
import com.vigor.forms.VigorForm;
import com.vigor.utils.VigorUtils;

@Service
public class DetermineMissingExonsService implements EvaluateModel {

	@Override
	public Model determine(Model model, VigorForm form) {
		String exon_size_string = form.getVigorParametersList().get("min_exon_size");
		String percent_coverage_string = form.getVigorParametersList().get("exon_percentage_coverage");
		int minExonSize = 50;
		int exonPercentCoverage = 65;
		if (VigorUtils.is_Integer(exon_size_string)) {
			minExonSize = Integer.parseInt(exon_size_string);
		}
		if (VigorUtils.is_Integer(percent_coverage_string)) {
			exonPercentCoverage = Integer.parseInt(percent_coverage_string);
		}
		
		Map<Range,Range> missingExons = findMissingExonRanges(model, exonPercentCoverage, minExonSize);
		NucleotideSequence NTSequence = model.getAlignment().getVirusGenome().getSequence();
		ProteinSequence AASequence = model.getAlignment().getViralProtein().getSequence();
		for (Map.Entry<Range, Range> entry : missingExons.entrySet()) {
			getPairWiseAlignment(entry.getKey(),entry.getValue(),NTSequence,AASequence);
			
			
			
			
		}
		return model;
	}
	
	
	public void getPairWiseAlignment(Range NTRange, Range AARange, NucleotideSequence NTSequence, ProteinSequence AASequence){
		
		List<NucleotideSequence> NTSequenceResidues = new ArrayList<NucleotideSequence>();
		NucleotideSequence NTSequenceResidue1 = NTSequence.toBuilder(NTRange).build();	
		NucleotideSequence NTSequenceResidue2 = NTSequence.toBuilder(Range.of(NTRange.getBegin()+1,NTRange.getEnd())).build();		
		NucleotideSequence NTSequenceResidue3 = NTSequence.toBuilder(Range.of(NTRange.getBegin()+2,NTRange.getEnd())).build();	
		NTSequenceResidues.add(NTSequenceResidue1);
		NTSequenceResidues.add(NTSequenceResidue2);
		NTSequenceResidues.add(NTSequenceResidue3);
		ProteinSequence AASequenceResidue = AASequence.toBuilder(AARange).build();
		for(NucleotideSequence NTSequenceResidue : NTSequenceResidues){
		performJillionPairWiseAlignment(NTSequenceResidue,AASequenceResidue);
		}
		
	}
	
    public void  performJillionPairWiseAlignment(NucleotideSequence queryNTSequence,ProteinSequence subjectAASequence){
    	
    	   	
        ProteinSequence queryAASequence = IupacTranslationTables.STANDARD.translate(queryNTSequence);
        AminoAcidSubstitutionMatrix blosom50 = BlosumMatrices.blosum50();
    	ProteinPairwiseSequenceAlignment actual = PairwiseAlignmentBuilder.createProtienAlignmentBuilder(queryAASequence, subjectAASequence, blosom50)
				.gapPenalty(-8, -8)	
				.useGlobalAlignment()
				.build();
    	
    	
    	System.out.println("alignment Length :"+actual.getAlignmentLength());
    	System.out.println("Query Range : "+actual.getQueryRange());
    	System.out.println("Subject Range : "+actual.getSubjectRange());
    	System.out.println("Gapped Query Alignment :"+actual.getGappedQueryAlignment());
    	System.out.println("Gapped Subject Alignment : "+actual.getGappedSubjectAlignment());
    	System.out.println("Score: "+actual.getScore());
       	System.out.println("Number of gap openings : "+actual.getNumberOfGapOpenings());
    	System.out.println("Number of mismatches : "+actual.getNumberOfMismatches());
    	
    	
    }

	public Map<Range,Range> findMissingExonRanges(Model model, int exonPercentCoverage, int minExonSize) {
		System.out.println("*********Model********Gene Symbol: " + model.getGeneSymbol() + "*******");
		Map<Range, Range> missingExons = new HashMap<Range, Range>();
		ViralProtein viralProtein = model.getAlignment().getViralProtein();
		List<Range> referenceExons = viralProtein.getNTfragments();
		System.out.println("Reference protein fragments : Spliceform : "
				+ model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform());
		referenceExons.forEach(System.out::println);
		List<Exon> exons = model.getExons();
		Range refExonNTRange;
		Range refExonAARange;
		Range modelExonNTRange;
		Range modelExonAARange;
		long preExonEnd = 0;
		List<Exon> tempExons = new ArrayList<Exon>();
		tempExons.addAll(model.getExons());
		boolean found = false;
		for (int i = 0; i < referenceExons.size(); i++) {
			found=false;
			refExonNTRange = referenceExons.get(i);
			if(i!=referenceExons.size()-1){
			refExonAARange = Range.of(refExonNTRange.getBegin() / 3, refExonNTRange.getEnd() / 3);
			}
			else{
				refExonAARange = Range.of(refExonNTRange.getBegin() / 3, (refExonNTRange.getEnd() / 3)-1);
			}
			for (int j = 0; j < exons.size(); j++) {
				Exon modelExon = exons.get(j);
				modelExonNTRange =modelExon.getRange();
				modelExonAARange = modelExon.getAlignmentFragment().getProteinSeqRange();
				long intersectionLength = modelExonAARange.intersection(refExonAARange).getLength();
				double intersectionPercent = ((100) * (intersectionLength)) / modelExonAARange.getLength();
				if (intersectionLength * 3 >= minExonSize || intersectionPercent >= exonPercentCoverage) {
					
					found = true;
					modelExon = adjustExonBoundaries(modelExon,refExonNTRange);
					tempExons.remove(exons.get(j));
					preExonEnd = modelExonNTRange.getEnd();
				}
			}
			if (!found) {
				if (i == 0) {
					missingExons.put(Range.of(0, exons.get(0).getRange().getBegin()-1), refExonAARange);
				} else {
					if (tempExons.size() > 0) {
						missingExons.put(Range.of(preExonEnd+1, tempExons.get(0).getRange().getBegin()-1), refExonAARange);
					} else {
						missingExons.put(
								Range.of(preExonEnd+1,
										model.getAlignment().getVirusGenome().getSequence().getLength() - 1),
								refExonAARange);
					}
				}
			}
		}
		return missingExons;
	}
	
	public Exon adjustExonBoundaries(Exon exon,Range refExonNTRange){
		
		Range modelExonNTRange = exon.getRange();
		Range modelExonProteinAlignRange = exon.getAlignmentFragment().getProteinSeqRange();
		
		return exon;
	}
	
	
	
}
