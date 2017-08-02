package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.stereotype.Service;

<<<<<<< HEAD:src/main/java/com/vigor/service/DetermineMissingExonsService.java
import com.vigor.component.AlignmentFragment;
import com.vigor.component.Exon;
import com.vigor.component.Model;
import com.vigor.component.ViralProtein;
import com.vigor.forms.VigorForm;
import com.vigor.utils.FormatVigorOutput;
import com.vigor.utils.VigorUtils;
=======

import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorUtils;
>>>>>>> 4aafe0f0632a193b0d31d36b3f8efd0b9439888a:src/main/java/org/jcvi/vigor/service/DetermineMissingExonsService.java

@Service
public class DetermineMissingExonsService implements EvaluateModel {
	private static final Logger LOGGER = LogManager.getLogger(DetermineMissingExonsService.class);
	private boolean isDebug = false;

	@Override
	public Model determine(Model inputModel, VigorForm form) {
		String minExonSizeParam = form.getVigorParametersList().get("min_referenceExon_size");
		String percentCoverageParam = form.getVigorParametersList().get("exon_percentage_coverage");
		Model model = new Model();
		try {
			model = (Model) inputModel.clone();
		} catch (CloneNotSupportedException e) {
			LOGGER.error(e.getMessage(), e);
		}
		List<Exon> exons = new ArrayList<Exon>();
		NucleotideSequence NTSequence = model.getAlignment().getVirusGenome().getSequence();
		ProteinSequence AASequence = model.getAlignment().getViralProtein().getSequence();

		Map<Range, Range> missingExons = findMissingExonRanges(model, minExonSizeParam, percentCoverageParam);
		for (Map.Entry<Range, Range> entry : missingExons.entrySet()) {
			Exon exon = performJillionPairWiseAlignment(entry.getKey(), entry.getValue(), NTSequence, AASequence,
					minExonSizeParam, percentCoverageParam);
			if (exon != null) {
				exons.add(exon);
			}
			
		}
		Collections.sort(exons, new Comparator<Exon>(){
			  public int compare(Exon e1, Exon e2){
				   if(e1.getRange().endsBefore(e2.getRange()))
					  return 1;
				  else if(e2.getRange().endsBefore(e1.getRange()))
					  return -1;
				  else return 0;
			   }
			});
		exons.addAll(model.getExons());
		model.setExons(exons);
		
		
		
		return model;
	}

	public Exon performJillionPairWiseAlignment(Range NTRange, Range AARange, NucleotideSequence NTSequence,
			ProteinSequence AASequence, String minExonSizeParam, String percentCoverageParam) {
		List<NucleotideSequence> NTSequenceResidues = new ArrayList<NucleotideSequence>();
		NucleotideSequence NTSequenceResidue1 = NTSequence.toBuilder(NTRange).build();
		NucleotideSequence NTSequenceResidue2 = NTSequence.toBuilder(Range.of(NTRange.getBegin() + 1, NTRange.getEnd()))
				.build();
		NucleotideSequence NTSequenceResidue3 = NTSequence.toBuilder(Range.of(NTRange.getBegin() + 2, NTRange.getEnd()))
				.build();
		NTSequenceResidues.add(NTSequenceResidue1);
		NTSequenceResidues.add(NTSequenceResidue2);
		NTSequenceResidues.add(NTSequenceResidue3);
		ProteinPairwiseSequenceAlignment actual = null;
		ProteinSequence subjectAASequence = AASequence.toBuilder(AARange).build();
		List<ProteinPairwiseSequenceAlignment> alignmentList = new ArrayList<ProteinPairwiseSequenceAlignment>();
		for (NucleotideSequence NTSequenceResidue : NTSequenceResidues) {
			ProteinSequence queryAASequence = IupacTranslationTables.STANDARD.translate(NTSequenceResidue);
			AminoAcidSubstitutionMatrix blosom50 = BlosumMatrices.blosum50();
			actual = PairwiseAlignmentBuilder
					.createProtienAlignmentBuilder(queryAASequence, subjectAASequence, blosom50).gapPenalty(-8, -8)
					.build();
			alignmentList.add(actual);
		}
		final Comparator<ProteinPairwiseSequenceAlignment> comp = (p1, p2) -> Float.compare(p1.getScore(),
				p2.getScore());
		ProteinPairwiseSequenceAlignment bestAlignment = alignmentList.stream().max(comp).get();

		/*
		 * System.out.println("alignment Length :"+bestAlignement.
		 * getAlignmentLength());
		 * System.out.println("Query Range : "+bestAlignement.getQueryRange());
		 * System.out.println("Subject Range : "+bestAlignement.getSubjectRange(
		 * )); System.out.println("Gapped Query Alignment :"+bestAlignement.
		 * getGappedQueryAlignment());
		 * System.out.println("Gapped Subject Alignment : "+bestAlignement.
		 * getGappedSubjectAlignment());
		 * System.out.println("Score: "+bestAlignement.getScore());
		 * System.out.println("Number of gap openings : "+bestAlignement.
		 * getNumberOfGapOpenings());
		 * System.out.println("Number of mismatches : "+bestAlignement.
		 * getNumberOfMismatches());
		 */

		boolean found = false;
		Range modelExonAARange = Range.of(bestAlignment.getSubjectRange().getRange().getBegin()+AARange.getBegin(),bestAlignment.getSubjectRange().getRange().getEnd()+AARange.getBegin());
		
		found = isExon(modelExonAARange, AARange, minExonSizeParam, percentCoverageParam);
		Exon exon = null;
		if (found) {
			exon = new Exon();
			Range range = bestAlignment.getQueryRange().getRange();
			Range modelExonNTRange = Range.of((range.getBegin() * 3)+NTRange.getBegin(), (range.getEnd() * 3)+NTRange.getBegin());
			exon.setRange(modelExonNTRange);
			// exon.setFrame(bestAlignment);
			AlignmentFragment alignmentFragment = new AlignmentFragment();
			alignmentFragment.setDirection(bestAlignment.getQueryRange().getDirection());
			// alignmentFragment.setFrame(frame);
			alignmentFragment.setNucleotideSeqRange(modelExonNTRange);
			alignmentFragment.setProteinSeqRange(modelExonAARange);
			alignmentFragment.setScore(bestAlignment.getScore());
			exon.setAlignmentFragment(alignmentFragment);

		}

		return exon;

	}

	public Map<Range, Range> findMissingExonRanges(Model model, String minExonSizeParam, String percentCoverageParam) {
		//System.out.println("*********Model********Gene Symbol: " + model.getGeneSymbol() + "*******");
		Map<Range, Range> missingExons = new HashMap<Range, Range>();
		ViralProtein viralProtein = model.getAlignment().getViralProtein();
		List<Range> referenceExons = viralProtein.getNTfragments();
		/*System.out.println("Reference protein fragments : Spliceform : "
				+ model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getSpliceform());
		referenceExons.forEach(System.out::println);*/
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
			found = false;
			refExonNTRange = referenceExons.get(i);
			if(i==0){
			refExonAARange = Range.of(0,refExonNTRange.getLength()/3 );
			}
			else{
				refExonAARange = Range.of(referenceExons.get(i-1).getEnd()+1,refExonNTRange.getLength()/3 );
			}
			
			for (int j = 0; j < exons.size(); j++) {
				Exon modelExon = exons.get(j);
				modelExonNTRange = modelExon.getRange();
				modelExonAARange = modelExon.getAlignmentFragment().getProteinSeqRange();
				
				if ( isExon(modelExonAARange, refExonAARange, minExonSizeParam, percentCoverageParam)) {
                    found = true;
					tempExons.remove(exons.get(j));
					preExonEnd = modelExonNTRange.getEnd();
				}
			}
			if (!found) {
				if (i == 0) {
					missingExons.put(Range.of(0, exons.get(0).getRange().getBegin() - 1), refExonAARange);
					List<String> status = new ArrayList<String>();
					status.addAll(model.getStatus());
					status.add("partial gene");
					model.setStatus(status);
					
				}
			   else {
					if (tempExons.size() > 0) {
						missingExons.put(Range.of(preExonEnd + 1, tempExons.get(0).getRange().getBegin() - 1),
								refExonAARange);
					} else {
						missingExons.put(
								Range.of(preExonEnd + 1,
										model.getAlignment().getVirusGenome().getSequence().getLength() - 1),
								refExonAARange);
					}
					
					if(i==referenceExons.size()-1){
						
						List<String> status = new ArrayList<String>();
						status.addAll(model.getStatus());
						status.add("partial gene");
						model.setStatus(status);
					}
				}
			}
		}
		return missingExons;
	}

	public boolean isExon(Range modelExonAARange, Range refExonAARange, String minExonSizeParam, String percentCoverageParam) {
		boolean found = false;
		int minExonSize = 50;
		int exonPercentCoverage = 65;
		if (VigorUtils.is_Integer(minExonSizeParam)) {
			minExonSize = Integer.parseInt(minExonSizeParam);
		}
		if (VigorUtils.is_Integer(percentCoverageParam)) {
			exonPercentCoverage = Integer.parseInt(percentCoverageParam);
		}
		long intersectionLength = modelExonAARange.intersection(refExonAARange).getLength();
		double intersectionPercent = ((100) * (intersectionLength)) / modelExonAARange.getLength();
		if (intersectionLength * 3 >= minExonSize || intersectionPercent >= exonPercentCoverage) {
			found = true;
		}
		return found;
	}
}
