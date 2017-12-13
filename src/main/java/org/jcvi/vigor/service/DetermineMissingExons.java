package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.*;
import org.jcvi.jillion.core.Direction;

@Service
public class DetermineMissingExons implements DetermineGeneFeatures {
	private static final Logger LOGGER = LogManager
			.getLogger(DetermineMissingExons.class);
	private int minExonSize=30;
	//private int exonPercentCoverage = 65;
	private int maxIntronSize=2500;
	private int min_missing_AA_size=10;

	@Override
	public List<Model> determine(Model model, VigorForm form) {
		String minExonSizeParam = form.getVigorParametersList().get(
				"min_exon_size");
		String maxIntronSizeparam = form.getVigorParametersList().get("max_intron_size");
		String minMissingAASizeParam = form.getVigorParametersList().get("min_missing_AA_size");
		if (VigorUtils.is_Integer(minExonSizeParam)) {
			minExonSize = Integer.parseInt(minExonSizeParam);
		}
		if (VigorUtils.is_Integer(maxIntronSizeparam)) {
			maxIntronSize = Integer.parseInt(maxIntronSizeparam);
		}
		if (VigorUtils.is_Integer(minMissingAASizeParam)) {
			min_missing_AA_size = Integer.parseInt(minMissingAASizeParam);
		}
		model.getExons().sort(Exon.Comparators.Ascending);	
		System.out.println("Model "+ model.getGeneSymbol());
		System.out.println("Before"+model.getExons().size());
		model = findMissingExons(model);	
		System.out.println("After"+model.getExons().size());
		model.getExons().sort(Exon.Comparators.Ascending);
		List<Model> models = new ArrayList<Model>();
		models.add(model);
		return models;
	}

	public Exon performJillionPairWiseAlignment(Range NTRange, Range AARange,
			NucleotideSequence NTSequence, ProteinSequence AASequence,Direction modelDirection) {
		Exon exon = null;			
		NucleotideSequence NTSubSequence = NTSequence.toBuilder(NTRange)
				.build();
		ProteinPairwiseSequenceAlignment actual = null;
		ProteinSequence subjectAASequence = AASequence.toBuilder(AARange)
				.build();
		Map<Frame, ProteinPairwiseSequenceAlignment> alignments = new HashMap<Frame, ProteinPairwiseSequenceAlignment>();
		ProteinPairwiseSequenceAlignment bestAlignment = null;		
		for(Frame frame: Frame.forwardFrames()){
		ProteinSequence queryAASequence = IupacTranslationTables.STANDARD.translate(NTSubSequence,frame);
		AminoAcidSubstitutionMatrix blosom50 = BlosumMatrices.blosum50();

		actual = PairwiseAlignmentBuilder
				.createProtienAlignmentBuilder(queryAASequence,
						subjectAASequence, blosom50).gapPenalty(-8, -8)
				.build();
		alignments.put(frame, actual);
		}		
		if(alignments!=null && alignments.size()>0){
		Set<Frame> frameSet = alignments.keySet();
		for (Frame myFrame : frameSet) {
			if (bestAlignment == null) {
				bestAlignment = alignments.get(myFrame);
		}
			if (bestAlignment.getScore() < alignments.get(myFrame).getScore()) {
				bestAlignment = alignments.get(myFrame);
				if(myFrame==Frame.TWO){
					NTRange = Range.of(NTRange.getBegin()+1,NTRange.getEnd());
				}else if (myFrame == Frame.THREE){
					NTRange = Range.of(NTRange.getBegin()+2,NTRange.getEnd());
				}
			   }			
		}
		
		if(bestAlignment.getQueryRange().getDirection().equals(modelDirection)){
   		Range modelExonAARange = Range.of(bestAlignment.getSubjectRange()
				.getRange().getBegin()
				+ AARange.getBegin(), bestAlignment.getSubjectRange()
				.getRange().getEnd()
				+ AARange.getBegin());
			exon = new Exon();
			Range range = bestAlignment.getQueryRange().getRange();
			Range modelExonNTRange = Range.of(
					(range.getBegin() * 3) + NTRange.getBegin(),
					(((range.getEnd()+1) * 3 )-1)+ NTRange.getBegin());
			exon.setRange(modelExonNTRange);
			AlignmentFragment alignmentFragment = new AlignmentFragment();
			alignmentFragment.setDirection(bestAlignment.getQueryRange()
					.getDirection());
			alignmentFragment.setNucleotideSeqRange(modelExonNTRange);
			alignmentFragment.setProteinSeqRange(modelExonAARange);
			alignmentFragment.setScore(bestAlignment.getScore());
			exon.setAlignmentFragment(alignmentFragment);
			exon.setFrame(Frame.ONE);
		}
		
		}
		return exon;
	}
/*
	public Model findMissingExonsWithSpliceFormAbsent(Model model) {
		List<Exon> exons = model.getExons();
		List<Exon> newExons = new ArrayList<Exon>();
		List<Range> sequenceGaps =model.getAlignment().getVirusGenome().getSequenceGaps();
		long refProtSeqLength = model.getAlignment().getViralProtein()
				.getSequence().getLength();
		long NTSeqLength = model.getAlignment().getVirusGenome().getSequence()
				.getLength();
		boolean temp = true;
		for (int i = 0; i < exons.size(); i++) {
			Range currentAARange = exons.get(i).getAlignmentFragment()
					.getProteinSeqRange();
			Range currentNTRange = exons.get(i).getAlignmentFragment()
					.getNucleotideSeqRange();
			Range nextAARange = null;
			Range nextNTRange = null;
			Range missingAAalignRange = null;
			Range missingNTalignRange = null;
			if (i != exons.size() - 1) {
				nextAARange = exons.get(i + 1).getAlignmentFragment()
						.getProteinSeqRange();
				nextNTRange = exons.get(i + 1).getAlignmentFragment()
						.getNucleotideSeqRange();
			}
			if (i == 0 && currentAARange.getBegin() > minMissingAAalign
					&& currentNTRange.getBegin() > minMissingAAalign*3 && temp) {
				missingAAalignRange = Range
						.of(0, currentAARange.getBegin() - 1);
				missingNTalignRange = Range
						.of(0, currentNTRange.getBegin() - 1);
			}
			else if (i == exons.size() - 1) {
				if (refProtSeqLength - minMissingAAalign > currentAARange.getEnd()
						&& NTSeqLength - minMissingAAalign*3 > currentNTRange.getEnd()) {
					missingAAalignRange = Range.of(currentAARange.getEnd() + 1,
							refProtSeqLength - 1);
					missingNTalignRange = Range.of(currentNTRange.getEnd() + 1,
							NTSeqLength);
				}
			} else {

				if (nextAARange.getBegin() - currentAARange.getEnd() > minMissingAAalign
						&& nextNTRange.getBegin() - currentNTRange.getEnd() > minMissingAAalign*3) {
					missingAAalignRange = Range.of(currentAARange.getEnd() + 1,
							nextAARange.getBegin() + 1);
					missingNTalignRange = Range.of(currentNTRange.getEnd() + 1,
							nextNTRange.getBegin() + 1);

				}
			}
           
			if (missingAAalignRange != null && missingNTalignRange != null) {
			long maxSearchLength = maxIntronSize+(missingAAalignRange.getLength()*3);
			if(missingNTalignRange.getLength()>maxSearchLength){
			   if(i==0 &&temp){
				missingNTalignRange = Range.of(missingNTalignRange.getEnd()-maxSearchLength,missingNTalignRange.getEnd());}
			   else{
				missingNTalignRange = Range.of(missingNTalignRange.getBegin(),missingNTalignRange.getBegin()+maxSearchLength);
			   }
			}
			 boolean sequenceGap=false;
			 Exon exon=null;
			if(i==0 || i==exons.size()-1){	
		  		if(sequenceGaps!=null){
				if(i==0&&temp){
					Collections.reverse(sequenceGaps);
				}
			   	for(Range range:sequenceGaps){
					if(range.intersects(missingNTalignRange)){
						Range intersection = range.intersection(missingNTalignRange);
						Range leftOver=null;
						if(i==0&&temp){
						leftOver = Range.of(intersection.getEnd()+1,missingNTalignRange.getEnd());
						}else{
						leftOver = Range.of(missingNTalignRange.getBegin(),intersection.getBegin()-1);
						}
						if(leftOver!=null && leftOver.getLength()>=20){
						missingNTalignRange=leftOver;
						}else{
						sequenceGap=true;
						}
						break;
					}
				}
				
				}
			}
	
				if(!sequenceGap){
				exon = performJillionPairWiseAlignment(
						missingNTalignRange, missingAAalignRange, model
								.getAlignment().getVirusGenome().getSequence(),
						model.getAlignment().getViralProtein().getSequence(),false,model.getDirection());
				}
				if(exon!=null){
					newExons.add(exon);
				}
			}
			if(temp){
			temp=false;
			Collections.reverse(sequenceGaps);
			i--;
			}
		}
		newExons.addAll(exons);
		model.setExons(newExons);
		return model;
	}
*/
	
	public Model findMissingExons(Model model) {
		if(model.getGeneSymbol().equals("Model gi|392976861|ref|YP_006491245.1|")){
			System.out.println("found");
		}
		//List<Model> models = new ArrayList<Model>();
		List<Exon> exons = model.getExons();
		List<Exon> missingExons = new ArrayList<Exon>();
		long proteinLength = model.getAlignment().getViralProtein().getSequence().getLength();
		NucleotideSequence NTSeq = model.getAlignment().getVirusGenome().getSequence();
		ProteinSequence AASeq = model.getAlignment().getViralProtein().getSequence();
		List<Range> sequenceGaps = model.getAlignment().getVirusGenome().getSequenceGaps();
		long NTSeqLength = NTSeq.getLength();
		Range preAARange = null;
		Range preNTRange = null;
		Exon prevExon=null;
		boolean afterLastExon=false;
		for(int i=0;i<=exons.size();i++){
			int j=0;
			if(i==exons.size()){
			   j=exons.size()-1;
			   preAARange=null;
			   preNTRange=null;
			   afterLastExon = true;
			}else{
				j=i;
			}
			Exon exon = exons.get(j);
			Range AARange = exon.getAlignmentFragment().getProteinSeqRange();
			Range NTRange = exon.getRange();
			Range missingAARange=null;
			Range missingNTRange = null;
			if(preAARange!=null){
				missingAARange = Range.of(preAARange.getEnd()+1,AARange.getBegin()-1);
				missingNTRange = Range.of(preNTRange.getEnd()+1,NTRange.getBegin()-1);
			}else{
				if(i!=exons.size()){
				missingAARange = Range.of(0,AARange.getBegin());
				missingNTRange = Range.of(0,NTRange.getBegin());	
				
				}else{
				missingAARange = Range.of(AARange.getEnd()+1,proteinLength-1);
				missingNTRange = Range.of(NTRange.getEnd()+1,NTSeqLength-1);
				}
			}
			long temp = maxIntronSize+(missingAARange.getLength()*3);
			if(missingNTRange.getLength()>temp){
				missingNTRange = Range.of(missingNTRange.getBegin(),missingNTRange.getBegin()+temp);
			}
			boolean sequenceGap=false;
			if(sequenceGaps!=null){
				for(Range range:sequenceGaps){
					if(range.intersects(missingNTRange)){
						Range intersection = range.intersection(missingNTRange);
						Range leftOver=null;
						if(i==0){
						leftOver = Range.of(intersection.getEnd()+1,missingNTRange.getEnd());
						}else{
						leftOver = Range.of(missingNTRange.getBegin(),intersection.getBegin()-1);
						}
						if(leftOver!=null && leftOver.getLength()>=20){
							missingNTRange=leftOver;
						}else{
						sequenceGap=true;
						}
						break;
					}
				}
				
			}	
			if(!sequenceGap){
			if(missingAARange.getLength()>=min_missing_AA_size){
			Exon determinedExon = performJillionPairWiseAlignment(missingNTRange,
						missingAARange, NTSeq,AASeq,model.getDirection());
			if(determinedExon.getRange().getLength()>=minExonSize){
			    missingExons.add(determinedExon);
			    if(prevExon!=null){
			    prevExon.set_3p_adjusted(false);
			    }
			    if(afterLastExon){
			    exon.set_3p_adjusted(false);
			    }else
			    {
			    exon.set_5p_adjusted(false);
			    }
			}
			}
			}
			
		preNTRange=NTRange;
		preAARange=AARange;
		prevExon = exon;
		}
		/*ViralProtein viralProtein = model.getAlignment().getViralProtein();
		List<Range> referenceExons = viralProtein.getNTfragments();
		List<Range> sequenceGaps = model.getAlignment().getVirusGenome().getSequenceGaps();
		NucleotideSequence NTSequence = model.getAlignment().getVirusGenome()
				.getSequence();
		List<Exon> foundMissingExons = new ArrayList<Exon>();
		ProteinSequence AASequence = model.getAlignment().getViralProtein()
				.getSequence();
		long temp=AASequence.getLength();
		List<Exon> exons = model.getExons();
		Range refExonNTRange;
		Range refExonAARange;
		Range modelExonNTRange;
		Range modelExonAARange;
		Range prevRefExonNTRange = Range.of(0, 0);
		long preExonEnd = 0;
		long preExonProteinEnd=0;
		List<Exon> tempExons = new ArrayList<Exon>();
		tempExons.addAll(model.getExons());
		Range prevRefExonAARange = Range.of(0, 0);
		boolean found = false;
		for (int i = 0; i < referenceExons.size(); i++) {
			found = false;
			Map<Range,Range> missingAlignFrags = new HashMap<Range,Range>();
			refExonNTRange = referenceExons.get(i);
		
			if (i == 0) {
				double a = Math.floor(refExonNTRange.getLength()/3);
				long b = Math.round(a);
				refExonAARange = Range.of(0,b);
			} else {
				refExonAARange = Range.of(
						prevRefExonAARange.getEnd() + 1,
						prevRefExonAARange.getEnd() + 1
								+ refExonNTRange.getLength() / 3 - 1);
			}

			for (int j = 0; j < exons.size(); j++) {
				Exon modelExon = exons.get(j);
				modelExonNTRange = modelExon.getRange();
				modelExonAARange = modelExon.getAlignmentFragment()
						.getProteinSeqRange();

				if (isExon(modelExonAARange, refExonAARange)) {
					found = true;
					tempExons.remove(exons.get(j));
					preExonEnd = modelExonNTRange.getEnd();
					preExonProteinEnd=modelExonAARange.getEnd();
				}				
			}

			if (!found) {
				Range missingNTRange = null;
				Range missingAARange = null;
				if (i == 0) {
					missingNTRange = Range.of(0, exons.get(0).getRange()
							.getBegin() - 1);
					missingAARange = Range.of(0,exons.get(0).getAlignmentFragment().getProteinSeqRange().getBegin()-1);

				} else {
					long intronSize = refExonNTRange.getBegin()
							- prevRefExonNTRange.getEnd() - 1;

					if (tempExons.size() > 0) {

						missingNTRange = Range.of(preExonEnd + 1 + intronSize,
								tempExons.get(0).getRange().getBegin() - 1);
						missingAARange = Range.of(preExonProteinEnd+1,tempExons.get(0).getAlignmentFragment().getProteinSeqRange().getBegin()-1);
					} else {
						
						long startTemp = preExonEnd+1+intronSize;
						long endTemp = startTemp+refExonNTRange.getLength();
						if(endTemp>NTSequence.getLength()-1){
							endTemp = NTSequence.getLength()-1;
						}
						missingNTRange = Range.of(startTemp,
								endTemp);
						missingAARange = Range.of(preExonProteinEnd+1,AASequence.getLength()-1);
					}
				}
				long maxSearchLength = maxIntronSize+(missingAARange.getLength()*3);
				if(missingNTRange.getLength()>maxSearchLength){
				   if(i==0){
					   missingNTRange = Range.of(missingNTRange.getEnd()-maxSearchLength,missingNTRange.getEnd());}
				   else{
					   missingNTRange = Range.of(missingNTRange.getBegin(),missingNTRange.getBegin()+maxSearchLength);
				   }
				}
				boolean sequenceGap=false;
					Exon exon=null;
					if(i==0 || i==referenceExons.size()-1){	
				  		if(sequenceGaps!=null){
						if(i==0){
							Collections.reverse(sequenceGaps);
						}
					   	for(Range range:sequenceGaps){
							if(range.intersects(missingNTRange)){
								Range intersection = range.intersection(missingNTRange);
								Range leftOver=null;
								if(i==0){
								leftOver = Range.of(intersection.getEnd()+1,missingNTRange.getEnd());
								}else{
								leftOver = Range.of(missingNTRange.getBegin(),intersection.getBegin()-1);
								}
								if(leftOver!=null && leftOver.getLength()>=20){
									missingNTRange=leftOver;
								}else{
								sequenceGap=true;
								}
								break;
							}
						}
						
						}
					}
				if(i==0){
						Collections.reverse(sequenceGaps);
					}
				boolean isAlignMissing=true;
				try{
				if(!sequenceGap){
					if(missingNTRange.getLength()>10 && missingAARange.getLength()>10)
					{
				exon = performJillionPairWiseAlignment(missingNTRange,
						missingAARange, NTSequence, AASequence,true,model.getDirection());}
				}else{
					isAlignMissing=false;
					}
				}
				catch(Exception e){
					System.out.println("Exception caught");
				}
				if(i==0&&exon==null&&isAlignMissing){
				model.setPartial5p(true);
				}else if(i==referenceExons.size()-1&&exon==null){
				model.setPartial3p(true);	
				}
			   if(exon!=null){
				    preExonEnd=exon.getRange().getEnd();
				    preExonProteinEnd=exon.getAlignmentFragment().getProteinSeqRange().getEnd();
					foundMissingExons.add(exon);
				}
			}
			prevRefExonNTRange = refExonNTRange;
			prevRefExonAARange = refExonAARange;
		}
		foundMissingExons.addAll(exons);
		model.setExons(foundMissingExons);
		}*/
		missingExons.addAll(model.getExons());
		model.setExons(missingExons);
		return model;
	}

	/*public boolean isExon(Range modelExonAARange, Range refExonAARange) {
		boolean found = false;
		long intersectionLength = modelExonAARange.intersection(refExonAARange)
				.getLength();
		double intersectionPercent = ((100) * (intersectionLength))
				/ refExonAARange.getLength();
		if (intersectionPercent >= exonPercentCoverage) {
			found = true;
		}		
		return found;
	}*/
		

}
