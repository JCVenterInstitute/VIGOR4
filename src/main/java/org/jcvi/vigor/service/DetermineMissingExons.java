package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.*;
import org.jcvi.jillion.core.Direction;

@Service
public class DetermineMissingExons implements DetermineGeneFeatures {
	private static final Logger LOGGER = LogManager
			.getLogger(DetermineMissingExons.class);
	private int minExonSize=30;
	private int maxIntronSize=2500;
	private int min_missing_AA_size=10;

	@Override
	public List<Model> determine(Model model, VigorForm form) {
		List<Model> outModels=new ArrayList<Model>();
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
		try{
		model.getExons().sort(Exon.Comparators.Ascending);	
		System.out.println("Model "+ model.getGeneSymbol());
		System.out.println("Before"+model.getExons().size());
		model = findMissingExons(model);	
		System.out.println("After"+model.getExons().size());
		model.getExons().sort(Exon.Comparators.Ascending);
		outModels.add(model);	
		}
		catch(Exception e){
		  LOGGER.error(e.getMessage(),e);
		}
		return outModels;
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
	
	public Model findMissingExons(Model model) {
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

		missingExons.addAll(model.getExons());
		model.setExons(missingExons);
		return model;
	}


}
