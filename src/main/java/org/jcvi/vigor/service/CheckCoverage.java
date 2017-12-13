package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.springframework.stereotype.Service;

@Service
public class CheckCoverage implements EvaluateModel {

	
	@Override
	public Model evaluate(Model model,VigorForm form) {
		//getInternalStops(model, cds);
	    determineHomology(model);
		
		
		return null;
	}
    public void determineHomology(Model model){
    	
    	ProteinSequence querySeq = determineCDS(model);
       	ProteinSequence subSeq = model.getAlignment().getViralProtein().getSequence();
    	AminoAcidSubstitutionMatrix blosom50 = BlosumMatrices.blosum50();
    	ProteinPairwiseSequenceAlignment actual = PairwiseAlignmentBuilder
				.createProtienAlignmentBuilder(querySeq,
						subSeq, blosom50).gapPenalty(-8, -8)
				.build();
    	double identity = actual.getPercentIdentity();
   	    Map<String,Double> scores = model.getScores();
   	    scores.put("percentIdentity",identity);
   	    model.setScores(scores);
   	    
    }
	public ProteinSequence determineCDS(Model model){
		model.getExons().sort(Exon.Comparators.Ascending);
		String insertionString = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing().getInsertionString();
		NucleotideSequenceBuilder virusGenomeSeqBuilder = new NucleotideSequenceBuilder(model.getAlignment().getVirusGenome().getSequence());
		if(model.getInsertRNAEditingRange()!=null && insertionString !=null){
		virusGenomeSeqBuilder.insert((int)model.getInsertRNAEditingRange().getBegin(), insertionString);
		}
		NucleotideSequence virusGenomeSeq = virusGenomeSeqBuilder.build();
		AminoAcid replacementAA = model.getAlignment().getViralProtein().getGeneAttributes().getStopTranslationException().getReplacementAA();
		List<Exon> exons = model.getExons();
		long replacementOffset=0;
		if(model.getReplaceStopCodonRange()!=null){
		replacementOffset = model.getReplaceStopCodonRange().getBegin();
		}
		if(replacementOffset != 0){
		replacementOffset = getTranslatedProteinCooridnate(model.getExons(),replacementOffset);
		}
		NucleotideSequenceBuilder NTSeqBuilder=new NucleotideSequenceBuilder("");
		NucleotideSequence NTSeq=null;				
		for(Exon exon : exons){
			NTSeqBuilder.append(virusGenomeSeq.toBuilder(exon.getRange()));
		}
		NTSeq = NTSeqBuilder.build();
		ProteinSequence translatedSeq = IupacTranslationTables.STANDARD.translate(NTSeq);
		ProteinSequenceBuilder proteinSeqBuilder = new ProteinSequenceBuilder(translatedSeq);
		if(replacementOffset !=0 && replacementAA !=null){
		proteinSeqBuilder.replace((int)replacementOffset,replacementAA);
		}
		return translatedSeq;		
	}
	
	public long getTranslatedProteinCooridnate(List<Exon> exons,long NTOffset){
		long translatedProteinLength=0;
		for(Exon exon:exons){
		    Range exonRange = exon.getRange();
		    long length = exonRange.getLength();
		    if(exonRange.intersects(Range.of(NTOffset))){
		        long startCoordinate = exonRange.getBegin()+exon.getFrame().getFrame()-1;
		        long difference = NTOffset-startCoordinate;
		        long proteinCoordinate = difference/3;
		        proteinCoordinate = proteinCoordinate+length;		        
		    }
		    translatedProteinLength = (long)(translatedProteinLength+Math.ceil(length/3));
		}	
		return translatedProteinLength;		
	}
	
	//find internal stops(remove range from the list of stop codon has to be ignored
	public List<Range> getInternalStops(Model model, NucleotideSequence cds){
		//dont have to translate, just find stops by calling a function
		List<Range> internalStops = new ArrayList<Range>();
		Map<Frame,List<Long>> stops = IupacTranslationTables.STANDARD.findStops(cds);
		for(Map.Entry<Frame,List<Long>> pair :stops.entrySet()){
			if(pair.getKey().equals(Frame.ONE)){
			List<Long> cdsStops = (List<Long>)pair.getValue();
			for(Long stop : cdsStops){
				Range NTStopRange = VigorFunctionalUtils.getNTRange(model.getExons(), Range.of(stop));
				internalStops.add(NTStopRange);			
			}
		}
		}
		return internalStops;
	}
		

}
