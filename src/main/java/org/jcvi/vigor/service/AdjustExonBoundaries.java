package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.service.EvaluateModel;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.Ribosomal_Slippage;
import org.jcvi.vigor.forms.VigorForm;

public class AdjustExonBoundaries implements EvaluateModel {

	public List<Model> determine(Model model,VigorForm form){
		List<Model> models=null;
		return models;
	}
	
	public Model adjustRibosomalSlippage(Model model){
	   Ribosomal_Slippage riboSlippage= model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
	   List<Model> models = new ArrayList<Model>();
	   long CDSStart =model.getExons().get(0).getRange().getBegin();
	   long CDSEnd = model.getExons().get(model.getExons().size()-1).getRange().getEnd();
	   NucleotideSequence cds = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(CDSStart,CDSEnd))
				.build();
	   List<Range> matches = cds.findMatches(riboSlippage.getSlippage_motif()).collect(Collectors.toList());
	  for(Range motif:matches){
		  Model newModel = new Model();
		  newModel = Model.deepClone(model);
		  Range slippagePoint = Range.of(motif.getEnd()+riboSlippage.getSlippage_offset(),motif.getEnd()+riboSlippage.getSlippage_offset());
		 for(int i=0;i<newModel.getExons().size();i++){
			 Range exonRange = newModel.getExons().get(i).getRange();
			 if(exonRange.intersects(slippagePoint)){
				String location = determineLocation(exonRange,slippagePoint);
				 if(location.equals("middle")){
					 Exon exon = new Exon();
					 exon = Exon.deepClone(exon);
					 exon.set_5p_adjusted(true);
					 exon.setRange(Range.of(slippagePoint.getBegin()+riboSlippage.getSlippage_offset(),exonRange.getEnd()));
					 exon.setFrame(newModel.getExons().get(i).getFrame().shift(riboSlippage.getSlippage_frameshift()));
					 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),slippagePoint.getBegin()));
					 newModel.getExons().add(exon);
					 
				 } else if(location.equals("start")){
					 newModel.getExons().get(i).setRange(Range.of(slippagePoint.getBegin(),exonRange.getEnd()));
				 }else{
					 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),slippagePoint.getBegin()));
				 }
			 }else if(i!=newModel.getExons().size()-1){
				 Range nextExonRange = newModel.getExons().get(i+1).getRange();
				 if(slippagePoint.intersects(Range.of(exonRange.getEnd()+1,nextExonRange.getBegin()-1))){
					 
				 }
			 }
		 }
		
		 models.add(newModel);
		 
	   }
		return model;
		
	}
	
	public Model adjustRNAEditing(Model model){
		
		return model;
	}
	
	public String determineLocation(Range searchRange,Range inputRange){
		long length = searchRange.getLength();
        Range start = Range.of(searchRange.getBegin(),searchRange.getBegin()+length/4);
        Range middle = Range.of(start.getEnd()+1,start.getEnd()+1+length/2);
        if(inputRange.isSubRangeOf(start))
        	return "start";
        else if(inputRange.isSubRangeOf(middle))
        	return "middle";
        else return "end";
	}
	
	
}
