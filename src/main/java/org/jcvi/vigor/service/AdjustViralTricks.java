package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.service.DetermineGeneFeatures;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.RNA_Editing;
import org.jcvi.vigor.component.Ribosomal_Slippage;
import org.jcvi.vigor.component.StopTranslationException;
import org.jcvi.vigor.forms.VigorForm;
import org.springframework.stereotype.Service;

@Service
public class AdjustViralTricks implements DetermineGeneFeatures {

	private static final Logger LOGGER = LogManager.getLogger(ModelGenerationService.class);
    private int leakyStopNotFoundScore=80;
	
	@Override
	public List<Model> determine(Model model,VigorForm form){
		List<Model> riboAdjustedmodels=null;
		List<Model> outputModels = new ArrayList<Model>();
		String leakyStopScoreParam = form.getVigorParametersList().get("LeakyStop_notFound_score");
		if (VigorUtils.is_Integer(leakyStopScoreParam)) {
			leakyStopNotFoundScore = Integer.parseInt(leakyStopScoreParam);
		}
		try{
		riboAdjustedmodels = adjustRibosomalSlippage(model);
		List<Model> rnaEditedModels= new ArrayList<Model>();
		for(Model riboAdjustedModel : riboAdjustedmodels){
			rnaEditedModels.addAll(adjustRNAEditing(riboAdjustedModel));
		}			
		for(Model rnaEditeddModel : rnaEditedModels){
			outputModels.add(checkForLeakyStop(rnaEditeddModel));
		}
		}catch(CloneNotSupportedException e){
			LOGGER.error(e.getMessage(),e);
            System.exit(0);
		}catch(Exception e){
			LOGGER.error(e.getMessage(),e);
            System.exit(0);
		}
		
		return outputModels;
	}
	
	//change findMatches method.While extending the exon boundaries, check that there should not be any internal stops.
	public List<Model> adjustRibosomalSlippage(Model model) throws CloneNotSupportedException{
	   Ribosomal_Slippage riboSlippage= model.getAlignment().getViralProtein().getGeneAttributes().getRibosomal_slippage();
	   List<Model> models = new ArrayList<Model>();
	   if(riboSlippage.isHas_ribosomal_slippage()){
	   long CDSStart =model.getExons().get(0).getRange().getBegin();
	   long CDSEnd = model.getExons().get(model.getExons().size()-1).getRange().getEnd();
	   NucleotideSequence cds = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(CDSStart,CDSEnd))
				.build();
	   List<Range> matches = new ArrayList<Range>();
	   matches = cds.findMatches(riboSlippage.getSlippage_motif()).collect(Collectors.toList());
	   matches=matches.stream().map(x->x=Range.of(x.getBegin()+CDSStart,x.getEnd()+CDSStart)).sequential().collect(Collectors.toList());
	   for(Range match:matches){
		  Model newModel = new Model();
		  newModel = model.clone();
		 // Range slippagePoint = Range.of(match.getEnd()+riboSlippage.getSlippage_offset(),match.getEnd()+riboSlippage.getSlippage_offset());
		 //test logic(vigor3 annotation)
		 Range slippagePoint = Range.of(match.getEnd()+riboSlippage.getSlippage_offset());
		 for(int i=0;i<newModel.getExons().size();i++){
			 Range exonRange = newModel.getExons().get(i).getRange();
			 if(exonRange.intersects(slippagePoint)){
				String pointOfOccurance = determineLocation(exonRange,slippagePoint,3);
				 if(pointOfOccurance.equals("middle")){
					 Exon exon = new Exon();
					 exon = newModel.getExons().get(i).clone();
					 exon.set_5p_adjusted(true);
					 exon.setRange(Range.of(slippagePoint.getBegin()+1+riboSlippage.getSlippage_frameshift(),exonRange.getEnd()));
					 exon.setFrame(Frame.ONE);
					 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),slippagePoint.getBegin()));
					 newModel.getExons().get(i).set_3p_adjusted(true);
					 newModel.getExons().add(exon);
										 
				 } else if(pointOfOccurance.equals("start")){
					 newModel.getExons().get(i).setRange(Range.of(slippagePoint.getBegin()+1+riboSlippage.getSlippage_frameshift(),exonRange.getEnd()));
					 newModel.getExons().get(i).set_5p_adjusted(true);
					 if(i!=0){
					 Range prevExonRange = newModel.getExons().get(i-1).getRange();
					 newModel.getExons().get(i-1).setRange(Range.of(prevExonRange.getBegin(),slippagePoint.getBegin()));
					 newModel.getExons().get(i-1).set_3p_adjusted(true);
					 }
				 }else{
					 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),slippagePoint.getBegin()));
					 newModel.getExons().get(i).set_3p_adjusted(true);
					 if(i!=newModel.getExons().size()-1){
					 Range nextExonRange = newModel.getExons().get(i+1).getRange();
					 newModel.getExons().get(i+1).setRange(Range.of(slippagePoint.getBegin()+1+riboSlippage.getSlippage_frameshift(),nextExonRange.getEnd()));
					 newModel.getExons().get(i+1).set_5p_adjusted(true);
					 }
				 }
			 }else if(i!=newModel.getExons().size()-1){
				 Range nextExonRange = newModel.getExons().get(i+1).getRange();
				 if(slippagePoint.intersects(Range.of(exonRange.getEnd()+1,nextExonRange.getBegin()-1))){
					 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),slippagePoint.getBegin()));
					 newModel.getExons().get(i).set_3p_adjusted(true);
                     newModel.getExons().get(i+1).setRange(Range.of(slippagePoint.getBegin()+1+riboSlippage.getSlippage_frameshift(),nextExonRange.getEnd()));
                     newModel.getExons().get(i+1).set_5p_adjusted(true);
				 }
			 }
		 }
		
		 models.add(newModel);
		 
	   }
	   }
	  if(models.size()==0){
		  models.add(model);
	  }
		return models;
		
	}
	
	public List<Model> adjustRNAEditing(Model model) throws CloneNotSupportedException{
		List<Model> models = new ArrayList<Model>();
		if(model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing().isHas_RNA_editing()){
		RNA_Editing rna_editing = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing();
		long CDSStart =model.getExons().get(0).getRange().getBegin();
		long CDSEnd = model.getExons().get(model.getExons().size()-1).getRange().getEnd();
		NucleotideSequence cds = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(CDSStart,CDSEnd))
					.build();
		List<Range> matches = cds.findMatches(rna_editing.getRegExp()).distinct().collect(Collectors.toList());
		//offset must be used to determine pointOfInsertion
		matches=matches.stream().map(x->x=Range.of(x.getBegin()+CDSStart,x.getEnd()+CDSStart)).sequential().collect(Collectors.toList());
		for(Range match:matches){
			  Model newModel = new Model();
			  newModel = model.clone();
			  Range pointOfInsertion=Range.of(match.getEnd()+rna_editing.getOffset(),match.getEnd()+rna_editing.getInsertionString().length()-1);
			  newModel.setInsertRNAEditingRange(pointOfInsertion);
			  for(int i=0;i<newModel.getExons().size();i++){
					 Range exonRange = newModel.getExons().get(i).getRange();
					 if(exonRange.intersects(pointOfInsertion)){
						String pointOfOccurance = determineLocation(exonRange,pointOfInsertion,2);
						 if(pointOfOccurance.equals("start")){
							 newModel.getExons().get(i).setRange(Range.of(pointOfInsertion.getBegin(),exonRange.getEnd()));
							 if(i!=0){
							 Range prevExonRange = newModel.getExons().get(i-1).getRange();
							 newModel.getExons().get(i-1).setRange(Range.of(prevExonRange.getBegin(),pointOfInsertion.getBegin()-1));
							 }
						 }else{
							 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),pointOfInsertion.getBegin()-1));
							 if(i!=newModel.getExons().size()-1){
							 Range nextExonRange = newModel.getExons().get(i+1).getRange();
							 newModel.getExons().get(i+1).setRange(Range.of(pointOfInsertion.getBegin(),nextExonRange.getEnd()));
							 }else{
							     Exon exon = new Exon();
							     exon.setRange(Range.of(pointOfInsertion.getEnd()+1,exonRange.getEnd()));
							     exon.setFrame(Frame.ONE);
                                 Frame sequenceFrame= VigorFunctionalUtils.getSequenceFrame(exon.getRange().getBegin()+exon.getFrame().getFrame()-1);
                                 exon.setSequenceFrame(sequenceFrame);
                                 exon.setAlignmentFragment(newModel.getExons().get(i).getAlignmentFragment());
                                 newModel.getExons().add(exon);
                             }
						 }
					 }else if(i!=newModel.getExons().size()-1){
						 Range nextExonRange = newModel.getExons().get(i+1).getRange();
						 if(pointOfInsertion.intersects(Range.of(exonRange.getEnd()+1,nextExonRange.getBegin()-1))){
							 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),pointOfInsertion.getBegin()-1));
		                     newModel.getExons().get(i+1).setRange(Range.of(pointOfInsertion.getBegin(),nextExonRange.getEnd()));
						 }
					 }
			}
			 models.add(newModel);
		}
		
		}
		if(models.size()==0){
			  models.add(model);
		  }
			return models;
	}

	public Model checkForLeakyStop(Model model){
		Range range=null;
		StopTranslationException stopTransExce = model.getAlignment().getViralProtein().getGeneAttributes().getStopTranslationException();
		Map<String,Double> scores = model.getScores();
		if(stopTransExce.isHasStopTranslationException()){
            long CDSStart =model.getExons().get(0).getRange().getBegin();
            long CDSEnd = model.getExons().get(model.getExons().size()-1).getRange().getEnd();
            NucleotideSequence cds = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(CDSStart,CDSEnd))
                    .build();
			Optional<Range> match = cds.findMatches(stopTransExce.getMotif()).distinct().findFirst();
			if(match.isPresent()){
			   Range leakyStopRange =  Range.of(match.get().getBegin()+CDSStart,match.get().getEnd()+CDSStart);
			   long start = leakyStopRange.getEnd()+stopTransExce.getOffset();
			   range = Range.of(start,start+2);
			   scores.put("leakyStopScore",100.00);
			   model.setReplaceStopCodonRange(range);			 			   
			 }else{
			   scores.put("leakyStopScore", (double)leakyStopNotFoundScore);				
			}
			if(model.getScores()!=null){
			    scores.putAll(model.getScores());
            }
			model.setScores(scores);
		}		
		return model;	
	}
	
	
	public String determineLocation(Range searchRange,Range inputRange, int noOfLocations){
		long length = searchRange.getLength();
		if(noOfLocations==3){
        Range start = Range.of(searchRange.getBegin(),searchRange.getBegin()+length/4);
        Range middle = Range.of(start.getEnd()+1,start.getEnd()+1+length/2);
        if(inputRange.isSubRangeOf(start))
        	return "start";
        else if(inputRange.isSubRangeOf(middle))
        	return "middle";
        else return "end";
		}
		else{
			Range start = Range.of(searchRange.getBegin(),searchRange.getBegin()+length/2);
			if(inputRange.isSubRangeOf(start)) return "start";
			else return "end";
		}
	}
	
    
	
	
}
