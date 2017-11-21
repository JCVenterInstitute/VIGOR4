package org.jcvi.vigor.service;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.service.EvaluateModel;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.RNA_Editing;
import org.jcvi.vigor.component.Ribosomal_Slippage;
import org.jcvi.vigor.component.Splicing.SpliceSite;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.forms.VigorForm;
import org.springframework.stereotype.Service;

@Service
public class AdjustExonBoundaries implements EvaluateModel {

	private static final Logger LOGGER = LogManager.getLogger(ModelGenerationService.class);
	private long defaultSearchWindow =50;
	private long minIntronLength = 20;
	
	public List<Model> determine(Model model,VigorForm form){
		List<Model> models=null;
		
		try{
		models = adjustRibosomalSlippage(model);
		for(Model riboAdjustedModel : models){
		models = adjustSpliceSites(riboAdjustedModel);
		}
		for(Model spliceAdjustedModel : models){
			models = adjustRNAEditing(spliceAdjustedModel);
		}
			
		}catch(CloneNotSupportedException e){
			LOGGER.error(e.getMessage(),e);
		}catch(Exception e){
			LOGGER.error(e.getMessage(),e);
		}
		return models;
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
		  Range slippagePoint = Range.of(match.getBegin()+riboSlippage.getSlippage_offset()-1);
		 for(int i=0;i<newModel.getExons().size();i++){
			 Range exonRange = newModel.getExons().get(i).getRange();
			 if(exonRange.intersects(slippagePoint)){
				String pointOfOccurance = determineLocation(exonRange,slippagePoint,3);
				 if(pointOfOccurance.equals("middle")){
					 Exon exon = new Exon();
					 exon = newModel.getExons().get(i).clone();
					 exon.set_5p_adjusted(true);
					 exon.setRange(Range.of(slippagePoint.getBegin()+1+riboSlippage.getSlippage_frameshift(),exonRange.getEnd()));
					 exon.setFrame(newModel.getExons().get(i).getFrame().shift(Math.abs(riboSlippage.getSlippage_frameshift())));
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
		List<Range> matches = cds.findMatches(rna_editing.getRegExp(),true).distinct().collect(Collectors.toList());
		//offset must be used to determine pointOfInsertion
		matches=matches.stream().map(x->x=Range.of(x.getBegin()+CDSStart,x.getEnd()+CDSStart)).sequential().collect(Collectors.toList());
		for(Range match:matches){
			  Model newModel = new Model();
			  newModel = model.clone();
			  Range pointOfInsertion=Range.of(match.getEnd(),match.getEnd()+rna_editing.getInsertionString().length()-1);
			  Exon newExon = new Exon();
			  newExon.setInsertionString(rna_editing.getInsertionString());
			  newExon.setRange(pointOfInsertion);
			  for(int i=0;i<newModel.getExons().size();i++){
					 Range exonRange = newModel.getExons().get(i).getRange();
					 if(exonRange.intersects(pointOfInsertion)){
						String pointOfOccurance = determineLocation(exonRange,pointOfInsertion,2);
						 if(pointOfOccurance.equals("start")){
							 newModel.getExons().get(i).setRange(Range.of(pointOfInsertion.getBegin()+1,exonRange.getEnd()));
							 if(i!=0){
							 Range prevExonRange = newModel.getExons().get(i-1).getRange();
							 newModel.getExons().get(i-1).setRange(Range.of(prevExonRange.getBegin(),pointOfInsertion.getBegin()));
							 }
						 }else{
							 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),pointOfInsertion.getBegin()));
							 if(i!=newModel.getExons().size()-1){
							 Range nextExonRange = newModel.getExons().get(i+1).getRange();
							 newModel.getExons().get(i+1).setRange(Range.of(pointOfInsertion.getBegin()+1,nextExonRange.getEnd()));
							 }
						 }
					 }else if(i!=newModel.getExons().size()-1){
						 Range nextExonRange = newModel.getExons().get(i+1).getRange();
						 if(pointOfInsertion.intersects(Range.of(exonRange.getEnd()+1,nextExonRange.getBegin()-1))){
							 newModel.getExons().get(i).setRange(Range.of(exonRange.getBegin(),pointOfInsertion.getBegin()));
		                     newModel.getExons().get(i+1).setRange(Range.of(pointOfInsertion.getBegin()+1,nextExonRange.getEnd()));
						 }
					 }
				 }
			  newModel.getExons().add(newExon);
			  models.add(newModel);
		}
		if(models.size()==0){
			  models.add(model);
		  }
		}
			return models;
	}
	
	/**
	 * 
	 * @param model
	 * @return : models with permutations and combinations of different splice sites found for each splice region. Exon boundaries are adjusted to the splice region.
	 * @throws CloneNotSupportedException 
	 */
	public List<Model> adjustSpliceSites(Model model) throws CloneNotSupportedException{
    	List<Model> models = new ArrayList<Model>();
    	if(model.getAlignment().getViralProtein().getIntrons().size()>=1){
    	List<Model> tempModels = new ArrayList<Model>();
    	VirusGenome virusGenome = model.getAlignment().getVirusGenome();
    	List<SpliceSite> splicePairs = new ArrayList<SpliceSite>();
    	if(model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getNonCanonical_spliceSites()!=null){
    	splicePairs.addAll(model.getAlignment().getViralProtein().getGeneAttributes().getSplicing().getNonCanonical_spliceSites());}
    	for(int i=0;i<model.getExons().size()-1;i++){
    			if(i!=model.getExons().size()-1){
    				Exon upExon = model.getExons().get(i);
	        		Exon downExon = model.getExons().get(i+1);
    			Range currentExon = upExon.getRange();
    			Range nextExon = downExon.getRange();
    			if(nextExon.getBegin()-currentExon.getEnd()<=minIntronLength){
    				upExon.set_3p_adjusted(true);
    				downExon.set_5p_adjusted(true);
    			}
    			//discuss with paolo and check about this condition
    			if(upExon.is_3p_adjusted()!=true && downExon.is_5p_adjusted()!=true)   {			
    			//Check if donor and acceptor are found at start and end of intron respectively
    			boolean foundSplicePair = false;
    			String expectedDonor = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(currentExon.getEnd()+1,currentExon.getEnd()+2)).build().toString();
    			String expectedAcceptor = model.getAlignment().getVirusGenome().getSequence().toBuilder(Range.of(nextExon.getBegin()-2,nextExon.getBegin()-1)).build().toString();
    			for(SpliceSite var:splicePairs){
    		    	if(var.donor.equals(expectedDonor) && var.acceptor.equals(expectedAcceptor)){
    		    		boolean isCompatible = checkSplicePairCompatibility(currentExon,nextExon,currentExon,nextExon,upExon.getFrame(),downExon.getFrame());
    		    		if(isCompatible){
    		    		foundSplicePair=true;
    		    		tempModels.add(model);
    		    		}
    		    	}
    		    }
    			boolean isBoundaryAdjusted=false;
    			boolean isPesudogene=false;
    			/*//*******************
    			if(!foundSplicePair){
    				Range intronRange = Range.of(currentExon.getEnd()+1,nextExon.getBegin()-1);
    			    if(intronRange.getLength()<=minIntronLength){
    				Map<Frame,List<Long>> intronStops = VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, intronRange);
    				List<Long> upStops=new ArrayList<Long>();
    				List<Long> downStops=new ArrayList<Long>();
    				if(intronStops.get(upExon.getSequenceFrame())!=null){
    					upStops = intronStops.get(upExon.getSequenceFrame());
    				}
    				if(intronStops.get(downExon.getSequenceFrame())!=null){
    					downStops = intronStops.get(downExon.getSequenceFrame());
    				}
    			    if(upStops.size()==0 && downStops.size()==0 && intronRange.getLength()%3==0 ){
    			    	Range adjustedrange = Range.of(currentExon.getBegin(),nextExon.getBegin()-1);
    			    	model.getExons().get(i).setRange(adjustedrange);
    			    	model.getExons().get(i).set_3p_adjusted(true);
    			    	model.getExons().get(i+1).set_5p_adjusted(true);
    			    	isBoundaryAdjusted=true;
    			    	final int exonIndex=i;
    			    	if(models.size()>0){
    			    		models.stream().forEach(x->x.getExons().get(exonIndex).setRange(adjustedrange));
    			    	}if(tempModels.size()>0){
    			    		tempModels.stream().forEach(x->x.getExons().get(exonIndex).setRange(adjustedrange));
    			    	}
    			    	    					
    			    } else{
    			    	isPesudogene=true;
    			    	model.setPseudogene(true);
    					return Arrays.asList(model);
    			    }
    				}
    			}    		
    			*/
    			//**********************
    			
    			if(!foundSplicePair&&!isBoundaryAdjusted&&!isPesudogene){
    			//determine Donor search window
    			long donorStart = currentExon.getEnd()-defaultSearchWindow;
    			if(donorStart<0){
    				donorStart = 0;
    			}
    			long donorEnd = currentExon.getEnd()+defaultSearchWindow;
    			if(donorEnd>model.getAlignment().getVirusGenome().getSequence().getLength()){
    				donorEnd = model.getAlignment().getVirusGenome().getSequence().getLength();
    			}
    			if(donorStart<currentExon.getBegin()){
    				donorStart=currentExon.getBegin();
    			}
    			Range donorSearchWindow = Range.of(donorStart,donorEnd); 
    			Map<Frame,List<Long>> donorInternalStopsMap=VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, donorSearchWindow);
    		    if(donorInternalStopsMap.get(upExon.getSequenceFrame())!=null){
    		    List<Long> donorInternalStops = donorInternalStopsMap.get(upExon.getSequenceFrame());
    		    if(donorInternalStops.size()>0){
    		    	Collections.sort(donorInternalStops);
        		    Long stopCoordinate = donorInternalStops.get(0);
        		    if(stopCoordinate>=currentExon.getBegin()){
        		    donorSearchWindow = Range.of(donorStart,stopCoordinate);
        		    }
    		    }
    		    }
    		    NucleotideSequence donorSearchSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(donorSearchWindow)
    					.build();
    		    
    		    //determine acceptor search window  		    		 
    		    long acceptorStart = nextExon.getBegin()-defaultSearchWindow;
    			if(acceptorStart<0){
    				acceptorStart = 0;
    			}
    			long acceptorEnd = nextExon.getBegin()+defaultSearchWindow;
    			if(acceptorEnd>model.getAlignment().getVirusGenome().getSequence().getLength()){
    				acceptorEnd = model.getAlignment().getVirusGenome().getSequence().getLength();
    			}
    			if(acceptorEnd>nextExon.getEnd()){
    				acceptorEnd = nextExon.getEnd();
    			}
    			Range acceptorSearchWindow = Range.of(acceptorStart,acceptorEnd); 
    			Map<Frame,List<Long>> acceptorInternalStopsMap=VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, acceptorSearchWindow);
    			if(acceptorInternalStopsMap.get(downExon.getSequenceFrame())!=null){
    		    List<Long> acceptorInternalStops = acceptorInternalStopsMap.get(downExon.getSequenceFrame());
    		    if(acceptorInternalStops.size()>0){
    		    	Collections.sort(acceptorInternalStops,Collections.reverseOrder());
        		    Long stopCoordinate = acceptorInternalStops.get(0);
        		    if(stopCoordinate<nextExon.getEnd()){
        		    acceptorSearchWindow = Range.of(stopCoordinate,acceptorEnd);
        		    }
    		    }
    			}    						
    		    NucleotideSequence acceptorSearchSeq = model.getAlignment().getVirusGenome().getSequence().toBuilder(acceptorSearchWindow)
    					.build();    		 
    	        for(SpliceSite splicePair : splicePairs){
    	        	String donor = splicePair.donor;
    	        	String acceptor = splicePair.acceptor;
    	        	List<Range> donorRanges=donorSearchSeq.findMatches(donor,true).distinct().collect(Collectors.toList());
    	           	List<Range> acceptorRanges=acceptorSearchSeq.findMatches(acceptor,true).distinct().collect(Collectors.toList());
    	           	final long acceptorStartTemp1 = acceptorSearchWindow.getBegin();
    	           	final long donorStartTemp1 = donorSearchWindow.getBegin();
    	          donorRanges=donorRanges.stream().map(range->Range.of(range.getBegin()+donorStartTemp1,range.getEnd()+donorStartTemp1)).collect(Collectors.toList());
    	          acceptorRanges = acceptorRanges.stream().map(range->Range.of(range.getBegin()+acceptorStartTemp1,range.getEnd()+acceptorStartTemp1)).collect(Collectors.toList());
    	          boolean isNewSpliceSite = true;
    	    	  tempModels.add(model);
    	        	for(Range donorRange : donorRanges){
    	        		for(Range acceptorRange : acceptorRanges){
    	        		Range foundUpExonRange = Range.of(currentExon.getBegin(),donorRange.getBegin()-1);
    	        		Range foundDownExonRange=Range.of(acceptorRange.getEnd()+1,nextExon.getEnd());
    	        		boolean isCompatible = checkSplicePairCompatibility(currentExon,nextExon,foundUpExonRange,foundDownExonRange,upExon.getFrame(),downExon.getFrame());
    	        		if(isCompatible){
    	        		   if(isNewSpliceSite&&models.size()>0){
    	        			   tempModels.clear();
    	        			   tempModels.addAll(models);
    	        			   models.clear();
    	        		   }    	        		
    	        			for(Model newModelprev : tempModels){
    	        			Model newModel = newModelprev.clone();
    	        			newModel.getExons().get(i).setRange(foundUpExonRange);
    	        			newModel.getExons().get(i+1).setRange(foundDownExonRange);
    	        			float donorScore = VigorFunctionalUtils.generateScore(currentExon.getEnd(), donorRange.getBegin());
    	        			float acceptorScore = VigorFunctionalUtils.generateScore(nextExon.getBegin(), acceptorRange.getBegin());
    	        			float spliceScore = donorScore + acceptorScore;
    	        			Map<String,Float> scores=null;
    	        			if(newModel.getScores()!=null){
    	        			 scores = newModel.getScores();
    	        			 if(scores.containsKey("spliceScore")){
    	        				float existingScore = scores.get("spliceScore");
    	        				float spliceScoreSum = existingScore+spliceScore;
    	        				scores.replace("spliceScore", spliceScoreSum);
    	        				newModel.setScores(scores);
    	        			 
    	        			}else{
    	        				scores = new HashMap<String,Float>();
    	        				scores.put("spliceScore", spliceScore);
    	        				newModel.setScores(scores);
    	        			}   } 	        		   	        		  	        		    
    	        			models.add(newModel);    			 			
    	                    
    	        			}
    	        			
    	        			isNewSpliceSite=false;
    	        		}
    	        		}
      	        	}
    	        	
       	        }
      			}
    			}
    		}
    	}
    	}
    	if(models.size()==0){
    		models.add(model);
    	}
    	return models;
    }
	
	public boolean checkSplicePairCompatibility(Range upExonRange,Range downExonRange,Range foundUpExonRange,Range foundDownExonRange,Frame upExonFrame, Frame downExonFrame){
		int upNucleotides = (int)(upExonRange.getLength() - (upExonFrame.getFrame()-1))%3;
		int downNucleotides = downExonFrame.getFrame()-1;
		int extendedUpNucleotides = Math.abs((int)(upExonRange.getEnd()-foundUpExonRange.getEnd()));
		int extendedDownNucleotides=Math.abs((int)(downExonRange.getBegin()-foundDownExonRange.getBegin()));
		int sum = upNucleotides+downNucleotides+extendedDownNucleotides+extendedUpNucleotides;
		sum = sum%3;
		if(sum==0){
			return true;
		}else{
			return false;
		}
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
