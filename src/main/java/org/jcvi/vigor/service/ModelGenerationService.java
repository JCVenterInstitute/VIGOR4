package org.jcvi.vigor.service;

import org.jcvi.vigor.component.*;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.FormatVigorOutput;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by snettem on 5/16/2017.
 */
@Service
public class ModelGenerationService {
	
	private long minCondensation;
	private long minIntronLength;
	private static boolean isDebug = false;
	private int AAOverlapOffset;
	private int NTOverlapOffset;

	@Autowired
	private GeneModelGenerationService geneModelGenerationService;
	private static final Logger LOGGER = LogManager.getLogger(ModelGenerationService.class);

	public void generateModels(List<Alignment> alignments, VigorForm form) {
		isDebug = form.isDebug();
		Map<String,String> params  = form.getVigorParametersList();
		if(VigorUtils.is_Integer(params.get("AAOverlap_offset"))){
			AAOverlapOffset=Integer.parseInt(params.get("AAOverlap_offset"));
		}
		if(VigorUtils.is_Integer(params.get("NTOverlap_offset"))){
			NTOverlapOffset=Integer.parseInt(params.get("NTOverlap_offset"));
		}
		if(VigorUtils.is_Integer(params.get("min_condensation"))){
			minCondensation= Integer.parseInt(params.get("min_condensation"));
		}
		minIntronLength=minCondensation*3;
		alignments =  mergeIdenticalProteinAlignments(alignments);
		List<Model> candidateModels = determineCandidateModels(alignments, form);
		geneModelGenerationService.generateGeneModel(candidateModels, form);
	}

	public List<Alignment> mergeIdenticalProteinAlignments(List<Alignment> alignments){
		
		List<Alignment> allOutAlignments = new ArrayList<Alignment>();
		Map<String,List<Alignment>> protein2AlignmentsMap = new HashMap<String,List<Alignment>>();
	    for(Alignment alignment : alignments){
	    	alignment.getAlignmentFragments().sort(AlignmentFragment.Comparators.Ascending);
	    	String proteinID = alignment.getViralProtein().getProteinID();
	    	if(protein2AlignmentsMap.containsKey(proteinID)){
	    		List<Alignment> existingAlignments = protein2AlignmentsMap.get(proteinID);
	    		existingAlignments.add(alignment);
	    		Collections.sort(existingAlignments,new Comparator<Alignment>(){
	    			@Override
	    			public int compare(Alignment alignment1,Alignment alignment2){
	    				return alignment1.getAlignmentFragments().get(alignment1.getAlignmentFragments().size()-1).compareTo(alignment2.getAlignmentFragments().get(0));
	    			}
	    		});
	    	    protein2AlignmentsMap.put(proteinID,existingAlignments);
	    	}else{
	    		protein2AlignmentsMap.put(proteinID,new ArrayList<>(Arrays.asList(alignment)));
	    	}
	    }
	    for(String proteinID : protein2AlignmentsMap.keySet()){
	    	List<Alignment> alignmentsTemp = protein2AlignmentsMap.get(proteinID);
	        Alignment mergedAlignment = alignmentsTemp.get(0);
	        alignmentsTemp.remove(mergedAlignment);
	    	for(Alignment alignment : alignmentsTemp){
	    		mergedAlignment.getAlignmentFragments().addAll(alignment.getAlignmentFragments());
	    		Map<String,Float> scores = mergedAlignment.getAlignmentScore();
	    		Float score = scores.get("ExonerateScore");
	    		score = score+alignment.getAlignmentScore().get("ExonerateScore");
	    		scores.put("ExonerateScore", score);
	    		mergedAlignment.setAlignmentScore(scores);
	    	}
	    	allOutAlignments.add(mergedAlignment);
	    }
	       
	  /*
	    boolean preMerge=false;
	   for(String proteinID : protein2AlignmentsMap.keySet()){
		   	List<Alignment> alignmentsTemp = protein2AlignmentsMap.get(proteinID);
	    	List<Alignment> outAlignments = new ArrayList<Alignment>();
	    	boolean isLastAdded=false;
	    	if(alignmentsTemp.size()>1){
	    		for(int i=0;i<alignmentsTemp.size()-1;i++){	    		    
	    			List<AlignmentFragment> fragments = alignmentsTemp.get(i).getAlignmentFragments();
	    			List<AlignmentFragment> nextFragments = alignmentsTemp.get(i+1).getAlignmentFragments();
	    			Range intersection = fragments.get(fragments.size()-1).getNucleotideSeqRange().intersection(nextFragments.get(0).getNucleotideSeqRange());
	    			boolean isMergeAllowed=true;
	    			if(intersection.getLength()==0){
	    				long end = fragments.get(fragments.size()-1).getNucleotideSeqRange().getEnd();
	    				long start = nextFragments.get(0).getNucleotideSeqRange().getBegin();
	    				if(start<end){
	    					isMergeAllowed=false;
	    				}
	    			}
	    			if(intersection.getLength()%3==0 && isMergeAllowed){
	    				if(preMerge&&outAlignments.size()>0){
	    					fragments.clear();
		    				fragments.addAll(outAlignments.get(outAlignments.size()-1).getAlignmentFragments());
		    				outAlignments.remove(outAlignments.get(outAlignments.size()-1));
		    			}
	    				if(intersection.getLength()>0){	    					
	    					AlignmentFragment mergedFragment = VigorFunctionalUtils.mergeTwoFragments(fragments.get(fragments.size()-1), nextFragments.get(0));
	    					nextFragments.remove(0);
	    					fragments.remove(fragments.size()-1);
	    					fragments.add(mergedFragment);	    										
	    				}
	    		    	
	    		    	fragments.addAll(nextFragments);
	    		    	alignmentsTemp.get(i).setAlignmentFragments(fragments);
	    		    	outAlignments.add(alignmentsTemp.get(i));
	    		    	preMerge=true;
	    		    }else {preMerge=false;	
	    		    if(i+1==alignmentsTemp.size()-1){
	    		    	isLastAdded=true;
	    		    }
	    		    outAlignments.add(alignmentsTemp.get(i));
	    		    } 		    	    
	    		}
	    		if(!preMerge&&isLastAdded){
	    			outAlignments.add(alignmentsTemp.get(alignmentsTemp.size()-1));
	    		}
	    		allOutAlignments.addAll(outAlignments);	    	
	    	}else{
	    		allOutAlignments.addAll(alignmentsTemp);
	    	}
	    }	    	*/	
		return allOutAlignments;
	}
	
	
	
	/**
	 * 
	 * @param alignments
	 * @param form
	 * @return all the possible models of each alignment and after splitting
	 *         models at the sequence gaps
	 */
	public List<Model> determineCandidateModels(List<Alignment> alignments, VigorForm form) {
		//System.out.println("Number of alignments are " + alignments.size());
		List<Model> initialModels = new ArrayList<Model>();
		List<Model> candidateModels = new ArrayList<Model>();
		try{
		for (int i = 0; i < alignments.size(); i++) {
			Alignment alignment = alignments.get(i);
			alignment.getAlignmentFragments().sort(AlignmentFragment.Comparators.Ascending);
			initialModels.addAll(alignmentToModels(alignment, alignment.getAlignmentTool_name()));
		}
		if (isDebug) {
			System.out.println("************Initial Models*************");
			FormatVigorOutput.printModels2(initialModels);
		}
        
		List<Range> sequenceGaps=new ArrayList<Range>();
		// get sequence gaps
		if(initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps() !=null){
		sequenceGaps = initialModels.get(0).getAlignment().getVirusGenome().getSequenceGaps();
		}
		List<Range> validSequenceGaps = new ArrayList<Range>();
		String minGapLenString = "";
		
		if (sequenceGaps !=null && sequenceGaps.size() > 0) {
			minGapLenString = form.getVigorParametersList().get("min_seq_gap_length");
			long minGapLength=0;
			if (VigorUtils.is_Integer(minGapLenString)) {
				minGapLength = Long.parseLong(minGapLenString);
			}

			for (Range gapRange : sequenceGaps) {
				if (gapRange.getLength() >= minGapLength) {
					validSequenceGaps.add(gapRange);
				}
			}
		}

		// split models at sequence gaps
        System.out.println("Count of initial models"+ initialModels.size());
		for (Model model : initialModels) {
			Range diffRange=null;
			for(int i=0;i<model.getExons().size()-1;i++){
			 diffRange = Range.of(model.getExons().get(i).getRange().getEnd(), model.getExons().get(i+1).getRange().getBegin());
			}
			if(diffRange!=null && diffRange.getLength()>=20){
			
			List<Model> newModelsreturned = splitModelAtSequenceGaps(model,validSequenceGaps);
			
			//candidateModels.addAll(splitModelAtSequenceGaps(model, validSequenceGaps));
			if(newModelsreturned.size()>1){
			System.out.println("Before Splitting");
			model.getExons().stream().forEach(System.out::println);
			System.out.println("After splitting at the sequence gaps");
			newModelsreturned.stream().forEach(System.out::println);
			}
			candidateModels.addAll(newModelsreturned);
			}else
			{
				candidateModels.add(model);
			}
		}
       System.out.println("Count after splitting"+candidateModels.size());
		if (isDebug) {
			System.out.println("********After splitting models at the Genome sequence gaps**********");
			FormatVigorOutput.printModels2(candidateModels);
		}}
		catch(CloneNotSupportedException e){
			LOGGER.error(e.getMessage(),e);
		}catch(Exception e){
			LOGGER.error(e.getMessage(),e);
		}

		return candidateModels;
	}

	/**
	 * 
	 * @param alignment
	 * @param alignmentTool
	 * @return Models of each alignment.
	 */
	public List<Model> alignmentToModels(Alignment alignment, String alignmentTool) {
		if(alignment.getViralProtein().getProteinID().equals("P21277.1")){
			System.out.println("I will break");
		}
		
		Map<Direction, List<AlignmentFragment>> alignmentFragsGroupedList = alignment.getAlignmentFragments().stream()
				.collect(Collectors.groupingBy(w -> w.getDirection()));
		Set<Direction> keyset = alignmentFragsGroupedList.keySet();
		Iterator<Direction> iter = keyset.iterator();
		List<Model> models = new ArrayList<Model>();
		for (Direction direction : keyset) {
			System.out.println("for protein " +alignment.getViralProtein().getProteinID());
			List<List<AlignmentFragment>> ListOfCompatibleFragsList = generateCompatibleFragsChains(
					alignmentFragsGroupedList.get(iter.next()), alignmentTool);
			Iterator<List<AlignmentFragment>> iter1 = ListOfCompatibleFragsList.iterator();
			while (iter1.hasNext()) {
				List<AlignmentFragment> compatibleFragsList = (List<AlignmentFragment>) iter1.next();
			    compatibleFragsList=mergeAlignmentFragments(compatibleFragsList, alignment.getVirusGenome());
			    alignment.setAlignmentFragments(compatibleFragsList);
				List<Exon> exons = determineVirusGenomeExons(compatibleFragsList);
				Model model = new Model();
				model.setExons(exons);
				model.setAlignment(alignment);
				model.setGeneSymbol(alignment.getViralProtein().getProteinID());
				model.setDirection(direction);
				model = generateScores(model, alignment);
				List<String> statusList = new ArrayList<String>();
				statusList.add("Initial Model");
				model.setStatus(statusList);
				models.add(model);
			}

		}

		return models;

	}
	/**
	 * 
	 * @param alignment
	 * @return merge two fragments if the intron length is <minIntronLength and if there are no stops in between and if intron length is divisible by 3. Also merge fragments if the missing protein alignment length between two
	 * fragments is less than the minimum_condensation;
	 */
	public List<AlignmentFragment> mergeAlignmentFragments(List<AlignmentFragment> fragments,VirusGenome virusGenome){
		fragments.sort(AlignmentFragment.Comparators.Ascending);
		List<AlignmentFragment> outFragments = new ArrayList<AlignmentFragment>();
		boolean isPreMerge=false;
		if(fragments.size()>1){
    		for(int i=0;i<fragments.size();i++){    			
    			if(i!=fragments.size()-1){
    				AlignmentFragment upFragment = fragments.get(i);
	        		AlignmentFragment  downFragment = fragments.get(i+1);
    			Range currentFragment = upFragment.getNucleotideSeqRange();
    			Range nextFragment = downFragment.getNucleotideSeqRange();
    			Range intronRange=null;
        try{
        intronRange = Range.of(currentFragment.getEnd()+1,nextFragment.getBegin()-1);
        }
        catch(Exception e){
        	System.out.println("Exception " +e.getMessage());
        }  Range missingAAalignRange=Range.of(0,0);
        try{
        if(intronRange.getEnd()-intronRange.getBegin()==-1){
			intronRange=Range.ofLength(0);
		}
      
       
        if(upFragment.getProteinSeqRange().getEnd()+1<=downFragment.getProteinSeqRange().getBegin()-1)  {	
		missingAAalignRange = Range.of(upFragment.getProteinSeqRange().getEnd()+1,downFragment.getProteinSeqRange().getBegin()-1);
        }
       
		if(missingAAalignRange.getBegin()-intronRange.getEnd()==-1 || missingAAalignRange==null){
			missingAAalignRange = Range.ofLength(0);
		}	
        }
        catch(Exception e){
        	System.out.println();
        }
		if((intronRange.getLength()<=minIntronLength && missingAAalignRange.getLength()<=minCondensation)){
		Map<Frame,List<Long>> intronStops = VigorFunctionalUtils.findStopsInSequenceFrame(virusGenome, intronRange);
		List<Long> upStops=new ArrayList<Long>();
		List<Long> downStops=new ArrayList<Long>();
		Frame upSeqFrame= VigorFunctionalUtils.getSequenceFrame(upFragment.getNucleotideSeqRange().getBegin()+upFragment.getFrame().getFrame()-1);
		if(intronStops.get(upSeqFrame)!=null){
			upStops = intronStops.get(upSeqFrame);
		}
		Frame downSeqFrame= VigorFunctionalUtils.getSequenceFrame(upFragment.getNucleotideSeqRange().getBegin()+upFragment.getFrame().getFrame()-1);
		if(intronStops.get(downSeqFrame)!=null){
			downStops = intronStops.get(downSeqFrame);
		}
	    if(upStops.size()==0 && downStops.size()==0 && (intronRange.getLength()%3==0)){
	    	
	    	if(isPreMerge){
	    		upFragment = outFragments.get(outFragments.size()-1);
	    		outFragments.remove(upFragment);
	    		currentFragment = upFragment.getNucleotideSeqRange();
	    	}
	    	Range adjustedNTrange = Range.of(currentFragment.getBegin(),nextFragment.getEnd());
	    	Range adjustedAArange = Range.of(upFragment.getProteinSeqRange().getBegin(),downFragment.getProteinSeqRange().getEnd());
	    	upFragment.setNucleotideSeqRange(adjustedNTrange);
	    	upFragment.setProteinSeqRange(adjustedAArange);
	    	outFragments.add(upFragment);
	    	isPreMerge=true;
	    		    	    					
	    } else{
	    	if(!isPreMerge){
	    	outFragments.add(fragments.get(i));
	    	}
	    	isPreMerge=false;
	    }
		}else{
			if(!isPreMerge){
				outFragments.add(fragments.get(i));
			}
			isPreMerge=false;
		}
				
   }else{
	   if(!isPreMerge){
		   outFragments.add(fragments.get(i));
	   }
	   }
   }
    }else{
    		 outFragments.addAll(fragments);
    		}
		if(outFragments.size()==0){
			outFragments.addAll(fragments);
		}
	
		return outFragments;
	}

	public Model generateScores(Model model, Alignment alignment) {

		Map<String, Float> alignmentScores = alignment.getAlignmentScore();
		model.setScores(alignmentScores);

		return model;
	}

	/**
	 * 
	 * @param initModels
	 * @param form
	 * @param genome
	 * @return Models are split at sequence gaps and new list of models are
	 *         returned
	 */
	public List<Model> splitModelAtSequenceGaps(Model initModel, List<Range> validSequenceGaps) throws CloneNotSupportedException {

		List<Model> newModels = new ArrayList<Model>();
		Model model = new Model();
		model = initModel.clone();
		if (validSequenceGaps.size() > 0) {
			Exon nextExon;
			Exon currentExon;
			Range diffRange;
			List<Exon> firstGroup = new ArrayList<Exon>();
			List<Exon> secondGroup = new ArrayList<Exon>();
			Model firstModel;
			Model secondModel=null;
				List<Exon> modelExons = new ArrayList<Exon>();
				modelExons = model.getExons();
				secondGroup.addAll(modelExons);
				boolean startExist = true;
				
				for (int j = 0; j < modelExons.size(); j++) {
										
					if (j != modelExons.size() - 1) {
						
						
						nextExon = modelExons.get(j + 1);
						currentExon = modelExons.get(j);
						diffRange = Range.of(currentExon.getRange().getEnd(), nextExon.getRange().getBegin());
						if (diffRange.getLength() >= 20) {
							firstGroup.add(modelExons.get(j));
							
							boolean temp = true;
							for (int k = 0; k < validSequenceGaps.size(); k++) {
								if (diffRange.intersects(validSequenceGaps.get(k)) && temp) {
									secondModel = new Model();
									secondModel =model.clone();
									firstModel = new Model();
									firstModel = model.clone();
									if(currentExon.getRange().getEnd()<validSequenceGaps.get(k).getBegin()){
									currentExon.set_3p_adjusted(true);
									}
								    if(nextExon.getRange().getBegin()>validSequenceGaps.get(k).getEnd()){
								    	nextExon.set_5p_adjusted(true);
								    }
									List<Exon> tempFirst = new ArrayList<>();
									tempFirst.addAll(firstGroup);
									List<Exon> tempSecond = new ArrayList<>();
									tempSecond.addAll(secondGroup);
									firstModel.setExons(tempFirst);
									firstModel.setPartial3p(true);
									firstModel.getStatus().add("Model splitted at sequence gaps");
									if (!startExist) {
										firstModel.setPartial5p(true);
									}
									secondGroup.removeAll(tempFirst);
									tempSecond.removeAll(tempFirst);
									secondModel.setExons(tempSecond);
									secondModel.setPartial5p(true);
									secondModel.getStatus().add("Model splitted at sequence gaps");
									newModels.add(firstModel);
									startExist = false;
									firstGroup.clear();
									temp = false;

								}
							}
						}
					}
				}
				if (secondModel !=null && secondModel.getExons().size() > 0) {
					newModels.add(secondModel);
				}
				if(newModels.size()==0){
					newModels.add(initModel);
				}

				return newModels;
			} 
			
		 else {
			newModels.add(initModel);
			return newModels;
		}

	}

	/**
	 *
	 * @param alignmentFragments
	 *            : alignment fragments that are grouped based on direction is
	 *            input to the function
	 * @return List<List<AlignmentFragment>>: compatible(check for overlap) list
	 *         of alignment fragments and their permutations and combinations is
	 *         grouped
	 */

	public List<List<AlignmentFragment>> generateCompatibleFragsChains(List<AlignmentFragment> alignmentFragments,
			String alignmentTool) {
		
				List<AlignmentFragment> compatibleFragsList = new ArrayList<AlignmentFragment>();
				List<AlignmentFragment> clonedCompatibleFragsList = null;
				List<List<AlignmentFragment>> ListOfCompatibleFragsList = new ArrayList<List<AlignmentFragment>>();
				compatibleFragsList.add(alignmentFragments.get(0).clone());
				ListOfCompatibleFragsList.add(compatibleFragsList);
				List<List<AlignmentFragment>> tempList = null;
				if (alignmentFragments.size() >1) {
					for (int j = 1; j < alignmentFragments.size(); j++) {
						boolean temp=true;
						tempList = new ArrayList<List<AlignmentFragment>>();
						for(int k=0;k<ListOfCompatibleFragsList.size();k++){
						List<AlignmentFragment> currentList = ListOfCompatibleFragsList.get(k);
						AlignmentFragment currentFrag = currentList.get(currentList.size()-1);
						clonedCompatibleFragsList=null;
						long NTEnd = currentFrag.getNucleotideSeqRange().getEnd();
						long AAEnd = currentFrag.getProteinSeqRange().getEnd();
						long nextNTStart = alignmentFragments.get(j).getNucleotideSeqRange().getBegin();
						long nextAAStart = alignmentFragments.get(j).getProteinSeqRange().getBegin();
						if (nextNTStart >= NTEnd - NTOverlapOffset && nextAAStart >= AAEnd - AAOverlapOffset) {
							currentList.add(alignmentFragments.get(j).clone());
						}else if(temp){
							clonedCompatibleFragsList = new ArrayList<AlignmentFragment>();
							for(int i=0;i<currentList.size()-1;i++){
								clonedCompatibleFragsList.add(currentList.get(i).clone());
							}
						    clonedCompatibleFragsList.add(alignmentFragments.get(j).clone());
							tempList.add(clonedCompatibleFragsList);
							if(clonedCompatibleFragsList.size()==1){
								temp=false;
							}
						}
					}						
						if(tempList.size()>0){
						ListOfCompatibleFragsList.addAll(tempList);
						}
						
					}
				}
				
			/*
				if (compatibleFragsList.size() > 2 && generateSubChains(alignmentTool)) {
					int temp = 1;
					for (int k = 0; k < compatibleFragsList.size() - 2; k++) {
						List<AlignmentFragment> subChain = new ArrayList<AlignmentFragment>();
						for(AlignmentFragment alignFrag : compatibleFragsList){
							subChain.add(alignFrag);
						}
						for (int j = 1; j <= temp; j++) {
							subChain.remove(1);
						}
						ListOfCompatibleFragsList.add(subChain);
						temp++;
					}
				}*/
			
			
		
		return ListOfCompatibleFragsList;
	}

	/**
	 *
	 * @return returns true if subChains has to be generated. eg: for exonerate
	 *         this function returns false as subChains are not required to be
	 *         generated
	 */

	/*public boolean generateSubChains(String alignmentTool) {

		if (alignmentTool.equals("blast")) {
			return true;
		}
		return false;
	}*/

	/**
	 *
	 * @param compatibleFragsList
	 * @return List<Exon>: nucleotideSequence range of alignment fragment will
	 *         be set as range of exon. 5' and 3' edges of exon are not
	 *         modified.
	 *
	 */
	public List<Exon> determineVirusGenomeExons(List<AlignmentFragment> compatibleFragsList) {
		Iterator<AlignmentFragment> iter = compatibleFragsList.iterator();
		List<Exon> exons = new ArrayList<Exon>();
		while (iter.hasNext()) {
			Exon exon = new Exon();
			AlignmentFragment alignmentFragment = (AlignmentFragment) iter.next();
			exon.setRange(alignmentFragment.getNucleotideSeqRange());
			exon.setAlignmentFragment(alignmentFragment);
			exon.setFrame(alignmentFragment.getFrame());
			Frame sequenceFrame= VigorFunctionalUtils.getSequenceFrame(exon.getRange().getBegin()+exon.getFrame().getFrame()-1);
			exon.setSequenceFrame(sequenceFrame);
			exons.add(exon);
		}
		return exons;
	}
	
    

}
