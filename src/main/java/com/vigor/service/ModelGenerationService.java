package com.vigor.service;

import com.vigor.component.*;
import com.vigor.forms.VigorForm;
import com.vigor.utils.FormatVigorOutput;
import com.vigor.utils.VigorUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by snettem on 5/16/2017.
 */
@Service
public class ModelGenerationService {

	//@Autowired
	//private GenerateGeneModelService generateGeneModelService;
	private static final Logger LOGGER = LogManager.getLogger ( ModelGenerationService.class );
	public void generateModels(List<Alignment> alignments,VigorForm form){
		
		List<Model> models = new ArrayList<Model>();
		List<Model> candidateModels = new ArrayList<Model>();
		for(int i=0;i<alignments.size();i++){
			models = alignmentToModels(alignments.get(i),form.getAlignmentTool());
			candidateModels.addAll(models);
	    }
		 if(form.isDebug ()) {
	            System.out.println ( "Initial List Of Models" );
	            FormatVigorOutput.printModels ( candidateModels );
	        }
	//	generateGeneModelService.generateGeneModel(models, form);
	}
	
	
    /**
     *
     * @param alignment : alignment results in multiple models.
     */
  
     public List<Model> alignmentToModels(Alignment alignment,String alignmentTool){
       Map<Direction, List<AlignmentFragment>> alignmentFragsGroupedList =
                alignment.getAlignmentFragments ().stream().collect(Collectors.groupingBy( w -> w.getDirection ()));
        Set<Direction> keyset = alignmentFragsGroupedList.keySet ();
        Iterator iter = keyset.iterator ();
        List<Model> models = new ArrayList<Model> (  );
        for(Direction direction : keyset){
         List<List<AlignmentFragment>> ListOfCompatibleFragsList =  generateCompatibleFragsChains ( alignmentFragsGroupedList.get ( iter.next () ), alignmentTool );
         Iterator iter1 = ListOfCompatibleFragsList.iterator ();
         while(iter1.hasNext ()){
             List<AlignmentFragment> compatibleFragsList = (List<AlignmentFragment>)iter1.next ();
             List<Exon> exons = determineVirusGenomeExons ( compatibleFragsList);
             Model model = new Model();
             model.setExons ( exons );
             model.setAlignment ( alignment );
             model.setGeneSymbol ( alignment.getViralProtein ().getProteinID () );
             model.setDirection (direction);
             model = generateScores ( model,alignment );
             model.setStatus ( Arrays.asList ( "Initial Model" ));
             models.add ( model );
         }
         }
   	 
       
	
		 /*
        * Split the model at gaps(Ns) in the Input virusGenome sequence
        */
       		
	//	List<Model> newModels = splitModelAtSequencingGaps(models, form,alignment.getVirusGenome());
	/*	if(form.isDebug ()) {
            System.out.println ( "Models after splitting each Model at sequence gaps" );
            FormatVigorOutput.printModels ( newModels );
        }    */
	             
        return models;

    }

    public Model generateScores(Model model,Alignment alignment){

     Map<String,Float> alignmentScores = alignment.getAlignmentScore ();
     model.setScores ( alignmentScores );

     return model;
    }
    
    public List<Model> splitModelAtSequencingGaps(List<Model> initModels, VigorForm form,VirusGenome genome){
    	List<Model> newModels = new ArrayList<Model>();
    	
    	//newModels.addAll(initModels);
    	List<Range> sequenceGaps = genome.getSequence().getRangesOfNs();
    	String minGapLenString="";
    	if(sequenceGaps.size()>0){
    	minGapLenString = form.getVigorParametersList().get("min_gap_length");
    	long minGapLength=20;
    	if(VigorUtils.is_Integer(minGapLenString)){
     		minGapLength = Long.parseLong(minGapLenString);
        }
    	List<Range> validSequenceGaps = new ArrayList<Range>();
    	for(Range gapRange : sequenceGaps){
    		if(gapRange.getLength()>=minGapLength){
    			validSequenceGaps.add(gapRange);
    		}
    	}
    	Model newModel;
    	Model tempModel;
    	try{
    	for(int i=0;i<initModels.size();i++){
    	Model model = initModels.get(i);
    	tempModel=(Model)model.clone();
    	for(int k=0; k<model.getExons().size();k++){
    		newModel = new Model();
    	    newModel = (Model) model.clone();
         	for(int j=0;j<validSequenceGaps.size();j++){
    		if(model.getExons().get(k).getRange().endsBefore(sequenceGaps.get(i))){
    			 
    			List<Exon> tempExons = new ArrayList<Exon>();
    			List<Exon> newExons = new ArrayList<Exon>();
    			for(int x=k;x>=0;x--){
    			  newExons.add(model.getExons().get(x));
    			}
    			for(int x=k+1;x<model.getExons().size();x++){
    				tempExons.add(model.getExons().get(x));
    			}
    			newModel.setExons(newExons);
    			tempModel.setExons(tempExons);
    			newModels.add(newModel);
      		}
    	}
    	if(k==model.getExons().size()-1){
    		newModels.add(tempModel);
    	}
    	}
    	}
    	}
    	catch(CloneNotSupportedException e){
    		LOGGER.error(e.getMessage(),e);
    	}
    	}
    	
    	return newModels;
    }



    /**
     *
     * @param alignmentFragments : alignment fragments that are grouped based on direction is input to the function
     * @return List<List<AlignmentFragment>>: compatible(check for overlap) list of alignment fragments and their permutations and combinations is grouped
     */

    public List<List<AlignmentFragment>> generateCompatibleFragsChains(List<AlignmentFragment> alignmentFragments,String alignmentTool) {
        List<List<AlignmentFragment>> ListOfCompatibleFragsList = new ArrayList<List<AlignmentFragment>> ();
        int AAOverlapOffset =10;
        int NTOverlapOffset =30;
        boolean tempFlag = true;
        for (int i = 0; i < alignmentFragments.size ();i++) {


            if (alignmentFragments.get ( i ).isSubChain ()) {
                tempFlag = generateSubChains (alignmentTool);
            }
            if (tempFlag) {
                long NTStart = alignmentFragments.get ( i ).getNucleotideSeqRange ().getBegin ();
                long NTEnd = alignmentFragments.get ( i ).getNucleotideSeqRange ().getEnd ();
                long AAStart = alignmentFragments.get ( i ).getProteinSeqRange ().getBegin ();
                long AAEnd = alignmentFragments.get ( i ).getProteinSeqRange ().getEnd ();
                List<AlignmentFragment> compatibleFragsList = new ArrayList<AlignmentFragment> ();
                compatibleFragsList.add ( alignmentFragments.get ( i ) );
                if (i != alignmentFragments.size () - 1) {
                    for (int j = i + 1; j < alignmentFragments.size (); j++) {
                        long nextNTStart = alignmentFragments.get ( j ).getNucleotideSeqRange ().getBegin ();
                        long nextNTEnd = alignmentFragments.get ( j ).getNucleotideSeqRange ().getEnd ();
                        long nextAAStart = alignmentFragments.get ( j ).getProteinSeqRange ().getBegin ();
                        long nextAAEnd = alignmentFragments.get ( j ).getProteinSeqRange ().getEnd ();
                        if (nextNTStart >= NTEnd - NTOverlapOffset && nextAAStart >= AAEnd - AAOverlapOffset) {
                            alignmentFragments.get ( j ).setSubChain ( true );
                            compatibleFragsList.add ( alignmentFragments.get ( j ) );
                        }
                    }
                }
                ListOfCompatibleFragsList.add ( compatibleFragsList );
                if (compatibleFragsList.size () > 2) {
                    int temp = 1;
                    for (int k = 0; k < compatibleFragsList.size () - 2; k++) {
                        List<AlignmentFragment> subChain = new ArrayList<AlignmentFragment> ();
                        subChain.addAll ( compatibleFragsList );
                        for (int j = 1; j <= temp; j++) {
                            subChain.remove ( 1 );
                        }
                        ListOfCompatibleFragsList.add ( subChain );
                        temp++;
                    }
                }
            }
        }
        return ListOfCompatibleFragsList;
    }


    /**
     *
     * @return returns true if subChains has to be generated.
     * eg: for exonerate this function returns false as subChains are not required to be generated
     */

    public boolean generateSubChains(String alignmentTool){
        alignmentTool="exonerate";
        if(alignmentTool.equals ( "blast" )){
            return true;
        }
        return false;
    }


    /**
     *
     * @param compatibleFragsList
     * @return List<Exon>: nucleotideSequence range of alignment fragment will be set as range of exon.
     * 5' and 3' edges of exon are not modified.
     *
     */
    public List<Exon> determineVirusGenomeExons(List<AlignmentFragment> compatibleFragsList){
        Iterator iter = compatibleFragsList.iterator ();
        List<Exon> exons = new ArrayList<Exon> (  );
        while(iter.hasNext ()){
            Exon exon = new Exon();
            AlignmentFragment alignmentFragment = (AlignmentFragment)iter.next ();
            exon.setRange (alignmentFragment.getNucleotideSeqRange ());
            exon.setAlignmentFragment ( alignmentFragment );
            exons.add ( exon );
        }
        return exons;
    }

 }



