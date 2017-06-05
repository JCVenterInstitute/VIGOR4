package com.vigor.service;

import com.vigor.component.*;
import com.vigor.forms.VigorForm;
import com.vigor.utils.FormatVigorOutput;
import com.vigor.utils.VigorUtils;

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

    /**
     *
     * @param alignment : alignment results in multiple models.
     */
  
	@Autowired
	private VirusGenomeService virusGenomeService;
    public void generateModel(Alignment alignment,VigorForm form){
       Map<Direction, List<AlignmentFragment>> alignmentFragsGroupedList =
                alignment.getAlignmentFragments ().stream().collect(Collectors.groupingBy( w -> w.getDirection ()));
        Set<Direction> keyset = alignmentFragsGroupedList.keySet ();
        Iterator iter = keyset.iterator ();
        List<Model> models = new ArrayList<Model> (  );
        for(Direction direction : keyset){
         List<List<AlignmentFragment>> ListOfCompatibleFragsList =  generateCompatibleFragsChains ( alignmentFragsGroupedList.get ( iter.next () ), form );
         System.out.println("Direction"+direction);
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
         
        if(form.isDebug ()) {
             System.out.println ( "Initial List Of Models" );
             FormatVigorOutput.printModels ( models );
         }
         
        /*
         * Split the model at gaps(Ns) in the Input virusGenome sequence
         */
        
        splitModelAtSequencingGaps(models, form,alignment.getVirusGenome());
         
        

    }

    public Model generateScores(Model model,Alignment alignment){

     Map<String,Float> alignmentScores = alignment.getAlignmentScore ();
     model.setScores ( alignmentScores );

     return model;
    }
    
    public List<Model> splitModelAtSequencingGaps(List<Model> initModels, VigorForm form,VirusGenome genome){
    	List<Model> newModels = new ArrayList<Model>();
    	newModels.addAll(initModels);
    	List<Range> sequenceGaps = genome.getSequence().getRangesOfNs();
    	String minGapLenString="";
    	minGapLenString = form.getVigorParametersList().get("min_gap_length");
    	long minGapLength=20;
    	if(VigorUtils.is_Integer(minGapLenString)){
     		minGapLength = Long.parseLong(minGapLenString);
        }
    	List<Range> filteredRangesOfNs = new ArrayList<Range>();
    	for(int i=0;i<initModels.size();i++){
    	Model model = initModels.get(i);
    	
    	}
    	
    	return newModels;
    }



    /**
     *
     * @param alignmentFragments : alignment fragments that are grouped based on direction is input to the function
     * @return List<List<AlignmentFragment>>: compatible(check for overlap) list of alignment fragments and their permutations and combinations is grouped
     */

    public List<List<AlignmentFragment>> generateCompatibleFragsChains(List<AlignmentFragment> alignmentFragments,VigorForm form) {
        List<List<AlignmentFragment>> ListOfCompatibleFragsList = new ArrayList<List<AlignmentFragment>> ();
        int AAOverlapOffset =10;
        int NTOverlapOffset =30;
        boolean tempFlag = true;
        for (int i = 0; i < alignmentFragments.size ();i++) {


            if (alignmentFragments.get ( i ).isSubChain ()) {
                tempFlag = generateSubChains (form.getAlignmentTool ());
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



