package com.vigor.service;

import com.vigor.component.*;
import org.jcvi.jillion.core.Direction;
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
    public void generateModel(Alignment alignment){
        Map<Direction, List<AlignmentFragment>> alignmentFragsGroupedList =
                alignment.getAlignmentFragments ().stream().collect(Collectors.groupingBy( w -> w.getDirection ()));
        Set<Direction> keyset = alignmentFragsGroupedList.keySet ();
        Iterator iter = keyset.iterator ();
        List<Model> models = new ArrayList<Model> (  );
        while(iter.hasNext ()){
         List<List<AlignmentFragment>> ListOfCompatibleFragsList =  generateCompatibleFragsChains ( alignmentFragsGroupedList.get ( iter.next () ) );
         Iterator iter1 = ListOfCompatibleFragsList.iterator ();
         while(iter1.hasNext ()){
             List<AlignmentFragment> compatibleFragsList = (List<AlignmentFragment>)iter1.next ();
             List<Exon> exons = determineVirusGenomeExons ( compatibleFragsList);
             Model model = new Model();
             model.setExons ( exons );
             model.setAlignment ( alignment );
             model.setGeneSymbol ( alignment.getViralProtein ().getProteinID () );
         }
        }
    }

    /**
     *
     * @param alignmentFragments : alignment fragments that are grouped based on direction is input to the function
     * @return List<List<AlignmentFragment>>: compatible(check for overlap) list of alignment fragments and their permutations and combinations is grouped
     */

    public List<List<AlignmentFragment>> generateCompatibleFragsChains(List<AlignmentFragment> alignmentFragments) {
        List<List<AlignmentFragment>> ListOfCompatibleFragsList = new ArrayList<List<AlignmentFragment>> ();
        int AAOverlapOffset =10;
        int NTOverlapOffset =30;
        for (int i = 0; i < alignmentFragments.size ();i++) {
            long NTStart = alignmentFragments.get ( i ).getNucleotideSeqRange ().getBegin ();
            long NTEnd = alignmentFragments.get ( i ).getNucleotideSeqRange ().getEnd ();
            long AAStart = alignmentFragments.get ( i ).getProteinSeqRange ().getBegin ();
            long AAEnd = alignmentFragments.get ( i ).getProteinSeqRange ().getEnd ();
            List<AlignmentFragment> compatibleFragsList = new ArrayList<AlignmentFragment> ();
            compatibleFragsList.add ( alignmentFragments.get ( i ) );
            if(i!=alignmentFragments.size ()-1) {
                for (int j = i + 1; j < alignmentFragments.size (); j++) {
                    long nextNTStart = alignmentFragments.get ( j ).getNucleotideSeqRange ().getBegin ();
                    long nextNTEnd = alignmentFragments.get ( j ).getNucleotideSeqRange ().getEnd ();
                    long nextAAStart = alignmentFragments.get ( j ).getProteinSeqRange ().getBegin ();
                    long nextAAEnd = alignmentFragments.get ( j ).getProteinSeqRange ().getEnd ();
                    if (nextNTStart >= NTEnd - NTOverlapOffset && nextAAStart >= AAEnd - AAOverlapOffset) {
                        compatibleFragsList.add ( alignmentFragments.get ( j ) );
                    }
                }
            }
            ListOfCompatibleFragsList.add ( compatibleFragsList );
            ListOfCompatibleFragsList.addAll (generateSubChains ( compatibleFragsList ));
            }
        return ListOfCompatibleFragsList;
    }

    /**
     *
     * @param chain: Compatible list of alignment fragments
     * @return List<List<AlignmentFragment>>: Returns List of valid combinations of alignmentFragments
     */

    public List<List<AlignmentFragment>> generateSubChains(List<AlignmentFragment> chain) {
        List<List<AlignmentFragment>> subChainsList = new ArrayList<List<AlignmentFragment>> ();
        int temp = 1;
        if (chain.size () > 2) {
            for (int i = 0; i < chain.size () - 2; i++) {
                List<AlignmentFragment> subChain = new ArrayList<AlignmentFragment> ();
                subChain.addAll ( chain );
                for (int j = 1; j <= temp; j++) {
                    subChain.remove ( 1 );
                }
                subChainsList.add ( subChain );
                temp++;
            }
        }
        return subChainsList;
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
            exon.setDirection ( alignmentFragment.getDirection () );
            exons.add ( exon );
        }
        return exons;
    }




 }



