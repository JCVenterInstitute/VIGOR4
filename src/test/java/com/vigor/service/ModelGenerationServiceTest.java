package com.vigor.service;

import com.vigor.component.AlignmentFragment;
import static org.junit.Assert.*;
import org.jcvi.jillion.core.Range;
import org.junit.Test;
import java.util.ArrayList;
import java.util.List;


/**
 * Created by snettem on 5/19/2017.
 */


public class ModelGenerationServiceTest {




    @Test
    public void generateModel() throws Exception {
    }

    @Test
    public void generateCompatibleFragsChains() throws Exception {
        AlignmentFragment alignmentFragment1 = new AlignmentFragment ();
        List<AlignmentFragment> alignmentFragmentList = new ArrayList<AlignmentFragment> ();

        Range range;
        Range.Builder builder= new Range.Builder ();
        builder.setBegin(200);
        builder.setEnd(302);
        range = builder.build();
        alignmentFragment1.setNucleotideSeqRange (range);
        builder.setBegin ( 1 );
        builder.setEnd ( 34 );
        range = builder.build();
        alignmentFragment1.setProteinSeqRange ( range );
        alignmentFragmentList.add ( alignmentFragment1 );
       AlignmentFragment alignmentFragment2 = new AlignmentFragment ();
        builder.setBegin ( 645);
        builder.setEnd ( 933 );
        range = builder.build();
        alignmentFragment2.setNucleotideSeqRange ( range );
        builder.setBegin ( 3 );
        builder.setEnd ( 99 );
        range = builder.build ();
        alignmentFragment2.setProteinSeqRange ( range );
        alignmentFragmentList.add ( alignmentFragment2 );
        AlignmentFragment alignmentFragment3 = new AlignmentFragment ();
        builder.setBegin ( 780 );
        builder.setEnd ( 1068 );
        range = builder.build ();
        alignmentFragment3.setNucleotideSeqRange ( range );
        builder.setBegin ( 32 );
        builder.setEnd ( 128 );
        range = builder.build ();
        alignmentFragment3.setProteinSeqRange ( range );
        alignmentFragmentList.add ( alignmentFragment3 );
        AlignmentFragment alignmentFragment4 = new AlignmentFragment ();
        builder.setBegin ( 1230 );
        builder.setEnd ( 1452 );
        range=builder.build ();
        alignmentFragment4.setNucleotideSeqRange ( range );
        builder.setBegin ( 124 );
        builder.setEnd ( 326 );
        range = builder.build ();
        alignmentFragment4.setProteinSeqRange ( range );
        alignmentFragmentList.add ( alignmentFragment4 );

        List<List<AlignmentFragment>> expListOfAlignmentFragsList = new ArrayList<List<AlignmentFragment>>(  );
        List<AlignmentFragment> expAlignmentFragsList = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList.add ( alignmentFragment1);
        expAlignmentFragsList.add ( alignmentFragment3 );
        expAlignmentFragsList.add ( alignmentFragment4 );
        expListOfAlignmentFragsList.add( expAlignmentFragsList );
        List<AlignmentFragment> expAlignmentFragsList1 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList1.add ( alignmentFragment1 );
        expAlignmentFragsList1.add ( alignmentFragment4 );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList1);
        List<AlignmentFragment> expAlignmentFragsList2 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList2.add ( alignmentFragment2 );
        expAlignmentFragsList2.add ( alignmentFragment4 );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList2 );
        List<AlignmentFragment> expAlignmentFragsList3 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList3.add ( alignmentFragment3 );
        expAlignmentFragsList3.add ( alignmentFragment4 );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList3 );
        List<AlignmentFragment> expAlignmentFragsList4 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList4.add ( alignmentFragment4 );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList4 );
        ModelGenerationService modelGenerationService = new ModelGenerationService ();
        List<List<AlignmentFragment>> result = modelGenerationService.generateCompatibleFragsChains ( alignmentFragmentList );
        assertEquals(expListOfAlignmentFragsList,result);



    }

}