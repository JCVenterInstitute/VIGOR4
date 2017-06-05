package com.vigor.service;

import com.vigor.component.Alignment;
import com.vigor.component.AlignmentFragment;
import static org.junit.Assert.*;

import com.vigor.component.ViralProtein;
import com.vigor.forms.VigorForm;
import com.vigor.utils.VigorUtils;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.FastaRecord;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.junit.Before;
import org.junit.Test;
import org.springframework.beans.factory.annotation.Autowired;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;


/**
 * Created by snettem on 5/19/2017.
 */


public class ModelGenerationServiceTest {


    @Autowired
    private ModelGenerationService modelGenerationService;
    private Alignment alignment;

    @Before
    public void getAlignment() throws Exception{
        AlignmentGenerationServiceTest alignmentGenerationServiceTest = new AlignmentGenerationServiceTest ();
        alignmentGenerationServiceTest.generateAlignment ();
        alignment = alignmentGenerationServiceTest.alignment;
    }


    @Test
    public void generateModelTest() throws Exception{
        ModelGenerationService modelGenerationService = new ModelGenerationService ();
        VigorForm form = new VigorForm ();
      //  modelGenerationService.generateModel ( alignment,form);

    }


    @Test
    public void generateCompatibleFragsChainsTest() throws Exception{


        List<AlignmentFragment> alignmentFragments = alignment.getAlignmentFragments ();
        List<List<AlignmentFragment>> expListOfAlignmentFragsList = new ArrayList<List<AlignmentFragment>>(  );
        List<AlignmentFragment> expAlignmentFragsList = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList.add ( alignmentFragments.get ( 0 ));
        expAlignmentFragsList.add ( alignmentFragments.get ( 2 ) );
        expAlignmentFragsList.add ( alignmentFragments.get ( 3 ) );
        expListOfAlignmentFragsList.add( expAlignmentFragsList );
        List<AlignmentFragment> expAlignmentFragsList1 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList1.add ( alignmentFragments.get ( 0 ) );
        expAlignmentFragsList1.add ( alignmentFragments.get ( 3 ) );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList1);
        List<AlignmentFragment> expAlignmentFragsList2 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList2.add ( alignmentFragments.get ( 1 ) );
        expAlignmentFragsList2.add ( alignmentFragments.get ( 3 ) );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList2 );
   /*     List<AlignmentFragment> expAlignmentFragsList3 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList3.add ( alignmentFragments.get ( 2 ) );
        expAlignmentFragsList3.add ( alignmentFragments.get ( 3 ) );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList3 );
        List<AlignmentFragment> expAlignmentFragsList4 = new ArrayList<AlignmentFragment> (  );
        expAlignmentFragsList4.add ( alignmentFragments.get ( 3 ) );
        expListOfAlignmentFragsList.add ( expAlignmentFragsList4 );*/


        VigorForm form = new VigorForm ();
        form.setAlignmentTool ( "Exonerate" );
     //   assertEquals (expListOfAlignmentFragsList,modelGenerationService.generateCompatibleFragsChains ( alignmentFragments,form ));

    }
    
    
    
    


}