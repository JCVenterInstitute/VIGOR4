package com.vigor.service;

import com.vigor.component.Alignment;
import com.vigor.component.AlignmentFragment;
import com.vigor.component.ViralProtein;
import com.vigor.utils.VigorUtils;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.junit.Before;
import org.junit.Test;
import org.springframework.core.io.ClassPathResource;
import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import static org.junit.Assert.*;

/**
 * Created by snettem on 5/19/2017.
 */
public class AlignmentGenerationServiceTest {

    public Alignment alignment;


    @Before
    public void generateAlignment() throws Exception {

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
        alignmentFragment1.setDirection ( Direction.FORWARD );
        alignmentFragmentList.add ( alignmentFragment1 );
        AlignmentFragment alignmentFragment2 = new AlignmentFragment ();
        builder.setBegin ( 645);
        builder.setEnd ( 933 );
        range = builder.build();
        alignmentFragment2.setNucleotideSeqRange ( range );
        alignmentFragment2.setDirection ( Direction.FORWARD );
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
        alignmentFragment3.setDirection ( Direction.FORWARD );
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
        alignmentFragment4.setDirection ( Direction.FORWARD );
        builder.setBegin ( 124 );
        builder.setEnd ( 326 );
        range = builder.build ();
        alignmentFragment4.setProteinSeqRange ( range );
        alignmentFragmentList.add ( alignmentFragment4 );

        Map<String,Float> alignmentcores = new HashMap<String,Float> (  );
        alignmentcores.put ( "exonerateScore",233f );
        alignment = new Alignment ();
        alignment.setAlignmentFragments ( alignmentFragmentList );

        ViralProtein viralProtein = new ViralProtein ();
               ProteinFastaDataStore datastore1 = new ProteinFastaFileDataStoreBuilder (new File ( VigorUtils.getVirusDatabasePath ()+File.separator+"flua_db" ))
                .hint( DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
                .build();
        ProteinFastaRecord record = datastore1.get ( "seg3prot2C" );
        viralProtein.setProteinID ( record.getId ());
        viralProtein.setSequence ( record.getSequence () );
        viralProtein.setDefline ( record.getComment () );
        ViralProteinService viralProteinService = new ViralProteinService ();
        viralProteinService.setViralProteinAttributes ( viralProtein );
        alignment.setViralProtein ( viralProtein );

        alignment.setAlignmentScore ( alignmentcores );



    }

    @Test
    public void chooseAlignmentTool() throws Exception {
    }

}