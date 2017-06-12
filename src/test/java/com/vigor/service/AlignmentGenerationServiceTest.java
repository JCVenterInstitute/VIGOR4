package com.vigor.service;

import com.vigor.component.Alignment;
import com.vigor.component.AlignmentEvidence;

import com.vigor.component.VirusGenome;

import org.jcvi.jillion.core.datastore.DataStoreProviderHint;

import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.junit.Assert;

import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import java.util.stream.Stream;


/**
 * Created by snettem on 5/19/2017.
 */
public class AlignmentGenerationServiceTest {

	
	public List<Alignment> alignments=new ArrayList<Alignment>();
	   
    @Test
    public void generateAlignmentTest() throws IOException{
		 NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder (new File("C:/Users/snettem/Downloads/sequence.txt") )
                .hint ( DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED )
                .build ();
        Stream<NucleotideFastaRecord> records = dataStore.records ();
        Iterator<NucleotideFastaRecord> i = records.iterator ();
       
            NucleotideFastaRecord record = i.next ();
            VirusGenome virusGenome = new VirusGenome ( record.getSequence (), record.getComment (), false, false );
            // Call referenceDBGenerationService methods to generate the alignmentEvidence object. 
           
            AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
            alignmentEvidence.setReference_db("flua_db");
            ExonerateService exonerateService = new ExonerateService();
           List<Alignment> outputAlignments = exonerateService.getAlignment(virusGenome, alignmentEvidence);
           alignments.addAll(outputAlignments);
           Assert.assertEquals(7, outputAlignments.size());
                      


    }


 
}