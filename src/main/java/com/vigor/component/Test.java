/*package com.vigor.component;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.stream.Stream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.SequenceBuilder;
import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.core.residue.nt.Nucleotide;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;

import com.vigor.forms.VigorForm;


public class Test {

	private static final Logger LOGGER = LogManager.getLogger(Test.class);
	VirusGenome genome1 = new VirusGenome(null, null, false, false);
	VigorForm form = new VigorForm();
	
	public void test()
	{
		try{
		NucleotideFastaDataStore datastore = new NucleotideFastaFileDataStoreBuilder(new File("C:/Users/snettem/Downloads/sequence.fasta"))
                .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
.build();
		 NucleotideFastaRecord fasta = datastore.get("KT388711.1");
		 System.out.println("Fasta comment is"+ fasta.getComment());
		 Sequence<Nucleotide> seq = fasta.getSequence();
		 genome1.setSequence(fasta.getSequence());
		Range range;
		 Builder builder= new Builder();
		 builder.setBegin(0);
		 builder.setEnd(genome1.getSequence().getLength()-1);
		 range = builder.build();
		 datastore.getNumberOfRecords();
		 
		 Stream<NucleotideFastaRecord> records = datastore.records();
		 Iterator i = records.iterator();
			
		 
		System.out.println("The sequence is "+genome1.getSequence());
		SequenceBuilder sq= genome1.getSequence().toBuilder();
		
		
		//System.out.println();
		ProteinFastaDataStore datastore1 = new ProteinFastaFileDataStoreBuilder(new File("C:/Users/snettem/Downloads/cdv_db"))
                .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
.build();
		ProteinFastaRecord pro = datastore1.get("AAT94550.1");
		System.out.println("protein sequence is "+pro.getSequence());
		
		System.out.println("protein sequence is "+pro.getComment());
		
		}
		catch(IOException e)
		{
			LOGGER.error(e.getMessage(),e);
			LOGGER.debug(e.getMessage(),e);
		}
		catch (DataStoreException e)
		{
			LOGGER.error(e.getMessage(),e);
			LOGGER.debug(e.getMessage(),e);
		}
		catch(Exception e)
		{
			LOGGER.error(e.getMessage(),e);
		}
	}
	
	public void test1()
	{
	//	genome.getSequence().getLength();
		//genome.getSequence().
	//	System.out.println("Length is "+genome.getSequence().getLength());
		
	}
	
	
}
*/