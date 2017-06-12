package com.vigor.service;


import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.vigor.forms.*;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.Rangeable;
import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.junit.Test;
import org.springframework.beans.factory.annotation.Autowired;

public class VirusGenomeServiceTest {
	
	
	private VirusGenomeService virusGenomeService= new VirusGenomeService();
	
	@Test
	public void findGapsTest(){
	
	/*		System.out.println("Entered test");
			File file = new File("C:"+File.separator+"Users"+File.separator+"snettem"+File.separator+"Desktop"+File.separator+"sequence.fasta");
		NucleotideFastaDataStore datastore;
		try {
			datastore = new NucleotideFastaFileDataStoreBuilder(file)
			.hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
			.build();
			 NucleotideFastaRecord fasta = datastore.get("KT388711.1");
			 VigorForm form = new VigorForm();
			 Map<String,String> list = new HashMap<String,String>();
			 list.put("min_gap_length", "20");
			 form.setVigorParametersList(list);
			 NucleotideSequence sequence = fasta.getSequence();
			 List<Range> ranges = virusGenomeService.findSequenceGapRanges(form,sequence);
			 for(int i=0;i<ranges.size();i++){
				 System.out.println(ranges.get(i));
			 }
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
	
	
		
	}

}