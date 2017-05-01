package com.vigor.service;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.stream.Stream;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import com.vigor.component.AlignmentEvidence;
import com.vigor.component.VirusGenome;
import com.vigor.forms.VigorForm;
import com.vigor.utils.LoadDefaultParameters;

@Service
public class VigorInitializationService {

	private static final Logger LOGGER = LogManager.getLogger(VigorInitializationService.class);

	public void initializeVigor(CommandLine inputs) {
		try {
            boolean isComplete = false;
            boolean isCircular = false;
            if(inputs.hasOption('C'))
            {
            	isComplete=true;
            	
            }
            if(inputs.hasOption('0'))
            {
            	isComplete=true;
            	isCircular=true;
            }
            VigorForm form = new VigorForm();
			form = loadVigorParameters(form,inputs);
			
			NucleotideFastaDataStore datastore = new NucleotideFastaFileDataStoreBuilder(new File(inputs.getOptionValue('i')))
	                                             .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
	                                             .build();
            Stream<NucleotideFastaRecord> records = datastore.records();
   		    Iterator<NucleotideFastaRecord> i = records.iterator();
   		    while(i.hasNext()){
   		    	NucleotideFastaRecord record = i.next();
   		    	VirusGenome genome = new VirusGenome(record.getSequence(),record.getComment(),isComplete,isCircular);
   		    	
   		        
   		    }

				

		} catch (Exception e) {
			e.printStackTrace();
			LOGGER.error(e.getMessage(), e);

		}

	}

	// load all the vigor parameters from Vigor.ini file
	public VigorForm loadVigorParameters(VigorForm form,CommandLine inputs) {

		HashMap<String, String> vigorParameterList = LoadDefaultParameters.loadVigorParameters();

		
		AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
        form.setAlignmentEvidence(alignmentEvidence);
		if (inputs.hasOption('A') && !(inputs.hasOption('d'))) {
			  
			  alignmentEvidence.setReference_db("virus_db");
			  form.setAlignmentEvidence(alignmentEvidence);
		} else if (inputs.hasOption('d')) {
			alignmentEvidence.setReference_db(inputs.getOptionValue('d'));
			form.setAlignmentEvidence(alignmentEvidence);
		}
		if (inputs.hasOption('s')) {
			vigorParameterList.put("min_gene_size", inputs.getOptionValue('s'));
		}
		if (inputs.hasOption('c')) {
			vigorParameterList.put("min_gene_coverage", inputs.getOptionValue('c'));
		}
		if (inputs.hasOption('f')) {
			vigorParameterList.put("frameshift_sensitivity", inputs.getOptionValue('f'));
		}
		if (inputs.hasOption('K')) {
			vigorParameterList.put("candidate_selection", inputs.getOptionValue('K'));
		}
		if (inputs.hasOption('l')) {
			vigorParameterList.put("use_locus_tags", "0");
		}
		if (inputs.hasOption('L')) {
			vigorParameterList.put("use_locus_tags", "1");
		}
		if (inputs.hasOption('m')) {
			vigorParameterList.put("min_candidate_pctsimilarity", "0");
			vigorParameterList.put("min_candidate_sbjcoverage", "0");
			vigorParameterList.put("mature_pep_mincoverage", "0");
			vigorParameterList.put("mature_pep_minsimilarity", "0");
			vigorParameterList.put("mature_pep_minidentity", "0");
			vigorParameterList.put("min_pseudogene_identity", "0");
			vigorParameterList.put("min_pseudogene_similarity", "0");
			vigorParameterList.put("min_pseudogene_coverage", "0");
		}
		if (inputs.hasOption('P')) {
			// Implement Parse parameters logic
		}
		if (inputs.hasOption('e')) {
			vigorParameterList.put("candidate_evalue", inputs.getOptionValue('e'));
		}
		if (inputs.hasOption('j')) {
			vigorParameterList.put("jcvi_rules", "0");
		}

		form.setVigorParametersList(vigorParameterList);
		return form;
	}

}
