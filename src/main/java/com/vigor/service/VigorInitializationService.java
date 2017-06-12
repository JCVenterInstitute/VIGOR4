package com.vigor.service;

import java.io.File;

import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.sun.javafx.collections.MappingChange;
import com.vigor.component.Alignment;
import com.vigor.utils.VigorUtils;
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

/** This class is under service layer. The methods in this class has functionality to
 * initialize the vigor application
 */
@Service
public class VigorInitializationService {

    private static final Logger LOGGER = LogManager.getLogger ( VigorInitializationService.class );

    @Autowired
    private ReferenceDBGenerationService referenceProteinGenerationService;
    @Autowired
    private AlignmentGenerationService alignmentGenerationService;

   
    /**
     * @param inputs: User provided command line inputs
     *                Retrieve each Genomic sequence from the input file and determine AlignmentEvidence
     */

    public void initializeVigor(CommandLine inputs) {
        try {
            boolean isComplete = false;
            boolean isCircular = false;
            if (inputs.hasOption ( 'C' )) {
                isComplete = true;

            }
            if (inputs.hasOption ( '0' )) {
                isComplete = true;
                isCircular = true;
            }
            VigorForm form = new VigorForm();
            form = loadDefaultParameters ( inputs, form );
            NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder ( new File ( inputs.getOptionValue ( 'i' ) ) )
                    .hint ( DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED )
                    .build ();
            Stream<NucleotideFastaRecord> records = dataStore.records ();
            Iterator<NucleotideFastaRecord> i = records.iterator ();
            while (i.hasNext ()) {
                NucleotideFastaRecord record = i.next ();
                VirusGenome virusGenome = new VirusGenome ( record.getSequence (), record.getComment (), isComplete, isCircular );
                // Call referenceDBGenerationService methods to generate the alignmentEvidence object. 
                alignmentGenerationService.GenerateAlignment (virusGenome,form);


            }

        } catch (Exception e) {
            e.printStackTrace ();
            LOGGER.error ( e.getMessage (), e );

        }

    }

    /**
     * load all the vigor parameters from Vigor.ini file
     *
     * @param inputs: Input parameters and values provided by user
     * @return form : output form object has the AlignmentEvidence object and the VigorParametersList.
     * Few vigor parameters will be overridden by the default parameters of vigor.ini file and saved to
     * VigorParametersList attribute of the form.
     */
    public VigorForm loadDefaultParameters(CommandLine inputs,VigorForm form) {

        Map<String, String> vigorParameterList = LoadDefaultParameters.loadVigorParameters ( VigorUtils.getVigorParametersPath ());
        AlignmentEvidence alignmentEvidence = new AlignmentEvidence ();
        form = new VigorForm ();
        if (inputs.hasOption ( 'A' ) && !(inputs.hasOption ( 'd' ))) {

            alignmentEvidence.setReference_db ( vigorParameterList.get ( "reference_db" ) );
             // here call the method from ReferenceDBGenerationService which determines the reference_db
            // and set the alignment evidence

        } else if (inputs.hasOption ( 'd' )) {
            System.out.println("Reference_db is "+inputs.getOptionValue ( "d"));
            alignmentEvidence.setReference_db ( inputs.getOptionValue ( 'd' ) );

        }

        vigorParameterList = loadVirusSpecificParameters (vigorParameterList,alignmentEvidence.getReference_db () );

        if (inputs.hasOption ( 's' )) {
            vigorParameterList.put ( "min_gene_size", inputs.getOptionValue ( 's' ) );
        }
        if (inputs.hasOption ( 'c' )) {
            vigorParameterList.put ( "min_gene_coverage", inputs.getOptionValue ( 'c' ) );
        }
        if (inputs.hasOption ( 'f' )) {
            vigorParameterList.put ( "frameshift_sensitivity", inputs.getOptionValue ( 'f' ) );
        }
        if (inputs.hasOption ( 'K' )) {
            vigorParameterList.put ( "candidate_selection", inputs.getOptionValue ( 'K' ) );
        }
        if (inputs.hasOption ( 'l' )) {
            vigorParameterList.put ( "use_locus_tags", "0" );
        }
        if (inputs.hasOption ( 'L' )) {
            vigorParameterList.put ( "use_locus_tags", "1" );
        }
        if (inputs.hasOption ( 'm' )) {
            vigorParameterList.put ( "min_candidate_pctsimilarity", "0" );
            vigorParameterList.put ( "min_candidate_sbjcoverage", "0" );
            vigorParameterList.put ( "mature_pep_mincoverage", "0" );
            vigorParameterList.put ( "mature_pep_minsimilarity", "0" );
            vigorParameterList.put ( "mature_pep_minidentity", "0" );
            vigorParameterList.put ( "min_pseudogene_identity", "0" );
            vigorParameterList.put ( "min_pseudogene_similarity", "0" );
            vigorParameterList.put ( "min_pseudogene_coverage", "0" );
        }

        if (inputs.hasOption ( 'e' )) {
            vigorParameterList.put ( "candidate_evalue", inputs.getOptionValue ( 'e' ) );
        }
        if (inputs.hasOption ( 'j' )) {
            vigorParameterList.put ( "jcvi_rules", "0" );
        }
        if (inputs.hasOption ( 'P' )) {
           Map<String,String> temp = Pattern.compile ("~~").splitAsStream ( inputs.getOptionValue ( 'P' ).trim () )
                                            .map( s -> s.split ( "=", 2 ) )
                                            .collect ( Collectors.toMap ( a -> a[0], a -> a.length > 1 ? a[1] : "" ) );
            for ( String key : temp.keySet() ) {
                if(vigorParameterList.containsKey ( key )){
                    vigorParameterList.put ( key,temp.get ( key ) );
                }
            }
        }
        form.setVigorParametersList ( vigorParameterList );
        form.setAlignmentEvidence(alignmentEvidence);
        return form;
    }

    /**
     *
     * @param vigorParametersList: Default vigor parameters from vigor.ini file
     * @param reference_db : Virus Specific reference_db
     * @return vigorParametersList : Default vigor Parameters will be overridden by virus specific parameters
     */
    public Map<String, String> loadVirusSpecificParameters(Map<String,String> vigorParametersList,String reference_db){
        Map<String, String> virusSpecificParameters = LoadDefaultParameters.loadVigorParameters
                (VigorUtils.getVirusSpecificParametersPath ()+ File.separator+ reference_db+".ini" );
        vigorParametersList.putAll ( virusSpecificParameters );
        return vigorParametersList;
    }
}