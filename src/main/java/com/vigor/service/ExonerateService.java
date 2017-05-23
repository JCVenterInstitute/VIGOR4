package com.vigor.service;

import com.vigor.component.AlignmentFragment;
import com.vigor.forms.VigorForm;
import com.vigor.utils.VigorUtils;
import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.datastore.DataStore;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.Nucleotide;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * Created by snettem on 5/17/2017.
 */
@Data
@Service
public class ExonerateService {

    private static final Logger LOGGER = LogManager.getLogger ( VigorInitializationService.class );
    private String proteinSequenceID;
    private ProteinSequence proteinSequence ;
    private float score;
    /*AlignmentFragment has below properties:
    private double score;
	private Direction direction;
	private Range proteinSeqRange;
	private Range nucleotideSeqRange;
    private Frame frame;
     */
    private List<AlignmentFragment> alignmentFragments;


    public void parseExonerateOutput(VigorForm form){
        try {
            ProteinFastaDataStore dataStore = new ProteinFastaFileDataStoreBuilder ( new File ( VigorUtils.getVirusDatabasePath()+File.separator+(form.getAlignmentEvidence ().getReference_db () ) ))
                    .hint ( DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED )
                    .build ();
            proteinSequence = dataStore.getSequence ( "seg7prot2ADR3_V27A&S31N" );
            proteinSequenceID = "seg7prot2ADR3_V27A&S31N";
            score = 458;
        }
        catch(IOException e){
            LOGGER.debug ( e.getMessage (),e );
        }
        catch(Exception e)
        {
            LOGGER.debug ( e.getMessage (),e );
        }
        score=233;

        AlignmentFragment alignmentFragment = new AlignmentFragment ();
       // alignmentFragment.setDirection ( );




    }

}
