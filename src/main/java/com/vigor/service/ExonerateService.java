package com.vigor.service;

import com.vigor.component.Alignment;
import com.vigor.component.AlignmentEvidence;
import com.vigor.component.AlignmentFragment;
import com.vigor.component.ViralProtein;
import com.vigor.component.VirusGenome;
import com.vigor.forms.VigorForm;
import com.vigor.utils.VigorUtils;
import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.exonerate.Exonerate;
import org.jcvi.jillion.align.exonerate.vulgar.VulgarProtein2Genome;
import org.jcvi.jillion.core.Direction;
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
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Created by snettem on 5/17/2017.
 */
@Data
@Service
public class ExonerateService {

    private static final Logger LOGGER = LogManager.getLogger ( VigorInitializationService.class );
   

    public List<Alignment> getAlignment(VirusGenome virusGenome,AlignmentEvidence alignmentEvidence)
    {
    
    	File file = queryExonerate(virusGenome,alignmentEvidence.getReference_db());
    	List<Alignment> alignments = parseExonerateOutput(file,alignmentEvidence,virusGenome);
    	
    	return alignments;
    }
       
    
    
    public File queryExonerate(VirusGenome virusGenome,String referenceDB)
    {
    	//logic to query exonerate. Just provide path to execute below command on Linux machine. 
    	//exonerate --model protein2genome -q /home/snettem/Exonerate_output/flua_db -t /home/snettem/Exonerate_output/sequence.fasta --showquerygff true --showtargetgff true
    	
    	File file = new File("C:/git/VIGOR4/VigorWorkSpace/Output_ExonGaps.txt");
    	
    	return file;
    }

    public List<Alignment> parseExonerateOutput(File file,AlignmentEvidence alignmentEvidence,VirusGenome virusGenome)
    {
    	List<Alignment> alignments = new ArrayList<Alignment>();
    	try
    	{
    		
    		List<VulgarProtein2Genome> output = Exonerate.parseVulgarOutput(file);
    		for(int i=0;i<output.size();i++){
    			VulgarProtein2Genome record = output.get(i);
    			int iter =0;
    		    int iter1=0;
    			List<AlignmentFragment> alignmentFragments = new ArrayList<AlignmentFragment>();
    			Alignment alignment = new Alignment();
    			for(int j=0;j<record.getQueryRanges().size();j++){
    			AlignmentFragment alignmentFragment = new AlignmentFragment();
    			alignmentFragment.setDirection(record.getTargetStrand().get());
    			if(record.getQueryRanges().get(iter1).getBegin()>record.getQueryRanges().get(iter1).getEnd()){
    				iter1++;
    				continue;
    				    				
    			}
    			else{
    			alignmentFragment.setNucleotideSeqRange(record.getTargetExons().get(iter));
    			
    			alignmentFragment.setProteinSeqRange(record.getQueryRanges().get(iter1));
    			//frame should be added to alignmentFragment
    			alignmentFragment.setDirection(record.getTargetStrand().get());
    			alignmentFragments.add(alignmentFragment);
    			iter++;
    			iter1++;
    			}
    			}
    			alignment.setAlignmentFragments(alignmentFragments);
    		    Map<String,Float> alignmentScores = new HashMap<String,Float>();
    		    alignmentScores.put("ExonerateScore",record.getScore());
    		    alignment.setAlignmentScore(alignmentScores);
    		    alignment.setAlignmentTool_name("Exonerate");
    		  
    		    ProteinFastaDataStore datastore = new ProteinFastaFileDataStoreBuilder(new File (VigorUtils.getVirusDatabasePath()+File.separator+alignmentEvidence.getReference_db()))
    		            .hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
    		            .build();
    		    ProteinFastaRecord fasta = datastore.get(record.getQueryId());
    		    ViralProtein viralProtein = new ViralProtein();
    		    viralProtein.setProteinID(fasta.getId());
    		    viralProtein.setDefline(fasta.getComment());
    		    viralProtein.setSequence(fasta.getSequence());
    		    alignment.setViralProtein(viralProtein);
    		    alignment.setVirusGenome(virusGenome);
    		    alignment.setAlignmentEvidence(alignmentEvidence);
    		    alignments.add(alignment);			
    		}
          
        }
       catch(Exception e)
        {
            LOGGER.debug ( e.getMessage (),e );
        }
    	return alignments;      
        
    }

}
