package org.jcvi.vigor.service;

import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorUtils;

import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.exonerate.Exonerate2;
import org.jcvi.jillion.align.exonerate.vulgar.VulgarProtein2Genome2;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import com.jcraft.jsch.Channel;
import com.jcraft.jsch.ChannelExec;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.Session;

/**
 * Created by snettem on 5/17/2017.
 */
@Data
@Service
public class ExonerateService {

	private static final Logger LOGGER = LogManager.getLogger(ExonerateService.class);

	@Autowired
	private ViralProteinService viralProteinService;
	
	

	public List<Alignment> getAlignment(VirusGenome virusGenome, AlignmentEvidence alignmentEvidence) {

		/*File file = queryExonerate(virusGenome, alignmentEvidence.getReference_db());
		System.out.println("File created in folder");
		if (file.exists()) {
			System.out.println("File exists");
		} else {
			System.out.println("File does not exist");
		}*/
		File outputFile = new File(VigorUtils.getVigorWorkSpace() + "/exonerate.txt");
		
		List<Alignment> alignments = parseExonerateOutput(outputFile, alignmentEvidence, virusGenome);
		return alignments;
	}

	

	public List<Alignment> parseExonerateOutput(File file, AlignmentEvidence alignmentEvidence,
			VirusGenome virusGenome) {
		List<Alignment> alignments = new ArrayList<Alignment>();
		try {

			List<VulgarProtein2Genome2> Jalignments = Exonerate2.parseVulgarOutput(file);
			for (int i = 0; i < Jalignments.size(); i++) {
				VulgarProtein2Genome2 Jalignment = Jalignments.get(i);
				Alignment alignment = new Alignment();
				Map<String, Float> alignmentScores = new HashMap<String, Float>();
				alignmentScores.put("ExonerateScore", Jalignment.getScore());
				alignment.setAlignmentScore(alignmentScores);
				alignment.setAlignmentTool_name("exonerate");
				alignment.setAlignmentEvidence(alignmentEvidence);
				List<AlignmentFragment> alignmentFragments = new ArrayList<AlignmentFragment>();
				
				for (int j = 0; j < Jalignment.getAlignmentFragments().size(); j++) {
					AlignmentFragment alignmentFragment = new AlignmentFragment();
					alignmentFragment.setDirection(Jalignment.getAlignmentFragments().get(j).getDirection());
					alignmentFragment.setFrame(Jalignment.getAlignmentFragments().get(j).getFrame());
					alignmentFragment
							.setNucleotideSeqRange(Jalignment.getAlignmentFragments().get(j).getNucleotideSeqRange());
					alignmentFragment
							.setProteinSeqRange(Jalignment.getAlignmentFragments().get(j).getProteinSeqRange());
					alignmentFragments.add(alignmentFragment);
				}
				ProteinFastaDataStore datastore = new ProteinFastaFileDataStoreBuilder(new File(
						VigorUtils.getVirusDatabasePath() + File.separator + alignmentEvidence.getReference_db()))
								.hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED).build();
				ProteinFastaRecord fasta = datastore.get(Jalignment.getQueryId());
				ViralProtein viralProtein = new ViralProtein();
				viralProtein.setProteinID(fasta.getId());
				viralProtein.setDefline(fasta.getComment());
				viralProtein.setSequence(fasta.getSequence());
				alignment.setAlignmentFragments(alignmentFragments);
				alignment.setViralProtein(viralProtein);
				alignment.setVirusGenome(virusGenome);
				alignment.setAlignmentEvidence(alignmentEvidence);
				alignments.add(alignment);
			}

		} catch (Exception e) {
			LOGGER.debug(e.getMessage(), e);
		}
		return alignments;

	}
	
	
	public void testing()
	{
		
	}

}
