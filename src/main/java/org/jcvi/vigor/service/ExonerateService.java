package org.jcvi.vigor.service;

import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.utils.GenerateExonerateOutput;
import org.jcvi.vigor.utils.VigorUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.align.exonerate.Exonerate2;
import org.jcvi.jillion.align.exonerate.vulgar.VulgarProtein2Genome2;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.aa.ProteinFastaDataStore;
import org.jcvi.jillion.fasta.aa.ProteinFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.aa.ProteinFastaRecord;
import org.springframework.stereotype.Service;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by snettem on 5/17/2017.
 */

@Service
public class ExonerateService {

	private static final Logger LOGGER = LogManager.getLogger(ExonerateService.class);

	

	public List<Alignment> getAlignment(VirusGenome virusGenome, AlignmentEvidence alignmentEvidence) {
       
		String outputFilePath = GenerateExonerateOutput.queryExonerate(virusGenome,alignmentEvidence.getReference_db(), VigorUtils.getVigorWorkSpace(),null);
		File outputFile = new File(outputFilePath);
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
				Map<String, Double> alignmentScores = new HashMap<String, Double>();
				alignmentScores.put("ExonerateScore",(double)Jalignment.getScore());
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
