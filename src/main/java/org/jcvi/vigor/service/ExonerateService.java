package org.jcvi.vigor.service;

import org.jcvi.jillion.core.datastore.DataStoreException;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.GenerateExonerateOutput;
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
import java.io.IOException;
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

	public List<Alignment> getAlignment(VirusGenome virusGenome, AlignmentEvidence alignmentEvidence, String workspace,String exoneratePath) throws ServiceException {

		String referenceDB= getClass().getClassLoader().getResource(alignmentEvidence.getReference_db()).getPath();
		String outputFilePath = GenerateExonerateOutput.queryExonerate(virusGenome,referenceDB,workspace,null,exoneratePath);
		File outputFile = new File(outputFilePath);
		return parseExonerateOutput(outputFile, alignmentEvidence, virusGenome);

	}

	public List<Alignment> parseExonerateOutput(File file, AlignmentEvidence alignmentEvidence,
			VirusGenome virusGenome) throws ServiceException{
		List<Alignment> alignments = new ArrayList<Alignment>();
		List<VulgarProtein2Genome2> Jalignments;
		try {
			Jalignments = Exonerate2.parseVulgarOutput(file);
		} catch (IOException e) {
			throw new ServiceException(String.format("Error parsing exonerate output %s", file.getName()));
		}
		String dbFilePath = this.getClass().getClassLoader().getResource(alignmentEvidence.getReference_db()).getPath();
		try (ProteinFastaDataStore datastore = new ProteinFastaFileDataStoreBuilder(new File(dbFilePath))
				.hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED).build();
		) {
			for (VulgarProtein2Genome2 Jalignment: Jalignments) {
				Alignment alignment = new Alignment();
				Map<String, Double> alignmentScores = new HashMap<String, Double>();
				alignmentScores.put("exonerateScore", (double) Jalignment.getScore());
				alignment.setAlignmentScore(alignmentScores);
				alignment.setAlignmentTool_name("exonerate");
				alignment.setAlignmentEvidence(alignmentEvidence);
				List<AlignmentFragment> alignmentFragments = new ArrayList<AlignmentFragment>();

				for (VulgarProtein2Genome2.AlignmentFragment fragment: Jalignment.getAlignmentFragments()) {
					AlignmentFragment alignmentFragment = new AlignmentFragment();
					alignmentFragment.setDirection(fragment.getDirection());
					alignmentFragment.setFrame(fragment.getFrame());
					alignmentFragment.setNucleotideSeqRange(fragment.getNucleotideSeqRange());
					alignmentFragment.setProteinSeqRange(fragment.getProteinSeqRange());
					alignmentFragments.add(alignmentFragment);
				}
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
		} catch (IOException e) {
			LOGGER.error(String.format("Problem reading virus database file %s", dbFilePath), e);
			throw new ServiceException(e);
		}
		return alignments;
	}
}
