package org.jcvi.vigor.service;

import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.exception.VigorException;
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
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by snettem on 5/17/2017.
 */

@Service
public class ExonerateService implements AlignmentService {

	private static final Logger LOGGER = LogManager.getLogger(ExonerateService.class);

	private Path exoneratePath;

	public void setExoneratePath(Path exoneratePath) throws VigorException {
		if (! (exoneratePath.toFile().exists() && exoneratePath.toFile().canExecute()) ) {
			LOGGER.warn("exonerate path {} does not exist or is not executable", exoneratePath);
			throw new VigorException(String.format("exonerate path %s does not exist or is not executable", exoneratePath));
		}
		this.exoneratePath = exoneratePath;
	}

	@Override
	public List<Alignment> getAlignment(AlignmentEvidence alignmentEvidence, VirusGenome virusGenome, String referenceDB, String workspace) throws ServiceException {

		try {
			String outputFilePath = GenerateExonerateOutput.queryExonerate(virusGenome, referenceDB, workspace, null, exoneratePath.toString());
			File outputFile = new File(outputFilePath);
			return parseExonerateOutput(outputFile, alignmentEvidence, virusGenome, referenceDB);
		} catch (VigorException e ) {
			throw new ServiceException(String.format("error getting alignment got %s: %s", e.getClass().getSimpleName(), e.getMessage()), e);
		}

	}

	public List<Alignment> parseExonerateOutput(File exonerateOutput, AlignmentEvidence alignmentEvidence,
			VirusGenome virusGenome, String referenceDB) throws ServiceException{
		List<Alignment> alignments = new ArrayList<Alignment>();
		List<VulgarProtein2Genome2> Jalignments;
		try {
			Jalignments = Exonerate2.parseVulgarOutput(exonerateOutput);
		} catch (IOException e) {
			throw new ServiceException(String.format("Error parsing exonerate output %s", exonerateOutput.getName()));
		}
		try (ProteinFastaDataStore datastore = new ProteinFastaFileDataStoreBuilder(new File(referenceDB))
				.hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED).build();
		) {
			for (VulgarProtein2Genome2 Jalignment: Jalignments) {
				Alignment alignment = new Alignment();
				Map<String, Double> alignmentScores = new HashMap<String, Double>();
				alignmentScores.put("exonerateScore", (double) Jalignment.getScore());
				alignment.setAlignmentScore(alignmentScores);
				alignment.setAlignmentTool_name("exonerate");
				List<AlignmentFragment> alignmentFragments = new ArrayList<>();

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
				alignment.setAlignmentEvidence(alignmentEvidence.copy());
				alignments.add(alignment);
			}
		} catch (IOException e) {
			LOGGER.error(String.format("Problem reading virus database file %s", referenceDB), e);
			throw new ServiceException(e);
		}
		return alignments;
	}
}
