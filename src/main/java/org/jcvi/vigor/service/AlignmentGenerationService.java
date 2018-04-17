package org.jcvi.vigor.service;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.FormatVigorOutput;

import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/**
 * Created by snettem on 5/9/2017.
 */

@Service
public class AlignmentGenerationService {

	private static final Logger LOGGER = LogManager.getLogger(AlignmentGenerationService.class);
	private boolean isDebug = false;

	@Autowired
	private ViralProteinService viralProteinService;
	
	@Autowired
	private ExonerateService exonerateService;

	public List<Alignment> generateAlignment(VirusGenome virusGenome, VigorForm form) throws ServiceException {
		isDebug = form.isDebug();
		AlignmentEvidence alignmentEvidence = form.getAlignmentEvidence();
		String alignmentTool = chooseAlignmentTool(alignmentEvidence);
		String min_gap_length = form.getVigorParametersList().get("min_seq_gap_length");
		String exoneratePath = form.getVigorParametersList().get("exonerate_path");
		String workspace = form.getVigorParametersList().get("output_directory");
		String referenceDB = Paths.get(form.getVigorParametersList().get("reference_database_path"), alignmentEvidence.getReference_db()).toString();

		List<Range> sequenceGaps = VirusGenomeService.findSequenceGapRanges(min_gap_length,
				virusGenome.getSequence());
		Map<Frame,List<Long>> internalStops =VirusGenomeService.findInternalStops(virusGenome.getSequence());
		virusGenome.setInternalStops(internalStops);
		virusGenome.setSequenceGaps(sequenceGaps);
		List<Alignment> alignments = Collections.EMPTY_LIST;
		AlignmentService alignmentService = getAlignmentService(alignmentTool, form);
		alignments = alignmentService.getAlignment(form.getAlignmentEvidence(), virusGenome, referenceDB, workspace);
		for (int i=0; i < alignments.size(); i++ ) {
			alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), form));
		}
		if (isDebug) {
			FormatVigorOutput.printAlignments(alignments);
		}

		return alignments;
	}
	
	// TODO why does this take alignmentEvidence rather than the configuration?
	public String chooseAlignmentTool(AlignmentEvidence alignmentEvidence) {
		return "exonerate";
	}

	private AlignmentService getAlignmentService(String alignmentTool, VigorForm form) throws ServiceException {
		if ("exonerate".equals(alignmentTool)) {
			try {
				exonerateService.setExoneratePath(Paths.get(form.getVigorParametersList().get("exonerate_path")));
				return exonerateService;
			} catch (VigorException e) {
				throw new ServiceException(e);
			}
		}
		throw new ServiceException(String.format("Unsupported alignment tool %s", alignmentTool));
	}

}
