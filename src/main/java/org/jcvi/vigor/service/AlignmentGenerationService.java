package org.jcvi.vigor.service;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.FormatVigorOutput;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
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
	private ExonerateService exonerateService;

	@Autowired
	private ViralProteinService viralProteinService;
	

	public List<Alignment> GenerateAlignment(VirusGenome virusGenome,VigorForm form) throws ServiceException {
		isDebug = form.isDebug();
		String alignmentTool = chooseAlignmentTool(form.getAlignmentEvidence());
		String min_gap_length = form.getVigorParametersList().get("min_seq_gap_length");
		String exoneratePath = form.getVigorParametersList().get("exonerate_path");
		String workspace = form.getVigorParametersList().get("output_directory");
		List<Range> sequenceGaps = VirusGenomeService.findSequenceGapRanges(min_gap_length,
				virusGenome.getSequence());
		Map<Frame,List<Long>> internalStops =VirusGenomeService.findInternalStops(virusGenome.getSequence());
		virusGenome.setInternalStops(internalStops);
		virusGenome.setSequenceGaps(sequenceGaps);
		List<Alignment> alignments = Collections.EMPTY_LIST;
		if (alignmentTool.equals("exonerate")) {
			alignments = generateExonerateAlignment(virusGenome, form.getAlignmentEvidence(),workspace,exoneratePath);
			for (int i=0; i < alignments.size(); i++ ) {
				alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), form));
			}
			if (isDebug) {
				FormatVigorOutput.printAlignments(alignments);
			}
		}
		// TODO if alignment service isn't available
		return alignments;
	}
	
	public List<Alignment> generateExonerateAlignment(VirusGenome virusGenome, AlignmentEvidence alignmentEvidence,String workspace,String exoneratePath) throws ServiceException {
	             return exonerateService.getAlignment(virusGenome, alignmentEvidence,workspace,exoneratePath);
	}

	public String chooseAlignmentTool(AlignmentEvidence alignmentEvidence) {
		return "exonerate";
	}

}
