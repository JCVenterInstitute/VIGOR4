package org.jcvi.vigor.service;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.FormatVigorOutput;
import java.util.List;
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
	private ModelGenerationService modelGenerationService;
	@Autowired
	private ViralProteinService viralProteinService;
	

	public void GenerateAlignment(VirusGenome virusGenome, VigorForm form) {
		isDebug = form.isDebug();
		try {
			String alignmentTool = chooseAlignmentTool(form.getAlignmentEvidence());
			String min_gap_length = form.getVigorParametersList().get("min_gap_length");
			List<Range> sequenceGaps = VirusGenomeService.findSequenceGapRanges(min_gap_length,
					virusGenome.getSequence());
			virusGenome.setSequenceGaps(sequenceGaps);
			if (alignmentTool.equals("exonerate")) {
				List<Alignment> alignments = exonerateService.getAlignment(virusGenome, form.getAlignmentEvidence());
				alignments = alignments.stream()
						.map(alignment -> viralProteinService.setViralProteinAttributes(alignment))
						.collect(Collectors.toList());
				if (isDebug) {
					System.out.println("*******Initial List of alignments*****");
					FormatVigorOutput.printAlignments(alignments);
					
				}
							
				modelGenerationService.generateModels(alignments, form);
			}
		} catch (Exception e) {
			LOGGER.error(e.getMessage(), e);
		}

	}

	public String chooseAlignmentTool(AlignmentEvidence alignmentEvidence) {

		return "exonerate";
	}

}
