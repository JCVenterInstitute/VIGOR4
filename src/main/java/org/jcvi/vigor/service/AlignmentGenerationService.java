package org.jcvi.vigor.service;

import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.FormatVigorOutput;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.utils.VigorConfiguration;
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

    public List<Alignment> generateAlignment ( VirusGenome virusGenome, VigorForm form ) throws VigorException {

        isDebug = form.getConfiguration().get(ConfigurationParameters.Verbose).equals("true") ? true : false;
        AlignmentEvidence alignmentEvidence = form.getAlignmentEvidence();
        String alignmentModule = form.getConfiguration().get(ConfigurationParameters.AlignmentModule);
        AlignmentTool alignmentTool = AlignmentToolFactory.getAlignmentTool(alignmentModule);
        form.setAlignmentTool(alignmentTool);
        VigorConfiguration vigorConfig = form.getConfiguration();
        String min_gap_length = vigorConfig.get(ConfigurationParameters.SequenceGapMinimumLength);
        String tempDir = vigorConfig.get(ConfigurationParameters.TemporaryDirectory);
        String referenceDB = alignmentEvidence.getReference_db();
        List<Range> sequenceGaps = VirusGenomeService.findSequenceGapRanges(min_gap_length,
                virusGenome.getSequence());
        Map<Frame, List<Long>> internalStops = VirusGenomeService.findInternalStops(virusGenome.getSequence());
        virusGenome.setInternalStops(internalStops);
        virusGenome.setSequenceGaps(sequenceGaps);
        List<Alignment> alignments;
        AlignmentService alignmentService = getAlignmentService(alignmentTool, form);
        Path workspace;
        try {
            workspace = Files.createTempDirectory(Paths.get(tempDir), "vigor4");
            form.setTempDirectoryPath(workspace.toString());
        } catch (IOException e) {
            throw new VigorException(String.format("Unable to create temporary directory under %s", tempDir));
        }
        alignments = alignmentService.getAlignment(form, virusGenome, referenceDB, workspace.toString());
        for (int i = 0; i < alignments.size(); i++) {
            alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), form));
        }
        if (isDebug) {
            FormatVigorOutput.printAlignments(alignments);
        }
        return alignments;
    }

    /**
     * @param alignmentTool
     * @param form
     * @return AlignmentService with alignment algorithm
     * @throws ServiceException
     */
    private AlignmentService getAlignmentService ( AlignmentTool alignmentTool, VigorForm form ) throws ServiceException {

        if (alignmentTool != null && "exonerate".equals(alignmentTool.getToolName())) {
            try {
                String exoneratePath = form.getConfiguration().get(ConfigurationParameters.ExoneratePath);
                LOGGER.debug("Using exonerate path {}", exoneratePath);
                if (exoneratePath == null || exoneratePath.isEmpty()) {
                    throw new VigorException("Exonerate path is not set");
                }
                exonerateService.setExoneratePath(Paths.get(exoneratePath));
                return exonerateService;
            } catch (VigorException e) {
                throw new ServiceException(e);
            }
        }
        throw new ServiceException(String.format("Unsupported alignment tool %s", alignmentTool));
    }
}
