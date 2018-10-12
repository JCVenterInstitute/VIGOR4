package org.jcvi.vigor.service;

import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.FormatVigorOutput;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/**
 * Created by snettem on 5/9/2017.
 */
@Service
public class AlignmentGenerationService {

    private static final Logger LOGGER = LogManager.getLogger(AlignmentGenerationService.class);
    @Autowired
    private ViralProteinService viralProteinService;
    @Autowired
    private ExonerateService exonerateService;

    public List<Alignment> generateAlignment ( VirusGenome virusGenome, String referenceDB, VigorConfiguration config ) throws VigorException {
        boolean isDebug = config.getOrDefault(ConfigurationParameters.Verbose, false);
        String alignmentModule = config.get(ConfigurationParameters.AlignmentModule);
        AlignmentTool alignmentTool = AlignmentToolFactory.getAlignmentTool(alignmentModule);
        String tempDir = config.get(ConfigurationParameters.TemporaryDirectory);
        AlignmentService alignmentService = getAlignmentService(alignmentTool);
        Path workspace;
        try {
            workspace = Files.createTempDirectory(Paths.get(tempDir), "vigor4");
        } catch (IOException e) {
            throw new VigorException(String.format("Unable to create temporary directory under %s", tempDir));
        }
        List<Alignment> alignments = alignmentService.getAlignment(config, virusGenome, referenceDB, workspace.toString());
        for (int i = 0; i < alignments.size(); i++) {
            alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), config));
        }
        if (isDebug) {
            FormatVigorOutput.printAlignments(alignments);
        }
        return alignments;
    }

    /**
     * @param alignmentTool
     * @return AlignmentService with alignment algorithm
     * @throws ServiceException
     */
    private AlignmentService getAlignmentService ( AlignmentTool alignmentTool) throws ServiceException {

        if (alignmentTool != null && "exonerate".equals(alignmentTool.getToolName())) {
                return exonerateService;
        }
        throw new ServiceException(String.format("Unsupported alignment tool %s", alignmentTool));
    }
}
