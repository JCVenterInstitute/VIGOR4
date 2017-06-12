package com.vigor.service;

import com.vigor.component.*;
import com.vigor.forms.VigorForm;
import com.vigor.utils.VigorUtils;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;



/**
 * Created by snettem on 5/9/2017.
 */
@Service
public class AlignmentGenerationService {

    private static final Logger LOGGER = LogManager.getLogger ( AlignmentGenerationService.class );

    @Autowired
    private ExonerateService exonerateService;
    @Autowired
    private ModelGenerationService modelGenerationService;

    public void GenerateAlignment(VirusGenome virusGenome,VigorForm form) {
        String alignmentTool = chooseAlignmentTool ( form.getAlignmentEvidence () );
        //create vigor workspace and generate the exonerate output file in that workspace
        // VigorUtils.getVigorWorkSpace();
        if (alignmentTool.equals ( "exonerate" )) {
        List<Alignment> alignments = exonerateService.getAlignment(virusGenome, form.getAlignmentEvidence());
        modelGenerationService.generateModels(alignments, form);
        }
        
    }

    public String chooseAlignmentTool(AlignmentEvidence alignmentEvidence) {

        return "exonerate";
    }



}
