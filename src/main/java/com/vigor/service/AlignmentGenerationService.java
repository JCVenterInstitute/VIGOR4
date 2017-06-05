package com.vigor.service;

import com.vigor.component.*;
import com.vigor.forms.VigorForm;
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
    private ViralProteinService viralProteinService;
    @Autowired
    private ModelGenerationService modelGenerationService;
    @Autowired
    private ExonerateService exonerateService;

    public void GenerateAlignment(VirusGenome virusGenome,VigorForm form) {
        String alignmentTool = chooseAlignmentTool ( form.getAlignmentEvidence () );
        Alignment alignment = new Alignment ();
        if (alignmentTool.equals ( "jillion" )) {
            //Query jillion for exonerate output object. Input to Jillion will be VirusGenome, reference Db and candidate evalue;
           ExonerateService exonerateService = new ExonerateService ();
           exonerateService.parseExonerateOutput (form);


        }
       //modelGenerationService.generateModel (alignment,viralProtein);






    }



    public String chooseAlignmentTool(AlignmentEvidence alignmentEvidence) {

        return "jillion";
    }



}
