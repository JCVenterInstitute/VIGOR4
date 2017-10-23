package org.jcvi.vigor.component;

import lombok.Data;

import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;


/**
 * Created by snettem on 5/8/2017.
 */

@Component
@Scope("prototype")
@Data

public class StopTranslationException{

    private boolean hasStopTranslationException=false;
    private AminoAcid replacementAA;

}
