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
public class StopTranslationException {

    private final boolean hasStopTranslationException;
    private final AminoAcid replacementAA;
    private final String motif;
    private final int offset;

    public StopTranslationException ( boolean hasStopTranslationException, AminoAcid replacementAA, String motif, int offset ) {

        this.hasStopTranslationException = hasStopTranslationException;
        this.replacementAA = replacementAA;
        this.motif = motif;
        this.offset = offset;
    }

    public static StopTranslationException NO_EXCEPTION = new StopTranslationException(false, AminoAcid.Unknown_Amino_Acid, "", 0);
}
