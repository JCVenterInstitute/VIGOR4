package org.jcvi.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import java.util.Collections;
import java.util.List;

/**
 * Created by snettem on 5/8/2017.
 */
@Component
@Scope("prototype")
@Data
public class StartTranslationException {

    private final boolean hasStartTranslationException;
    private final List<String> alternateStartCodons;

    public StartTranslationException ( boolean hasStartTranslationException, List<String> alternateStartCodons ) {

        this.hasStartTranslationException = hasStartTranslationException;
        this.alternateStartCodons = alternateStartCodons;
    }

    public static StartTranslationException NO_EXCEPTION = new StartTranslationException(false, Collections.EMPTY_LIST);
}
