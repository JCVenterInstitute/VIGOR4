package org.jcvi.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

@Component
@Scope("prototype")
@Data
@SuppressWarnings("serial")
public class MaturePeptide {
    private Alignment alignment;
    private boolean partial5p=false;
    private boolean partial3p=false;
}
