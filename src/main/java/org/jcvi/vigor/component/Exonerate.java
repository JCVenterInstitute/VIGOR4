package org.jcvi.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

@Component
@Scope("prototype")
@Data
public class Exonerate extends AlignmentTool {

    private String name;

    public Exonerate ( String name ) {

        this.name = name;
    }

    @Override
    public String getToolName () {

        return this.name;
    }
}
