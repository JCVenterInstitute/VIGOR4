package org.jcvi.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

@Component
@Scope("prototype")
@Data
public abstract class AlignmentTool {

    public abstract String getToolName();

    @Override
    public String toString() {
        return this.getToolName();
    }

}
