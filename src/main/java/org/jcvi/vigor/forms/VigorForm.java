package org.jcvi.vigor.forms;

import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.utils.VigorConfiguration;
import lombok.Data;

@Data
public class VigorForm {

    private VigorConfiguration configuration;
    private AlignmentEvidence alignmentEvidence;

    // TODO remove these
    private String AlignmentOutputTempFile;
    private String tempDirectoryPath;

    public VigorForm ( VigorConfiguration configuration ) {
        this.configuration = new VigorConfiguration(configuration);
    }
}
