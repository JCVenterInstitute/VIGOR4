package org.jcvi.vigor.service;

import org.jcvi.vigor.component.AlignmentTool;
import org.jcvi.vigor.component.Exonerate;

public class AlignmentToolFactory {

    public static AlignmentTool getAlignmentTool ( String alignmentModule ) {

        if (alignmentModule.equalsIgnoreCase("exonerate")) {
            return new Exonerate("exonerate");
        }
        return null;
    }
}
