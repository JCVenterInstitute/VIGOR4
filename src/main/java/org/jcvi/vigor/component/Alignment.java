package org.jcvi.vigor.component;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.core.Direction;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

@Component
@Scope("prototype")
@Data
@SuppressWarnings("serial")
public class Alignment implements Serializable {

    private transient AlignmentTool alignmentTool;
    private List<AlignmentFragment> alignmentFragments = Collections.EMPTY_LIST;
    private Map<String, Double> alignmentScore = Collections.EMPTY_MAP;
    private transient VirusGenome virusGenome;
    private transient ViralProtein viralProtein;
    private transient AlignmentEvidence alignmentEvidence;

    // TODO fix this
    public Direction getDirection() {
        return alignmentFragments.stream()
                                 .anyMatch(af -> af.getDirection() == Direction.REVERSE) ? Direction.REVERSE : Direction.FORWARD;
    }
}
