package org.jcvi.vigor.component;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

import java.util.List;
import java.util.Map;

@Component
@Scope("prototype")
@Data
public class VirusGenome {

    // sequence, defline and id :-> single object from jillion
    private NucleotideSequence sequence;
    private String defline;
    private String id;
    private Boolean isComplete = false;
    private Boolean isCircular = false;
    private List<Range> sequenceGaps;
    private Map<Frame, List<Long>> internalStops;

    public VirusGenome ( NucleotideSequence sequence, String defline, String id, boolean isComplete, boolean isCircular ) {

        this.sequence = sequence;
        this.defline = defline;
        this.id = id;
        this.isComplete = isComplete;
        this.isCircular = isCircular;
    }

    public VirusGenome () {

    }
}