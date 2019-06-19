package org.jcvi.vigor.component;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

import java.util.*;
import java.util.stream.Collectors;

@Component
@Scope("prototype")
@Data
public class VirusGenome {

    // sequence, defline and id :-> single object from jillion
    private NucleotideSequence sequence;
    private String defline;
    private String id;
    private Boolean isCircular = false;
    private List<Range> sequenceGaps = Collections.EMPTY_LIST;
    private Map<Frame, List<Long>> internalStops = Collections.EMPTY_MAP;

    public VirusGenome ( NucleotideSequence sequence, String defline, String id, boolean isCircular ) {

        this.sequence = sequence;
        this.defline = defline;
        this.id = id;
        this.isCircular = isCircular;
    }

    public VirusGenome () {

    }

    public VirusGenome(VirusGenome copyFrom) {
        this(copyFrom.getSequence(), copyFrom.getDefline(), copyFrom.getId(), copyFrom.getIsCircular());
        setSequenceGaps(new ArrayList<>(copyFrom.getSequenceGaps()));
        setInternalStops(copyFrom.getInternalStops()
                                 .entrySet()
                                 .stream()
                                 .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));
    }
}