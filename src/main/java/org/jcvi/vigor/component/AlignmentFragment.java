package org.jcvi.vigor.component;

import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.internal.core.util.JillionUtil;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

import java.util.Comparator;

@Component
@Scope("prototype")
@Data
public class AlignmentFragment implements Comparable<AlignmentFragment>, Cloneable {

    private final Direction direction;
    private final Range proteinSeqRange;
    private final Range nucleotideSeqRange;
    private final Frame frame;

    public AlignmentFragment ( Range proteinSeqRange, Range nucleotideRange, Direction direction, Frame frame ) {

        this.proteinSeqRange = proteinSeqRange;
        this.nucleotideSeqRange = nucleotideRange;
        this.direction = direction;
        this.frame = frame;
    }
    
    public int compareTo ( AlignmentFragment compareFragment ) {

        return Range.Comparators.ARRIVAL.compare(getNucleotideSeqRange(),
                                                 compareFragment.getNucleotideSeqRange());
    }

    public enum Comparators implements Comparator<AlignmentFragment> {
        Descending {
            @Override
            public int compare ( AlignmentFragment e1, AlignmentFragment e2 ) {

                return -1 * JillionUtil.compare(e1.getProteinSeqRange().getBegin(), e2.getProteinSeqRange().getBegin());
            }
        },
        Ascending {
            @Override
            public int compare ( AlignmentFragment e1, AlignmentFragment e2 ) {

                return JillionUtil.compare(e1.getProteinSeqRange().getBegin(), e2.getProteinSeqRange().getBegin());
            }
        };
    }
}
