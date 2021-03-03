package org.jcvi.vigor.component;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.utils.VigorFunctionalUtils;

import java.util.ArrayList;
import java.util.List;

public class MappedNucleotideSequence  {
    private static final Logger LOGGER = LogManager.getLogger(MappedNucleotideSequence.class);

    private final NucleotideSequence sequence;
    private final List<Range> originalRanges;

    public MappedNucleotideSequence(NucleotideSequence sequence, List<Range> originalRanges) {
        this.sequence = sequence;
        this.originalRanges = new ArrayList(originalRanges);
    }

    public NucleotideSequence getSequence() {
        return sequence;
    }

    public long getOriginalCoordinate(long coordinate) {
        try {
            List<Range> ranges = getOriginalRanges(Range.of(coordinate));
            return ranges.get(0).getBegin();
        } catch (RuntimeException e) {
            throw new RuntimeException("Unable to get sequence coordinate for " + coordinate);
        }
    }

    /**
     * Given a range in the CDS sequence, return corresponding range or ranges in the original sequence
     * @param range
     * @return
     */
    public List<Range> getOriginalRanges(Range range) {
        LOGGER.debug("For sequence {} getting originalRanges for CDS range {}", sequence, range);
        List<Range> ranges = new ArrayList<>();
        long index = 0;
        Range addedRange;
        for (Range originalRange: originalRanges) {
            // negative ranges are insertion ranges which don't have original coordinates
            if (originalRange.getBegin() >= 0) {
                // this is the range in the CDS
                Range cdsRange = Range.ofLength(originalRange.getLength()).toBuilder().shift(index).build();
                Range intersection = cdsRange.intersection(range);
                if (intersection.getLength() > 0) {
                    addedRange = Range.ofLength(intersection.getLength()).toBuilder().shift(originalRange.getBegin() + (intersection.getBegin() - cdsRange.getBegin())).build();
                    LOGGER.debug("Given CDS range {}, original range {} adding range {}", cdsRange, originalRange, addedRange);
                    ranges.add(addedRange);
                }
            }
            index += originalRange.getLength();
        }
        if (ranges.isEmpty()) {
            throw new RuntimeException("Unable to get NT range(s) for AA range " + range.toString());
        }
        return VigorFunctionalUtils.mergeAdjacentRanges(ranges);
    }
}

