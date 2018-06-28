package org.jcvi.vigor.utils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;

public class VigorFunctionalUtils {

    public static NucleotideSequence getCDS ( Model model ) {

        NucleotideSequence virusGenomeSeq = model.getAlignment().getVirusGenome().getSequence();
        NucleotideSequenceBuilder NTSeqBuilder = new NucleotideSequenceBuilder("");
        for (Exon exon : model.getExons()) {
            NTSeqBuilder.append(virusGenomeSeq.toBuilder(exon.getRange()));
        }
        NucleotideSequence cds = NTSeqBuilder.build();
        return cds;
    }

    private static Frame[] FRAMES = { Frame.ONE, Frame.TWO, Frame.THREE };

    public static Frame getSequenceFrame ( long coordinate ) {

        return FRAMES[ (int) coordinate % 3 ];
    }

    public static double generateScore ( long referenceCoordinate, long pointOfOccurance ) {

        long distance;
        double score;
        if (pointOfOccurance - referenceCoordinate < 0) {
            distance = referenceCoordinate - pointOfOccurance;
        } else {
            distance = pointOfOccurance - referenceCoordinate;
        }
        score = 100f / ( 1f + distance );
        return score;
    }

    public static boolean isInFrameWithExon ( List<Exon> exons, long match ) {

        boolean isInFrame = false;
        for (Exon exon : exons) {
            if (exon.getRange().intersects(Range.of(match))) {
                Frame exonFrame = getSequenceFrame(exon.getRange().getBegin() + exon.getFrame().getFrame() - 1);
                Frame matchFrame = getSequenceFrame(match);
                if (exonFrame.equals(matchFrame)) isInFrame = true;
            }
        }
        return isInFrame;
    }

    public static boolean intheSequenceGap ( List<Range> sequenceGaps, Range match ) {

        if (sequenceGaps != null) {
            for (Range range : sequenceGaps) {
                if (range.intersects(match)) {
                    return true;
                }
            }
        }
        return false;
    }

    public static Map<Frame, List<Long>> frameToSequenceFrame ( Map<Frame, List<Long>> rangeFrameMaP ) {

        Map<Frame, List<Long>> outRangeFrameMap = new HashMap<Frame, List<Long>>();
        for (Frame frame : rangeFrameMaP.keySet()) {
            if (rangeFrameMaP.get(frame).size() > 0) {
                List<Long> ranges = rangeFrameMaP.get(frame);
                Long coordinate = ranges.get(0);
                Frame seqFrame = getSequenceFrame(coordinate);
                outRangeFrameMap.put(seqFrame, ranges);
            }
        }
        return outRangeFrameMap;
    }

    public static Range get5pNearestSequenceGap ( List<Range> sequenceGaps, Range exonRange ) {

        Range nearestRange = null;
        for (Range range : sequenceGaps) {
            if (range.getEnd() <= exonRange.getBegin())
                nearestRange = range;
        }
        return nearestRange;
    }

    public static long getNTRange ( List<Exon> exons, long CDSNTCoordinate ) {

        exons.sort(Exon.Comparators.Ascending);
        long outputStart = 0;
        long refCurLength;
        long refBases = 0;
        for (int i = 0; i < exons.size(); i++) {
            long bases;
            if (i == 0) {
                bases = exons.get(i).getRange().getBegin();
            } else {
                bases = exons.get(i).getRange().getBegin() - exons.get(i - 1).getRange().getEnd() - 1;
            }
            refBases = refBases + bases;
            refCurLength = exons.get(i).getRange().getEnd() - refBases;
            if (refCurLength >= CDSNTCoordinate) {
                outputStart = CDSNTCoordinate + refBases;
                break;
            }
        }
        return outputStart;
    }

    public static Map<Frame, List<Long>> findStopsInSequenceFrame ( VirusGenome virusGenome, Range searchRange ) {

        Map<Frame, List<Long>> stops = IupacTranslationTables.STANDARD.findStops(virusGenome.getSequence().toBuilder(searchRange).build());
        Map<Frame, List<Long>> stopsTemp = new HashMap<Frame, List<Long>>();
        final long startCoordinate = searchRange.getBegin();
        for (Frame frame : stops.keySet()) {
            List<Long> temp = stops.get(frame).stream().map(x -> x + startCoordinate).collect(Collectors.toList());
            stopsTemp.put(frame, temp);
        }
        Map<Frame, List<Long>> outStopsMap = VigorFunctionalUtils.frameToSequenceFrame(stopsTemp);
        return outStopsMap;
    }
}
