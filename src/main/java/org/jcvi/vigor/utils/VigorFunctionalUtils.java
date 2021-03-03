package org.jcvi.vigor.utils;

import java.util.*;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.component.MappedNucleotideSequence;

public class VigorFunctionalUtils {

    private static final Logger LOGGER = LogManager.getLogger(VigorFunctionalUtils.class);

    public static MappedNucleotideSequence getCDS( Model model) {
        return getCDS(model, true);
    }

    /**
     * @param model
     * @return coding sequence of model
     */
    public static MappedNucleotideSequence getCDS (Model model, boolean trimStop) {

        List<Range> originalRanges = new ArrayList<>();
        model.getExons().sort(Exon.Comparators.Ascending);
        List<Exon> exons = model.getExons();
        NucleotideSequence virusGenomeSeq = model.getAlignment().getVirusGenome().getSequence();
        boolean inserted = false;
        NucleotideSequenceBuilder ntSequenceBuilder = new NucleotideSequenceBuilder("");
        for (int i = 0; i < exons.size(); i++) {
            Range exonRange = exons.get(i).getRange();
            // trim off stop codon unless the model is partial
            if (i == exons.size() - 1 && trimStop && !model.isPartial3p()) {
                exonRange = Range.of(exonRange.getBegin(), exonRange.getEnd() - 3);
            }
            ntSequenceBuilder.append(virusGenomeSeq.toBuilder(exonRange));
            originalRanges.add(exonRange);
            if (!inserted && model.getInsertRNAEditingRange() != null && model.getInsertRNAEditingRange().getBegin() == (exonRange.getEnd() + 1)) {
                String insertionString = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing().getInsertionString();
                ntSequenceBuilder.append(insertionString);
                inserted = true;
                originalRanges.add(Range.of(-(insertionString.length()), -1));
            }
        }
        return new MappedNucleotideSequence(ntSequenceBuilder.build(), originalRanges);
    }

    private static Frame[] FRAMES = { Frame.ONE, Frame.TWO, Frame.THREE };

    public static Frame getSequenceFrame ( long coordinate ) {

        return FRAMES[ (int) coordinate % 3 ];
    }

    /**
     * TODO This does not take into account any viral tricks except insertion
     * @param model
     * @param proteinRange
     * @return
     */
    public static List<Range> proteinRangeToCDSRanges ( Model model, Range proteinRange ) {

        LOGGER.debug("For model {} getting CDS ranges for AA Range {}", model, proteinRange);
        Range ntRange = Range.ofLength(proteinRange.getLength() * 3 ).toBuilder().shift(proteinRange.getBegin() * 3).build();
        List<Range> cdsRanges =  model.getCDS().getOriginalRanges(ntRange);
        LOGGER.debug("For model {} for AA range {} returning NT Range {}", model, proteinRange, String.join(",", cdsRanges.stream().map(Range::toString).collect(Collectors.toList())));
        return cdsRanges;
    }

    public static double generateProximityScore ( long referenceCoordinate, long pointOfOccurance ) {

        return 100d / ( 1d + Math.abs(pointOfOccurance - referenceCoordinate) );
    }

    public static boolean isInFrameWithExon ( List<Exon> exons, long match ) {

        boolean isInFrame = false;
        Range matchRange = Range.of(match);
        for (Exon exon : exons) {
            if (exon.getRange().intersects(matchRange)) {
                Frame exonFrame = getSequenceFrame(exon.getRange().getBegin() + exon.getFrame().getFrame() - 1);
                Frame matchFrame = getSequenceFrame(match);
                if (exonFrame.equals(matchFrame)) {
                    isInFrame = true;
                    break;
                }
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

        Map<Frame, List<Long>> outRangeFrameMap = new HashMap<>();
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
        Map<Frame, List<Long>> stopsTemp = new HashMap<>();
        final long startCoordinate = searchRange.getBegin();
        for (Frame frame : stops.keySet()) {
            List<Long> temp = stops.get(frame).stream().map(x -> x + startCoordinate).collect(Collectors.toList());
            stopsTemp.put(frame, temp);
        }
        Map<Frame, List<Long>> outStopsMap = VigorFunctionalUtils.frameToSequenceFrame(stopsTemp);
        return outStopsMap;
    }


    public static long getDirectionBasedCoordinate(long coordinate,long seqLength,Direction direction){
        return ((direction==Direction.REVERSE) ? seqLength+1-coordinate : coordinate);
    }

    public static Range getDirectionBasedRange(Range range, long seqLength, Direction direction) {
        long begin = getDirectionBasedCoordinate(range.getBegin(Range.CoordinateSystem.RESIDUE_BASED), seqLength, direction);
        long end = getDirectionBasedCoordinate(range.getEnd(Range.CoordinateSystem.RESIDUE_BASED), seqLength, direction);
        return Range.of(Range.CoordinateSystem.RESIDUE_BASED, Math.min(begin, end), Math.max(begin, end));
    }

    public static List<Range> mergeAdjacentRanges(List<Range> ranges) {
        List<Range> returnedRanges = new ArrayList<>(ranges.size());
        Range r,s;
        if (! ranges.isEmpty()) {
            returnedRanges.add(ranges.get(0));
            for (int i=1; i < ranges.size(); i++) {
                r = ranges.get(i);
                s = returnedRanges.get(returnedRanges.size() -1);
                if (r.getBegin() == s.getEnd() + 1) {
                    returnedRanges.set(returnedRanges.size() -1,
                                       s.toBuilder().expandEnd(r.getLength()).build());
                } else {
                    returnedRanges.add(r);
                }
            }
        }
        return returnedRanges;
    }

}
