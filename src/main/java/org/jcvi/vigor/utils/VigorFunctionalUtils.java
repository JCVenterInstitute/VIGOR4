package org.jcvi.vigor.utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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

public class VigorFunctionalUtils {

    private static final Logger LOGGER = LogManager.getLogger(VigorFunctionalUtils.class);

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

    public static List<Range> proteinRangeToCDSRanges ( Model model, Range proteinRange ) {

        List<Range> ranges = new ArrayList<>();
        long pBegin = proteinRange.getBegin() * 3;
        long proteinNTLength = proteinRange.getLength() * 3;
        LOGGER.trace("getting ranges for proteinRange {}, ntBegin {}, ntLength {}", proteinRange, pBegin, proteinNTLength);
        List<Exon> exons = model.getExons();
        exons.sort(Exon.Comparators.Ascending);
        boolean inserted = false;
        long proteinBases = 0;
        int insertedLength;
        long exonLength;
        Range addedRange;
        long rangeStart;
        long frameAdjustment;
        Exon exon;
        for (int i = 0; i < exons.size() && proteinNTLength > 0; i++) {
            exon = exons.get(i);
            frameAdjustment = exon.getFrame().getFrame() - 1;
            Range exonRange = exon.getRange();
            long adjustedBegin = exonRange.getBegin() + frameAdjustment;
            long adjustedEnd = exonRange.getEnd() + frameAdjustment;
            // TODO should ranges include stop codons?
            // trim off stop codon unless the model is already partial
            if (i == exons.size() - 1 && !model.isPartial3p()) {
                adjustedEnd -= 3;
            }
            exonLength = adjustedEnd - adjustedBegin;
            LOGGER.trace("checking range {}-{}, ntBegin {} against ntcount {}", adjustedBegin, adjustedEnd, pBegin, proteinBases);
            if (proteinBases + exonLength > pBegin) {
                rangeStart = adjustedBegin + ( pBegin - proteinBases );
                if (proteinNTLength > exonLength) {
                    addedRange = Range.of(rangeStart, adjustedEnd);
                } else {
                    addedRange = Range.of(rangeStart, rangeStart + proteinNTLength - 1);
                }
                ranges.add(addedRange);
                LOGGER.trace("added range {}", addedRange);
                proteinNTLength -= addedRange.getLength();
                LOGGER.trace("{} bases left to to account for", proteinNTLength > 0 ? proteinNTLength : 0);
            }
            proteinBases += exonLength;
            if (!inserted && model.getInsertRNAEditingRange() != null && model.getInsertRNAEditingRange().getBegin() == adjustedEnd + 1) {
                insertedLength = model.getAlignment().getViralProtein().getGeneAttributes().getRna_editing().getInsertionString().length();
                proteinBases += insertedLength;
                if (proteinBases > pBegin) {
                    LOGGER.trace("accounting for RNA editing of legth {}", insertedLength);
                    proteinNTLength -= ( proteinBases - pBegin );
                }
                inserted = true;
            }
        }
        return ranges;
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

}
