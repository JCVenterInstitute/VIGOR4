package org.jcvi.vigor.utils;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;



public class VigorFunctionalUtils {
	
	/*public static AlignmentFragment mergeTwoFragments(AlignmentFragment frag1 , AlignmentFragment frag2){
		Range prevExonNTRange = frag1.getNucleotideSeqRange();
		Range nextExonNTRange = frag2.getNucleotideSeqRange();
		Range prevExonAARange = frag1.getProteinSeqRange();
		Range nextExonAARange = frag2.getProteinSeqRange();
		Range mergedExonNTRange = Range.of(prevExonNTRange.getBegin(),nextExonNTRange.getEnd());
		Range mergedExonAARange = Range.of(prevExonAARange.getBegin(),nextExonAARange.getEnd());
		frag1.setProteinSeqRange(mergedExonAARange);
		frag1.setNucleotideSeqRange(mergedExonNTRange);
		return frag1;
	}*/

	public static NucleotideSequence getCDS(Model model){
        NucleotideSequence virusGenomeSeq = model.getAlignment().getVirusGenome().getSequence();
        NucleotideSequenceBuilder NTSeqBuilder=new NucleotideSequenceBuilder("");
        for(Exon exon : model.getExons()){
            NTSeqBuilder.append(virusGenomeSeq.toBuilder(exon.getRange()));
        }
        NucleotideSequence cds = NTSeqBuilder.build();
	    return cds;
    }
	private static Frame[] FRAMES = { Frame.ONE, Frame.TWO, Frame.THREE};
	public static Frame getSequenceFrame(long coordinate){
		return FRAMES[(int) coordinate % 3];
	}

	// TODO more descriptive name?
	public static double generateScore(long referenceCoordinate,long pointOfOccurance){
		return 100d/(1d + Math.abs(pointOfOccurance - referenceCoordinate));
	}

	public static boolean isInFrameWithExon(List<Exon> exons,long match){
	    boolean isInFrame=false;
	    for(Exon exon: exons){
	        if(exon.getRange().intersects(Range.of(match))){
	            Frame exonFrame = getSequenceFrame(exon.getRange().getBegin()+exon.getFrame().getFrame()-1);
	            Frame matchFrame = getSequenceFrame(match);
	            if(exonFrame.equals(matchFrame)) isInFrame=true;
            }
        }

	    return isInFrame;
    }

    public static boolean intheSequenceGap(List<Range> sequenceGaps,Range match){
        if (sequenceGaps != null) {
            for (Range range : sequenceGaps) {
                if (range.intersects(match)) {
                    return true;
                }
            }
        }
        return false;
    }

    public static Map<Frame,List<Long>> frameToSequenceFrame(Map<Frame,List<Long>> rangeFrameMaP){
		Map<Frame,List<Long>> outRangeFrameMap = new HashMap<Frame,List<Long>>();
		for(Frame frame : rangeFrameMaP.keySet()){
			if(rangeFrameMaP.get(frame).size()>0){
			List<Long> ranges = rangeFrameMaP.get(frame);
			Long coordinate = ranges.get(0);
			Frame seqFrame = getSequenceFrame(coordinate);
			outRangeFrameMap.put(seqFrame, ranges);	
			}
		}
		return outRangeFrameMap;
	}
	/*public static List<Range> rangesMatchingSequenceFrame(List<Range> ranges, Frame frame){
		List<Range> outRanges = new ArrayList<Range>();
		for(Range range: ranges){
			Frame outFrame = getSequenceFrame(range.getBegin());
			if(outFrame.equals(frame)){
				outRanges.add(range);
			}
		}
		return outRanges;

	}*/

	public static Range get5pNearestSequenceGap(List<Range> sequenceGaps , Range exonRange){
	    Range nearestRange = null;
	    for(Range range : sequenceGaps){
	       if(range.getEnd()<exonRange.getBegin())
             nearestRange=range;
        }
        return nearestRange;
    }
	public static long getNTRange(List<Exon> exons,long CDSNTCoordinate){
		exons.sort(Exon.Comparators.Ascending);
		long outputStart=0;
		long refCurLength=0;
		long refBases=0;
		for(int i=0;i<exons.size();i++){
			long bases;
			if(i==0) {
				bases = exons.get(i).getRange().getBegin();
			}else{
				bases = exons.get(i).getRange().getBegin()-exons.get(i-1).getRange().getEnd()-1;
			}
			refBases = refBases+bases;
			refCurLength=exons.get(i).getRange().getEnd()-refBases;
		   if(refCurLength>=CDSNTCoordinate){
			    outputStart=CDSNTCoordinate+refBases;
			    break;
		   }

		}
		return outputStart;
	}

	public static Map<Frame,List<Long>> findStopsInSequenceFrame(VirusGenome virusGenome,Range searchRange){
		Map<Frame,List<Long>> stops = IupacTranslationTables.STANDARD.findStops(virusGenome.getSequence().toBuilder(searchRange).build());
		Map<Frame,List<Long>> stopsTemp = new HashMap<Frame,List<Long>>();
		final long startCoordinate = searchRange.getBegin();
		for(Frame frame : stops.keySet()){
			List<Long> temp = stops.get(frame).stream().map(x->x+startCoordinate).collect(Collectors.toList());
			stopsTemp.put(frame, temp);
		}
		Map<Frame,List<Long>> outStopsMap = VigorFunctionalUtils.frameToSequenceFrame(stopsTemp);
		return outStopsMap;
		
	}
   /* public static long get5partialStart(){
	    long start=0;

    }*/
	/*public Range getGenomeSequenceCoordinates(List<Exon> exons, Range insertRNAEditingRange,Range proteinRange){

    }*/
	
	
			
}
