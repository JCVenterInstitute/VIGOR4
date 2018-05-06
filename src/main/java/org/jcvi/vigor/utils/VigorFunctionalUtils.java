package org.jcvi.vigor.utils;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.VirusGenome;



public class VigorFunctionalUtils {
	
	public static AlignmentFragment mergeTwoFragments(AlignmentFragment frag1 , AlignmentFragment frag2){
		Range prevExonNTRange = frag1.getNucleotideSeqRange();
		Range nextExonNTRange = frag2.getNucleotideSeqRange();
		Range prevExonAARange = frag1.getProteinSeqRange();
		Range nextExonAARange = frag2.getProteinSeqRange();
		Range mergedExonNTRange = Range.of(prevExonNTRange.getBegin(),nextExonNTRange.getEnd());
		Range mergedExonAARange = Range.of(prevExonAARange.getBegin(),nextExonAARange.getEnd());
		frag1.setProteinSeqRange(mergedExonAARange);
		frag1.setNucleotideSeqRange(mergedExonNTRange);
		return frag1;
	}
	
	public static Frame getSequenceFrame(long coordinate){
		Frame frame=Frame.ONE;
		coordinate=coordinate+1;
		if(coordinate>2){
		long remin = coordinate%3;
		if(remin==0){
			frame=Frame.THREE;
		}
		else if(remin==1){
			frame=Frame.ONE;
		}else
		{
			frame = Frame.TWO;
		}
		}else if(coordinate==1)frame=Frame.ONE;
		else if (coordinate==2)frame=Frame.TWO;
		else frame=Frame.THREE;
		return frame;
	}
	public static double generateScore(long referenceCoordinate,long pointOfOccurance){
		long distance = 0;
		double score =0;
			
			if (pointOfOccurance - referenceCoordinate < 0) {
				distance = referenceCoordinate - pointOfOccurance;
				
			} else {
				distance = pointOfOccurance - referenceCoordinate;
			}
			
			  
			  score = 100f/(1f+distance);
			 				
		return score;
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
	public static List<Range> rangesMatchingSequenceFrame(List<Range> ranges, Frame frame){
		List<Range> outRanges = new ArrayList<Range>();
		for(Range range: ranges){
			Frame outFrame = getSequenceFrame(range.getBegin());
			if(outFrame.equals(frame)){
				outRanges.add(range);
			}
		}
		return outRanges;
	}
		
	public static long getNTRange(List<Exon> exons,long CDSNTCoordinate){
		exons.sort(Exon.Comparators.Ascending);
		long outputStart=0;
		long outputEnd=0;
		long refCurLength=0;
		long refPreLength=0;
		long difference=0;
		long refBases=0;
		for(int i=0;i<exons.size();i++){
			long bases=0;
			if(i==0) {
				bases = exons.get(i).getRange().getBegin();
			}else{
				bases = exons.get(i).getRange().getBegin()-exons.get(i-1).getRange().getEnd()-1;
			}
		   refCurLength=exons.get(i).getRange().getLength()-bases+refCurLength;
		   if(refCurLength>=CDSNTCoordinate){
			    outputStart=CDSNTCoordinate+refBases;
			    break;
		   }
		   refBases = refBases+bases;
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

	/*public Range getGenomeSequenceCoordinates(List<Exon> exons, Range insertRNAEditingRange,Range proteinRange){

    }*/
	
	
			
}
