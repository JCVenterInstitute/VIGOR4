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
		if(coordinate>2){
		long remin = coordinate%3;
		if(remin==0){
			frame=Frame.ONE;
		}
		else if(remin==1){
			frame=Frame.TWO;
		}else
		{
			frame = Frame.THREE;
		}
		}else if(coordinate==0)frame=Frame.ONE;
		else if (coordinate==1)frame=Frame.TWO;
		else frame=Frame.THREE;
		return frame;
	}
	public static float generateScore(long referenceCoordinate,long pointOfOccurance){
		long distance = 0;
		float score =0;
			
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
		
	public static Range getNTRange(List<Exon> exons,Range CDSNTRange){
		exons.sort(Exon.Comparators.Ascending);
		long outputStart=0;
		long outputEnd=0;
		long refCurLength=0;
		long refPreLength=0;
		long difference=0;
		long CDSNtStart=CDSNTRange.getBegin();
		for(int i=0;i<exons.size();i++){
		   refCurLength=exons.get(i).getRange().getLength()+refCurLength;
		   if(refPreLength==0){
			   refPreLength=refCurLength;
		   }
		   if(refCurLength>=CDSNtStart){
			    difference = CDSNtStart-refPreLength;
			    outputStart=exons.get(i).getRange().getBegin()+difference;
			    long leftOver = exons.get(i).getRange().getLength()-difference;
			    if(leftOver==1){
			    	outputEnd = exons.get(i+1).getRange().getBegin()+2;
			    }else if(leftOver==2){
			    	outputEnd = exons.get(i+1).getRange().getBegin()+1;
			    }else{
			    	outputEnd = outputStart+2;
			    }
			    break;
		   }
		
		  refPreLength=refCurLength;
		}
		return Range.of(outputStart,outputEnd);
        }
	
	public static Long getNTCoordinate(List<Exon> exons, Long coordinate){
		Range outputRange = getNTRange(exons,Range.of(coordinate));
		return outputRange.getBegin();
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
	
	
			
}
