package org.jcvi.vigor.utils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;

public class VigorFunctionalUtils {
	
	public static Frame getFrame(Range inputRange){
		Range range = Range.of(inputRange.getBegin(),inputRange.getEnd());
		Frame frame=Frame.ONE;
		if(range.getBegin()>2){
		long remin = range.getBegin()%3;
		if(remin==0){
			frame=Frame.ONE;
		}
		else if(remin==1){
			frame=Frame.TWO;
		}else
		{
			frame = Frame.THREE;
		}
		}else if(range.getBegin()==0)frame=Frame.ONE;
		else if (range.getBegin()==1)frame=Frame.TWO;
		else frame=Frame.THREE;
		return frame;
	}
	public static Map<Range,Float> generateScore(long centroid,List<Range> ranges){
		Map<Range,Float> rangeScoreMap = new HashMap<Range,Float>();
		float difference = 0;
		float score =0;
		for (Range range : ranges) {
			
			if (range.getBegin() - centroid < 0) {
				difference = centroid - range.getBegin();
				
			} else if(range.getBegin() -centroid ==0){
			    difference = 0.1f;
			}
			else{
				difference = range.getBegin() - centroid;
				}
			  
			  score = 100/difference;
			  rangeScoreMap.put(Range.of(range.getBegin(),range.getEnd()), score);				
		}		
		return rangeScoreMap;
	}
	
	
}
