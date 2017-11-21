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
public class AlignmentFragment implements Comparable<AlignmentFragment>,Cloneable {

	private double score;
	private Direction direction;
	private Range proteinSeqRange;
	private Range nucleotideSeqRange;
    private Frame frame;
	
	public AlignmentFragment(Range proteinSeqRange,Range nucleotideRange,double score,Direction direction,Frame frame){
		this.score=score;
		this.proteinSeqRange=proteinSeqRange;
		this.nucleotideSeqRange=nucleotideRange;
		this.direction=direction;
		this.frame=frame;
	}
	public AlignmentFragment(){
		
	}
	
  public AlignmentFragment clone(){
		
		AlignmentFragment frag=null;
		try {
			frag = (AlignmentFragment)(super.clone());
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return frag;
	}
	
	public int compareTo(AlignmentFragment compareFragment){
		Range compareRange = compareFragment.getNucleotideSeqRange();
		if(compareRange.endsBefore(this.nucleotideSeqRange)){
			return 1;
		}else if(this.nucleotideSeqRange.endsBefore(compareRange)){
			return -1;
		}else return 0;
	}
	
	public enum Comparators implements Comparator<AlignmentFragment>{
	Descending{
			@Override
			  public int compare(AlignmentFragment e1, AlignmentFragment e2) {
                return -1 * JillionUtil.compare(e1.getProteinSeqRange().getBegin(), e2.getProteinSeqRange().getBegin());
            }
		}
	,
	
	Ascending{

        @Override
        public int compare(AlignmentFragment e1, AlignmentFragment e2) {
            return JillionUtil.compare(e1.getProteinSeqRange().getBegin(),e2.getProteinSeqRange().getBegin());
        }
        
    };
	 
	}

	
}
