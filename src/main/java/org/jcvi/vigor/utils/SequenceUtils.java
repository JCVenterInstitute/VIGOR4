package org.jcvi.vigor.utils;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.Nucleotide;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.Triplet;

import java.util.Iterator;
public class SequenceUtils {


	
	public static Triplet getNextTriplet(Iterator<Nucleotide> iter){
		
		Nucleotide first= getNextNucleotide(iter);
		Nucleotide second= getNextNucleotide(iter);
		Nucleotide third= getNextNucleotide(iter);
		if(first==null || second ==null || third ==null){
			return null;
		}
		return Triplet.create(first, second, third);
	}
	
	private static Nucleotide getNextNucleotide(Iterator<Nucleotide> iter){
		if(!iter.hasNext()){
			return null;
		}
		Nucleotide n = iter.next();
		return n;
	}
	
	
	public static Iterator<Nucleotide> handleFrame(NucleotideSequence sequence, Frame frame) {
	    Iterator<Nucleotide> iter;
	    if(frame.onReverseStrand()){
	        iter = sequence.toBuilder().reverseComplement().iterator();
	          switch(frame){
                        case NEGATIVE_THREE:
                                        if(iter.hasNext()){
                                                iter.next();
                                        }
                        case NEGATIVE_TWO:
                                        if(iter.hasNext()){
                                                iter.next();
                                        }
                                        break;
                        default:
                                       
                                break;
                }
	    }else{
	        iter = sequence.iterator();
		switch(frame){
			case THREE:
					if(iter.hasNext()){
						iter.next();
					}
			case TWO:
					if(iter.hasNext()){
						iter.next();
					}
					break;
			default:
					break;
		}
	    }
	    return iter;
	}

	
}
