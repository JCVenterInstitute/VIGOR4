package org.jcvi.vigor.component;

import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;


@Component
@Scope("prototype")
@Data
@SuppressWarnings("serial")
public class AlignmentFragment implements Serializable {

	private double score;
	private Direction direction;
	private Range proteinSeqRange;
	private Range nucleotideSeqRange;
    private Frame frame;
	private boolean isSubChain=false;

	
	 public static AlignmentFragment deepClone(AlignmentFragment alignmentFragment) {
		   try {
		     ByteArrayOutputStream baos = new ByteArrayOutputStream();
		     ObjectOutputStream oos = new ObjectOutputStream(baos);
		     oos.writeObject(alignmentFragment);
		     ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		     ObjectInputStream ois = new ObjectInputStream(bais);
		     return (AlignmentFragment)(ois.readObject());
		   }
		   catch (Exception e) {
		     e.printStackTrace();
		     return null;
		   }
		 }

	
	
}
