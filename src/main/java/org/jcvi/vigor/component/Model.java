package org.jcvi.vigor.component;
import java.util.*;
import java.util.stream.Collectors;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Model implements Cloneable{
	
	private List<Exon> exons;
	private Alignment alignment;
	private Map<String,Double> scores;
	private String geneSymbol;
	private List<String> status;
	private Direction direction;
    private boolean partial5p=false;
    private boolean partial3p=false;
    private boolean isPseudogene=false;
    private Range replaceStopCodonRange;
    private Range insertRNAEditingRange;
    private NucleotideSequence cds;
    private ProteinSequence tanslatedSeq;
    private String proteinID;
   public Model clone() throws CloneNotSupportedException {
	   Model model = (Model) super.clone();
	   model.exons = exons.stream().map(x->x.clone()).collect(Collectors.toList());
	   model.setScores(null);
	   Map<String,Double> scoresCopy = new HashMap<String,Double>();
	   String key;
	   if(this.scores!=null) {
		   Iterator<String> it = this.scores.keySet().iterator();
		   while (it.hasNext()) {
			   key = it.next();
			   scoresCopy.put(key, this.scores.get(key));
		   }
	   }
       List<String> statusCopy = new ArrayList<String>();
	   if(this.status!=null) {
		   statusCopy.addAll(this.status);

	   }
	   model.setStatus(statusCopy);
	   model.setScores(scoresCopy);

	/*  
	   
	   try {
	     ByteArrayOutputStream baos = new ByteArrayOutputStream();
	     ObjectOutputStream oos = new ObjectOutputStream(baos);
	     oos.writeObject(model);
	     ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
	     ObjectInputStream ois = new ObjectInputStream(bais);
	     return (Model)(ois.readObject());
	   }
	   catch (Exception e) {
	     e.printStackTrace();
	     return null;
	   }
	 }*/
return model;

}
}
