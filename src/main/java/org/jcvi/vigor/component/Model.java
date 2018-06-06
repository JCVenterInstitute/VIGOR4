package org.jcvi.vigor.component;
import java.util.*;
import java.util.stream.Collectors;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.vigor.utils.NoteType;
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
    private Range ribosomalSlippageRange;
    private Range insertRNAEditingRange;
    private NucleotideSequence cds;
    private ProteinSequence tanslatedSeq;
    private String geneID;
    private List<MaturePeptideMatch> maturePeptides;
    private List<NoteType> notes;

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
	   return model;
   }

   public Range getRange() {
	   List<Exon> exons = getExons();
	   if (! exons.isEmpty()) {
		   return Range.of(exons.get(0).getRange().getBegin(),
				   exons.get(exons.size() - 1).getRange().getEnd());
	   }
	   return Range.ofLength(0);
   }
}
