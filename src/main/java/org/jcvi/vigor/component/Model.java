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
	
	private List<Exon> exons = Collections.EMPTY_LIST;
	private Alignment alignment;
	private Map<String,Double> scores = Collections.EMPTY_MAP;
	private String geneSymbol;
	private List<String> status = Collections.EMPTY_LIST;
	private Direction direction;
    private boolean partial5p=false;
    private boolean partial3p=false;
    private boolean isPseudogene=false;
    private Range replaceStopCodonRange;
    private Range ribosomalSlippageRange;
    private Range insertRNAEditingRange;
    private NucleotideSequence cds;
    private ProteinSequence translatedSeq;
    private String geneID;
    private List<MaturePeptideMatch> maturePeptides = Collections.EMPTY_LIST;
    private EnumMap<NoteType,String> notes = new EnumMap(NoteType.class);

   public Model clone() throws CloneNotSupportedException {
	   Model model = (Model) super.clone();
	   model.setExons(this.getExons().stream().map(x->x.clone()).collect(Collectors.toList()));
	   model.setScores(new HashMap<>(this.scores));
	   model.setStatus(new ArrayList<>(this.status));

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
