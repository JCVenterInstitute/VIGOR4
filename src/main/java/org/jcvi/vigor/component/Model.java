package org.jcvi.vigor.component;

import java.util.*;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.utils.NoteType;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Model implements Cloneable {

    private List<Exon> exons = Collections.EMPTY_LIST;
    private Alignment alignment;
    private Map<String, Double> scores = Collections.EMPTY_MAP;
    private String geneSymbol;
    private List<String> status = Collections.EMPTY_LIST;
    private Direction direction;
    private boolean partial5p = false;
    private boolean partial3p = false;
    private boolean isPseudogene = false;
    private Range replaceStopCodonRange;
    private Range ribosomalSlippageRange;
    private Range insertRNAEditingRange;
    private ProteinSequence translatedSeq;
    private String geneID;
    private List<NoteType> notes = Collections.EMPTY_LIST;
    private List<MaturePeptideMatch> maturePeptides = Collections.EMPTY_LIST;

    public Model clone () throws CloneNotSupportedException {

        Model model = (Model) super.clone();
        model.setExons(this.getExons().stream().map(x -> x.clone()).collect(Collectors.toList()));
        model.setScores(new HashMap<>(this.scores));
        model.setStatus(new ArrayList<>(this.status));
        model.setNotes(new ArrayList<>(this.notes));
        return model;
    }

    public Range getRange () {

        List<Exon> exons = getExons();
        if (!exons.isEmpty()) {
            return Range.of(exons.get(0).getRange().getBegin(),
                    exons.get(exons.size() - 1).getRange().getEnd());
        }
        return Range.ofLength(0);
    }

    public String getProteinID() {
        return getAlignment().getViralProtein().getProteinID();
    }

    @Override
    public String toString(){
        return "\n Gene Symbol : "+geneSymbol+"\n isPseudogene : "+isPseudogene
                +"\n isPartial3p : "+partial3p
                +"\n isPartial5p : "+partial5p
                +"\n Exons : "+exons.stream().map(Object::toString).collect(Collectors.joining("\n "));
    }

}
