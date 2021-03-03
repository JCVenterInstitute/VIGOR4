package org.jcvi.vigor.component;

import java.util.*;
import java.util.stream.Collectors;

import com.google.common.base.MoreObjects;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.vigor.component.MappedNucleotideSequence;
import org.jcvi.vigor.utils.NoteType;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Model implements Cloneable {

    private List<Exon> exons = new ArrayList<>();
    private Alignment alignment;
    private Map<String, Double> scores = new HashMap<>();
    private String geneSymbol;
    private List<String> status = new ArrayList<>();
    private Direction direction;
    private boolean partial5p = false;
    private boolean partial3p = false;
    private boolean isPseudogene = false;
    private Range replaceStopCodonRange;
    private Range ribosomalSlippageRange;
    private Range insertRNAEditingRange;
    private ProteinSequence translatedSeq;
    private String geneID;
    private List<String> notes = new ArrayList<>();
    private List<MaturePeptideMatch> maturePeptides = new ArrayList<>();

    public Model clone () throws CloneNotSupportedException {

        Model model = (Model) super.clone();
        model.setExons(this.getExons().stream().map(Exon::clone).collect(Collectors.toList()));
        model.setScores(new HashMap<>(this.scores));
        model.setStatus(new ArrayList<>(this.status));
        model.setNotes(new ArrayList<>(this.notes));
        return model;
    }

    /**
     * Get the nucleotide range from the beginning of the first exon to the end of the last exon
     *
     * @return
     */
    public Range getRange () {

        List<Exon> exons = getExons();
        if (!exons.isEmpty()) {
            long start = exons.get(0).getRange().getBegin();
            long end = exons.get(exons.size() - 1).getRange().getEnd();
            return Range.of(Math.min(start, end), Math.max(start, end));
        }
        return Range.ofLength(0);
    }

    public String getProteinID() {
        return getAlignment().getViralProtein().getProteinID();
    }

    @Override
    public String toString() {
        return MoreObjects.toStringHelper(this)
                          .omitNullValues()
                          .add("ID", hashCode())
                          .add("Gene", geneSymbol)
                          .add("protein", getProteinID())
                          .add("exons", exons.stream().map(Object::toString).collect(Collectors.joining(",")))
                          .toString();
    }

    public void addNote(NoteType note) {
        getNotes().add(note.toString());
    }

    public void addNote(String note) {
        getNotes().add(note);
    }

    public MappedNucleotideSequence getCDS() {
        return getCDS(true);
    }

    public MappedNucleotideSequence getCDS(boolean trimStops) {
        // can't cache while exons are changing underneath us
        return VigorFunctionalUtils.getCDS(this, trimStops);
    }

}
