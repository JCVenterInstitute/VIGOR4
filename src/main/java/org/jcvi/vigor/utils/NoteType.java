package org.jcvi.vigor.utils;

public enum NoteType {
    RNA_Editing("non-templated G's inserted during transcription"),
    Sequence_Gap("coding region disrupted by sequencing gap"),
    Gene(""),
    StopCodonReadThrough("Translation Exception"),
    StopCodonInterruption("CDS interrupted by stop codon"),
    StopCodonsInterruption("CDS interrupted by many stop codons");
    private final String text;

    NoteType ( String text ) {

        this.text = text;
    }

    @Override
    public String toString () {

        return text;
    }
}
