package org.jcvi.vigor.utils;


public enum NoteType {
    RNA_Editing("non-templated G's inserted during transcription"),
    RNA_Editing_Location("location of RNA editing"),
    Sequence_Gap("coding region disrupted by sequencing gap"),
    Gene(""),
    Remark(""),
    StopCodonReadThrough("");



    private final String text;

    NoteType(String text){
        this.text = text;
    }
    @Override
    public String toString(){
        return text;
    }


}
