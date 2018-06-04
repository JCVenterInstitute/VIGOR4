package org.jcvi.vigor.component;


import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

/**
 * Created by snettem on 5/8/2017.
 */

@Component
@Scope("prototype")
@Data
public class RNA_Editing  {

    private final boolean has_RNA_editing;
    private final int offset;
    private final String regExp;
    private final String insertionString;
    private final String note;

    public RNA_Editing (boolean has_editing, int offset, String regExp, String insertionString, String note) {
        this.has_RNA_editing = has_editing;
        this.offset = offset;
        this.regExp = regExp;
        this.insertionString = insertionString;
        this.note = note;
    }

    public static final RNA_Editing NO_EDITING = new RNA_Editing(false, 0, "", "", "");
}