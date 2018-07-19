package org.jcvi.vigor.component;

import lombok.Data;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

/**
 * Created by snettem on 5/8/2017.
 */
@Component
@Scope("prototype")
@Data
public class RNA_Editing {

    private static Logger LOGGER = LogManager.getLogger(RNA_Editing.class);
    private final boolean has_RNA_editing;
    private final int offset;
    private final String regExp;
    private final String insertionString;
    private final String note;

    public RNA_Editing ( boolean has_editing, int offset, String regExp, String insertionString, String note ) {

        this.has_RNA_editing = has_editing;
        this.offset = offset;
        this.regExp = regExp;
        this.insertionString = insertionString;
        this.note = note;
    }

    /**
     *
     * @param rnaEditingString
     * @return
     * @throws IllegalArgumentException
     *
     * TODO what is the first element if it's not the offset?
     */
    public static RNA_Editing parseFromString(String rnaEditingString) throws IllegalArgumentException {
         String[] temp = rnaEditingString.split("/");
         if (temp.length !=4) {
             throw new IllegalArgumentException(
                     String.format("Bad format for rna editing \"%s\". Format is OFFSET/INSERTION/MOTIF/NOTE", rnaEditingString)
             );
         }
        int rna_editing_offset = 0;
        if (VigorUtils.is_Integer(temp[0])) {
            rna_editing_offset = Integer.parseInt(temp[0]);
        } else if (! temp[0].trim().isEmpty()) {
            LOGGER.warn("Bad offset value {} for RNA editing. Full string {}", temp[0], rnaEditingString);
        }
        return new RNA_Editing(true, rna_editing_offset, temp[2], temp[1], temp[3]);
    }

    public static final RNA_Editing NO_EDITING = new RNA_Editing(false, 0, "", "", "");
}