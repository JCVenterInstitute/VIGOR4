package com.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

/**
 * Created by snettem on 5/8/2017.
 */

@Component
@Scope("prototype")
@Data
public class RNA_Editing {

    private boolean has_RNA_editing=false;
    private int size;
    private String regExp;
    private String replacementString;
    private String note;

}
