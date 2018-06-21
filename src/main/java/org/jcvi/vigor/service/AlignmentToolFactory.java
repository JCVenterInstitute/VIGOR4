package org.jcvi.vigor.service;

import org.jcvi.vigor.component.AlignmentTool;
import org.jcvi.vigor.component.Exonerate;

import java.io.File;

public class AlignmentToolFactory {

    public static AlignmentTool getAlignmentTool(String refDB){
        refDB= new File(refDB).getName();
        if(refDB.matches("flua_db|veev_db|rsv_db|flub_db|fluc_db")){
            return new Exonerate("exonerate");
        }
        return null;
    }

}
