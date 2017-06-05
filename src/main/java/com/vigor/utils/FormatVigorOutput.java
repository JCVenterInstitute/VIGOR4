package com.vigor.utils;

import com.sun.javafx.binding.StringFormatter;
import com.vigor.component.Alignment;
import com.vigor.component.Exon;
import com.vigor.component.Model;
import com.vigor.forms.VigorForm;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by snettem on 5/24/2017.
 */
public class FormatVigorOutput {

    public static void printModels(List<Model> models){
        System.out.println("********************************************************************************************************************");
        System.out.println(String.format ( "%-32s%-20s%-20s%-20s", "Gene_Symbol", "Direction", "Status","Range"));

        for(Model model : models) {

            String status = model.getStatus ().stream ().collect ( Collectors.joining(System.lineSeparator ()+"%30s"));
            String ranges = model.getExons ().stream ().map(e->e.getRange ().toString ()).collect ( Collectors.joining(System.lineSeparator ()+"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t") );
            System.out.println(String.format ( "%-32s%-20s%-20s%-20s", model.getGeneSymbol (),model.getDirection (),status,ranges));
            System.out.println(System.lineSeparator ());
        }
    }


}
