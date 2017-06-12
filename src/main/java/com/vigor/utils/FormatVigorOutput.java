package com.vigor.utils;

import com.vigor.component.Model;
import java.util.List;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;

/**
 * Created by snettem on 5/24/2017.
 */
public class FormatVigorOutput {

    public static void printModels(List<Model> models){
        System.out.println("********************************************************************************************************************");
        System.out.println(String.format ( "%-32s%-20s%-20s%-20s", "Gene_Symbol", "Direction","NTSeqRange","AASeqRange"));

        for(Model model : models) {
         
        System.out.print(String.format("%-32s",model.getGeneSymbol()));
        System.out.print(String.format("%-20s", model.getDirection()));
        List<Range> NTranges = model.getExons ().stream ().map(e->e.getRange()).collect(Collectors.toList());
        List<Range> AAranges = model.getExons().stream().map(e->e.getAlignmentFragment().getProteinSeqRange()).collect(Collectors.toList()) ;
        for(int i=0;i<NTranges.size();i++){
        	System.out.print(String.format("%-20s", NTranges.get(i).getBegin()+"-"+NTranges.get(i).getEnd()));
            System.out.println(String.format("%-20s", AAranges.get(i).getBegin()+"-"+AAranges.get(i).getEnd()));
            System.out.print(String.format("%-52s",""));
        }
         System.out.println("");           
        }
    }
}
