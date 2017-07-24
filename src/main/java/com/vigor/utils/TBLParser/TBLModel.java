package com.vigor.utils.TBLParser;

import java.util.List;

import com.vigor.component.Exon;

import lombok.Data;

@Data
public class TBLModel {
	
	private String virusGenomeID;
	private List<TBLFragment> TBLFragments;
	   

}
