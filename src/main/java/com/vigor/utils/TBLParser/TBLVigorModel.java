package com.vigor.utils.TBLParser;

import java.util.List;

import com.vigor.component.Exon;

import lombok.Data;

@Data
public class TBLVigorModel {

	private String virusGenomeID;
	private String viralProteinID;
	private List<Exon> exons;
}
