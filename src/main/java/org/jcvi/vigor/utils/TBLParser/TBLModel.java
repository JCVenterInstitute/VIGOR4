package org.jcvi.vigor.utils.TBLParser;

import java.util.List;

import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.component.Exon;

import lombok.Data;

@Data
public class TBLModel {
	
	private String virusGenomeID;
	private List<Exon> exons;
	private String viralProteinID;
	private String geneID;
	private String product;
	private String note;
	private String gene;
	private boolean isPseudoGene=false;
	private boolean is5Partial=false;
	private boolean is3Partial=false;
	private boolean isRiboSlippage=false;
	private Range stopCodonReadThrough;

}
