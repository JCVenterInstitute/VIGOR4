package com.vigor.utils.TBLParser;

import org.jcvi.jillion.core.Range;
import lombok.Data;

@Data
public class TBLFragment {
	private String proteinID;
	private String product;
	private String note;
	private String gene;
	private Range range;

}
