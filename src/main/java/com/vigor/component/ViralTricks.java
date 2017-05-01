package com.vigor.component;

import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Data
public class ViralTricks {

	private int geneVariation;
	private String StopCodon_ReadThru;
	
	//splicing
	private boolean Splicing;
	private String SpliceForm;
	private String nonCanonical_splicing;
	
	//ribosomal Slippage
	private boolean ribosomal_slippage;
	private String slippage_motif;
	private int slippageOffset;
	private int slippage_frameshift;
	
	
	
	
	
	
	
}
