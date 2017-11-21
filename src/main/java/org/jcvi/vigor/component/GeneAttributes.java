package org.jcvi.vigor.component;


import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class GeneAttributes{

	private Ribosomal_Slippage ribosomal_slippage;
	private Splicing splicing;
	private StartTranslationException startTranslationException;
	private StopTranslationException stopTranslationException;
	private RNA_Editing rna_editing;
	private StructuralSpecifications structuralSpecifications;
	
		
}
