package com.vigor.component;

import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import org.jcvi.jillion.core.Sequence;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class ViralProtein {
	
	private GeneStructure geneStructure;
	private ProteinSequence sequence;
    private GeneAttributes geneAttributes;
    private String proteinID;
    private String defline;



    
}
