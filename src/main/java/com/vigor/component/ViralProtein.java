package com.vigor.component;

import org.springframework.stereotype.Component;
import org.jcvi.jillion.core.Sequence;

import lombok.Data;

@Component
@Data
public class ViralProtein {
	
	private GeneStructure geneStructure;
	private Sequence<String> sequence;
    private ViralTricks viralTricks;
    
    public ViralProtein(Sequence<String> sequence)
    {
    	this.sequence=sequence;
    }
    
}
