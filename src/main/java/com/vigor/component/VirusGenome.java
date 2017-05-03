package com.vigor.component;

import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.residue.nt.Nucleotide;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class VirusGenome {
	
	private Sequence<Nucleotide> sequence;
	private String defline;
	private Boolean isComplete = false;
	private Boolean isCircular = false;

	public VirusGenome(Sequence<Nucleotide> sequence, String defline, boolean isComplete, boolean isCircular){

		this.sequence=sequence;
		this.defline=defline;
		this.isComplete=isComplete;
		this.isCircular=isCircular;
	}


	
	
}