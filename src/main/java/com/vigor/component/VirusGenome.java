package com.vigor.component;

import org.jcvi.jillion.core.Sequence;
import org.jcvi.jillion.core.residue.nt.Nucleotide;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class VirusGenome {

    // sequence, defline and id :-> single object from jillion
	private NucleotideSequence sequence;
	private String defline;
	private String id;
	private Boolean isComplete = false;
	private Boolean isCircular = false;

	public VirusGenome(NucleotideSequence sequence, String defline, boolean isComplete, boolean isCircular){

		this.sequence=sequence;
		this.defline=defline;
		this.isComplete=isComplete;
		this.isCircular=isCircular;
	}


	
	
}