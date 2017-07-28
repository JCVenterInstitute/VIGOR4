package org.jcvi.vigor.component;
import java.util.List;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class GeneStructure {

	private List<Exon> exons;
	private List<Intron> introns;
	
	
	
	
}
