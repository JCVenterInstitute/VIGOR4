package com.vigor.component;

import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Exon {

	private Range range;
	private Frame frame;
	private boolean is_5p_adjusted=false;
	private boolean is_3p_adjusted=false;

	
}
