package com.vigor.component;

import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

import java.util.HashMap;

@Component
@Scope("prototype")
@Data
public class AlignmentFragment {

	private double score;
	private Direction direction;
	private Range proteinSeqRange;
	private Range nucleotideSeqRange;
    private Frame frame;
	private boolean isSubChain=false;

	
	
}
