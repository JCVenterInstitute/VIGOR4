package com.vigor.component;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


import org.jcvi.jillion.core.Direction;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Model {
	
	//These exons are the Hsps converted to exons
	private List<Exon> Exons;
	private Alignment alignment;
	private Map<String,Float> scores;
	private String geneSymbol;
	private List<String> status;
	private Direction direction;



}
