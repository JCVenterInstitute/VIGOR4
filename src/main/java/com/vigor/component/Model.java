package com.vigor.component;

import java.util.ArrayList;
import java.util.HashMap;


import org.jcvi.jillion.core.Direction;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Scope("prototype")
@Data
public class Model {
	
	//These exons are the Hsps converted to exons
	private ArrayList<Exon> Exons;
	private Direction direction;
	private ViralProtein viralProtein;
	private Alignment alignment;
	private HashMap<String,Float> scores;
	private String geneSymbol;	
	
}
