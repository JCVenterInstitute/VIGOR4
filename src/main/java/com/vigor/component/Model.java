package com.vigor.component;

import java.util.ArrayList;
import java.util.HashMap;

import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Data
public class Model {
	
	//These exons are the Hsps converted to exons
	private ArrayList<Exon> Exons;
	private int strand;
	private ViralProtein viralProtein;
	private Alignment alignment;
	private HashMap<String,Float> scores;
	private String geneSymbol;	
	
}
