package com.vigor.component;
import java.util.ArrayList;

import org.springframework.stereotype.Component;

import lombok.Data;

@Component
@Data
public class GeneStructure {

	private ArrayList<Exon> exons;
	private ArrayList<Intron> introns;
	
	
	
	
}
