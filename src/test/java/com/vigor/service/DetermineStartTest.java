package com.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Before;
import org.junit.Test;

import com.vigor.component.Alignment;
import com.vigor.component.Model;
import com.vigor.utils.VigorTestUtils;

public class DetermineStartTest {
	
	private List<Alignment> alignments;
	private List<Model> models=new ArrayList<Model>();
	private ModelGenerationService modelGenerationService = new ModelGenerationService();
	private ViralProteinService viralProteinService = new ViralProteinService();
	private DetermineStart determineStart = new DetermineStart();
	
	@Before
	public void getModel() {
			alignments = VigorTestUtils.getAlignments();
			alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment))
					.collect(Collectors.toList());
			models.addAll(modelGenerationService.alignmentToModels(alignments.get(0), "exonerate"));
	}
	
		
	@Test
	public void findStart() throws CloneNotSupportedException{
		List<String> startCodons = new ArrayList<String>();
		startCodons.add("ATG");
		Model model = models.get(0);
		determineStart.findStart(startCodons, model, "5");	    
		
		
	}
	
	
	
	

}
