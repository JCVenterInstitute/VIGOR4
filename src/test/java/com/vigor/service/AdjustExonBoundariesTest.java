package com.vigor.service;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Before;

import com.vigor.component.Alignment;
import com.vigor.component.Model;
import com.vigor.utils.VigorTestUtils;

public class AdjustExonBoundariesTest {

	private List<Alignment> alignments;
	private List<Model> models=new ArrayList<Model>();
	private ModelGenerationService modelGenerationService = new ModelGenerationService();
	private ViralProteinService viralProteinService = new ViralProteinService();
	private AdjustExonBoundaries adjustExonBoundaries = new AdjustExonBoundaries();
	
	@Before
	public void getModel() {
			alignments = VigorTestUtils.getAlignments();
			alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment))
					.collect(Collectors.toList());
			models.addAll(modelGenerationService.alignmentToModels(alignments.get(0), "exonerate"));
	}
	
	
	
}
