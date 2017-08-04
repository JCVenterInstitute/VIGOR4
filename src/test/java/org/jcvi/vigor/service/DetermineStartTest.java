package org.jcvi.vigor.service;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.Before;
import org.junit.Test;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;

public class DetermineStartTest {
	
	private List<Alignment> alignments;
	private List<Model> models=new ArrayList<Model>();
	private ModelGenerationService modelGenerationService = new ModelGenerationService();
	private ViralProteinService viralProteinService = new ViralProteinService();
	private DetermineStart determineStart = new DetermineStart();
	private static ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	private static File file = new File(classLoader.getResource("vigorUnitTestInput/sequence.fasta"). getFile());
	
	@Before
	public void getModel() {
			alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"flua_db",VigorUtils.getVigorWorkSpace());
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
