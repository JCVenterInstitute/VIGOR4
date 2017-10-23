package org.jcvi.vigor.service;


import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.junit.Before;
import org.junit.Test;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;

public class AdjustExonBoundariesTest {

	private List<Alignment> alignments;
	private List<Model> models=new ArrayList<Model>();
	private ModelGenerationService modelGenerationService = new ModelGenerationService();
	private ViralProteinService viralProteinService = new ViralProteinService();
	private AdjustExonBoundaries adjustExonBoundaries = new AdjustExonBoundaries();
	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	private File file = new File(classLoader.getResource("vigorUnitTestInput/sequence_veev.fasta"). getFile());
	
	@Before
	public void getModel() {
			alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"veev_db",VigorUtils.getVigorWorkSpace());
			alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment))
					.collect(Collectors.toList());
			models.addAll(modelGenerationService.alignmentToModels(alignments.get(0), "exonerate"));
	}
	
 /*   @Test
    public void adjustRibosomalSlippageTest(){
    	adjustExonBoundaries.adjustRibosomalSlippage(models.get(0));
    	
    	
    }*/
    /*@Test
    public void adjustSpliceSitesTest(){
    	adjustExonBoundaries.adjustSpliceSites(models.get(0));
    	
    	
    }
	*/
	
	@Test
    public void SpliceSitesTest(){
    	List<Model> outModels = adjustExonBoundaries.adjustSpliceSites(models.get(0));
    	Comparator<Model> bySpliceScore = (Model m1,Model m2)->m1.getScores().get("spliceScore").compareTo(m2.getScores().get("spliceScore"));
    	Optional<Model> outModel = outModels.stream().sorted(bySpliceScore.reversed()).findFirst();
    	assertEquals(13,outModels.size());
    	 	
    }
	
}
