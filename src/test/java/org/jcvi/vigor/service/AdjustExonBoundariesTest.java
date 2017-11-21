package org.jcvi.vigor.service;


import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.AppConfig;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;


@RunWith(SpringRunner.class)
@ContextConfiguration(classes = AppConfig.class)
public class AdjustExonBoundariesTest {

	private List<Alignment> alignments;
	private List<Model> models=new ArrayList<Model>();
	@Autowired
	private ModelGenerationService modelGenerationService;
	@Autowired
	private ViralProteinService viralProteinService;
	@Autowired
	private AdjustExonBoundaries adjustExonBoundaries;
	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	private File file = new File(classLoader.getResource("vigorUnitTestInput/chikv.ungapped.fasta.JF274082.1.ref.fasta"). getFile());
	
	@Before
	public void getModel() {
			alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"chikv_db",VigorUtils.getVigorWorkSpace(),"gi|392976865|ref|YP_006491243.1|");
			alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment,new VigorForm()))
					.collect(Collectors.toList());
			alignments.stream().forEach(x -> {
				if(x.getViralProtein().getProteinID().equals("gi|392976865|ref|YP_006491243.1|")){
					models.addAll(modelGenerationService.alignmentToModels(x, "exonerate"));
				}
			});
		
	}
	
    @Test
    public void adjustRibosomalSlippageTest() throws CloneNotSupportedException{
    	Optional<Model> testModel = models.stream().filter(x->x.getAlignment().getViralProtein().getProteinID().equals("gi|392976865|ref|YP_006491243.1|")).findFirst();
    	List<Model> models=null;
    	if(testModel.isPresent()){
    	models = adjustExonBoundaries.adjustRibosomalSlippage(testModel.get());
    	}
    	Range actual = models.get(0).getExons().get(0).getRange();    
    	assertEquals(Range.of(7553,9943),actual);
    }
    
    @Test
    public void adjustSpliceSitesTest() throws CloneNotSupportedException{
    	
    	adjustExonBoundaries.adjustSpliceSites(models.get(0));
    	
    	
    }
	
	/*
	@Test
    public void SpliceSitesTest() throws CloneNotSupportedException{
    	List<Model> outModels = adjustExonBoundaries.adjustSpliceSites(models.get(0));
    	Comparator<Model> bySpliceScore = (Model m1,Model m2)->m1.getScores().get("spliceScore").compareTo(m2.getScores().get("spliceScore"));
    	Optional<Model> outModel = outModels.stream().sorted(bySpliceScore.reversed()).findFirst();
    	assertEquals(13,outModels.size());
    	 	
    }
	*/
}
