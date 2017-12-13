package org.jcvi.vigor.service;


import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.AppConfig;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

@RunWith(SpringRunner.class)
@ContextConfiguration(classes = AppConfig.class)
public class AdjustUneditedExonBoundariesTest {

	private List<Alignment> alignments;
	private List<Model> models=new ArrayList<Model>();
	@Autowired
	private ModelGenerationService modelGenerationService;
	@Autowired
	private ViralProteinService viralProteinService;
	@Autowired
	private AdjustUneditedExonBoundaries adjustUneditedExonBoundaries;
	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	private File file = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua_1.fasta"). getFile());

	@Before
	public void getModel() {
			alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"flua_db",VigorUtils.getVigorWorkSpace(),"seg7prot2ADR2_S31N");
			alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment,new VigorForm()))
					.collect(Collectors.toList());
			alignments.stream().forEach(x -> {
					models.addAll(modelGenerationService.alignmentToModels(x, "exonerate"));
			});
		
	}
	
	@Test
	public void adjustSpliceSitesTest() throws CloneNotSupportedException{
		   Model testModel = models.get(0);
		   testModel.getExons().get(0).setRange(Range.of(8,31));
		   List<Model> outModels = adjustUneditedExonBoundaries.adjustSpliceSites(testModel);
	       Comparator<Model> bySpliceScore = (Model m1,Model m2)->m1.getScores().get("spliceScore").compareTo(m2.getScores().get("spliceScore"));
	       Optional<Model> outModel = outModels.stream().sorted(bySpliceScore.reversed()).findFirst();
	       assertEquals(Range.of(8,33),outModel.get().getExons().get(0).getRange());
	}
	
	
}
