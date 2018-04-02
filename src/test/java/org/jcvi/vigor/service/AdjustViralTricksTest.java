package org.jcvi.vigor.service;


import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
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
public class AdjustViralTricksTest {


	@Autowired
	private ModelGenerationService modelGenerationService;
	@Autowired
	private ViralProteinService viralProteinService;
	@Autowired
	private AdjustViralTricks adjustViralTricks;
	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader();

    @Test
    public void adjustRibosomalSlippageTest() throws CloneNotSupportedException{
<<<<<<< HEAD
	    List<Alignment> alignments;
	    List<Model> models=new ArrayList<Model>();
        File file = new File(classLoader.getResource("vigorUnitTestInput/Flua_RiboSlippage_Test.fasta"). getFile());
        alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"flua_db",VigorUtils.getVigorWorkSpace(),"seg3prot2A");
        alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment,new VigorForm()))
                .collect(Collectors.toList());
        alignments = modelGenerationService.mergeIdenticalProteinAlignments(alignments);
        alignments.stream().forEach(x -> {
            models.addAll(modelGenerationService.alignmentToModels(x, "exonerate"));
        });
    	Model testModel = models.get(0);
       	List<Model> outputModels = adjustViralTricks.adjustRibosomalSlippage(testModel);
       	Range actual = outputModels.get(0).getExons().get(0).getRange();
    	assertEquals(Range.of(9,581),actual);
    }

    @Test
    public void checkForLeakyStopTest(){
        List<Alignment> alignments;
        List<Model> models=new ArrayList<Model>();
        File file = new File(classLoader.getResource("vigorUnitTestInput/Veev_StopTranslationEx_Test.fasta"). getFile());
        alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"veev_db",VigorUtils.getVigorWorkSpace(),"399240871_NSP");
        alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment,new VigorForm()))
                .collect(Collectors.toList());
        alignments = modelGenerationService.mergeIdenticalProteinAlignments(alignments);
        alignments.stream().forEach(x -> {
            models.addAll(modelGenerationService.alignmentToModels(x, "exonerate"));
        });
        Model testModel = models.get(0);
        Model outputModel = adjustViralTricks.checkForLeakyStop(testModel);
        Range actual = outputModel.getExons().get(0).getRange();
        assertEquals(Range.of(9,581),actual);
=======
		assertTrue(String.format("No models created for file %s and chikv_db", file), models.size() > 0);
    	Model testModel = models.get(0);
    	List<Model> adjustedModels= adjustViralTricks.adjustRibosomalSlippage(testModel);
    	assertTrue(String.format("No adjusted models returned for test model %s", testModel), adjustedModels.size() > 0);
		assertTrue(String.format("No exons for adjusted model %s", adjustedModels.get(0)), adjustedModels.get(0).getExons().size() > 0);
       	Range actual = adjustedModels.get(0).getExons().get(0).getRange();
    	assertEquals(Range.of(7553,9943),actual);
>>>>>>> 66dc8b68e7543c0bd7eca9a56916a7ad74428382
    }



	@Test
    public void adjustRNAEditingTest() throws CloneNotSupportedException{
        List<Alignment> alignments;
        List<Model> models=new ArrayList<Model>();
        File file = new File(classLoader.getResource("vigorUnitTestInput/mmp_rna_editing_Test.fasta"). getFile());
        alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"mmp_db",VigorUtils.getVigorWorkSpace(),"AEY76114.1");
        alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment,new VigorForm()))
                .collect(Collectors.toList());
        alignments = modelGenerationService.mergeIdenticalProteinAlignments(alignments);
        alignments.stream().forEach(x -> {
            models.addAll(modelGenerationService.alignmentToModels(x, "exonerate"));
        });
        Model testModel = models.get(0);
        List<Model> outputModels = adjustViralTricks.adjustRNAEditing(testModel);
        Range actual = outputModels.get(0).getExons().get(0).getRange();
        assertEquals(Range.of(1861,2322),actual);

    }
}
