package org.jcvi.vigor.service;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jcvi.jillion.core.Range;
import org.junit.Before;
import org.junit.Test;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;

/**
 * Created by snettem on 5/19/2017.
 */

public class ModelGenerationServiceTest {

	private List<Alignment> alignments;
	private List<Model> models;
	private ModelGenerationService modelGenerationService = new ModelGenerationService();
	private static ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	private static File file = new File(classLoader.getResource("vigorUnitTestInput/sequence.fasta"). getFile());

	@Before
	public void getModels() throws IOException {
		alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"flua_db",VigorUtils.getVigorWorkSpace());
		for(Alignment alignment:alignments){
		models = modelGenerationService.alignmentToModels(alignment, "exonerate");
		}
	}

	@Test
	public void alignmentToModelsTest() throws IOException {
		
		Alignment alignment = alignments.get(0);
		List<Model> outputModels = modelGenerationService.alignmentToModels(alignment, "exonerate");
		assertEquals(1, outputModels.size());

	}

	@Test
	public void splitModelAtSequencingGapsTest() {
		
        Model model = models.get(0);
        List<Range> sequenceGaps = new ArrayList<Range>();
        sequenceGaps.add(Range.of(593,769));
        sequenceGaps.add(Range.of(1400,1554));
        sequenceGaps.add(Range.of(2938,3010));
        sequenceGaps.add(Range.of(4550,4904));
        sequenceGaps.add(Range.of(8508,8698));
        sequenceGaps.add(Range.of(9906,9965));
        sequenceGaps.add(Range.of(11619,11759));
        
        List<Model> models = modelGenerationService.splitModelAtSequencingGaps(model, sequenceGaps);
        assertEquals(1,models.size());
        		
	}

}
