package com.vigor.service;

import static org.junit.Assert.assertEquals;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.jcvi.jillion.core.Range;
import org.junit.Before;
import org.junit.Test;
import com.vigor.component.Alignment;
import com.vigor.component.Model;
import com.vigor.utils.VigorTestUtils;

/**
 * Created by snettem on 5/19/2017.
 */

public class ModelGenerationServiceTest {

	private List<Alignment> alignments;
	private List<Model> models;
	private ModelGenerationService modelGenerationService = new ModelGenerationService();

	@Before
	public void getModels() throws IOException {
		alignments = VigorTestUtils.getAlignments();
		for(Alignment alignment:alignments){
		models = modelGenerationService.alignmentToModels(alignment, "exonerate");
		}
	}

	@Test
	public void alignmentToModelsTest() throws IOException {
		
		Alignment alignment = alignments.get(0);
		List<Model> outputModels = modelGenerationService.alignmentToModels(alignment, "exonerate");
		assertEquals(2, outputModels.size());

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