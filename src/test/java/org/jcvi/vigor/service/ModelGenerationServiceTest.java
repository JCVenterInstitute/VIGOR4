package org.jcvi.vigor.service;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.junit.Test;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.utils.VigorTestUtils;

/**
 * Created by snettem on 5/19/2017.
 */

public class ModelGenerationServiceTest {

	private List<Alignment> alignments;
	private ModelGenerationService modelGenerationService = new ModelGenerationService();
    private ExonerateService exonerateService = new ExonerateService();
	private ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	private File file = new File(classLoader.getResource("vigorUnitTestInput/exonerate_flua.txt"). getFile());
	@Test
	public void alignmentToModelsTest() throws IOException {
		
		alignments = exonerateService.parseExonerateOutput(file, new AlignmentEvidence("flua_db"), new VirusGenome());
		Alignment alignment = alignments.get(0);
		List<Model> outputModels = modelGenerationService.alignmentToModels(alignment, "exonerate");
		assertEquals(1, outputModels.size());

	}

	@Test
	public void splitModelAtSequenceGapsTest() {
		
        Model model = new Model();
        List<Exon> exons = new ArrayList<Exon>();
        exons.add(new Exon(Range.of(20, 500),Frame.ONE));
        exons.add(new Exon(Range.of(600,1300),Frame.TWO));
        exons.add(new Exon(Range.of(1600,2800),Frame.ONE));
        exons.add(new Exon(Range.of(3011,4500),Frame.ONE));
        model.setExons(exons);
        model.setStatus(new ArrayList<String>());
        List<Range> sequenceGaps = new ArrayList<Range>();
        sequenceGaps.add(Range.of(1400,1554));
        sequenceGaps.add(Range.of(2938,3010));
        sequenceGaps.add(Range.of(4550,4904));
        sequenceGaps.add(Range.of(8508,8698));
        sequenceGaps.add(Range.of(9906,9965));
        sequenceGaps.add(Range.of(11619,11759));
        List<Model> models = modelGenerationService.splitModelAtSequenceGaps(model, sequenceGaps);
        assertEquals(3,models.size());
        	
	}

}
