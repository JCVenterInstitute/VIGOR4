package org.jcvi.vigor.service;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.junit.Before;
import org.junit.Test;

import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.VigorTestUtils;

public class DetermineMissingExonsServiceTest {
	private List<Alignment> alignments;
	private List<Model> models;
	private ModelGenerationService modelGenerationService = new ModelGenerationService();
	private DetermineMissingExonsService determineMissingExonsService = new DetermineMissingExonsService();
	private ViralProteinService viralProteinService = new ViralProteinService();
	
	@Before
	public void getModel(){
		models = new ArrayList<Model>();
		alignments = VigorTestUtils.getAlignments();
		alignments = alignments.stream()
				.map(alignment -> viralProteinService.setViralProteinAttributes(alignment))
				.collect(Collectors.toList());
		for(Alignment alignment:alignments){
		models.addAll(modelGenerationService.alignmentToModels(alignment, "exonerate"));
		models.stream().forEach(System.out::println);;
		}
	}
	
		
	@Test
	public void findMissingExonsTest(){
		Model model = models.get(0);
		model.getExons().remove(1);
		int missingExons = determineMissingExonsService.findMissingExonRanges(model, 65, 50).size();
		assertEquals(1,missingExons);
		
	}
	
	@Test 
	public void performPairWiseAlignment(){
		Model model = models.get(0);
		model.getExons().remove(1);
		Map<Range,Range> missingExons = determineMissingExonsService.findMissingExonRanges(model, 65, 50);
		NucleotideSequence NTSequence = model.getAlignment().getVirusGenome().getSequence();
		ProteinSequence AASequence = model.getAlignment().getViralProtein().getSequence();
		for (Map.Entry<Range, Range> entry : missingExons.entrySet()) {
			
			determineMissingExonsService.getPairWiseAlignment(entry.getKey(),entry.getValue(),NTSequence,AASequence);
			
			
			
			
		}
			
			
		
	}

}
