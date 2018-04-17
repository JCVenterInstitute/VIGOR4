package org.jcvi.vigor.service;

import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.DirectedRange;
import org.jcvi.jillion.core.Direction;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.ServiceException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.jcvi.vigor.forms.VigorForm;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;


@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class DetermineMissingExonsTest {
	@Autowired
	private ModelGenerationService modelGenerationService;
	@Autowired
	private DetermineMissingExons determineMissingExons ;
	@Autowired
	private ViralProteinService viralProteinService ;


	@Test
	public void findMissingExonsWithSpliceFormPresent() throws VigorException {
		ClassLoader classLoader = VigorTestUtils.class.getClassLoader();
		File file = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua.fasta"). getFile());
		List<Alignment> alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"flua_db",
				VigorUtils.getVigorWorkSpace(),null);
		List<Model> models = new ArrayList<>();
		for (int i=0;i<alignments.size();i++) {
			alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), new VigorForm()));
		}
		models.addAll(modelGenerationService.alignmentToModels(alignments.get(0), "exonerate"));
		assertTrue(String.format("Expected at least 1 model, got %s", models.size()), 1 >= models.size());
		Model model = models.get(0);
		assertTrue(String.format("Expected models %s to have at least 2 exons, got %s", model, model.getExons().size()),
				2 >= model.getExons().size());
		model.getExons().remove(1);
		int exons = determineMissingExons.findMissingExons(model).getExons().size();
		assertEquals(2, exons);
	}

	@Test
	public void performPairWiseAlignment() {

		Range NTRange = Range.of(579, 2199);
		Range AARange = Range.of(190, 231);

		NucleotideSequence NTSequence = new NucleotideSequenceBuilder(
				"TGATCCAAAATGGAAGATTTTGTGCGACAATGCTTCAATCCAATGATTGTCGAGCTTGCGGAAAAGGCAATGAAAGAATATGGGGAAGATCCGAAAATCGAAACGAACAAATTTGCCGCAATAT"
						+ "GCACACACTTAGAGGTCTGTTTCATGTATTCGGATTTCCACTTTATTGATGAACGGGGCGAATCAATAATTGTAGAATCTGGCGATCCAAATGCATTATTGAAACACCGATTTGAGATAATTGAAGGGAGAGACCGAA"
						+ "CGATGGCCTGGACAGTGGTGAATAGTATCTGCAACACCACAGGAGTCGAGAAACCTAAATTTCTCCCAGATTTGTATGACTACAAAGAGAATCGATTCATTGAAATTGGAGTAACACGGAGGGAAGTTCATATATAC"
						+ "TATCTAGAAAAGGCCAACAAGATAAAATCAGAGAAGACACACATTCACATATTCTCATTCACTGGAGAGGAAATGGCCACCAAAGCGGACTACACTCTTGACGAAGAGAGTAGGGCAAGAATCAAAACCAGGCTGTTC"
						+ "ACTATAAGGCAGGAAATGGCCAGTAGGGGTCTATGGGATTCCTTTCGTCAGTCCGAGAGAGGCGAAGAGACAGTTGAAGAAAGATTTGAAATCACAGGAACCATGCGCAGGCTTGCCGACCAAAGTCTCCCACCGAACT"
						+ "TCTCCAGCCTTGAAAACTTTAGAGCCTATGTGGATGGATTCGAACCGAACGGCTGCATTGAGGGCAAGCTTTCTCAAATGTCAAAAGAAGTGAACGCCCGAATTGAGCCATTTCTGAAGACAACACCACGCCCTCTCA"
						+ "AACTACCTGACGGGCCTCCCTGCTCTCAACGGTCGAAGTTCCTGCTGATGGATGCCCTTAAATTAAGCATCGAAGACCCGAGTCATGAGGGGGAGGGTATACCGCTATATGATGCAATCAAATGCATGAAGACATTTTT"
						+ "CGGCTGGAAAGAGCCCAACATTGTAAAACCACATGAAAAGGGCATAAACCCCAATTACCTCCTGGCTTGGAAGCAAGTGCTGGCAGAACTCCAAGATATTGAAAATGAGGAGAAAATCCCAAAAACAAAGAACATGAAGAA"
						+ "AACGAGCCAGTTGAAGTGGGCACTTGGTGAGAATATGGCACCGGAGAAGGTAGACTTTGAGGATTGCAAGGATGTTAGCGATCTGAGACAGTATGACAGTGATGAACCAGAGTCTAGATCGCTAGCAAGCTGGATCCAGAGT"
						+ "GAATTCAACAAGGCATGTGAATTGACAGATTCAAGTTGGATTGAGCTTGATGAAATAGGGGAAGACATTGCTCCAATTGAGCACATTGCGAGTATGAGAAGAAACTACTTCACAGCGGAAGTATCCCATTGCAGGGCTACTGAA"
						+ "TACATAATGAAAGGAGTGTACATAAACACAGCCTTGTTGAATGCATCCTGTGCAGCCATGGATGACTTCCAACTGATTCCAATGATAAGCAAATGCAGGACCAAAGAAGGGAGGCGGAAGACTAATCTGTATGGATTCATTATA"
						+ "AAAGGAAGATCCCATTTGAGAAATGACACCGATGTAGTAAACTTTGTGAGCATGGAATTCTCTCTTACTGACCCGAGGCTGGAGCCACACAAGTGGGAAAAGTACTGTGTTCTCGAGATAGGAGACATGCTCCTACGGACTGC"
						+ "AATAGGCCAAGTGTCAAGGCCCATGTTCCTGTATGTGAGAACCAATGGGACTTCCAAGATCAAGATGAAGTGGGGCATGGAAATGAGGCGATGCCTTCTTCAATCCCTTCAACAAATTGAGAGCATGATTGAAGCCGAGTCTTC"
						+ "TGTCAAAGAGAAGGACATGACCAAAGAATTCTTTGAAAACAAATCAGAAACATGGCCAATTGGAGAGTCACCCAAAGGGGTGGAGGAAGGCTCCATTGGGAAGGTGTGCAGAACCTTACTGGCAAAATCTGTATTCAACAGCCTATA"
						+ "TGCATCTCCACAACTCGAGGGATTTTCAGCTGAATCAAGAAAGTTGCTTCTCATTGTCCAGGCACTTAGGGACAACCTGGAACCTGGGACCTTCGATCTTGGGGGGCTATATGAAGCAATTGAGGAGTGCCTGATTAATGATCCCTGGG"
						+ "TTTTGCTTAATGCGTCTTGGTTCAACTCCTTCCTCACACATGCACTGAAATAGTTGTGGCAATGCTACTATTTGCTATCCATACTGTCCAAAA")
				.build();

		ProteinSequence AASequence = new ProteinSequenceBuilder(
				"MEDFVRQCFNPMIVELAEKTMKEYGEDLKIETNKFAAICTHLEVCFMYSDFHFINEQGESIIVELGDPNALLKHRFEIIEGRDRTMAWTVVNSICNTTGAEKPKFLPDLYDYKENRFIEIGVTRREVHIYYLEKANKI"
						+ "KSEKTHIHIFSFTGEEMATKADYTLDEESRARIKTRLFTIRQEMASRGLWDSFVSPREEKRQLKKGLKSQEQCASLPTKVSRRTSPALKILEPM")
				.build();

		Exon exon = determineMissingExons.performJillionPairWiseAlignment(NTRange, AARange, NTSequence,
				AASequence,Direction.FORWARD);
		assertEquals(exon.getAlignmentFragment().getProteinSeqRange(), AARange);

	}

	@Test
	public void testJillionPairwiseAlignment(){
		ProteinSequence querySequence = new ProteinSequenceBuilder("MEDFVRQCFNPMIVELAEKTMKEYGEDLKIETNKFAAICTHLEVCFMYSDFHFI").build();
		ProteinSequence subjectSequence = new ProteinSequenceBuilder("MEDFVRQCFNPMIVELAEKTMKEYGEDLKIETNKFAAICTHLEVCFMYSDFHFINEQGESIIVELGDPNALLKHRFEIIEGRDRTMAWTVVNSICNTTGAEKPKF").build();
		AminoAcidSubstitutionMatrix blosom50 = BlosumMatrices.blosum50();
		ProteinPairwiseSequenceAlignment actual = PairwiseAlignmentBuilder
				.createProtienAlignmentBuilder(querySequence, subjectSequence, blosom50).gapPenalty(-8, -8)
				.build();
		DirectedRange expected;
		expected =DirectedRange.create((Range.of(0,53)),Direction.FORWARD );
		DirectedRange queryRange = actual.getQueryRange();
		assertEquals(expected,queryRange);

	}
	

}
