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
import org.junit.Test;
import org.junit.runner.RunWith;
import org.jcvi.vigor.AppConfig;
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
@ContextConfiguration(classes = AppConfig.class)
public class DetermineMissingExonsTest {
	private List<Alignment> alignments;
	private List<Model> models= new ArrayList<Model>();
	@Autowired
	private ModelGenerationService modelGenerationService;
	@Autowired
	private DetermineMissingExons determineMissingExons ;
	@Autowired
	private ViralProteinService viralProteinService ;
	
	
	@Test
	public void findMissingExonsWithSpliceFormPresent() {
	    ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	    File file = new File(classLoader.getResource("vigorUnitTestInput/sequence_flua.fasta"). getFile());
		alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"flua_db",VigorUtils.getVigorWorkSpace(),null);
		alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment,new VigorForm()))
				.collect(Collectors.toList());
		models.addAll(modelGenerationService.alignmentToModels(alignments.get(0), "exonerate"));
		Model model = new Model();
		model = models.get(0);
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
	
	/*@Test
	public void findMissingExonsWithSpliceFormAbsent(){
		ClassLoader classLoader = VigorTestUtils.class.getClassLoader(); 
	    File file = new File(classLoader.getResource("vigorUnitTestInput/sequence_veev.fasta"). getFile());
		alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"veev_db",VigorUtils.getVigorWorkSpace());
	    alignments = alignments.stream().map(alignment -> viralProteinService.setViralProteinAttributes(alignment))
				.collect(Collectors.toList());
		models.addAll(modelGenerationService.alignmentToModels(alignments.get(0), "exonerate"));
		Model model = new Model();
		model = models.get(0);
		model.getExons().remove(0);
	    Model outputModel = determineMissingExons.findMissingExonsWithSpliceFormAbsent(model);
	    assertEquals(12,outputModel.getExons().size());
	}
	*/
	
	
}
