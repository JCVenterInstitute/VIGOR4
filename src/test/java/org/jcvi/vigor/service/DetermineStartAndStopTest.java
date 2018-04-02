package org.jcvi.vigor.service;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.Triplet;
import org.jcvi.vigor.AppConfig;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.AlignmentFragment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;


@RunWith(SpringRunner.class)
@ContextConfiguration(classes = AppConfig.class)
public class DetermineStartAndStopTest {
	
	private Model model = new Model();
	@Autowired
	private DetermineStart determineStart ;
	@Autowired
	private DetermineStop determineStop;
	
	@Before
	public void getModel() {
	    VirusGenome virusGenome = new VirusGenome();
		NucleotideSequence seq = new NucleotideSequenceBuilder("CGAAGGCTGGCCGATAGAAAACAGAAACTAAGCCAAGCAAGCAACAAACGAGACATCAGCAGTGATGC"
				+ "TGATTATGAAAATGATGATGATGCTACAGCGGCTGCAGGGATAGGAGGAATTTAACAGGATAATTGGACA"
				+ "GTAGAAACCAGATCAAAAGTAAGAAAAACTTAGGGTGAATGGCAATTCACAGATCAGCTCAACCAGACAT"
				+ "CATCAGCATACACGAAACCAACCTTCACATGGGATACCTCAGCATCCAAAACTCTCCTTCCCGAATGGAT"
				+ "CAGGATGCCTTCTTTTTTGAGAGGGATCCTGAGGCCGAAGGAGAGGCACCACGAAAACAAGAATCACTCT"
				+ "CAGATGTCATCGGACTCCTTGACGTCGTCCTATCCTACAAGCCCACCGAAATTGGAGAAGACAGAAGCTG"
				+ "GCTCCATAGTATCATCGACAACCCAAAAGAAAACAAGTCATCATGCAAATCTGACGATAACGATAAAGAC"
				+ "AGAGCAATCTCGACGTCGACCCAAGATCATAGATCAAGTGAGGAGAGTAGAGTCTCTAGGAGAACAGGTG"
				+ "AGTCAAAAACAGAGACACATGCTAGAATCCTTGATCAACAAGGTGTACACAGGGCCTCTAGGCGAGGAAC"
				+ "TAGTCCAAACCCTCTACCTGAGAATATGGGCAATGAAAGAAACACCAGAATAGAGGAAGATCCTTCAAAT"
				+ "GAGAGAAGACATCAGAGATCAGTATCTACGGNNNNNNNNNNNNNNNNNNNNNNNNNNTTTAATAAGAGGG"
				+ "AAGAAGACCAAGTTGAGGGATTTCCAGAAGAGGTACGAGGAAGTACATCCTTATCTGATGATGGAGAGAG"
				+ "TAGAACAAATAATAATGGAAGAAGCATGGAAACTAGCAGCACACATAGTACAAGAATAACTGATGTCATT"
				+ "ACCAACCCAAGTCCAGAGCTTGAAGATGCCGTTCTACAAAGGAATAAAAGACGGCCGACGACCATCAAGC").build();
		virusGenome.setSequence(seq);
		System.out.println(seq.getLength());
	    ViralProtein viralProtein = new ViralProtein();
	    //Dummy protein sequence
	    ProteinSequence proteinSeq = new ProteinSequenceBuilder("MNRPFFQNFGRRPFPAPSIAWRPRRRRSAAPAFQAPRFGLANQIQQLTSAVSALVIGQSA"
                                                                +"RTQPPRARQQPRRPPPKKKQPPPPKKEKPTKKQPKKPTKPKPGKRQRMVLKLEADRLFDV"
                                                                +"KNEQGDVVGHALAMEGKVMKPLHVKGTIDHPVLSKLKFTKSSAYDMEFAPLPVNMKSEAF"
                                                                +"NYTSEHPEGFYNWHHGAVQYSGGRFTVPRGVGGKGDSGRPIMDNTGKVVAIVLGGADEGA"
                                                                +"RTALSVVTWNAKGKTIKTTPEGTEEWSAAAAITTLCLIGNMTFPCDRPPTCYNVNPATTL"
                                                                +"DILEQNVDHPLYDTLLTSITRCSSRRHKRSITDDFTLTSPYLGTCSYCHHTEPCFSPIKI"
                                                                +"EQVWDDADDGTIRIQTSAQFGYNQNGAADN").build();
	    viralProtein.setSequence(proteinSeq);
	    List<Exon> exons = new ArrayList<Exon>();
		Exon exon = new Exon();
		exon.setRange(Range.of(236,732));
		exon.setFrame(Frame.ONE);
		Exon exon1 = new Exon();
		exon1.setFrame(Frame.TWO);
		exon1.setRange(Range.of(756,977));
		AlignmentFragment frag = new AlignmentFragment();
		AlignmentFragment frag1 = new AlignmentFragment();
		frag.setProteinSeqRange(Range.of(0,166));
		frag1.setProteinSeqRange(Range.of(167,389));
		exon.setAlignmentFragment(frag);
		exon1.setAlignmentFragment(frag1);
		exons.add(exon);
		exons.add(exon1);
		model.setExons(exons);
		Alignment alignment = new Alignment();
		alignment.setVirusGenome(virusGenome);
		alignment.setViralProtein(viralProtein);
		model.setAlignment(alignment);
	}	
		
	@Test
	public void findStart() throws CloneNotSupportedException{
		List<Triplet> startCodons = new ArrayList<Triplet>();
	    Triplet triplet4 = Triplet.create('A','T','G');
	    Triplet triplet1 = Triplet.create('A', 'C', 'G');
	    Triplet triplet2 = Triplet.create('C', 'C', 'G');
	    Triplet triplet3 = Triplet.create('G', 'T', 'G');
	    startCodons.add(triplet1);
	    startCodons.add(triplet2);
	    startCodons.add(triplet3);
	    startCodons.add(triplet4);
	    List<Model> outputModels = determineStart.findStart(startCodons, model, "50");		
	    assertEquals(2,outputModels.size());
	}
	
	@Test
	public void findStop() throws CloneNotSupportedException{
	  List<Model> models = determineStop.findStop(model);
	  assertEquals(2,models.size());
	  //stop codons in frame with the last exon are at 927,951;
	  assertEquals(951,models.get(0).getExons().get(model.getExons().size()-1).getRange().getEnd());
	  assertEquals(927,models.get(1).getExons().get(model.getExons().size()-1).getRange().getEnd());
	}
	
}
