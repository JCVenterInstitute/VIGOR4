package org.jcvi.vigor.service;
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.Triplet;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Exon;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.component.ViralProtein;
import org.jcvi.vigor.component.VirusGenome;


public class DetermineStartTest {
	
	private Model model = new Model();
	private DetermineStart determineStart = new DetermineStart();
	
	
	@Before
	public void getModel() {
	    VirusGenome virusGenome = new VirusGenome();
		NucleotideSequence seq = new NucleotideSequenceBuilder("CGAAGGCTGGCCGATAGAAAACAGAAACTAAGCCAAGCAAGCAACAAACGAGACATCAGCAGTGATGC"
				+ "TGATTATGAAAATGATGATGATGCTACAGCGGCTGCAGGGATAGGAGGAATTTAACAGGATAATTGGACA"
				+ "GTAGAAACCAGATCAAAAGTAAGAAAAACTTAGGGTGAATGGCAATTCACAGATCAGCTCAACCAGACAT"
				+ "CATCAGCATACACGAAACCAACCTTCACAGTGGATACCTCAGCATCCAAAACTCTCCTTCCCGAATGGAT"
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
	    ViralProtein viralProtein = new ViralProtein();
	    List<Exon> exons = new ArrayList<Exon>();
		Exon exon = new Exon();
		exon.setRange(Range.of(236,882));
		exon.setFrame(Frame.ONE);
		exons.add(exon);
		model.setExons(exons);
		Alignment alignment = new Alignment();
		alignment.setVirusGenome(virusGenome);
		alignment.setViralProtein(viralProtein);
		model.setAlignment(alignment);
	}	
	/*	
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
		assertEquals(1,outputModels.size());
	}*/
}
