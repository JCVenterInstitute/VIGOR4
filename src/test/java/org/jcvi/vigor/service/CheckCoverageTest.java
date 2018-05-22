package org.jcvi.vigor.service;

import org.jcvi.jillion.align.AminoAcidSubstitutionMatrix;
import org.jcvi.jillion.align.BlosumMatrices;
import org.jcvi.jillion.align.pairwise.PairwiseAlignmentBuilder;
import org.jcvi.jillion.align.pairwise.ProteinPairwiseSequenceAlignment;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.IupacTranslationTables;
import org.jcvi.jillion.core.residue.aa.ProteinSequence;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.Nucleotide;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.*;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class CheckCoverageTest {

    @Autowired
    private CheckCoverage checkCoverage;


    /*@Test
    public void test(){
        ProteinSequence querySeq = new ProteinSequenceBuilder("MSLLTEVETYTLSIIPSGPLKAEIARRRRRQRLEDVFAGKNADLEALMEWIKTRPILSPLTKGILGFVFTLTVPSERGLQRRRFVQNALNGNGDPNNMDKAVKLYKKLKREMTFHGAKEVALSYSTGALASCMGLIYNRMGTVTAEGALGLVCATCEQIADAQHRSHRQMATTTNPLIRHENRMVLASTTAKAMEQMAGSSEQAAEAMEVASQAR").build();
        System.out.println("Length of the string is "+querySeq.getLength());
        ProteinSequence subSeq = new ProteinSequenceBuilder("MSLLTEVETYVSSIIPSGPLKAEIAQRLEDVFAGKNADLEALMEWIKTRPILSPLTKGILGFVFTLTVPSERGLQRRRFVQNALNGNGDPNNMDKAVQQQQKLYKKLKREMTFHGARKEVALSYSTGALASCMGLIYNRMGTVTAEGALGLVCATCEQIADAQHRSHRQMATTTNPLIRHENRMVLASTTAKAMEQMAGSSEQAAEAMEVASQAR").build();
        AminoAcidSubstitutionMatrix blosom90 = BlosumMatrices.blosum90();
        ProteinPairwiseSequenceAlignment actual = PairwiseAlignmentBuilder
                .createProtienAlignmentBuilder(querySeq,
                        subSeq, blosom90).gapPenalty(-8, -8)
                .build();
        System.out.println("Alignment score :"+actual.getScore());
        System.out.println("Number of mismatches "+actual.getNumberOfMismatches());
        System.out.println("Number of gap openings "+ actual.getNumberOfGapOpenings());
        int mismatches = actual.getNumberOfMismatches()+actual.getNumberOfGapOpenings();
        System.out.println("Alignment Length "+actual.getAlignmentLength());
        int matches = actual.getAlignmentLength()-mismatches;
        System.out.println("Number of matches "+ matches);
        //max of query or subject sequence length
        long maxSeqLength = Long.max(querySeq.getLength(),subSeq.getLength());
        double similarity = ((double)matches/maxSeqLength)*100;
        System.out.println("Similarity " +similarity);
        System.out.println("identity " + actual.getPercentIdentity());
        long coverage = Long.max(actual.getQueryRange().getLength(),actual.getSubjectRange().getLength());
        System.out.println(coverage);
        double percentCoverage = (coverage/querySeq.getLength())*100;
        System.out.println("%coverage "+percentCoverage);
        // coverage : max of sub , query alignment length;
    }*/

    @Test
    public void determineHomologyTest(){
        Model model = getTestModel();
        NucleotideSequence cds = checkCoverage.determineCDS(model);
        model = checkCoverage.determineHomology(model,cds);
        assertEquals("VT*KS*R",model.getTanslatedSeq().toString());
    }

    @Test
    public void getTranslatedProteinCoordinateTest(){
        List<Exon> exons = new ArrayList<Exon>();
        Exon exon1 = new Exon();
        exon1.setRange(Range.of(1,10));
        exon1.setFrame(Frame.ONE);
        Exon exon2 = new Exon();
        exon2.setRange(Range.of(11,23));
        exon2.setFrame(Frame.THREE);
        Exon exon3 = new Exon();
        exon3.setRange(Range.of(24,30));
        exon3.setFrame(Frame.TWO);
        exons.add(exon1);
        exons.add(exon2);
        exons.add(exon3);
        // NTOffset after considering the insertion string length;
        long PCoordinate = checkCoverage.getTranslatedProteinCooridnate(exons,11,Range.of(11,12));
        assertEquals(4,PCoordinate);
    }
    @Test
    public void getInternalStopsTest(){
        Model model = new Model();
        List<Exon> exons =  new ArrayList<Exon>();
        Exon exon1 = new Exon();
        exon1.setRange(Range.of(1,11));
        exon1.setFrame(Frame.ONE);
        Exon exon2 = new Exon();
        exon2.setRange(Range.of(12,22));
        exon2.setFrame(Frame.TWO);
        Exon exon3 = new Exon();
        exon3.setRange(Range.of(23,27));
        exon3.setFrame(Frame.THREE);
        exons.add(exon1);
        exons.add(exon2);
        exons.add(exon3);
        Alignment alignment = new Alignment();
        VirusGenome genome = new VirusGenome();
        genome.setSequence(new NucleotideSequenceBuilder("ATGAGTCTTCTAACCGAGGTCGTGACGTA").build());
        alignment.setVirusGenome(genome);
        model.setAlignment(alignment);
        model.setExons(exons);
        model.setReplaceStopCodonRange(Range.of(1,3));
        List<Range> internalStops = checkCoverage.getInternalStops(model);
        assertEquals(2,internalStops.size());
    }
    public Model getTestModel(){
        Model model = new Model();
        List<Exon> exons =  new ArrayList<Exon>();
        Exon exon1 = new Exon();
        exon1.setRange(Range.of(4,12));
        exon1.setFrame(Frame.ONE);
        Exon exon2 = new Exon();
        exon2.setRange(Range.of(19,22));
        exon2.setFrame(Frame.ONE);
        Exon exon3 = new Exon();
        exon3.setRange(Range.of(23,27));
        exon3.setFrame(Frame.THREE);
        exons.add(exon1);
        exons.add(exon2);
        exons.add(exon3);
        Alignment alignment = new Alignment();
        VirusGenome genome = new VirusGenome();
        genome.setSequence(new NucleotideSequenceBuilder("ATGAGTCTTCTAACCGAGGTCGTGACGTA").build());
        alignment.setVirusGenome(genome);
        model.setExons(exons);
        ViralProtein vp = new ViralProtein();
        vp.setSequence(new ProteinSequenceBuilder("VFTKSRR").build());
        GeneAttributes attributes = new GeneAttributes();
        StopTranslationException translationEx = new StopTranslationException();
        translationEx.setReplacementAA(AminoAcid.Threonine);
        attributes.setStopTranslationException(translationEx);
        RNA_Editing rna_editing = new RNA_Editing();
        rna_editing.setInsertionString("AAA");
        attributes.setRna_editing(rna_editing);
        vp.setGeneAttributes(attributes);
        alignment.setViralProtein(vp);
        model.setAlignment(alignment);
        model.setInsertRNAEditingRange(Range.of(13,15));
        model.setReplaceStopCodonRange(Range.of(7,9));
        return model;
    }

    @Test
    public void determineCDS(){
        Model model = getTestModel();
        NucleotideSequence cds = checkCoverage.determineCDS(model);
        //  ProteinSequence translatedSeq = IupacTranslationTables.STANDARD.translate(cds);
        //  System.out.println(translatedSeq);
        assertEquals("GTCTTCTAAAAATCGTGACGT",cds.toString());
    }

}
