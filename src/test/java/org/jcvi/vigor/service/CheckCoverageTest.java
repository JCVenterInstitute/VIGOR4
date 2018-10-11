package org.jcvi.vigor.service;

import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.Frame;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.jillion.core.residue.aa.ProteinSequenceBuilder;
import org.jcvi.jillion.core.residue.nt.NucleotideSequence;
import org.jcvi.jillion.core.residue.nt.NucleotideSequenceBuilder;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.testing.category.Fast;
import org.jcvi.vigor.testing.category.Isolated;
import org.jcvi.vigor.component.*;
import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

@Category({Fast.class, Isolated.class})
@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class CheckCoverageTest {

    @Autowired
    private CheckCoverage checkCoverage;

    @Test
    public void determineHomologyTest () {

        Model model = getTestModel();
        NucleotideSequence cds = checkCoverage.determineCDS(model);
        model = checkCoverage.determineHomology(model, cds);
        assertEquals("VT*KS*", model.getTranslatedSeq().toString());
    }

    @Test
    public void getTranslatedProteinCoordinateTest () {

        List<Exon> exons = new ArrayList<Exon>();
        Exon exon1 = new Exon();
        exon1.setRange(Range.of(1, 10));
        exon1.setFrame(Frame.ONE);
        Exon exon2 = new Exon();
        exon2.setRange(Range.of(11, 23));
        exon2.setFrame(Frame.THREE);
        Exon exon3 = new Exon();
        exon3.setRange(Range.of(24, 30));
        exon3.setFrame(Frame.TWO);
        exons.add(exon1);
        exons.add(exon2);
        exons.add(exon3);
        // NTOffset after considering the insertion string length;
        long PCoordinate = checkCoverage.getTranslatedProteinCoordinate(exons, 11, Range.of(11, 12));
        assertEquals(4, PCoordinate);
    }

    @Test
    public void getInternalStopsTest () {

        Model model = new Model();
        List<Exon> exons = new ArrayList<Exon>();
        Exon exon1 = new Exon();
        exon1.setRange(Range.of(1, 11));
        exon1.setFrame(Frame.ONE);
        Exon exon2 = new Exon();
        exon2.setRange(Range.of(12, 22));
        exon2.setFrame(Frame.TWO);
        Exon exon3 = new Exon();
        exon3.setRange(Range.of(23, 27));
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
        model.setReplaceStopCodonRange(Range.of(1, 3));
        List<Range> internalStops = checkCoverage.getInternalStops(model);
        assertEquals(2, internalStops.size());
    }

    public Model getTestModel () {

        Model model = new Model();
        List<Exon> exons = new ArrayList<Exon>();
        Exon exon1 = new Exon();
        exon1.setRange(Range.of(4, 12));
        exon1.setFrame(Frame.ONE);
        Exon exon2 = new Exon();
        exon2.setRange(Range.of(19, 22));
        exon2.setFrame(Frame.ONE);
        Exon exon3 = new Exon();
        exon3.setRange(Range.of(23, 27));
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
        StopTranslationException translationEx = new StopTranslationException(true, AminoAcid.Threonine, "", 0);
        attributes.setStopTranslationException(translationEx);
        RNA_Editing rna_editing = new RNA_Editing(true, 0, "", "AAA", "");
        attributes.setRna_editing(rna_editing);
        vp.setGeneAttributes(attributes);
        alignment.setViralProtein(vp);
        model.setAlignment(alignment);
        model.setInsertRNAEditingRange(Range.of(13, 15));
        model.setReplaceStopCodonRange(Range.of(7, 9));
        return model;
    }

    @Test
    public void determineCDS () {

        Model model = getTestModel();
        NucleotideSequence cds = checkCoverage.determineCDS(model);
        assertEquals("GTCTTCTAAAAATCGTGA", cds.toString());
    }
}
