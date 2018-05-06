package org.jcvi.vigor.service;


import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.jcvi.vigor.utils.VigorUtils;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class AdjustUneditedExonBoundariesTest {

    final static Logger logger = LogManager.getLogger(AdjustUneditedExonBoundariesTest.class);

    @Autowired
    private ModelGenerationService modelGenerationService;
    @Autowired
    private ViralProteinService viralProteinService;
    @Autowired
    private AdjustUneditedExonBoundaries adjustUneditedExonBoundaries;




    @Test
    public void adjustSpliceSitesTest() throws CloneNotSupportedException, VigorException {
        ClassLoader classLoader = VigorTestUtils.class.getClassLoader();
        File virusGenomeSeqFile = new File(classLoader.getResource("vigorUnitTestInput/Flua_SpliceSites_Test.fasta"). getFile());
        File alignmentOutput = new File(classLoader.getResource("vigorUnitTestInput/Flua_SpliceSites_Test.txt"). getFile());
        String referenceDB = classLoader.getResource("vigorResources/data3/flua_db").getFile().toString();
        assertTrue("couldn't find reference DB", referenceDB != null);
        logger.info("using alignmentOutput file {} and reference database {}", alignmentOutput, referenceDB);
        List<Alignment> alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile,
                referenceDB,alignmentOutput);
        for (int i=0; i<alignments.size(); i++) {
            alignments.set(i,viralProteinService.setViralProteinAttributes(alignments.get(i), new VigorForm()));
        }
        logger.info("found {} alignments for protein {}", alignments.size());
        List<Model> models = alignments.stream()
                                       .flatMap(x -> modelGenerationService.alignmentToModels(x, "exonerate").stream())
                                       .collect(Collectors.toList());
        logger.info("{} models for {} alignments", models.size(), alignments.size());
        assertTrue("no models found for alignments", models.size() > 0);
        Model testModel = models.get(0);
        assertTrue(String.format("no exons found for test model %s", testModel) , testModel.getExons().size() > 0);
        testModel.getExons().get(0).setRange(Range.of(11,30));
        List<Model> outModels = adjustUneditedExonBoundaries.adjustSpliceSites(testModel);
        Comparator<Model> bySpliceScore = Comparator.comparing( (m)->m.getScores().get("spliceScore"));
        Optional<Model> outModel = outModels.stream().sorted(bySpliceScore.reversed()).findFirst();
        assertTrue("No adjusted model found", outModel.isPresent());
        assertEquals(Range.of(11,40),outModel.get().getExons().get(0).getRange());
    }


}
