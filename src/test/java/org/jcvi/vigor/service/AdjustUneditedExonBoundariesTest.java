package org.jcvi.vigor.service;


import static junit.framework.TestCase.assertTrue;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.notNullValue;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.file.Paths;
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
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
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
    @Autowired
    private VigorInitializationService initializationService;

    @Test
    public void adjustSpliceSitesTest() throws CloneNotSupportedException, VigorException {
        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        File resources = new File("src/test/resources");
        File virusGenomeSeqFile = new File(resources.getAbsolutePath()+File.separator+"vigorUnitTestInput/Flua_SpliceSites_Test.fasta");
         File alignmentOutput = new File(resources.getAbsolutePath()+File.separator+"vigorUnitTestInput/Flua_SpliceSites_Test.txt");
        String referenceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(referenceDBPath);
        assertThat("reference database path must be set", referenceDBPath, is(notNullValue()));
        String referenceDB = Paths.get(referenceDBPath, "flua_db").toString();
        assertTrue("couldn't find reference DB", referenceDB != null);
        logger.info("using alignmentOutput file {} and reference database {}", alignmentOutput, referenceDB);
        List<Alignment> alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile,
                referenceDB,alignmentOutput, config);
        for (int i=0; i<alignments.size(); i++) {
            alignments.set(i,viralProteinService.setViralProteinAttributes(alignments.get(i), new VigorForm(config)));
        }
        logger.info("found {} alignments", alignments.size());
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
