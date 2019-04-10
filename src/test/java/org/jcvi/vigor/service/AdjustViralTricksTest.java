package org.jcvi.vigor.service;

import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.notNullValue;
import static org.junit.Assert.*;
import static org.junit.Assume.assumeTrue;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.jcvi.vigor.testing.category.Fast;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.testing.category.ReferenceDatabase;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.junit.Test;
import org.junit.experimental.categories.Category;
import org.junit.runner.RunWith;
import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.utils.VigorTestUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringRunner;

@Category({Fast.class, ReferenceDatabase.class})
@RunWith(SpringRunner.class)
@ContextConfiguration(classes = Application.class)
public class AdjustViralTricksTest {

    @Autowired
    private ModelGenerationService modelGenerationService;
    @Autowired
    private ViralProteinService viralProteinService;
    @Autowired
    private AdjustViralTricks adjustViralTricks;
    @Autowired
    private VigorInitializationService initializationService;

    @Test
    public void adjustRibosomalSlippageTest () throws VigorException, CloneNotSupportedException {

        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String referenceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(referenceDBPath);
        assertThat("reference database path must be set", referenceDBPath, is(notNullValue()));
        String referenceDB = Paths.get(referenceDBPath, "flua_db").toString();
        List<Alignment> alignments;
        List<Model> models = new ArrayList<Model>();
        File virusGenomeSeqFile = new File(getClass().getResource("/vigorUnitTestInput/Flua_RiboSlippage_Test.fasta").getFile());
        File alignmentOutput = new File(getClass().getResource("/vigorUnitTestInput/Flua_RiboSlippage_Test.txt").getFile());
        alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile, referenceDB, alignmentOutput, config);
        for (int i = 0; i < alignments.size(); i++) {
            alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), config));
        }
        alignments = modelGenerationService.mergeIdenticalProteinAlignments(alignments);
        alignments.stream().forEach(x -> {
            models.addAll(modelGenerationService.alignmentToModels(x, config));
        });
        Model testModel = models.get(0);
        List<Model> outputModels = adjustViralTricks.adjustRibosomalSlippage(testModel);
        Range actual = outputModels.get(0).getExons().get(0).getRange();
        assertEquals(Range.of(9, 578), actual);
    }

    @Test
    public void checkForLeakyStopTest () throws VigorException, CloneNotSupportedException {

        VigorConfiguration defaultConfig = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String referenceDBPath = defaultConfig.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(referenceDBPath);
        assertThat("reference database path must be set", referenceDBPath, is(notNullValue()));
        List<VigorConfiguration> configurations = new ArrayList<>();
        configurations.add(defaultConfig);
        configurations.addAll(initializationService.loadVirusConfiguration(new File(referenceDBPath, "veev_db.ini") ));
        VigorConfiguration config = initializationService.mergeConfigurations(configurations);
        String referenceDB = Paths.get(referenceDBPath, "veev_db").toString();
        List<Alignment> alignments;
        List<Model> models = new ArrayList<Model>();
        File virusGenomeSeqFile = new File(getClass().getResource("/vigorUnitTestInput/Veev_StopTranslationEx_Test.fasta").getFile());
        File alignmentOutout = new File(getClass().getResource("/vigorUnitTestInput/Veev_StopTranslationEx_Test.txt").getFile());
        alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile, referenceDB, alignmentOutout, config);
        for (int i = 0; i < alignments.size(); i++) {
            alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), config));
        }
        alignments = modelGenerationService.mergeIdenticalProteinAlignments(alignments);
        alignments.stream().forEach(x -> {
            models.addAll(modelGenerationService.alignmentToModels(x, config));
        });
        Model testModel = models.get(0);
        double leakStopNotFoundScore = 1;
        List<Model> checkedModels = adjustViralTricks.checkForLeakyStop(testModel, leakStopNotFoundScore);
        Model outputModel = checkedModels.get(0);
        Range actual = outputModel.getReplaceStopCodonRange();
        assertEquals("TGA", outputModel.getAlignment().getVirusGenome().getSequence().toBuilder().trim(Range.of(5627, 5629)).build().toString());
        assertEquals(Range.of(5627, 5629), actual);
    }

    @Test
    public void adjustRNAEditingTest () throws VigorException, CloneNotSupportedException {

        VigorConfiguration config = initializationService.mergeConfigurations(initializationService.getDefaultConfigurations());
        String referenceDBPath = config.get(ConfigurationParameters.ReferenceDatabasePath);
        VigorTestUtils.assumeReferenceDB(referenceDBPath);
        assertThat("reference database path must be set", referenceDBPath, is(notNullValue()));
        Path referenceDB = Paths.get(referenceDBPath, "mmp_db");
        assumeTrue("This test requires the mmp_db", referenceDB.toFile().exists());
        List<Alignment> alignments;
        List<Model> models = new ArrayList<Model>();
        File virusGenomeSeqFile = new File(getClass().getResource( "/vigorUnitTestInput/mmp_rna_editing_Test.fasta").getFile());
        File alignmentOutput = new File(getClass().getResource("/vigorUnitTestInput/mmp_rna_editing_Test.txt").getFile());
        alignments = VigorTestUtils.getAlignments(virusGenomeSeqFile, referenceDB.toString(), alignmentOutput, config);
        for (int i = 0; i < alignments.size(); i++) {
            alignments.set(i, viralProteinService.setViralProteinAttributes(alignments.get(i), config));
        }
        alignments = modelGenerationService.mergeIdenticalProteinAlignments(alignments);
        alignments.stream().forEach(x -> {
            models.addAll(modelGenerationService.alignmentToModels(x, config));
        });
        Model testModel = models.get(0);
        List<Model> outputModels = adjustViralTricks.adjustRNAEditing(testModel);
        Range actual = outputModels.get(0).getExons().get(0).getRange();
        assertEquals(Range.of(1861, 2322), actual);
    }
}
