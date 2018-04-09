package org.jcvi.vigor.service;


import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.jcvi.jillion.core.Range;
import org.jcvi.vigor.Application;
import org.jcvi.vigor.component.Alignment;
import org.jcvi.vigor.component.Model;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.service.exception.ServiceException;
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

    private List<Alignment> alignments;
    private List<Model> models=new ArrayList<Model>();
    @Autowired
    private ModelGenerationService modelGenerationService;
    @Autowired
    private ViralProteinService viralProteinService;
    @Autowired
    private AdjustUneditedExonBoundaries adjustUneditedExonBoundaries;
    private ClassLoader classLoader = VigorTestUtils.class.getClassLoader();
    private File file = new File(classLoader.getResource("vigorUnitTestInput/Flua_SpliceSites_Test.fasta"). getFile());

    @Before
    public void getModel() throws ServiceException {
        alignments = VigorTestUtils.getAlignments(file.getAbsolutePath(),"flua_db",VigorUtils.getVigorWorkSpace(),"seg8prot2A");
        for (int i=0; i<alignments.size(); i++) {
            alignments.set(i,viralProteinService.setViralProteinAttributes(alignments.get(i), new VigorForm()));
        }
        alignments.stream().forEach(x -> {
            models.addAll(modelGenerationService.alignmentToModels(x, "exonerate"));
        });

    }

    @Test
    public void adjustSpliceSitesTest() throws CloneNotSupportedException{
        Model testModel = models.get(0);
        testModel.getExons().get(0).setRange(Range.of(11,30));
        List<Model> outModels = adjustUneditedExonBoundaries.adjustSpliceSites(testModel);
        Comparator<Model> bySpliceScore = Comparator.comparing( (m)->m.getScores().get("spliceScore"));
        Optional<Model> outModel = outModels.stream().sorted(bySpliceScore.reversed()).findFirst();
        assertEquals(Range.of(11,40),outModel.get().getExons().get(0).getRange());
    }


}
