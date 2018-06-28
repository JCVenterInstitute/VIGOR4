package org.jcvi.vigor.service;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

@RunWith(Suite.class)
@Suite.SuiteClasses( {
        AlignmentGenerationServiceTest.class,
        ModelGenerationServiceTest.class,
        CheckCoverageTest.class,
        DetermineMissingExonsTest.class,
        PeptideServiceTest.class,
        DetermineStartAndStopTest.class,
        AdjustViralTricksTest.class,
        AdjustUneditedExonBoundariesTest.class
})
public class TestSuite {

}

