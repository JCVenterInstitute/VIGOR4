package org.jcvi.vigor.service;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

@RunWith(Suite.class)

@Suite.SuiteClasses({
        AlignmentGenerationServiceTest.class,
        ModelGenerationServiceTest.class,
        DetermineMissingExonsTest.class,
        DetermineStartAndStopTest.class,
        AdjustViralTricksTest.class,
        AdjustUneditedExonBoundariesTest.class,
        PeptideServiceTest.class
})
public class TestSuit {
}

