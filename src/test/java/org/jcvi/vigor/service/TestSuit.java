package org.jcvi.vigor.service;

/**
 * Created by snettem on 5/25/2017.
 */


import org.junit.runner.RunWith;
import org.junit.runners.Suite;

@RunWith(Suite.class)

@Suite.SuiteClasses({
        AlignmentGenerationServiceTest.class,
        ModelGenerationServiceTest.class,
        DetermineMissingExonsTest.class,
        DetermineStartAndStopTest.class,
        AdjustViralTricksTest.class,
        AdjustUneditedExonBoundariesTest.class        
})
public class TestSuit {
}

