package org.jcvi.vigor.RegressionTest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;
import org.junit.runner.JUnitCore;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class Vigor4RegressionTestRunner {

    private static Logger LOGGER = LogManager.getLogger(Vigor4RegressionTestRunner.class);

    public static void main ( String[] args ) {

        try {
            Map<String, String> optsList = parseArguments(args);
            String outputDir = getOutputDirectory(optsList);
            JUnitCore jUnitCore = new JUnitCore();
            ValidateVigor4Models.prepare(outputDir);
            jUnitCore.run(ValidateVigor4Models.class);
        } catch (VigorException e) {
            LOGGER.error(e);
            System.exit(1);
        }
    }

    public static Map<String, String> parseArguments ( String[] args ) {

        Map<String, String> optsList = new HashMap<String, String>();
        for (int i = 0; i < args.length; i++) {
            switch (args[ i ].charAt(0)) {
                case '-':
                    if (args[ i ].length() < 2)
                        throw new IllegalArgumentException("Not a valid argument: " + args[ i ]);
                    else {
                        if (args.length - 1 == i)
                            if (args[ i ] == "-h") {
                                printHelp();
                            } else
                                throw new IllegalArgumentException("Expected arg after: " + args[ i ]);
                        // -opt
                        optsList.put(args[ i ], args[ i + 1 ]);
                        i++;
                    }
                    break;
            }
        }
        return optsList;
    }

    public static String getOutputDirectory ( Map<String, String> optsList ) throws VigorException {

        String outputDirectory = "";
        if (optsList.get("-o") != null) {
            outputDirectory = optsList.get("-o");
        } else throw new VigorException("Please provide output directory \"(-o <outputDirectory> )\"");
        File outputDir = new File(outputDirectory);
        if (!outputDir.exists()) {
            throw new VigorException(String.format("Output directory %s does not exist", outputDirectory));
        }
        return outputDir.getAbsolutePath();
    }

    private static void printHelp () {

        LOGGER.info("-o <outputDirectory>");
    }
}
