package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.Vigor;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.regex.Pattern;

public class VigorUtils {

    private static Logger LOGGER = LogManager.getLogger(Vigor.class);
    private static Pattern hypPattern = Pattern.compile("^HYP\\b");

    public static String getVigorParametersPath () {

        return Paths.get("vigorResources", "config", "vigor.ini").toString();
    }

    public static boolean is_Integer ( String value ) {

        try {
            Integer.parseInt(value);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    public static String removeQuotes ( String in ) {

        return in.replaceAll("^\"|\"$", "");
    }

    public static String nameToLocus ( String name, String prefix, boolean isPseudoGene ) {

        if (prefix == null || prefix.isEmpty()) {
            return null;
        }
        String symbol = name;
        if (isPseudoGene || hypPattern.matcher(symbol).find()) {
            symbol = symbol.replace("-", "p");
        } else {
            symbol = symbol.replace("-[0-9].*$", "");
        }
        symbol = symbol.replace("[^a-zA-Z0-9-]", "");
        return prefix + symbol;
    }

    public static String putativeName ( String name, boolean truncated3p, boolean truncated5p ) {

        String putativeName = name;
        if (!name.contains("putative")) {
            putativeName = "putative " + name;
        }
        if (truncated3p && truncated5p) {
            return putativeName + ", fragment";
        } else if (truncated3p) {
            return putativeName + ", N-terminal";
        } else if (truncated5p) {
            return putativeName + ", C-terminal";
        } else {
            return name;
        }
    }

    public static void deleteTempFiles ( String path ) {

        Path workspace = new File(path).toPath();
        try {
            // clean up
            if (workspace != null) {
                Files.walk(workspace)
                        .map(Path:: toFile)
                        .sorted(Comparator.reverseOrder())
                        .forEach(File:: delete);
            }
        } catch (IOException e) {
            LOGGER.warn("Error deleting temporary working directory {}", workspace);
        }
    }
}
