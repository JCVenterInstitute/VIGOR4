package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.*;
import java.util.*;
import java.util.regex.Pattern;

public class VigorUtils {

    private static Logger LOGGER = LogManager.getLogger(VigorUtils.class);
    private static Pattern hypPattern = Pattern.compile("^HYP\\b");

    public static boolean is_Integer ( String value ) {

        try {
            Integer.parseInt(value);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    public static String getDefaultConfigurationPath () {

        return Paths.get("vigorResources", "config", "defaults.ini").toString();
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

        if (! deleteDirectory(new File(path).toPath())) {
            LOGGER.warn("Error deleting temporary working directory {}", path);
        }
    }

    public static boolean deleteDirectory (Path path) {
        try {
            // clean up
            if (path != null) {
                Files.walk(path)
                     .map(Path::toFile)
                     .sorted(Comparator.reverseOrder())
                     .forEach(File::delete);
            }
            return true;
        } catch (NoSuchFileException e) {
            return true;
        } catch (IOException e) {
            return false;
        }

    }

    // TODO this is a security issue but it requires the user's assistance in that they pass the file path
    public static String expandTilde(String path) throws IOException {
        return new BufferedReader(new InputStreamReader(Runtime.getRuntime().exec(new String[] {"/bin/bash", "-c", "echo " + path}).getInputStream())).readLine();
    }

    public enum FileCheck {
        READ, WRITE, EXISTS, NOT_EXISTS, EXECUTE, ABSOLUTE, RELATIVE, DIRECTORY, FILE
    }

    public static String checkFilePath(String description, String path, FileCheck ... modes) throws VigorException {
        if (path == null) {
            throw new VigorException(String.format("%s not set", description));
        }
        File testFile = new File(path);

        EnumSet<FileCheck> checks = EnumSet.noneOf(FileCheck.class);
        if (modes.length > 0) {
            checks = EnumSet.copyOf(Arrays.asList(modes));
        }

        if (checks.contains(FileCheck.EXISTS) && ! testFile.exists()) {
            throw new VigorException(String.format("%s %s does not exist", description, path));
        }
        if (checks.contains(FileCheck.NOT_EXISTS) && testFile.exists()) {
            throw new VigorException(String.format("%s %s exists", description, path));
        }

        List<String> errors = new ArrayList<>();

        if (checks.contains(FileCheck.READ) && ! testFile.canRead()) {
            errors.add("not readable");
        }

        if (checks.contains(FileCheck.WRITE) && ! testFile.canWrite()) {
            errors.add("not writable");
        }

        if (checks.contains(FileCheck.EXECUTE) && ! testFile.canExecute()) {
            errors.add("not executable");
        }

        if (checks.contains(FileCheck.ABSOLUTE) && ! testFile.isAbsolute()) {
            errors.add("not an absolute path");
        }

        if (checks.contains(FileCheck.RELATIVE) && testFile.isAbsolute()) {
            errors.add("not a relative path");
        }

        if (checks.contains(FileCheck.DIRECTORY) && ! testFile.isDirectory()) {
            errors.add("not a directory");
        }
        if (checks.contains(FileCheck.FILE) && ! testFile.isFile()) {
            errors.add("not a file");
        }

        if (! errors.isEmpty()) {
            throw new VigorException(String.format("%s %s is %s", description, path, String.join(",", errors)));
        }

        return path;
    }
}
