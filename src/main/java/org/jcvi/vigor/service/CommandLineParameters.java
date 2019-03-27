package org.jcvi.vigor.service;

public final class CommandLineParameters {

    public final static String outputPrefix = "output_prefix";
    public final static String inputFile = "input_fasta";
    public final static String referenceDB = "reference_database";
    public final static String minCoverage = "min_coverage";
    public final static String circularGenome = "circular_genome";
    public final static String frameshiftSensitivity = "frameshift_sensitivity";
    public final static String locusTag = "locus_tag";
    public final static String ignoreReferenceRequirements = "ignore_reference_requirements";
    public final static String parameters = "parameters";
    public final static String verbose = "verbose";
    public final static String ignoreRefID = "ignore_refID";
    public final static String configFile = "config_file";
    public final static String referenceDB_Path = "reference_database_path";
    public final static String overwriteOutputFiles = "overwrite_output_files";
    public final static String virusSpecificConfig = "virusSpecificConfig";
    public final static String virusSpecificConfigPath = "virusSpecificConfigPath";
    public final static String temporaryDirectory = "temporaryDirectory";
    public final static String listDatabases = "listDatabases";

    /**
     * Not to be instantiated
     */
    private CommandLineParameters () {

    }
}
