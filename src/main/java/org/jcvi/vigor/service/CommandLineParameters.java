package org.jcvi.vigor.service;

public final class CommandLineParameters {
	public final static String outputPrefix = "output_prefix";
	public final static String inputFile = "input_fasta";
	public final static String referenceDB = "reference_database";
	public final static String genbankDB = "genbank_reference";
	public final static String eValue = "evalue";
	public final static String minCoverage = "min_coverage";

	public final static String circularGene = "circular_gene";
	public final static String completeGene = "complete_gene";
	public final static String frameshiftSensitivity = "frameshift_sensitivity";
	public final static String skipSelection = "skip_selection";
	public final static String locusTag = "locus_tag";
	public final static String jcviRules = "jcvi_rules";
	public final static String ignoreReferenceRequirements = "ignore_reference_requirements";
	public final static String parameters = "parameters";
	public final static String minGeneSize = "min_gene_size";
	public final static String verbose = "verbose";
	public final static String ignoreRefID = "ignore_refID";
	public final static String listConfigParameters = "list_config_parameters";
	public final static String configFile = "config_file";
	public final static String referenceDB_Path = "reference_database_path";
    public final static String overwriteOutputFiles = "overwrite_output_files";
    public final static String virusSpecificConfig = "virusSpecificConfig";
    public final static String virusSpecificConfigPath = "virusSpecificConfigPath";

    /**
     * Not to be instantiated
     */
	private CommandLineParameters() {

    }
}