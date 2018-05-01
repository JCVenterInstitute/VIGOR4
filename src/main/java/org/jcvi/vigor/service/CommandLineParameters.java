package org.jcvi.vigor.service;

public final class CommandLineParameters {
	public final static String outputPrefix = "output_prefix";
	public final static String inputFile = "input_fasta";
	public final static String referenceDB = "reference_database";
	public final static String genbankDB = "genbank_reference";
	public final static String eValue = "evalue";
	public final static String minCoverage = "min-converage";

	public final static String circularGene = "circular_gene";
	public final static String completeGene = "completeGene";
	public final static String frameshiftSensitivity = "frameshift_sensitivity";
	public final static String skipSelection = "skip_selection";
	public final static String useLocusTags = "use_locus_tags";
	public final static String jcviRules = "jcvi_rules";
	public final static String ignoreReferenceRequirements = "ignore_reference_requirements";
	public final static String parameters = "parameters";
	public final static String minGeneSize = "mine_gene_size";
	public final static String verbose = "verbose";
	public final static String ignoreRefID = "ignore_refID";
	public final static String listConfigParameters = "list_config_parameters";
    /**
     * Not to be instantiated
     */
	private CommandLineParameters() {

    }
}
