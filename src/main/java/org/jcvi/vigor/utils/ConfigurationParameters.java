package org.jcvi.vigor.utils;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

public enum ConfigurationParameters {

    CandidateEvalue("candidate_evalue", ""),
    CandidateSelection("candidate_selection", "TODO"), // TODO
    CandidateMinimumSimilarity("min_candidate_pctsimilarity", ""),
    CandidateMinimumSubjectCoverage("min_candidate_sbjcoverage", ""),

    CircularGene("circular_gene", "Gene is circular"),
    CompleteGene("complete_gene", "Gene is complete"),

    CondensationMinimum("min_condensation", ""),

    ExonMaximumSize("max_exon_size", "Maximum sequence length of an exon"),
    ExonMinimumSize("min_exon_size", "Minimum sequence length of an exon"),

    ExoneratePath("exonerate_path", "Path to exonerate tool"),

    FrameShiftSensitivity("frameshift_sensitivity", "How to handle frameshifts"), // TODO

    GeneMinimumSize("min_gene_size", "Minimum sequence length to be considered as a gene"),
    GeneMinimumCoverage("min_gene_coverage", "Minimum coverage of genes"), // TODO elaborate

    IntronMaximumSize("max_intron_size", "Maximum sequence length of an intron"),
    IntronMinimumSize("min_intron_size", "Minimum sequence length of an intron"),

    JCVIRules("jcvi_rules", ""),

    MaturePeptideMinimumCoverage("mature_pep_mincoverage", ""),
    MaturePeptideMinimumSimilarity("mature_pep_minsimilarity", ""),
    MaturePeptideMinimumIdentity("mature_pep_minidentity", ""),

    MinimumMissingAASize("min_missing_AA_size", "TODO"), // TODO

    OutputDirectory("output_directory", "Write output to this directory"),
    OutputPrefix("output_prefix", "Use this prefix output files"),

    PseudoGeneMinimumIdentity("min_pseudogene_identity", ""),
    PseudoGeneMinimumSimilarity("min_pseudogene_similarity", ""),
    PseudoGeneMinimumCoverage("min_pseudogene_coverage", ""),

    ReferenceDatabasePath("reference_database_path", "Directory containing reference database files"),

    ScoreFactorExonerate("exonerate_score_factor", "weight of exonerate score in evaluating alignment"),
    ScoreFactorStart("start_score_factor", ""),
    ScoreFactorSplicing("splicing_score_factor", ""),
    ScoreFactorStop("stop_score_factor", ""),
    ScoreFactorLeakyStop("leakystop_score_factor", ""),
    ScoreFactorLeakyStopNotFound("leakystop_notFound_score", ""),

    SequenceGapMinimumLength("min_seq_gap_length", "Minimum length of a sequence gap to ...."), // TODO
    StartCodons("StartCodons", "Comma separated list of expected start codons"),
    StartCodonSearchWindow("start_codon_search_window", "Number of nucleotides before and after a candidate site to check for a start codon"),
    StopCodonSearchWindow("stop_codon_search_window", "Number of nucleotides before and after a candidate site to check for a stop codon"),

    UseLocustags("use_locus_tags", "Include locus tags in output"),

    AAOverlap_offset("AAOverlap_offset", ""),
    NTOverlap_offset("NTOverlap_offset", ""),

    VirusSpecificConfigurationPath("virusSpecific_parameters","Directory containing virus specific configurations. TODO change configuration to path");

    static final Map<String, ConfigurationParameters> byConfigKey;
    static {
        byConfigKey = Arrays.stream(ConfigurationParameters.values()).collect(Collectors.toMap(cp -> cp.configKey, cp -> cp));
    }

    public final String description;
    public final String configKey;
    ConfigurationParameters(String configKey, String description) {
        this.configKey = configKey;
        this.description = description;
    }

    public static ConfigurationParameters getParameterByConfigKey(String configKey) {
        return byConfigKey.get(configKey);
    }

}
