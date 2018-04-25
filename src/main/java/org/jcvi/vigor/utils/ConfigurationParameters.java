package org.jcvi.vigor.utils;

import java.util.*;
import java.util.stream.Collectors;

public enum ConfigurationParameters {

    CandidateEvalue("candidate_evalue", "Evalue for identifying potential genes", Flags.VERSION_3),
    CandidateSelection("candidate_selection", "TODO"), // TODO
    CandidateMinimumSimilarity("min_candidate_pctsimilarity", ""),
    CandidateMinimumSubjectCoverage("min_candidate_sbjcoverage", ""),

    CircularGene("circular_gene", "Gene is circular", Flags.UNIMPLEMENTED, Flags.VERSION_3, Flags.VERSION_4),
    CompleteGene("complete_gene", "Gene is complete", Flags.UNIMPLEMENTED, Flags.VERSION_3, Flags.VERSION_4),

    CondensationMinimum("min_condensation", ""),

    ExonMaximumSize("max_exon_size", "Maximum sequence length of an exon"),
    ExonMinimumSize("min_exon_size", "Minimum sequence length of an exon"),

    ExoneratePath("exonerate_path", "Path to exonerate tool", Flags.VERSION_4, Flags.REQUIRED, Flags.PLATFORM_DEPENDENT),

    FrameShiftSensitivity("frameshift_sensitivity", "How to handle frameshifts", Flags.VERSION_3, Flags.VERSION_4, Flags.UNIMPLEMENTED), // TODO

    GeneMinimumSize("min_gene_size", "Minimum sequence length to be considered as a gene", Flags.VERSION_3, Flags.VERSION_4),
    GeneMinimumCoverage("min_gene_coverage", "Minimum coverage of genes",  Flags.VERSION_3, Flags.VERSION_4), // TODO elaborate

    IntronMaximumSize("max_intron_size", "Maximum sequence length of an intron", Flags.VERSION_3, Flags.VERSION_4),
    IntronMinimumSize("min_intron_size", "Minimum sequence length of an intron", Flags.VERSION_3, Flags.VERSION_4),

    JCVIRules("jcvi_rules", ""),

    MaturePeptideMinimumCoverage("mature_pep_mincoverage", "",  Flags.VERSION_3, Flags.VERSION_4),
    MaturePeptideMinimumSimilarity("mature_pep_minsimilarity", "",  Flags.VERSION_3, Flags.VERSION_4),
    MaturePeptideMinimumIdentity("mature_pep_minidentity", "",  Flags.VERSION_3, Flags.VERSION_4),

    MinimumMissingAASize("min_missing_AA_size", "TODO"), // TODO

    OutputDirectory("output_directory", "Write output to this directory",  Flags.VERSION_3, Flags.VERSION_4, Flags.COMMANDLINE),
    OutputPrefix("output_prefix", "Use this prefix output files", Flags.VERSION_3, Flags.VERSION_4, Flags.COMMANDLINE),

    PseudoGeneMinimumIdentity("min_pseudogene_identity", ""),
    PseudoGeneMinimumSimilarity("min_pseudogene_similarity", ""),
    PseudoGeneMinimumCoverage("min_pseudogene_coverage", ""),

    ReferenceDatabasePath("reference_database_path", "Directory containing reference database files", Flags.VERSION_4, Flags.REQUIRED),

    ScoreFactorExonerate("exonerate_score_factor", "weight of exonerate score in evaluating alignment", Flags.VERSION_4),
    ScoreFactorStart("start_score_factor", ""),
    ScoreFactorSplicing("splicing_score_factor", ""),
    ScoreFactorStop("stop_score_factor", ""),
    ScoreFactorLeakyStop("leakystop_score_factor", ""),
    ScoreFactorLeakyStopNotFound("leakystop_notFound_score", ""),

    SequenceGapMinimumLength("min_seq_gap_length", "Minimum length of a sequence gap to ...."), // TODO
    StartCodons("StartCodons", "Comma separated list of expected start codons"),
    StartCodonSearchWindow("start_codon_search_window", "Number of nucleotides before and after a candidate site to check for a start codon"),
    StopCodonSearchWindow("stop_codon_search_window", "Number of nucleotides before and after a candidate site to check for a stop codon"),

    UseLocustags("use_locus_tags", "Include locus tags in output", Flags.VERSION_3, Flags.VERSION_4),

    AAOverlap_offset("AAOverlap_offset", ""),
    NTOverlap_offset("NTOverlap_offset", ""),

    VirusSpecificConfigurationPath("virusSpecific_parameters","Directory containing virus specific configurations. TODO change configuration to path",Flags.VERSION_4);

    static final Map<String, ConfigurationParameters> byConfigKey;
    static {
        byConfigKey = Arrays.stream(ConfigurationParameters.values()).collect(Collectors.toMap(cp -> cp.configKey, cp -> cp));
    }

    public final String configKey;
    public final String description;

    private final EnumSet<Flags> flags = EnumSet.noneOf(Flags.class);
    public enum Flags { VERSION_3, VERSION_3_5, VERSION_4, UNIMPLEMENTED, PLATFORM_DEPENDENT, REQUIRED, COMMANDLINE }

    ConfigurationParameters(String configKey, String description, Flags... flags) {
        this.configKey = configKey;
        this.description = description;

        if (flags.length == 0) {
           flags = new Flags[] { Flags.VERSION_3, Flags.VERSION_4};
        }

        this.flags.addAll(Arrays.asList(flags));
    }

    public static ConfigurationParameters getParameterByConfigKey(String configKey) {
        return byConfigKey.get(configKey);
    }

    public boolean hasFlag(Flags configFlag) {
        return flags.contains(configFlag);
    }

}
