package org.jcvi.vigor.utils;

import com.google.common.collect.Sets;
import org.jcvi.vigor.component.RNA_Editing;
import org.jcvi.vigor.component.Ribosomal_Slippage;

import java.util.*;
import java.util.stream.Collectors;

import static org.jcvi.vigor.utils.ConfigurationParameterFunctions.*;


public enum ConfigurationParameters {

    CandidateEvalue("candidate_evalue", "Evalue for identifying potential genes", Flags.VERSION_3),
    CandidateSelection("candidate_selection", "TODO"), // TODO
    CandidateMinimumSimilarity("min_candidate_pctsimilarity", "",  toPercent),
    CandidateMinimumSubjectCoverage("min_candidate_sbjcoverage", "",  toPercent),

    CircularGene("circular_gene", "Gene is circular", toBoolean,  Flags.UNIMPLEMENTED, Flags.VERSION_3, Flags.VERSION_4),
    CompleteGene("complete_gene", "Gene is complete", toBoolean, Flags.UNIMPLEMENTED, Flags.VERSION_3, Flags.VERSION_4),

    CondensationMinimum("min_condensation", "", toInteger, Flags.VERSION_4),
    RelaxCondensationMinimum("relax_min_condenstaion", "", toInteger, Flags.VERSION_4),

    ExonMaximumSize("max_exon_size", "Maximum sequence length of an exon", toInteger),
    ExonMinimumSize("min_exon_size", "Minimum sequence length of an exon", toInteger),

    ExoneratePath("exonerate_path", "Path to exonerate tool",
                  Flags.VERSION_4,
                  Flags.REQUIRED,
                  Flags.PLATFORM_DEPENDENT,
                  Flags.COMMANDLINE_SET,
                  Flags.PROGRAM_CONFIG_SET),

    FrameShiftSensitivity("frameshift_sensitivity", "How to handle frameshifts", Flags.VERSION_3, Flags.VERSION_4, Flags.UNIMPLEMENTED), // TODO

    GeneMinimumSize("min_gene_size", "Minimum sequence length to be considered as a gene", Flags.VERSION_3, Flags.VERSION_4),
    GeneMinimumCoverage("min_gene_coverage", "Minimum coverage of genes",  Flags.VERSION_3, Flags.VERSION_4), // TODO elaborate

    IntronMaximumSize("max_intron_size", "Maximum sequence length of an intron", toInteger, Flags.VERSION_3, Flags.VERSION_4),
    IntronMinimumSize("min_intron_size", "Minimum sequence length of an intron", toInteger, Flags.VERSION_3, Flags.VERSION_4),

    JCVIRules("jcvi_rules", "",  toBoolean),
    Verbose("verbose", "Make console and .rpt file output more detailed", toBoolean,
            Flags.VERSION_3, Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET),

    MaturePeptideMinimumCoverage("mature_pep_mincoverage", "",  ConfigurationParameterFunctions.toPercent,
                                 Flags.VERSION_3, Flags.VERSION_4),
    MaturePeptideMinimumSimilarity("mature_pep_minsimilarity", "",
                                   toPercent,
                                   Flags.VERSION_3, Flags.VERSION_4),
    MaturePeptideMinimumIdentity("mature_pep_minidentity", "",
                                 toPercent,
                                 Flags.VERSION_3, Flags.VERSION_4),

    MinimumMissingAASize("min_missing_AA_size", "TODO", toInteger), // TODO

    OutputDirectory("output_directory", "Write output to this directory",  Flags.VERSION_3, Flags.VERSION_4,
                    Flags.COMMANDLINE_SET, Flags.REQUIRED),
    OutputPrefix("output_prefix", "Use this prefix output files", Flags.VERSION_3, Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.REQUIRED),
    OverwriteOutputFiles("overwrite_output_files", "Overwrite output files if they exist",
                         toBoolean,
                         Flags.VERSION_4,
                         Flags.COMMANDLINE_SET,
                         Flags.PROGRAM_CONFIG_SET),

    PseudoGeneMinimumIdentity("min_pseudogene_identity", "", toPercent),
    PseudoGeneMinimumSimilarity("min_pseudogene_similarity", "", toPercent),
    PseudoGeneMinimumCoverage("min_pseudogene_coverage", "", toPercent),

    ReferenceDatabasePath("reference_database_path", "Directory containing reference database files",
                          Flags.VERSION_4, Flags.REQUIRED, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET),

    ScoreFactorAlignment("alignment_score_factor", "weight of alignment score in evaluating alignment", toDouble, Flags.VERSION_4),
    ScoreFactorStart("start_score_factor", "" ,  toDouble, Flags.VERSION_4),
    ScoreFactorSplicing("splicing_score_factor", "",  toDouble, Flags.VERSION_4),
    ScoreFactorStop("stop_score_factor", "",  toDouble, Flags.VERSION_4),
    ScoreFactorLeakyStop("leakystop_score_factor", "", toDouble, Flags.VERSION_4),
    ScoreFactorLeakyStopNotFound("leakystop_notFound_score", "", toDouble, Flags.VERSION_4),

    SequenceGapMinimumLength("min_seq_gap_length", "Minimum length of a sequence gap to ....",  toInteger, Flags.VERSION_4), // TODO
    StartCodons("StartCodons", "Comma separated list of expected start codons", toListOfStrings),
    StartCodonSearchWindow("start_codon_search_window", "Number of nucleotides before and after a candidate site to check for a start codon",
                           toInteger, Flags.VERSION_4),
    StopCodonSearchWindow("stop_codon_search_window", "Number of nucleotides before and after a candidate site to check for a stop codon",
                          toInteger, Flags.VERSION_4),

    TemporaryDirectory("temporary_directory", "Directory under which Vigor creates temporary files and directories",
                       Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET, Flags.REQUIRED),
    Locustag("locus_tag", "Locus tag prefix to use in output", Flags.VERSION_3, Flags.VERSION_4),

    AAOverlap_offset("AAOverlap_offset", "", toInteger, Flags.VERSION_4),
    NTOverlap_offset("NTOverlap_offset", "", toInteger, Flags.VERSION_4),

    VirusSpecificConfiguration("virusSpecific_config", "Path to virus specific configuration file."
            , Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET),
    VirusSpecificConfigurationPath("virusSpecific_config_path","Directory containing virus specific configurations.",
                                   Flags.VERSION_4,
                                   Flags.COMMANDLINE_SET,
                                   Flags.PROGRAM_CONFIG_SET),

    MaxGeneOverlap("max_gene_overlap"," In reporting gene models, maximum overlap of genes allowed.",toInteger, Flags.VERSION_4),

    AlignmentModule("alignment_module", "Alignment Module", isMemberOfSet("exonerate"),
                    Flags.COMMANDLINE_SET,
                    Flags.PROGRAM_CONFIG_SET,
                    Flags.VERSION_4,
                    Flags.REQUIRED),

    RibosomalSlippage("ribosomal_slippage", "V4_Ribosomal_Slippage","Ribosomal slippage. Format is offset/frameshift/regex",
                      ConfigurationParameterFunctions.of(Ribosomal_Slippage.class, Ribosomal_Slippage::parseFromString),
                      Flags.VERSION_4, Flags.GENE_SET),
    RNAEditing("rna_editing",  "V4_rna_editing", "RNA editing. Format is offset/regex/insertion string/note",
               ConfigurationParameterFunctions.of(RNA_Editing.class, RNA_Editing::parseFromString),
               Flags.VERSION_4, Flags.GENE_SET ),
    SpliceForm("splice_form", "Splice form. Format is ([ie]\\\\d+)+", toSpliceForms, Flags.GENE_SET, Flags.VIRUS_SET),

    StopCodonReadthrough("stop_codon_readthrough","V4_stop_codon_readthrough", "format is amino acid/offset/regex",
                         toStopException,
                         Flags.VERSION_4, Flags.GENE_SET),
    TinyExon3("tiny_exon3", "Tiny exon 3. Format is regex:[offset]", Flags.VERSION_4, Flags.GENE_SET),
    TinyExon5("tiny_exon5", "Tiny exon 5. Format is regex:[offset]", Flags.VERSION_4, Flags.GENE_SET),
    SharedCDS("shared_cds", "Shared CDS. Format is CDS[,CDS]",
              ConfigurationParameterFunctions.toListOfStrings,
              Flags.VERSION_4, Flags.GENE_SET),
    MinFunctionalLength("min_functional_length", "Minimum functional length", toInteger, Flags.VERSION_4, Flags.GENE_SET),
    MaturePeptideDB("matpepdb", "Mature peptide database to use for gene", Flags.VERSION_4, Flags.GENE_SET),
    NonCanonicalSplicing("nancanonical_splicing", "Alternate splice sites", Flags.VERSION_4, Flags.VIRUS_SET, Flags.GENE_SET),
    ExcludesGene("excludes_gene", "Excludes gene. TODO",
                 ConfigurationParameterFunctions.toListOfStrings,
                 Flags.VERSION_4, Flags.GENE_SET),
    AlternateStartCodons("alternate_startcodon", "Alternate start codons for gene. Format is CODON[,CODON,..]",
                         ConfigurationParameterFunctions.toListOfStrings,
                         Flags.VERSION_4, Flags.GENE_SET),

    // Version 3.5 parameters
    CandidateBlastOpts("candidate_blastopts", "Blast options when generating candidate models", Flags.VERSION_3_5),
    SlippageFrameShift("slippage_frameshift", "Frame shift for ribosomal slippage", Flags.VERSION_3_5),
    SlippageMotif("slippage_motif", "Motif for ribosomal slippage", Flags.VERSION_3_5),
    SlippageOffset("slippage_offset", "Offset from motif for ribosomal slippage", Flags.VERSION_3_5),
    Variation("variation", "Variation. TODO", Flags.VERSION_3_5),

    // Common DB defline attributes that can be ignored
    DBProduct("product", "product", Flags.GENE_SET, Flags.METADATA),
    DBLength("length", "length", Flags.VERSION_3_5, Flags.GENE_SET, Flags.METADATA, Flags.IGNORE),
    DBGene("gene", "gene name", Flags.GENE_SET, Flags.METADATA),
    DBDB("db", "gene database file?", Flags.GENE_SET, Flags.VERSION_3_5, Flags.METADATA, Flags.IGNORE),
    DBGeneVariation("gene_variation", "?", Flags.GENE_SET, Flags.VERSION_3_5, Flags.IGNORE),
    DBStopCodonReadThru("stopcodon_readthru", "", Flags.GENE_SET, Flags.VERSION_3_5, Flags.IGNORE),
    DBOrganism("organism", "", Flags.GENE_SET,  Flags.IGNORE),
    DBCluser("cluster", "", Flags.GENE_SET, Flags.IGNORE),
    DBGeneSynonym("gene_synonym", "gene_synonym", "", Flags.GENE_SET, Flags.METADATA)

    ;


    static final Map<String, ConfigurationParameters> byConfigKey;
    static {
        byConfigKey = Arrays.stream(ConfigurationParameters.values()).collect(Collectors.toMap(cp -> cp.configKey, cp -> cp));
    }

    public final String configKey;
    public final String deflineConfigKey;
    public final String description;

    public final ValueFunction valueFunction;

    private final EnumSet<Flags> flags = EnumSet.noneOf(Flags.class);

    public enum Flags {
        VERSION_3,
        VERSION_3_5,
        VERSION_4,
        UNIMPLEMENTED,
        PLATFORM_DEPENDENT,
        // not config, just informational
        METADATA,
        // don't complain if present, but don't set
        IGNORE,
        REQUIRED,
        // settable only via the commandline
        COMMANDLINE_SET,
        // settable only in the program level config, ie not in virus specific config
        PROGRAM_CONFIG_SET,
        // settable at the virus level
        VIRUS_SET,
        // settable at the gene level
        GENE_SET;

        public static final EnumSet<Flags> SET_FLAGS = EnumSet.of( Flags.PROGRAM_CONFIG_SET,
                                                                   Flags.COMMANDLINE_SET,
                                                                   Flags.VIRUS_SET,
                                                                   Flags.GENE_SET);

        public static final EnumSet<Flags> ALL_VERSIONS = EnumSet.of(Flags.VERSION_3,
                                                                     Flags.VERSION_3_5,
                                                                     Flags.VERSION_4);
    }

    ConfigurationParameters(String configKey, String description, Flags ... flags) {
        this(configKey, configKey, description,  flags);
    }

    ConfigurationParameters(String configKey, String description, ValueFunction valueFunction, Flags ... flags) {
        this(configKey, configKey, description, valueFunction, flags);
    }

    ConfigurationParameters(String configKey, String deflineConfigKey, String description, Flags ... flags) {
        this(configKey, deflineConfigKey, description, null, flags);
    }

    ConfigurationParameters(String configKey, String deflineConfigKey, String description, ValueFunction valueFunction, Flags... flags) {
        this.configKey = configKey;
        this.deflineConfigKey = deflineConfigKey;
        this.description = description;
        this.valueFunction = valueFunction == null ? ConfigurationParameterFunctions.of(String.class, s -> s): valueFunction;

        if (flags.length == 0) {
           flags = new Flags[] { Flags.VERSION_3, Flags.VERSION_4};
        }

        this.flags.addAll(Arrays.asList(flags));
        // If not specified then the parameter may be set at any level
        if (Sets.intersection(this.flags, Flags.SET_FLAGS).isEmpty()) {
            this.flags.addAll(Flags.SET_FLAGS);
        }
        // If not specified then we assume it's for all know versions
        if (Sets.intersection(this.flags, Flags.ALL_VERSIONS).isEmpty()) {
            this.flags.addAll(Flags.ALL_VERSIONS);
        }
    }

    public static ConfigurationParameters getParameterByConfigKey(String configKey) {
        return byConfigKey.get(configKey);
    }

    public boolean hasFlag(Flags checkFlag) {
        return flags.contains(checkFlag);
    }

    public EnumSet<Flags> hasFlags(EnumSet<Flags> checkFlags) {
        EnumSet<Flags> result = EnumSet.copyOf(flags);
        result.removeIf(e -> ! checkFlags.contains(e));
        return result;
    }

    public boolean hasOneOrMoreFlags(Flags ... checkFlags) {
        return Arrays.stream(checkFlags).anyMatch(flags::contains);
    }

    public boolean hasAllFlags(Flags ... checkFlags) {
        return Arrays.stream(checkFlags).allMatch(flags::contains);
    }

    public String getEnvVarName() {
        return "VIGOR_" + configKey.toUpperCase();
    }

    public String getSystemPropertyName() {
        return "vigor." + configKey;
    }

    public Object stringToValue(String stringValue) {
        try {
            return this.valueFunction.valueClass.cast(valueFunction.valueFunction.apply(stringValue));
        } catch (InvalidValue e) {
            throw new InvalidValue(String.format("bad value for %s: \"%s\":  %s",
                                                 this.configKey,
                                                 stringValue,
                                                 e.getMessage()));
        } catch (RuntimeException e ) {
            throw new ConfigurationParameterFunctions.InvalidValue(String.format("bad value for %s: \"%s\": got %s %s",
                                                 this.configKey,
                                                 stringValue,
                                                 e.getClass().getSimpleName(),
                                                 e.getMessage())
            );
        }
    }

}
