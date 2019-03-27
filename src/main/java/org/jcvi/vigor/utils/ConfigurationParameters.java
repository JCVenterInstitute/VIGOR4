package org.jcvi.vigor.utils;

import com.google.common.collect.Sets;
import org.jcvi.vigor.component.RNA_Editing;
import org.jcvi.vigor.component.Ribosomal_Slippage;

import java.util.*;
import java.util.stream.Collectors;

import static org.jcvi.vigor.utils.ConfigurationParameterFunctions.*;


public enum ConfigurationParameters {
    Alias("alias", "One or more aliases used for the database", toListOfStrings, Flags.VERSION_4, Flags.METADATA_SET),
    AAOverlapMaximum("max_aa_overlap", "Maximum number of proteins that may overlap for alignment fragments to be considered compatible when generating a gene model", toInteger, Flags.VERSION_4),
    AlignmentModule("alignment_module", "Alignment Module", isMemberOfSet("exonerate"),
                    Flags.COMMANDLINE_SET,
                    Flags.PROGRAM_CONFIG_SET,
                    Flags.VERSION_4,
                    Flags.REQUIRED),
    AlternateStartCodons("alternate_startcodon", "Alternate start codons for gene. Format is CODON[,CODON,..]",
                         ConfigurationParameterFunctions.toListOfStrings,
                         Flags.VERSION_4, Flags.GENE_SET),
    CandidateBlastOpts("candidate_blastopts", "Blast options when generating candidate models", Flags.VERSION_3_5),

    CandidateEvalue("candidate_evalue", "Evalue for identifying potential genes", Flags.VERSION_3),
    CandidateMinimumSimilarity("min_candidate_pctsimilarity", "Minimum percent similarity between the candidate model and the reference protein (for definition of similarity see mature_pep_minsimilarity)", toPercent, Flags.VERSION_3),

    CandidateMinimumSubjectCoverage("min_candidate_sbjcoverage", "Minimum percentage of coverage in the alignment between candidate model and reference protein, based on the longest of the two.", toPercent, Flags.VERSION_3),
    CandidateSelection("candidate_selection", "", Flags.VERSION_3),

    CircularGene("circular_genome", "When this parameter is set to TRUE, VIGOR consider the genome as circular, enabling annotating genes spanning both ends of the sequence (which would be continuous when circularized).", toBoolean, Flags.UNIMPLEMENTED, Flags.VERSION_3, Flags.VERSION_4),

    Description("description", "Description of virus database", Flags.METADATA_SET),

    DBCluster("cluster", "", Flags.GENE_SET, Flags.IGNORE),

    DBDB("db", "gene database file backreference used for autoselecting reference database", Flags.GENE_SET, Flags.VERSION_3_5, Flags.METADATA),

    DBGene("gene", "Mandatory identifier for each sequence in the reference database. The value corresponds to the value found in the [gene:Gene_ID] section of the config file (e.g. gene=\"NS1\")", Flags.GENE_SET, Flags.METADATA),
    DBGeneSynonym("gene_synonym",  "", Flags.GENE_SET, Flags.METADATA),

    DBGeneVariation("gene_variation", "?", Flags.GENE_SET, Flags.VERSION_3_5),
    DBLength("length", "length", Flags.VERSION_3_5, Flags.GENE_SET, Flags.METADATA, Flags.IGNORE),

    DBOrganism("organism", "organism", Flags.GENE_SET, Flags.IGNORE),
    DBProduct("product", "product", Flags.GENE_SET, Flags.METADATA),

    DBStopCodonReadThru("stopcodon_readthru", "", Flags.GENE_SET, Flags.VERSION_3_5),
    ExcludesGene("excludes_gene", "Gene/sequence-level attribute describing the incompatibility of presence of two or more alternate genes (generally used only in large DNA viruses).",
                 ConfigurationParameterFunctions.toListOfStrings,
                 Flags.VERSION_4, Flags.GENE_SET),

    ExoneratePath("exonerate_path", "Path to exonerate tool",
                  Flags.VERSION_4,
                  Flags.REQUIRED,
                  Flags.PLATFORM_DEPENDENT,
                  Flags.COMMANDLINE_SET,
                  Flags.PROGRAM_CONFIG_SET),
    FrameShiftSensitivity("frameshift_sensitivity", "Dictates the sensitivity VIGOR should use in handling frame-shifts. Accepted values: 0, 1, 2, with 2 being the strictest, forcing VIGOR to create pseudogenes and raising error messages whenever a perfect model for a given gene cannot be created (this is used most for validating assemblies).", Flags.VERSION_3, Flags.VERSION_4, Flags.UNIMPLEMENTED),
    GeneMinimumCoverage("min_gene_coverage", "Minimum coverage of genes",  toPercent, Flags.VERSION_3, Flags.VERSION_4), // TODO elaborate

    GeneOptional("is_optional", "This parameter works in combination with complete_genome and frameshift_sensitivity: if it is set to TRUE, the absence of that particular gene for which is set suppresses the stricter behavior set by the use of the other parameters.", ConfigurationParameterFunctions.isPresentOrBoolean, Flags.GENE_SET),
    GeneRequired("is_required", "By setting this parameter to TRUE, VIGOR will produce an error message in the case it cannot create a valid gene model for that given gene.", ConfigurationParameterFunctions.isPresentOrBoolean, Flags.GENE_SET),

    IntronMaximumSize("max_intron_size", "Maximum sequence length of an intron", toPositiveInteger, Flags.VERSION_3, Flags.VERSION_4),
    IntronMinimumSize("min_intron_size", "Minimum sequence length of an intron", toPositiveInteger, Flags.VERSION_3, Flags.VERSION_4),

    Locustag("locus_tag", "Prefix to be used for the locus_tag attribute in the feature table (.tbl) and in the GFF 3 outputs", Flags.VERSION_3, Flags.VERSION_4, Flags.COMMANDLINE_SET),
    MaturePeptideDB("matpepdb", "Location and name of the mature peptides reference database (e.g. matpepdb=<vigordata>/flua_ha_mp)", Flags.VERSION_4, Flags.GENE_SET),
    MaturePeptideMinimumCoverage("mature_pep_mincoverage", "Minimum percent coverage of the sequence of the reference mature peptide to consider the prediction valid.", ConfigurationParameterFunctions.toPercent,
                                 Flags.VERSION_3, Flags.VERSION_4),
    MaturePeptideMinimumIdentity("mature_pep_minidentity", "Minimum percent identity between the sequence of the reference mature peptide and the sequence of the candidate mature peptide.",
                                 toPercent,
                                 Flags.VERSION_3, Flags.VERSION_4),
    MaturePeptideMinimumSimilarity("mature_pep_minsimilarity", "Minimum percent similarity between the sequence of the reference mature peptide and the sequence of the candidate mature peptide. Note: similarity is calculated using compatibility values found in Blosum40 distance matrix",
                                   toPercent,
                                   Flags.VERSION_3, Flags.VERSION_4),

    MaxAlignMergeAAGap("max_align_merge_aa_gap", "Maximum number of proteins in a gap between two alignments to consider them for merging.", toPositiveInteger, Flags.VERSION_4),
    MaxGeneOverlap("max_gene_overlap", " In reporting gene models, maximum overlap of genes allowed.", toPositiveInteger, Flags.VERSION_4),
    MinFunctionalLength("min_functional_len" , "Minimum functional length for a protein (expressed in aa) to be functional: if a premature stop codon makes it shorter than that, it should be annotated as pseudogene.", toPositiveInteger, Flags.VERSION_4, Flags.GENE_SET),
    MinimumMissingAASize("min_missing_AA_size", "Minimum number of proteins missing in a given alignment to search for missing exons.", toPositiveInteger),

    NTOverlapMaximum("max_nt_overlap", "Maximum number of nucleotides that may overlap for alignment fragments to be considered compatible when generating a gene model", toInteger, Flags.VERSION_4),
    NonCanonicalSplicing("noncanonical_splicing", "List of alternative splicing donor and acceptor sequence pairs. Format: noncanonical_splicing=donor+acceptor,donor+acceptor,... (e.g. noncanonical_splicing=AA+GT)", Flags.VERSION_4, Flags.VIRUS_SET, Flags.GENE_SET),
    Note("note","Information associated to a gene or protein variant (to be reported in the .tbl and GFF 3 outputs).", toListOfStrings, Flags.VIRUS_SET, Flags.GENE_SET, Flags.METADATA_SET),
    OutputDirectory("output_directory", "Write output to this directory", Flags.VERSION_3, Flags.VERSION_4,
                    Flags.COMMANDLINE_SET, Flags.REQUIRED),

    OutputPrefix("output_prefix", "Use this prefix output files", Flags.VERSION_3, Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.REQUIRED),
    OverwriteOutputFiles("overwrite_output_files", "Overwrite output files if they exist",
                         toBoolean,
                         Flags.VERSION_4,
                         Flags.COMMANDLINE_SET,
                         Flags.PROGRAM_CONFIG_SET),
    PseudoGeneMinimumCoverage("min_pseudogene_coverage", "Minimum percentage of coverage in the alignment between candidate pseudogene and reference protein, based on the longest of the two.", toPercent, Flags.UNIMPLEMENTED),

    PseudoGeneMinimumIdentity("min_pseudogene_identity", "Minimum percentage of identity in the alignment between candidate pseudogene and reference protein, based on the longest of the two.", toPercent, Flags.UNIMPLEMENTED),

    PseudoGeneMinimumSimilarity("min_pseudogene_similarity", "Minimum percentage of similarity in the alignment between candidate pseudogene and reference protein, based on the longest of the two.", toPercent, Flags.UNIMPLEMENTED),

    RNAEditing("rna_editing", "RNA editing. Format is rna_editing=offset/regex/insertion string/note (e.g. rna_editing=0/GGGG/[ARWMDHVN][ARWMDHVN][ARWMDHVN][ARWMDHVN][GRSKBDVN][GRSKBDVN][GRSKBDVN]/four non-templated G's inserted during transcription)",
               ConfigurationParameterFunctions.of(RNA_Editing.class, RNA_Editing::parseFromString),
               Flags.VERSION_4, Flags.GENE_SET),
    ReferenceDatabaseFile("reference_database_file", "full path the reference database file",
                          Flags.REQUIRED, Flags.VERSION_4, Flags.COMMANDLINE_SET),

    ReferenceDatabasePath("reference_database_path", "Directory containing reference database files",
                          Flags.VERSION_4, Flags.REQUIRED, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET),

    RelaxAlignMergeAAGap("relax_align_merge_aa_gap", "", toPositiveInteger, Flags.VERSION_4),
    RibosomalSlippage("ribosomal_slippage",  "Ribosomal slippage. Format is ribosomal_slippage=offset/frameshift/regex (e.g. ribosomal_slippage=-7/+1/[BDHKNTWY][BCHMNSVY][BCHMNSVY][BDHKNTWY][BDHKNTWY][BDHKNTWY][BCHMNSVY][BDGKNRSV][BCDHKMNSTVWY][BCHMNSVY])",
                      ConfigurationParameterFunctions.of(Ribosomal_Slippage.class, Ribosomal_Slippage::parseFromString),
                      Flags.VERSION_4, Flags.GENE_SET),
    ScoreFactorAlignment("alignment_score_factor", "Weight of alignment score in evaluating alignment", toDouble, Flags.VERSION_4),
    ScoreFactorLeakyStop("leakystop_score_factor", "Weight for scoring leaky stops for gene model selection.", toDouble, Flags.VERSION_4),
    ScoreFactorLeakyStopNotFound("leakystop_notFound_score", "Weight to be assigned to the gene model candidates scoring when an expected leaky stop codon is not found (note: in some viruses the leaky stop codon is not always present).", toDouble, Flags.VERSION_4),
    ScoreFactorSplicing("splicing_score_factor", "Weight to apply on splice sites (e.g. canonical v.s. non-canonical, how far from expected location, etc.) in the scoring of gene models.", toDouble, Flags.VERSION_4),
    ScoreFactorStart("start_score_factor", "Weight to apply on start codons (e.g. canonical v.s. non-canonical, how far from expected location, etc.) in the scoring of gene models.", toDouble, Flags.VERSION_4),
    ScoreFactorStop("stop_score_factor", "Weight to apply on stop codons (e.g. how far from expected location, etc.) in the scoring of gene models.", toDouble, Flags.VERSION_4),
    SequenceGapMinimumLength("min_seq_gap_length", "Minimum number of undefined nucleotides (i.e. Ns) to consider it a sequencing gap.", toPositiveInteger, Flags.VERSION_4),

    SharedCDS("shared_cds", "List of other genes (gene symbols) sharing the same region of the viral genome. Format: shared_cds=Gene1_ID,Gene_2ID (e.g. shared_cds=NSP1-2; shared_cds=NSP1-1,NSP1-3)",
              ConfigurationParameterFunctions.toListOfStrings,
              Flags.VERSION_4, Flags.GENE_SET),
    SlippageFrameShift("slippage_frameshift", "Frame shift for ribosomal slippage", Flags.VERSION_3_5),
    SlippageMotif("slippage_motif", "Motif for ribosomal slippage", Flags.VERSION_3_5),
    SlippageOffset("slippage_offset", "Offset from motif for ribosomal slippage", Flags.VERSION_3_5),
    SpliceForm("splice_form", "Describes the expected gene structure with actual lengths of exons and introns of the reference protein aligned back to its originating genome (or the closest possible genome, if the original one is not available). Format: e\\d+(i-?\\d+e\\d+)* The string always starts and ends with exons. (e.g. splice_form=\"e26i686e268\"; splice_form=\"e1656\") Note: given that this parameter is specific for each single isoform, it is allowed only on the sequence header.", toSpliceForms, Flags.GENE_SET, Flags.VIRUS_SET),

    StartCodonSearchWindow("start_codon_search_window", "Number of nucleotides before and after a candidate site to check for a start codon",
                           toInteger, Flags.VERSION_4),
    StartCodons("start_codons", "Comma separated list of expected start codons", toListOfStrings),
    StopCodonReadthrough("stop_codon_readthrough",  "Format is amino acid/offset/regex (e.g.stop_codon_readthrough=-11/R/[NATC][NATC][NT][NG][NA][NC][NTG][NAG][NATG][NATCG][NTC][NAG][NAG])",
                         toStopException,
                         Flags.VERSION_4, Flags.GENE_SET),
    StopCodonSearchWindow("stop_codon_search_window", "Number of nucleotides before and after a candidate site to check for a stop codon",
                          toPositiveInteger, Flags.VERSION_4),
    TemporaryDirectory("temporary_directory", "Directory under which Vigor creates temporary files and directories",
                       Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET, Flags.REQUIRED),
    TinyExon3("tiny_exon3", "Tiny exon 3. Format is regex:[offset]", toTinyExonMap, Flags.VERSION_4, Flags.GENE_SET),
    TinyExon5("tiny_exon5", "Tiny exon 5. Format is regex:[offset]", toTinyExonMap, Flags.VERSION_4, Flags.GENE_SET),
    Variation("variation", "Variation. TODO", Flags.VERSION_3_5),
    Verbose("verbose", "Make console and .rpt file output more detailed", toBoolean,
            Flags.VERSION_3, Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET),
    Version("version", "Database version", Flags.METADATA_SET, Flags.VERSION_4),
    VirusName("virus_name", "Common name of the virus", Flags.METADATA_SET, Flags.VERSION_4),
    VirusSpecificConfiguration("virusSpecific_config", "Path to virus specific configuration file."
            , Flags.VERSION_4, Flags.COMMANDLINE_SET, Flags.PROGRAM_CONFIG_SET),

    VirusSpecificConfigurationPath("virusSpecific_config_path", "Directory containing virus specific configurations.",
                                   Flags.VERSION_4,
                                   Flags.COMMANDLINE_SET,
                                   Flags.PROGRAM_CONFIG_SET);


    static final Map<String, ConfigurationParameters> byConfigKey;

    static {
        byConfigKey = Arrays.stream(ConfigurationParameters.values()).collect(Collectors.toMap(cp -> cp.configKey, cp -> cp));
    }

    public final String configKey;
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
        // settable in metadata section
        METADATA_SET,
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

        public static final EnumSet<Flags> SET_FLAGS = EnumSet.of(Flags.PROGRAM_CONFIG_SET,
                                                                  Flags.COMMANDLINE_SET,
                                                                  Flags.VIRUS_SET,
                                                                  Flags.GENE_SET);

        public static final EnumSet<Flags> ALL_VERSIONS = EnumSet.of(Flags.VERSION_3,
                                                                     Flags.VERSION_3_5,
                                                                     Flags.VERSION_4);
    }

    ConfigurationParameters(String configKey, String description, Flags... flags) {
        this(configKey, description, null,  flags);
    }



    ConfigurationParameters(String configKey, String description, ValueFunction valueFunction, Flags... flags) {
        this.configKey = configKey;
        this.description = description;
        this.valueFunction = valueFunction == null ? ConfigurationParameterFunctions.of(String.class, s -> s) : valueFunction;

        if (flags.length == 0) {
            flags = new Flags[]{Flags.VERSION_3, Flags.VERSION_4};
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
        result.removeIf(e -> !checkFlags.contains(e));
        return result;
    }

    public boolean hasOneOrMoreFlags(Flags... checkFlags) {
        return Arrays.stream(checkFlags).anyMatch(flags::contains);
    }

    public boolean hasAllFlags(Flags... checkFlags) {
        return Arrays.stream(checkFlags).allMatch(flags::contains);
    }

    public String getEnvVarName() {
        return "VIGOR_" + configKey.toUpperCase();
    }

    public String getSystemPropertyName() {
        return "vigor." + configKey;
    }

    public String valueToString(Object o) {
        try {
            return valueFunction.stringFunction.apply(o);
        } catch (InvalidValue | ClassCastException e) {
            throw new InvalidValue(String.format("bad value for %s: \"%s\":  %s",
                                                 this.configKey,
                                                 o,
                                                 e.getMessage()));
        } catch (RuntimeException e) {
            throw new ConfigurationParameterFunctions.InvalidValue(String.format("bad value for %s: \"%s\": got %s %s",
                                                                                 this.configKey,
                                                                                 o,
                                                                                 e.getClass().getSimpleName(),
                                                                                 e.getMessage())
            );
        }

    }

    public Object stringToValue(String stringValue) {
        try {
            return this.valueFunction.valueClass.cast(valueFunction.valueFunction.apply(stringValue));
        } catch (InvalidValue | ClassCastException e) {
            throw new InvalidValue(String.format("bad value for %s: \"%s\":  %s",
                                                 this.configKey,
                                                 stringValue,
                                                 e.getMessage()));
        } catch (RuntimeException e) {
            throw new ConfigurationParameterFunctions.InvalidValue(String.format("bad value for %s: \"%s\": got %s %s",
                                                                                 this.configKey,
                                                                                 stringValue,
                                                                                 e.getClass().getSimpleName(),
                                                                                 e.getMessage())
            );
        }
    }

}
