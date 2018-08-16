package org.jcvi.vigor.service;

import com.google.common.collect.ImmutableList;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.ConfigurationUtils;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorUtils;
import org.apache.commons.lang3.StringUtils;
import org.jcvi.vigor.component.SpliceSite;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.springframework.stereotype.Service;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by snettem on 5/16/2017.
 */
@Service
public class ViralProteinService {

    private static final Logger LOGGER = LogManager.getLogger(ViralProteinService.class);

    /**
     * @param alignment
     * @return viralProtein: For the given protein ID ViralProtein object is
     * generated and all the properties are defined;
     */
    public Alignment setViralProteinAttributes ( Alignment alignment, VigorConfiguration config ) throws VigorException {
        ViralProtein viralProtein = alignment.getViralProtein();
        Map<String, String> attributes = parseDeflineAttributes(StringUtils.normalizeSpace(viralProtein.getDefline()),
                                                                viralProtein.getProteinID());


        VigorConfiguration defaultConfig = config;
        viralProtein = setProteinAttributes(viralProtein, attributes);
        viralProtein.setConfiguration(getGeneConfiguration(viralProtein, defaultConfig, attributes));
        /* set geneStructure property of viralProtein */
        viralProtein = setAttributesFromConfig(alignment, viralProtein.getConfiguration());
        int min_intron_length = viralProtein.getConfiguration().getOrDefault(ConfigurationParameters.IntronMinimumSize, 0);
        viralProtein = determineGeneStructure(viralProtein, min_intron_length);

        alignment.setViralProtein(viralProtein);
        return alignment;
    }

    // TODO make sure configuration values are appropriate
    private ViralProtein setAttributesFromConfig(Alignment alignment, VigorConfiguration geneConfig) {
        ViralProtein viralProtein = alignment.getViralProtein();
        String product = viralProtein.getProduct();
        viralProtein.setProduct(product != null ? product : geneConfig.getOrDefault(ConfigurationParameters.DBProduct,""));
        String geneSynonym = viralProtein.getGeneSynonym();
        viralProtein.setGeneSynonym(geneSynonym != null ? geneSynonym : geneConfig.getOrDefault(ConfigurationParameters.DBGeneSynonym,""));

        String matpepDB = alignment.getAlignmentEvidence().getMatpep_db();
        matpepDB = ! VigorUtils.nullElse(matpepDB,"").isEmpty() ? matpepDB: geneConfig.getOrDefault(ConfigurationParameters.MaturePeptideDB, "");
        matpepDB = matpepDB.replace("<vigordata>", geneConfig.get(ConfigurationParameters.ReferenceDatabasePath));
        LOGGER.trace("setting mature peptide db for {} (gene {}) to {}", viralProtein.getProteinID(), viralProtein.getGeneSymbol(), matpepDB);
        alignment.getAlignmentEvidence().setMatpep_db(matpepDB);

        GeneAttributes attributes = viralProtein.getGeneAttributes();
        attributes.setRibosomal_slippage(geneConfig.getOrDefault(ConfigurationParameters.RibosomalSlippage, Ribosomal_Slippage.NO_SLIPPAGE));
        attributes.setRna_editing(geneConfig.getOrDefault(ConfigurationParameters.RNAEditing, RNA_Editing.NO_EDITING));
        List<SpliceSite> nonCanonicalSpliceSites = geneConfig.getOrDefault(ConfigurationParameters.NonCanonicalSplicing, Collections.EMPTY_LIST);
        if (nonCanonicalSpliceSites.isEmpty()) {
            attributes.setSpliceSites(SpliceSite.DEFAULT_SPLICE_SITES);
        } else {
            attributes.setSpliceSites(Stream.concat(nonCanonicalSpliceSites.stream(),
                                                    SpliceSite.DEFAULT_SPLICE_SITES.stream())
                                            .collect(ImmutableList.toImmutableList()));
        }
        attributes.setStopTranslationException(geneConfig.getOrDefault(ConfigurationParameters.StopCodonReadthrough, StopTranslationException.NO_EXCEPTION));
        attributes.setSpliceForms(geneConfig.getOrDefault(ConfigurationParameters.SpliceForm, Collections.EMPTY_LIST));
        // TODO unify start codons specified via config or command line and these.
        List<String> alternateStarts = geneConfig.getOrDefault(ConfigurationParameters.AlternateStartCodons, Collections.EMPTY_LIST);
        if (! alternateStarts.isEmpty()) {
            attributes.setStartTranslationException(new StartTranslationException(true, alternateStarts));
        }
        StructuralSpecifications specs = attributes.getStructuralSpecifications();
        specs.setShared_cds(geneConfig.getOrDefault(ConfigurationParameters.SharedCDS, Collections.EMPTY_LIST));
        specs.setMinFunctionalLength(geneConfig.getOrDefault(ConfigurationParameters.MinFunctionalLength, 0));
        specs.setExcludes_gene(geneConfig.getOrDefault(ConfigurationParameters.ExcludesGene, Collections.EMPTY_LIST));
        specs.setTiny_exon3(geneConfig.getOrDefault(ConfigurationParameters.TinyExon3, Collections.EMPTY_MAP));
        specs.setTiny_exon5(geneConfig.getOrDefault(ConfigurationParameters.TinyExon5, Collections.EMPTY_MAP));

        Boolean required = geneConfig.getOrDefault(ConfigurationParameters.GeneRequired, null);
        Boolean optional = geneConfig.getOrDefault(ConfigurationParameters.GeneOptional, null);

        if (optional != null && required != null) {
            LOGGER.warn("optional and required set for {}", viralProtein.getProteinID());
        } else {
            StructuralSpecifications structuralSpecifications = viralProtein.getGeneAttributes().getStructuralSpecifications();
            if (required != null) {
                structuralSpecifications.set_required(required);
            } else if (required != null) {
                structuralSpecifications.set_required(! optional);
            }
        }
        // TODO more

        return viralProtein;
    }


    public ViralProtein setProteinAttributes(ViralProtein viralProtein, Map<String,String> attributes) {

        /* Set product attribute */
        if (attributes.containsKey("product")) {
            viralProtein.setProduct(VigorUtils.removeQuotes(attributes.get("product")));
        }

        /* Set gene symbol */
        if (attributes.containsKey("gene")) {
            viralProtein.setGeneSymbol(attributes.get("gene"));
        }
        /* Set gene symbol */
        if (attributes.containsKey("gene_synonym")) {
            viralProtein.setGeneSynonym(attributes.get("gene_synonym"));
        }
        return viralProtein;
    }


    /*
     *
     * @param viralProtein
     * @return GeneStructure: has list of exons and introns of the viralProtein.
     *         These are determined from spliceform annotated in the defline*/

    public ViralProtein determineGeneStructure(ViralProtein viralProtein, int min_intron_length ) throws ServiceException {

        boolean is_ribosomal_slippage = viralProtein.getGeneAttributes().getRibosomal_slippage()
                .isHas_ribosomal_slippage();
        List<SpliceForm> splices = viralProtein.getGeneAttributes().getSpliceForms();
        List<Range> NTFragments = new ArrayList<Range>();
        List<Range> introns = new ArrayList<Range>();
        if (! is_ribosomal_slippage && splices.isEmpty()) {
            Range range = Range.of(0, 3 * ( viralProtein.getSequence().getLength() - 1 ));
            NTFragments.add(range);
        } else {
            if (! splices.isEmpty()) {
                Long dnaOffset = 0l;
                for (SpliceForm spliceForm : splices) {
                    if (spliceForm.length > 0) {
                        Range range = Range.of(dnaOffset, (dnaOffset + spliceForm.length) - 1);
                        if (spliceForm.type == SpliceForm.SpliceType.EXON) {
                            NTFragments.add(range);
                        } else if (spliceForm.length >= min_intron_length) {
                            introns.add(range);
                        }
                        dnaOffset = range.getEnd() + 1;
                    } else {
                        dnaOffset = dnaOffset - spliceForm.length;
                    }
                }
            } else {
                String exception = "Spliced reference missing for " + viralProtein.getProteinID();
                LOGGER.debug(exception);
            }
        }
        viralProtein.setNTfragments(NTFragments);
        viralProtein.setIntrons(introns);
        return viralProtein;
    }

    /**
     * @param defline of the protein
     * @return List of attributes in the defline
     */
    public static Map<String,String> parseDeflineAttributes ( String defline, String id ) throws VigorException {

        Pattern attributePattern = Pattern.compile("(?<key>\\b(?<!>)\\w+\\b)(?:\\s*=\\s*(?<value>\"[^\"]*\"|\\S+))?");
        Matcher matcher = attributePattern.matcher(defline);
        Map<String,String> attributes = new HashMap<>();
        String key,value,error;
        List<String> errors = new ArrayList<>();
        List<String> warnings = new ArrayList<>();

        while (matcher.find()) {
            key = matcher.group("key");
            value = VigorUtils.removeQuotes(VigorUtils.nullElse(matcher.group("value"),"").trim());
            if (attributes.containsKey(key)) {
                error = String.format("duplicate key \"%s\" previous value \"%s\" new value \"%s\"", key, attributes.get(key), value);
                if (value.equals(attributes.get(key))) {
                    warnings.add(error);
                } else {
                    errors.add(error);
                }
            } else {
                attributes.put(key, value);
            }
        }
        if (! warnings.isEmpty()) {
            LOGGER.warn(warnings.stream()
                                .collect(Collectors.joining("\n",
                                                            String.format("Problem parsing defline attributes for %s:\n", id),
                                                            "\nattributes line: " + defline)));
        }
        if (! errors.isEmpty()) {
            throw new VigorException(errors.stream()
                                           .collect(Collectors.joining("\n",
                                                                       String.format("Problem parsing defline attributes for %s:\n", id),
                                                                       "\nattributes line: " + defline))
            );
        }
        return attributes;
    }

    private VigorConfiguration getGeneConfiguration(ViralProtein viralProtein,
                                                    VigorConfiguration defaultConfig,
                                                    Map<String, String> attributes) throws VigorException {
        VigorConfiguration deflineConfig = new VigorConfiguration("defline: " + viralProtein.getProteinID());
        if (! attributes.isEmpty()) {
            // use defline config entries where they exist
            deflineConfig = ConfigurationUtils.configurationFromMap("defline: " + viralProtein.getProteinID(), p->p.deflineConfigKey, attributes, ConfigurationParameters.Flags.GENE_SET);
        }

        String geneSection = ConfigurationUtils.getGeneSectionName(viralProtein.getGeneSymbol());
        if (defaultConfig.hasSection(geneSection)) {
            LOGGER.trace("found gene section for {}", viralProtein.getGeneSymbol());
            Map<ConfigurationParameters,Object> geneConfigMap = defaultConfig.getSectionConfig(geneSection);
            if (! geneConfigMap.isEmpty()) {
                VigorConfiguration geneConfig = new VigorConfiguration("config: " + viralProtein.getGeneSymbol(), defaultConfig);
                geneConfig.putAll(geneConfigMap);
                defaultConfig = geneConfig;
            }
        } else {
            LOGGER.trace("No gene section for {}", viralProtein.getGeneSymbol());
        }
        deflineConfig.setDefaults(defaultConfig);
        return deflineConfig;
    }
}
