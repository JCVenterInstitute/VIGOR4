package org.jcvi.vigor.service;

import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.ConfigurationUtils;
import org.jcvi.vigor.utils.VigorConfiguration;
import org.jcvi.vigor.utils.VigorUtils;
import org.apache.commons.lang3.StringUtils;
import org.jcvi.vigor.component.Splicing.SpliceSite;
import org.jcvi.vigor.forms.VigorForm;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.springframework.stereotype.Service;
import sun.security.krb5.Config;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

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
    public Alignment setViralProteinAttributes ( Alignment alignment, VigorForm form ) throws VigorException {
        ViralProtein viralProtein = alignment.getViralProtein();
        Map<String, String> attributes = parseDeflineAttributes(StringUtils.normalizeSpace(viralProtein.getDefline()),
                                                                viralProtein.getProteinID());


        VigorConfiguration defaultConfig = form.getConfiguration();
        viralProtein = setProteinAttributes(viralProtein, attributes);
        viralProtein.setConfiguration(getGeneConfiguration(viralProtein, defaultConfig, attributes));
        /* set geneStructure property of viralProtein */
        viralProtein = setGeneAttributes(alignment, defaultConfig);
        viralProtein = setGeneAttributesFromConfig(alignment, viralProtein.getConfiguration());
        int min_intron_length = viralProtein.getConfiguration().getOrDefault(ConfigurationParameters.IntronMinimumSize, 0);
        viralProtein = determineGeneStructure(viralProtein, min_intron_length);

        alignment.setViralProtein(viralProtein);
        return alignment;
    }

    // TODO make sure configuration values are appropriate
    private ViralProtein setGeneAttributesFromConfig(Alignment alignment, VigorConfiguration geneConfig) {
        ViralProtein viralProtein = alignment.getViralProtein();
        GeneAttributes attributes = viralProtein.getGeneAttributes();
        attributes.setRibosomal_slippage(geneConfig.get(ConfigurationParameters.RibosomalSlippage));
        attributes.setRna_editing(geneConfig.get(ConfigurationParameters.RNAEditing));
        attributes.setSplicing(geneConfig.get(ConfigurationParameters.Splicing));
        attributes.setStopTranslationException(geneConfig.get(ConfigurationParameters.StopCodonReadthrough));

        List<String> alternateStarts = geneConfig.getOrDefault(ConfigurationParameters.AlternateStartCodons, Collections.EMPTY_MAP);
        if (! alternateStarts.isEmpty()) {
            attributes.setStartTranslationException(new StartTranslationException(true, alternateStarts));
        }
        StructuralSpecifications specs = attributes.getStructuralSpecifications();
        specs.setShared_cds(geneConfig.get(ConfigurationParameters.SharedCDS));
        specs.setMinFunctionalLength(geneConfig.get(ConfigurationParameters.MinFunctionalLength));
        specs.setExcludes_gene(geneConfig.get(ConfigurationParameters.ExcludesGene));
        specs.setTiny_exon3(geneConfig.get(ConfigurationParameters.TinyExon3));
        specs.setTiny_exon5(geneConfig.get(ConfigurationParameters.TinyExon5));
        // TODO more

        return viralProtein;
    }

    public ViralProtein setGeneAttributes ( Alignment alignment, VigorConfiguration configuration ) throws VigorException {
        ViralProtein viralProtein = alignment.getViralProtein();
        Map<String, String> attributes = parseDeflineAttributes(StringUtils.normalizeSpace(viralProtein.getDefline()),
                                                                viralProtein.getProteinID());
        return setGeneAttributes(alignment, attributes, configuration);
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

    /**
     * @param alignment
     * @return viralProtein object which has all the gene attributes set.
     */
    public ViralProtein setGeneAttributes ( Alignment alignment,
                                            Map<String,String> attributes,
                                            VigorConfiguration configuration ) throws VigorException {

        ViralProtein viralProtein = alignment.getViralProtein();
        GeneAttributes geneAttributes = new GeneAttributes();
        StructuralSpecifications structuralSpecifications = geneAttributes.getStructuralSpecifications();


        /* Set Splicing attributes */
        if (! attributes.getOrDefault("splice_form","").isEmpty()) {
            List<SpliceSite> splicePairs = new ArrayList<SpliceSite>();
            if (! attributes.getOrDefault("noncanonical_splicing","N").equalsIgnoreCase("N")) {
                Pattern.compile(";").splitAsStream(attributes.get("noncanonical_splicing").trim())
                       .map(s -> s.split(Pattern.quote("+")))
                       .map(a -> new Splicing.SpliceSite(a[ 0 ], a[ 1 ]))
                       .forEach(splicePairs::add);
            }
            splicePairs.add(new Splicing.SpliceSite("GT", "AG"));
            geneAttributes.setSplicing(new Splicing(true, splicePairs, attributes.get("splice_form")));

            //TODO this only set if splicing?
            if (attributes.containsKey("tiny_exon3")) {
                String attribute = attributes.get("tiny_exon3");
                int offset = 0;
                String[] temp = attribute.split(":",2);
                if (temp.length == 2) {
                    String regex = temp[ 0 ];
                    if (VigorUtils.is_Integer(temp[ 1 ])) {
                        offset = Integer.parseInt(temp[ 1 ]);
                        // TODO can not be specified for 0?
                    } else if (! temp[1].trim().isEmpty()) {
                        LOGGER.warn("For protein {} bad value for offset in tiny_exon3 {}, full value {}",
                                    viralProtein.getProteinID(), temp[1], attribute);
                    }
                    Map<String, Integer> temp1 = new HashMap<String, Integer>();
                    temp1.put(regex, offset);
                    structuralSpecifications.setTiny_exon3(temp1);
                } else {
                    LOGGER.warn("For protein {} bad value for tiny_exon3 {}", viralProtein.getProteinID(), attribute);
                }
            }

            //TODO this only set if splicing?
            if (attributes.containsKey("tiny_exon5")) {
                String attribute = attributes.get("tiny_exon5");
                String regex = null;
                int offset = 0;
                if (attribute.matches(".*?:.*")) {
                    String[] temp = attribute.split(":");
                    regex = temp[ 0 ];
                    if (VigorUtils.is_Integer(temp[ 1 ])) {
                        offset = Integer.parseInt(temp[ 1 ]);
                    } else {
                        regex = attribute;
                    }
                    Map<String, Integer> tempMap = new HashMap<String, Integer>();
                    tempMap.put(regex, offset);
                    structuralSpecifications.setTiny_exon5(tempMap);
                } else {
                    LOGGER.warn("For protein {} unexpected value for tiny_exon5 {}", viralProtein.getProteinID(), attribute);
                }
            }
        }


        /* Set RibosomalSlippage attributes */
        if (attributes.containsKey("V4_Ribosomal_Slippage")) {
            String attribute = attributes.get("V4_Ribosomal_Slippage");
            if (! (attribute == null || attribute.isEmpty() )) {
                geneAttributes.setRibosomal_slippage(Ribosomal_Slippage.parseFromString(attribute));
            }
        }

        /* Set StopTranslationException attributes */
        if (attributes.containsKey("V4_stop_codon_readthrough")) {
            String attribute = attributes.get("V4_stop_codon_readthrough");
            String[] temp = attribute.split("/");
            geneAttributes.setStopTranslationException(new StopTranslationException(true, AminoAcid.parse(temp[ 1 ]), temp[ 2 ], Integer.parseInt(temp[ 0 ])));
        }

        /* Set StartTranslationException attributes */
        if (attributes.containsKey("alternate_startcodon")) {
            String attribute = attributes.get("alternate_startcodon");
            geneAttributes.setStartTranslationException(new StartTranslationException(true, Arrays.asList(attribute.split(","))));
        }

        /* Set RNA_Editing attributes */
        if (attributes.containsKey("V4_rna_editing")) {
            String attribute = attributes.get("V4_rna_editing");
            if (! (attribute == null || attribute.isEmpty())) {
                geneAttributes.setRna_editing(RNA_Editing.parseFromString(attribute));
            }
        }

        /* Set StructuralSpecifications */
        if (attributes.containsKey("shared_cds")) {
            String attribute = attributes.get("shared_cds");
            structuralSpecifications.setShared_cds(Arrays.asList(attribute.split(",")));
        }
        if (attributes.containsKey("is_optional")) {
            structuralSpecifications.set_required(false);
        }
        if (attributes.containsKey("is_required")) {
            structuralSpecifications.set_required(true);
        }
        if (attributes.containsKey("excludes_gene")) {
            String attribute = attributes.get("excludes_gene");
            structuralSpecifications.setExcludes_gene(Arrays.asList(attribute.split(",")));
        }
        if (attributes.containsKey("min_functional_len")) {
            String attribute = attributes.get("min_functional_len");
            if (VigorUtils.is_Integer(attribute)) {
                structuralSpecifications.setMinFunctionalLength(Integer.parseInt(attribute));
            } else {
                // TODO fatal?
                LOGGER.warn("For protein {} bad value for min_functional_len {}", viralProtein.getProteinID(), attribute);
            }
        }

        /* Set maturepeptide DB attribute */
        String matPepDB = attributes.getOrDefault("matpepdb", "");
        if (matPepDB != null && matPepDB.contains("<vigordata>")) {
            matPepDB = matPepDB.replace("<vigordata>", configuration.get(ConfigurationParameters.ReferenceDatabasePath));
        }

        AlignmentEvidence alignmentEvidence = alignment.getAlignmentEvidence();
        alignmentEvidence.setMatpep_db(matPepDB);


        /* set geneAttributes property of viralProtein */
        viralProtein.setGeneAttributes(geneAttributes);
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
        boolean is_spliced = viralProtein.getGeneAttributes().getSplicing().isSpliced();
        List<Range> NTFragments = new ArrayList<Range>();
        List<Range> introns = new ArrayList<Range>();
        if (!( is_ribosomal_slippage ) && !( is_spliced )) {
            Range range = Range.of(0, 3 * ( viralProtein.getSequence().getLength() - 1 ));
            NTFragments.add(range);
        } else {
            String spliceform = viralProtein.getGeneAttributes().getSplicing().getSpliceform();
            if (spliceform != null) {
                if (spliceform.equals("") || !( spliceform.matches("([i,e]-?[0-9]*)+") )) {
                    String exception = String.format("For protein %s spliced reference missing/malformed splice_form tag: %s",
                            viralProtein.getProteinID(), spliceform);
                    LOGGER.warn(exception);
                    throw new ServiceException(exception);
                } else {
                    List<String> splices = new ArrayList<String>();
                    Matcher m = Pattern.compile("[e,i]-?[0-9]*").matcher(spliceform);
                    while (m.find()) {
                        splices.add(m.group(0));
                    }
                    Long dnaOffset = 0l;
                    for (int i = 0; i < splices.size(); i++) {
                        long nucleotides = Long.parseLong(splices.get(i).substring(1));
                        if (nucleotides > 0) {
                            Range range = Range.of(dnaOffset, ( dnaOffset + nucleotides ) - 1);
                            if (splices.get(i).matches("e-?[0-9]*")) {
                                dnaOffset = range.getEnd() + 1;
                                NTFragments.add(range);
                            } else {
                                if (nucleotides >= min_intron_length) {
                                    introns.add(range);
                                }
                                dnaOffset = range.getEnd() + 1;
                            }
                        } else {
                            dnaOffset = dnaOffset - nucleotides;
                        }
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

        Pattern attributePattern = Pattern.compile("(?<key>\\b(?<!>)\\w+\\b)(?:\\s*=\\s*(?<value>\"[^\"]*\"|\\w+))?");
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
            deflineConfig = ConfigurationUtils.configurationFromMap("defline: " + viralProtein.getProteinID(),attributes);
        }

        String geneSection = ConfigurationUtils.getGeneSectionName(viralProtein.getGeneSymbol());
        if (defaultConfig.hasSection(geneSection)) {
            Map<ConfigurationParameters,Object> geneConfigMap = defaultConfig.getSectionConfig(geneSection);
            if (! geneConfigMap.isEmpty()) {
                VigorConfiguration geneConfig = new VigorConfiguration("config: " + viralProtein.getGeneSymbol(), defaultConfig);
                geneConfig.putAll(geneConfigMap);
                defaultConfig = geneConfig;
            }
        }
        deflineConfig.setDefaults(defaultConfig);
        return deflineConfig;
    }
}
