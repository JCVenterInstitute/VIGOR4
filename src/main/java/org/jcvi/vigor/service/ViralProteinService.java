package org.jcvi.vigor.service;

import org.jcvi.vigor.component.*;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.ServiceException;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.VigorUtils;
import org.apache.commons.lang3.StringUtils;
import org.jcvi.vigor.component.Splicing.SpliceSite;
import org.jcvi.vigor.forms.VigorForm;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.Range;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.springframework.stereotype.Service;

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
    private String matPepDB = "";

    /**
     * @param alignment
     * @return viralProtein: For the given protein ID ViralProtein object is
     * generated and all the properties are defined;
     */
    public Alignment setViralProteinAttributes ( Alignment alignment, VigorForm form ) throws VigorException {

        ViralProtein viralProtein = setGeneAttributes(alignment.getViralProtein(), form);
        AlignmentEvidence alignmentEvidence = alignment.getAlignmentEvidence();
        alignmentEvidence.setMatpep_db(matPepDB);
        alignment.setAlignmentEvidence(alignmentEvidence);
        alignment.setViralProtein(viralProtein);
        return alignment;
    }

    /**
     * @param viralProtein
     * @return viralProtein object which has all the gene attributes set.
     */
    public ViralProtein setGeneAttributes ( ViralProtein viralProtein, VigorForm form ) throws VigorException {

        GeneAttributes geneAttributes = new GeneAttributes();
        String defline = viralProtein.getDefline();
        defline = StringUtils.normalizeSpace(defline);
        Map<String, String> attributes = parseDeflineAttributes(defline, viralProtein.getProteinID());
        StructuralSpecifications structuralSpecifications = new StructuralSpecifications();
        Splicing splicing = Splicing.NO_SPLICING;
        Ribosomal_Slippage ribosomal_slippage = Ribosomal_Slippage.NO_SLIPPAGE;
        RNA_Editing rna_editing = RNA_Editing.NO_EDITING;
        StopTranslationException stopTranslationException = StopTranslationException.NO_EXCEPTION;
        StartTranslationException startTranslationException = StartTranslationException.NO_EXCEPTION;

        /* Set product attribute */
        if (attributes.containsKey("product")) {
            String product = attributes.get("product").replaceAll("\"", "");
            viralProtein.setProduct(product);
        }

        /* Set gene symbol */
        if (attributes.containsKey("gene")) {
            String geneSymbol = attributes.get("gene");
            viralProtein.setGeneSymbol(geneSymbol);
        }
        /* Set gene symbol */
        if (attributes.containsKey("gene_synonym")) {
            String geneSynonym = attributes.get("gene_synonym");
            viralProtein.setGeneSynonym(geneSynonym);
        }

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
            splicing = new Splicing(true, splicePairs, attributes.get("splice_form"));
            //TODO this only set if splicing?
            if (attributes.containsKey("tiny_exon3")) {
                String attribute = attributes.get("tiny_exon3");
                int offset = 0;
                String[] temp = attribute.split(":");
                String regex = temp[ 0 ];
                if (VigorUtils.is_Integer(temp[ 1 ])) {
                    offset = Integer.parseInt(temp[ 1 ]);
                }
                Map<String, Integer> temp1 = new HashMap<String, Integer>();
                temp1.put(regex, offset);
                structuralSpecifications.setTiny_exon3(temp1);
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
                }
                Map<String, Integer> temp = new HashMap<String, Integer>();
                temp.put(regex, offset);
                structuralSpecifications.setTiny_exon5(temp);
            }
        }


        /* Set RibosomalSlippage attributes */
        if (attributes.containsKey("V4_Ribosomal_Slippage")) {
            String attribute = attributes.get("V4_Ribosomal_Slippage");
            if (! (attribute == null || attribute.isEmpty() )) {
                ribosomal_slippage = Ribosomal_Slippage.parseFromString(attribute);
            }
        }

        /* Set StopTranslationException attributes */
        if (attributes.containsKey("V4_stop_codon_readthrough")) {
            String attribute = attributes.get("V4_stop_codon_readthrough");
            String[] temp = attribute.split("/");
            stopTranslationException = new StopTranslationException(true, AminoAcid.parse(temp[ 1 ]), temp[ 2 ], Integer.parseInt(temp[ 0 ]));
        }

        /* Set StartTranslationException attributes */
        if (attributes.containsKey("alternate_startcodon")) {
            String attribute = attributes.get("alternate_startcodon");
            startTranslationException = new StartTranslationException(true, Arrays.asList(attribute.split(",")));
        }

        /* Set RNA_Editing attributes */
        if (attributes.containsKey("V4_rna_editing")) {
            String attribute = attributes.get("V4_rna_editing");
            if (! (attribute == null || attribute.isEmpty())) {
                rna_editing = RNA_Editing.parseFromString(attribute);
            }
        }

        /* Set StructuralSpecifications */
        if (attributes.containsKey("shared_cds")) {
            String attribute = attributes.get("shared_cds");
            structuralSpecifications.setShared_cds(Arrays.asList(attribute.split(",")));
        }
        if (attributes.containsKey("is_optional")) {
            structuralSpecifications.set_required(true);
        }
        if (attributes.containsKey("is_required")) {
            structuralSpecifications.set_required(false);
        }
        if (attributes.containsKey("excludes_gene")) {
            String attribute = attributes.get("excludes_gene");
            structuralSpecifications.setExcludes_gene(Arrays.asList(attribute.split(",")));
        }
        if (attributes.containsKey("min_functional_len")) {
            String attribute = attributes.get("min_functional_len");
            if (VigorUtils.is_Integer(attribute)) {
                structuralSpecifications.setMinFunctionalLength(Integer.parseInt(attribute));
            }
        }

        /* Set maturepeptide DB attribute */
        matPepDB = attributes.getOrDefault("matpepdb", "");
        if (matPepDB != null && matPepDB.contains("<vigordata>")) {
            matPepDB = matPepDB.replace("<vigordata>", form.getConfiguration().get(ConfigurationParameters.ReferenceDatabasePath));
        }
        //LOGGER.debug("For protein {} setting mapPepDB to {}", viralProtein.getProteinID(), matPepDB);
        /* Move all the different attribute objects to geneAttributes */
        geneAttributes.setRibosomal_slippage(ribosomal_slippage);
        geneAttributes.setRna_editing(rna_editing);
        geneAttributes.setSplicing(splicing);
        geneAttributes.setStartTranslationException(startTranslationException);
        geneAttributes.setStopTranslationException(stopTranslationException);
        geneAttributes.setStructuralSpecifications(structuralSpecifications);

        /* set geneAttributes property of viralProtein */
        viralProtein.setGeneAttributes(geneAttributes);

        /* set geneStructure property of viralProtein */
        int min_intron_length = form.getConfiguration().getOrDefault(ConfigurationParameters.IntronMinimumSize, 0);
        viralProtein = DetermineGeneStructure(viralProtein, min_intron_length);
        return viralProtein;
    }

    /*
     *
     * @param viralProtein
     * @return GeneStructure: has list of exons and introns of the viralProtein.
     *         These are determined from spliceform annotated in the defline*/

    public ViralProtein DetermineGeneStructure ( ViralProtein viralProtein, int min_intron_length ) throws ServiceException {

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
}
