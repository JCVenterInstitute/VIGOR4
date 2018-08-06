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
    private int min_intron_length;

    /**
     * @param alignment
     * @return viralProtein: For the given protein ID ViralProtein object is
     * generated and all the properties are defined;
     */
    public Alignment setViralProteinAttributes ( Alignment alignment, VigorForm form ) throws VigorException {

        min_intron_length = form.getConfiguration().getOrDefault(ConfigurationParameters.IntronMinimumSize, 0);
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
        List<String> deflineAttributes = parseDeflineAttributes(defline);
        Map<String, String> attributes = deflineAttributes.stream()
                .map(s -> s.split("=", 2))
                .collect(Collectors.toMap(a -> a[ 0 ].trim(),
                        a -> a.length > 1 ? VigorUtils.removeQuotes(a[ 1 ]) : "",
                        ( s1, s2 ) -> s1));
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
        if (attributes.containsKey("splice_form")) {
            if (!( attributes.get("splice_form").equals("") )) {
                List<SpliceSite> splicePairs = new ArrayList<SpliceSite>();
                if (attributes.containsKey("noncanonical_splicing")) {
                    if (!( attributes.get("noncanonical_splicing").equalsIgnoreCase("N") )) {
                        List<String> spliceSites = Pattern.compile(";").splitAsStream(attributes.get("noncanonical_splicing").trim()).collect(Collectors.toList());
                        for (String spliceSite : spliceSites) {
                            String[] temp = spliceSite.split(Pattern.quote("+"));
                            splicePairs.add(splicing.new SpliceSite(temp[ 0 ], temp[ 1 ]));
                        }
                    }
                }
                SpliceSite defaultSpliceSite = splicing.new SpliceSite("GT", "AG");
                splicePairs.add(defaultSpliceSite);
                splicing = new Splicing(true, splicePairs, attributes.get("splice_form"));
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
        }


        /* Set RibosomalSlippage attributes */
        if (attributes.containsKey("V4_Ribosomal_Slippage")) {
            String attribute = attributes.get("V4_Ribosomal_Slippage");
            String[] temp = attribute.split("/");
            if (temp.length == 3) {
                ribosomal_slippage = new Ribosomal_Slippage(true, temp[ 2 ], Integer.parseInt(temp[ 0 ]), Integer.parseInt(temp[ 1 ]));
            } else {
                LOGGER.warn("bad ribosomal slippage attribute {}", attribute);
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
            String[] temp = attribute.split("/");
            int rna_editing_offset = 0;
            if (VigorUtils.is_Integer(temp[ 0 ])) {
                rna_editing_offset = Integer.parseInt(temp[ 0 ]);
            }
            rna_editing = new RNA_Editing(true, rna_editing_offset, temp[ 2 ], temp[ 1 ], temp[ 3 ]);
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
        viralProtein = DetermineGeneStructure(viralProtein);
        return viralProtein;
    }

    /*
     *
     * @param viralProtein
     * @return GeneStructure: has list of exons and introns of the viralProtein.
     *         These are determined from spliceform annotated in the defline*/

    public ViralProtein DetermineGeneStructure ( ViralProtein viralProtein ) throws ServiceException {

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
    public List<String> parseDeflineAttributes ( String defline ) throws VigorException {

        Pattern pattern;
        Matcher matcher;
        List<String> deflineAttributes = new ArrayList<String>();
        try {
            /* parsing splicing attributes */
            pattern = Pattern.compile("splice_form=\"?(-?[a-zA-Z_0-9])*\"?");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("spliced=[YyNn]");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("(noncanonical_splicing=([a-zA-Z]*\\+[a-zA-Z]*,?)+)");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("product=\"([^\"]*\")");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("gene=\"\\S*\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("gene_synonym=\"\\S*\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }

            /* parsing ribosomal slippage attributes */
            pattern = Pattern.compile("V4_Ribosomal_Slippage=\\\"(-?\\+?\\d*)/(-?\\+?\\d*)/(\\S*)\\\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("V4_stop_codon_readthrough=\\\"(-?\\+?\\d*)/[A-Z]/(\\S*)\\\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("alternate_startcodon=\"[a-zA-Z,]*\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }

            /* parsing rna_editing attribute */
            pattern = Pattern.compile("V4_rna_editing=\\\"(-?\\+?\\d*)/([A-Z]*)/(\\S*)/(.*?)\\\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }

            /* parsing Polyprotein mature peptide attributes */
            pattern = Pattern.compile("matpepdb=\"\\S*\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }

            /* Parsing other structural tags */
            pattern = Pattern.compile("shared_cds=\"\\S*\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("\\bis_optional\\b");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("\\bis_required\\b");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("excludes_gene=\"[a-zA-Z,]*\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("tiny_exon3=\"\\w*(:\\w*)?\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("tiny_exon5=\"\\w*(:\\w*)?\"");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
            pattern = Pattern.compile("min_functional_len=\\d*");
            matcher = pattern.matcher(defline);
            if (matcher.find()) {
                deflineAttributes.add(matcher.group(0));
            }
        } catch (Exception e) {
            String message = String.format("Problem parsing defline %s", defline);
            LOGGER.error(message, e);
            throw new VigorException(message);
        }
        return deflineAttributes;
    }
}
