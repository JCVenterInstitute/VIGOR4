package org.jcvi.vigor.service;
import org.jcvi.vigor.component.*;
import org.jcvi.vigor.utils.VigorFunctionalUtils;
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
	 *
	 * @param alignment
	 * @return viralProtein: For the given protein ID ViralProtein object is
	 *         generated and all the properties are defined;
	 */
	public Alignment setViralProteinAttributes(Alignment alignment,VigorForm form) {
		String min_intronSize_param = form.getVigorParametersList().get("min_intron_size");
		if(VigorUtils.is_Integer(min_intronSize_param)){
			min_intron_length=Integer.parseInt(min_intronSize_param);
		}		
		ViralProtein viralProtein = setGeneAttributes(alignment.getViralProtein());
		AlignmentEvidence alignmentEvidence = alignment.getAlignmentEvidence();
		alignmentEvidence.setMatpep_db(matPepDB);
		alignment.setAlignmentEvidence(alignmentEvidence);
		alignment.setViralProtein(viralProtein);
		return alignment;
	}

	/**
	 *
	 * @param viralProtein
	 * @return viralProtein object which has all the gene attributes set.
	 */

	public ViralProtein setGeneAttributes(ViralProtein viralProtein) {

		try {

			GeneAttributes geneAttributes = new GeneAttributes();
			String defline = viralProtein.getDefline();
			defline = StringUtils.normalizeSpace(defline);
			List<String> deflineAttributes = parseDeflineAttributes(defline);
			Map<String, String> attributes = deflineAttributes.stream().map(s -> s.split("=", 2))
					.collect(Collectors.toMap(a -> a[0], a -> a.length > 1 ? a[1] : "", (s1, s2) -> s1));

			StructuralSpecifications structuralSpecifications = new StructuralSpecifications();
			Splicing splicing = new Splicing();
			Ribosomal_Slippage ribosomal_slippage = new Ribosomal_Slippage();
			StopTranslationException stopTranslationException = new StopTranslationException();
			StartTranslationException startTranslationException = new StartTranslationException();
			RNA_Editing rna_editing = new RNA_Editing();

			/* Set Splicing attributes */
			if (attributes.containsKey("splice_form")) {
				if (!(attributes.get("splice_form").equals(""))) {
					splicing.setSpliced(true);
					splicing.setSpliceform((attributes.get("splice_form")).replaceAll("^\"|\"$", ""));
					List<SpliceSite> splicePairs = new ArrayList<SpliceSite>();
					if (attributes.containsKey("noncanonical_splicing")) {
						if (!(attributes.get("noncanonical_splicing").equalsIgnoreCase("N"))) {
							
								List<String> spliceSites = Pattern.compile(";").splitAsStream(attributes.get("noncanonical_splicing").trim()).collect(Collectors.toList());

								for(String spliceSite : spliceSites){
									String[] temp = spliceSite.split("+");
									SpliceSite spliceSiteObj = splicing.new SpliceSite();
									spliceSiteObj.donor=temp[0];
									spliceSiteObj.acceptor=temp[1];
									splicePairs.add(spliceSiteObj);
									
								}
						}
					}
					    SpliceSite defaultSpliceSite = splicing.new SpliceSite();
					    defaultSpliceSite.donor="GT";
					    defaultSpliceSite.acceptor="AG";
					    splicePairs.add(defaultSpliceSite);
						splicing.setNonCanonical_spliceSites(splicePairs);			
						
					
					if (attributes.containsKey("tiny_exon3")) {
						String attribute = attributes.get("tiny_exon3").replaceAll("^\"|\"$", "");
						String regex = null;
						int offset = 0;
						String[] temp = attribute.split(":");
						regex = temp[0];
						if (VigorUtils.is_Integer(temp[1])) {
							offset = Integer.parseInt(temp[1]);

						}
						Map<String, Integer> temp1 = new HashMap<String, Integer>();
						temp1.put(regex, offset);
						structuralSpecifications.setTiny_exon3(temp1);
					}
					if (attributes.containsKey("tiny_exon5")) {
						String attribute = attributes.get("tiny_exon5").replaceAll("^\"|\"$", "");
						String regex = null;
						int offset = 0;
						if (attribute.matches(".*?:.*")) {
							String[] temp = attribute.split(":");
							regex = temp[0];
							if (VigorUtils.is_Integer(temp[1])) {
								offset = Integer.parseInt(temp[1]);
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
			if (attributes.containsKey("ribosomal_slippage")) {

				if (attributes.get("ribosomal_slippage").equalsIgnoreCase("Y")) {
					ribosomal_slippage.setHas_ribosomal_slippage(true);
					if (attributes.containsKey("slippage_frameshift")) {
						if (VigorUtils.is_Integer(attributes.get("slippage_frameshift"))) {
							int fs = Integer.parseInt(attributes.get("slippage_frameshift"));
							if (fs >= -3 && fs <= 3)
								ribosomal_slippage.setSlippage_frameshift(fs);
						}
					}
					if (attributes.containsKey("slippage_offset")) {
						if (VigorUtils.is_Integer(attributes.get("slippage_offset")))
							ribosomal_slippage.setSlippage_offset(
									Integer.parseInt(attributes.get("slippage_offset").replaceAll("^\"|\"$", "")));
					}
					if (attributes.containsKey("slippage_motif"))
						ribosomal_slippage
								.setSlippage_motif(attributes.get("slippage_motif").replaceAll("^\"|\"$", ""));
				}
			}

			/* Set StopTranslationException attributes */
			if (attributes.containsKey("stop_codon_readthru")) {
				String attribute = attributes.get("stop_codon_readthru").replaceAll("^\"|\"$", "");
				stopTranslationException.setHasStopTranslationException(true);
				String[] temp = attribute.split(":");
				AminoAcid replacementAA = AminoAcid.valueOf(temp[1]);
				stopTranslationException.setReplacementAA(replacementAA);
			}

			/* Set StartTranslationException attributes */
			if (attributes.containsKey("alternate_startcodon")) {
				String attribute = attributes.get("alternate_startcodon").replaceAll("^\"|\"$", "");
				startTranslationException.setAlternateStartCodons(Arrays.asList(attribute.split(",")));
				startTranslationException.setHasStartTranslationException(true);

			}

			/* Set RNA_Editing attributes */
			if (attributes.containsKey("rna_editing")) {
				String attribute = attributes.get("rna_editing").replaceAll("^\"|\"$", "");
				String[] temp = attribute.split("/");
				if (VigorUtils.is_Integer(temp[0])) {
					rna_editing.setSize(Integer.parseInt(temp[0]));
				}
				rna_editing.setRegExp(temp[1]);
				rna_editing.setHas_RNA_editing(true);
				String insertionString = temp[2].substring(2);
				rna_editing.setInsertionString(insertionString);
				rna_editing.setNote(temp[3]);
			}

			/* Set StructuralSpecifications */
			if (attributes.containsKey("shared_cds")) {
				String attribute = attributes.get("shared_cds").replaceAll("^\"|\"$", "");
				structuralSpecifications.setShared_cds(Arrays.asList(attribute.split(",")));

			}
			if (attributes.containsKey("is_optional ")) {
				structuralSpecifications.set_required(true);
			}
			if (attributes.containsKey("is_required ")) {
				structuralSpecifications.set_required(false);
			}
			if (attributes.containsKey("excludes_gene")) {
				String attribute = attributes.get("excludes_gene").replaceAll("^\"|\"$", "");
				structuralSpecifications.setExcludes_gene(Arrays.asList(attribute.split(",")));
			}
			if (attributes.containsKey("min_functional_len")) {
				String attribute = attributes.get("min_functional_len");
				if (VigorUtils.is_Integer(attribute)) {
					structuralSpecifications.setMinFunctionalLength(Integer.parseInt(attribute));
				}
			}

			/* Set maturepeptide DB attribute */
			if (attributes.containsKey("matpepdb")) {
				matPepDB = attributes.get("matpepdb");
			}

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

		} catch (Exception e) {
			LOGGER.error(e.getMessage(), e);
		}
		return viralProtein;

	}

	/**
	 *
	 * @param viralProtein
	 * @return GeneStructure: has list of exons and introns of the viralProtein.
	 *         These are determined from spliceform annotated in the defline
	 */
	public ViralProtein DetermineGeneStructure(ViralProtein viralProtein) {
        
		boolean is_ribosomal_slippage = viralProtein.getGeneAttributes().getRibosomal_slippage()
				.isHas_ribosomal_slippage();
		boolean is_spliced = viralProtein.getGeneAttributes().getSplicing().isSpliced();
		List<Range> NTFragments = new ArrayList<Range>();
		List<Range> introns = new ArrayList<Range>();
		if (!(is_ribosomal_slippage) && !(is_spliced)) {
			Range range = Range.of(0, 3 * (viralProtein.getSequence().getLength()-1));
			NTFragments.add(range);
			
		} else {
			String spliceform = viralProtein.getGeneAttributes().getSplicing().getSpliceform();
			if (spliceform != null) {
				if (spliceform.equals("") || !(spliceform.matches("([i,e]-?[0-9]*)+"))) {
					String exception = "Spliced reference missing/malformed splice_form tag:"
							+ viralProtein.getProteinID();
					System.out.println(exception);
					LOGGER.debug(exception);
				} else {

					List<String> splices = new ArrayList<String>();
					Matcher m = Pattern.compile("[e,i]-?[0-9]*").matcher(spliceform);
					while (m.find()) {
						splices.add(m.group(0));
					}
					Long dnaOffset = 0l;
					for (int i = 0; i < splices.size(); i++) {
						long nucleotides = Long.parseLong(splices.get(i).substring(1));
						if(nucleotides>0){
						Range range = Range.of(dnaOffset, (dnaOffset + nucleotides) - 1);
						if (splices.get(i).matches("e-?[0-9]*")) {
							dnaOffset = range.getEnd() + 1;
							NTFragments.add(range);
						} else {
							if(nucleotides>=min_intron_length){
								introns.add(range);
							}
							dnaOffset = range.getEnd() + 1;
						}}
						else{
							dnaOffset=dnaOffset-nucleotides;
						}
					}
				}

			} else {

				String exception = "Spliced reference missing for " + viralProtein.getProteinID();
				System.out.println(exception);
				LOGGER.debug(exception);
			}
		}
        viralProtein.setNTfragments(NTFragments);
        viralProtein.setIntrons(introns);
		return viralProtein;
	}

	/**
	 * @param defline
	 *            of the protein
	 * @return List of attributes in the defline
	 */

	public List<String> parseDeflineAttributes(String defline) {
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

			/* parsing ribosomal slippage attributes */
			pattern = Pattern.compile("(ribosomal_slippage=[Yy])");
			matcher = pattern.matcher(defline);
			if (matcher.find()) {
				deflineAttributes.add(matcher.group(0));
			}
			pattern = Pattern.compile("(slippage_motif=\"\\S*\")");
			matcher = pattern.matcher(defline);
			if (matcher.find()) {
				deflineAttributes.add(matcher.group(0));
			}
			pattern = Pattern.compile("(slippage_offset=(-?\\d*))");
			matcher = pattern.matcher(defline);
			if (matcher.find()) {
				deflineAttributes.add(matcher.group(0));
			}
			pattern = Pattern.compile("(slippage_frameshift=(-?\\d*))");
			matcher = pattern.matcher(defline);
			if (matcher.find()) {
				deflineAttributes.add(matcher.group(0));
			}

			/* Parsing Translation exception attributes */
			pattern = Pattern.compile("(stop_codon_readthru=[Yy](:[a-zA-Z])?)");
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
			pattern = Pattern.compile("rna_editing=\\d*/\\(.*\\)/\\$1[a-zA-Z]*/.*?/");
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
			pattern = Pattern.compile("shared_cds=\"[a-zA-Z,]*\"");
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
			LOGGER.error(e.getMessage(), e);
		}
		return deflineAttributes;
	}
}
