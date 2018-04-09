package org.jcvi.vigor.service;

import java.io.File;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Iterator;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.sourceforge.argparse4j.inf.Namespace;
import org.jcvi.vigor.utils.VigorLogging;
import org.jcvi.vigor.utils.VigorUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.fasta.nt.NucleotideFastaDataStore;
import org.jcvi.jillion.fasta.nt.NucleotideFastaFileDataStoreBuilder;
import org.jcvi.jillion.fasta.nt.NucleotideFastaRecord;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.component.VirusGenome;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.LoadDefaultParameters;

/**
 * This class is under service layer. The methods in this class has
 * functionality to initialize the vigor application
 */
@Service
public class VigorInitializationService {

	private static final Logger LOGGER = LogManager.getLogger(VigorInitializationService.class);

	@Autowired
	private AlignmentGenerationService alignmentGenerationService;

	/**
	 * @param inputs:
	 *            User provided command line inputs Retrieve each Genomic
	 *            sequence from the input file and determine AlignmentEvidence
	 */

	public void initializeVigor(Namespace inputs) {
		try {
			boolean isComplete = false;
			boolean isCircular = false;
			Boolean complete_gene = inputs.getBoolean("complete_gene");
			if (complete_gene != null && complete_gene) {
				isComplete = true;
			}
			Boolean circular_gene = inputs.getBoolean("circular_gene");
			if (circular_gene != null && circular_gene) {
				isComplete = true;
				isCircular = true;
			}
			VigorForm form = loadParameters(inputs);
			String inputFileName = inputs.getString("input_fasta");
			// TODO check
			NucleotideFastaDataStore dataStore = new NucleotideFastaFileDataStoreBuilder(
					new File(inputFileName)).hint(DataStoreProviderHint.RANDOM_ACCESS_OPTIMIZE_SPEED)
							.build();
			Stream<NucleotideFastaRecord> records = dataStore.records();
			Iterator<NucleotideFastaRecord> i = records.iterator();
			while (i.hasNext()) {
				NucleotideFastaRecord record = i.next();
				VirusGenome virusGenome = new VirusGenome(record.getSequence(), record.getComment(), record.getId(),
						isComplete, isCircular);
				// Call referenceDBGenerationService methods to generate alignmentEvidence.
				alignmentGenerationService.GenerateAlignment(virusGenome, form);
                
			}

		} catch (Exception e) {
			e.printStackTrace();
			LOGGER.error(e.getMessage(), e);

		}

	}

	/**
	 * load all the vigor parameters from Vigor.ini file
	 *
	 * @param inputs:
	 *            Input parameters and values provided by user
	 * @return form : output form object has the AlignmentEvidence object and
	 *         the VigorParametersList. Few vigor parameters will be overridden
	 *         by the default parameters of vigor.ini file and saved to
	 *         VigorParametersList attribute of the form.
	 */
	public VigorForm loadParameters(Namespace inputs) {

		Map<String, String> vigorParameterList = LoadDefaultParameters
				.loadVigorParameters(VigorUtils.getVigorParametersPath());
		AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
		VigorForm form = new VigorForm();
		String reference_db = inputs.getString("reference_database");
		if ("any".equals(reference_db)) {
			reference_db = vigorParameterList.get("reference_db");
			LOGGER.debug("autoselecting reference database; setting value to {}", reference_db);
		}

		LOGGER.debug("Reference_db is " + reference_db);
		// TODO check if reference db is a path, if so use that, else construct path using database path
		alignmentEvidence.setReference_db(reference_db);


		vigorParameterList = loadVirusSpecificParameters(vigorParameterList, alignmentEvidence.getReference_db());

		Integer min_gene_size = inputs.getInt("min_gene_size");
		if (min_gene_size != null) {
			vigorParameterList.put("min_gene_size", min_gene_size.toString());
		}

		String min_gene_coverage = inputs.getString("min_gene_coverage");
		if (min_gene_coverage != null ) {
			vigorParameterList.put("min_gene_coverage", min_gene_coverage);
		}
		String frameshift_sensitivity = inputs.getString("frameshift_sensitity");
		if (frameshift_sensitivity != null ) {
			vigorParameterList.put("frameshift_sensitivity", frameshift_sensitivity);
		}
		String candidate_selection = inputs.getString("skip_selection");
		if (candidate_selection != null ) {
			vigorParameterList.put("candidate_selection", candidate_selection);
		}
		Boolean use_locus_tags = inputs.getBoolean("use_locus_tags");
		if (use_locus_tags != null) {
			vigorParameterList.put("use_locus_tags", use_locus_tags ? "1": "0");
		}
		Boolean ignore_reference_requirements = inputs.getBoolean("ignore_reference_requirements");
		if (ignore_reference_requirements != null && ignore_reference_requirements) {
			vigorParameterList.put("min_candidate_pctsimilarity", "0");
			vigorParameterList.put("min_candidate_sbjcoverage", "0");
			vigorParameterList.put("mature_pep_mincoverage", "0");
			vigorParameterList.put("mature_pep_minsimilarity", "0");
			vigorParameterList.put("mature_pep_minidentity", "0");
			vigorParameterList.put("min_pseudogene_identity", "0");
			vigorParameterList.put("min_pseudogene_similarity", "0");
			vigorParameterList.put("min_pseudogene_coverage", "0");
		}

		String evalue = inputs.getString("evalue");
		if (evalue != null) {
			vigorParameterList.put("candidate_evalue", evalue);
		}
		Boolean jcvi_rules = inputs.getBoolean("jcvi_rules");
		if (jcvi_rules != null) {
			vigorParameterList.put("jcvi_rules", jcvi_rules ? "1": "0");
		}
		List<String> parameters = inputs.getList("parameters");
		if (parameters != null) {
			final Pattern splitter = Pattern.compile("~~");
			Map<String, String> temp = parameters.stream()
												 .flatMap(p -> splitter.splitAsStream(p.trim()))
												 .map(s -> s.split("=", 2))
												 .collect(Collectors.toMap(a -> a[0], a -> a.length > 1 ? a[1] : ""));
			for (String key : temp.keySet()) {
				if (vigorParameterList.containsKey(key)) {
					vigorParameterList.put(key, temp.get(key));
				} else {
					LOGGER.warn(VigorLogging.VIGOR4_USER_MESSAGE, "skipping unknown parameter \"{}\" with value \"{}\"", key, temp.get(key));
				}
			}
		}
		form.setVigorParametersList(vigorParameterList);
		form.setAlignmentEvidence(alignmentEvidence);
		return form;
	}

	/**
	 *
	 * @param vigorParametersList:
	 *            Default vigor parameters from vigor.ini file
	 * @param reference_db
	 *            : Virus Specific reference_db
	 * @return vigorParametersList : Default vigor Parameters will be overridden
	 *         by virus specific parameters
	 */
	public Map<String, String> loadVirusSpecificParameters(Map<String, String> vigorParametersList,
			String reference_db) {
		Map<String, String> virusSpecificParameters = LoadDefaultParameters.loadVigorParameters(
				VigorUtils.getVirusSpecificParametersPath() + File.separator + reference_db + ".ini");
		vigorParametersList.putAll(virusSpecificParameters);
		return vigorParametersList;
	}


}
