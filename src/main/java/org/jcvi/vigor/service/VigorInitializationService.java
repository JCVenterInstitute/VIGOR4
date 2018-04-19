package org.jcvi.vigor.service;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.LoadDefaultParameters;
import org.jcvi.vigor.utils.VigorLogging;
import org.jcvi.vigor.utils.VigorUtils;
import org.springframework.stereotype.Service;

import java.io.File;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * This class is under service layer. The methods in this class has
 * functionality to initialize the vigor application
 */
@Service
public class VigorInitializationService {

	private static final Logger LOGGER = LogManager.getLogger(VigorInitializationService.class);


	/**
	 * @param inputs:
	 *            User provided command line inputs Retrieve each Genomic
	 *            sequence from the input file and determine AlignmentEvidence
	 */

	public VigorForm initializeVigor(Namespace inputs) throws VigorException {

            boolean isComplete = false;
            boolean isCircular = false;
            Boolean complete_gene = inputs.getBoolean(CommandLineParameters.completeGene);
            if (complete_gene != null && complete_gene) {
                isComplete = true;
            }
            Boolean circular_gene = inputs.getBoolean(CommandLineParameters.circularGene);
            if (circular_gene != null && circular_gene) {
                isComplete = true;
                isCircular = true;
            }
            VigorForm form = loadParameters(inputs);
            String outputDir = form.getVigorParametersList().get("output_directory");
            String outputPrefix = form.getVigorParametersList().get("output_prefix");
            initiateReportFile(outputDir,outputPrefix);
            form.getVigorParametersList().put("circular_gene", isCircular ? "1" : "0");
            form.getVigorParametersList().put("complete_gene", isComplete ? "1" : "0");
            return form;

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
	public VigorForm loadParameters(Namespace inputs) throws VigorException{

		Map<String, String> vigorParameterList = LoadDefaultParameters
				.loadVigorParameters(VigorUtils.getVigorParametersPath());
		AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
		VigorForm form = new VigorForm();
		String reference_db_dir=vigorParameterList.get("reference_db_path");
		String reference_db= inputs.getString(CommandLineParameters.referenceDB);
		if ("any".equals(reference_db)) {
			//TODO initiate referenceDB generation service
			LOGGER.debug("autoselecting reference database from {}. setting value to {}", reference_db_dir, "TODO");
		}else{
		    File file = new File(reference_db);
		    if(file.exists() && file.isFile() ){
		        reference_db=file.getAbsolutePath();
            }else{
                reference_db=Paths.get(reference_db_dir,reference_db).toString();
            }
        }

		LOGGER.debug("Reference_db is {}", reference_db);
		alignmentEvidence.setReference_db(reference_db);

		vigorParameterList = loadVirusSpecificParameters(vigorParameterList, reference_db);

		String outputPath = inputs.getString(CommandLineParameters.outputPrefix);
		File outputFile= new File(outputPath);
		if(outputFile.getParentFile().exists() &&outputFile.getParentFile().isDirectory()){

        }else{
		    throw new VigorException("Invalid -o parameter.Please provide valid output directory followed by prefix");
        }
		vigorParameterList.put("output_prefix",outputFile.getName());
        vigorParameterList.put("output_directory",outputFile.getParentFile().getAbsolutePath());

		Integer min_gene_size = inputs.getInt(CommandLineParameters.minGeneSize);
		if (min_gene_size != null) {
			vigorParameterList.put("min_gene_size", min_gene_size.toString());
		}

		String min_gene_coverage = inputs.getString(CommandLineParameters.minCoverage);
		if (min_gene_coverage != null ) {
			vigorParameterList.put("min_gene_coverage", min_gene_coverage);
		}
		String frameshift_sensitivity = inputs.getString(CommandLineParameters.frameshiftSensitivity);
		if (frameshift_sensitivity != null ) {
			vigorParameterList.put("frameshift_sensitivity", frameshift_sensitivity);
		}
		String candidate_selection = inputs.getString(CommandLineParameters.skipSelection);
		if (candidate_selection != null ) {
			vigorParameterList.put("candidate_selection", candidate_selection);
		}
		Boolean use_locus_tags = inputs.getBoolean(CommandLineParameters.useLocusTags);
		if (use_locus_tags != null) {
			vigorParameterList.put("use_locus_tags", use_locus_tags ? "1": "0");
		}
		Boolean ignore_reference_requirements = inputs.getBoolean(CommandLineParameters.ignoreReferenceRequirements);
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

		String evalue = inputs.getString(CommandLineParameters.eValue);
		if (evalue != null) {
			vigorParameterList.put("candidate_evalue", evalue);
		}
		Boolean jcvi_rules = inputs.getBoolean(CommandLineParameters.jcviRules);
		if (jcvi_rules != null) {
			vigorParameterList.put("jcvi_rules", jcvi_rules ? "1": "0");
		}
		List<String> parameters = inputs.getList(CommandLineParameters.parameters);
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
	    String virusSpecificParametersPath = vigorParametersList.get("virusSpecific_parameters");
		Map<String, String> virusSpecificParameters = LoadDefaultParameters.loadVigorParameters(
				virusSpecificParametersPath + File.separator + reference_db + ".ini");
		vigorParametersList.putAll(virusSpecificParameters);
		return vigorParametersList;
	}

	public void initiateReportFile(String outputDir, String outputPrefix ){
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        FileAppender fa = FileAppender.newBuilder().withName("mylogger").withAppend(false).withFileName(new File(outputDir, outputPrefix+".rpt").toString())
                .withLayout(PatternLayout.newBuilder().withPattern("%-5p %d  [%t] %C{2} (%F:%L) - %m%n").build())
                .setConfiguration(lc.getConfiguration()).build();
        fa.start();
        lc.getConfiguration().addAppender(fa);
        lc.getRootLogger().addAppender(lc.getConfiguration().getAppender(fa.getName()));
        lc.updateLoggers();
    }
}
