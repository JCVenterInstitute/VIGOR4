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
import org.jcvi.vigor.utils.*;
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
            String outputDir = form.getConfiguration().get(ConfigurationParameters.OutputDirectory);
            String outputPrefix = form.getConfiguration().get(ConfigurationParameters.OutputPrefix);
            initiateReportFile(outputDir,outputPrefix);
            form.getConfiguration().put(ConfigurationParameters.CircularGene, isCircular ? "1" : "0");
            form.getConfiguration().put(ConfigurationParameters.CompleteGene, isComplete ? "1" : "0");
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

		VigorConfiguration propConfiguration = new VigorConfiguration("system-properties");
		String val;
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			val = System.getProperty("vigor."+ param.configKey);
			if (val != null) {
				propConfiguration.put(param, val);
			}
		}

		VigorConfiguration envConfiguration = new VigorConfiguration("environment");
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			val = System.getenv("VIGOR_" + param.configKey.toUpperCase());
			if (val != null) {
				envConfiguration.put(param, val);
			}
		}

		VigorConfiguration defaultConfiguration = LoadDefaultParameters
				.loadVigorConfiguration("defaults",Thread.currentThread().getContextClassLoader().getResource(VigorUtils.getVigorParametersPath()));
		VigorConfiguration commandLineConfig = new VigorConfiguration("commandline");

		AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
		VigorForm form = new VigorForm(defaultConfiguration);
		String reference_db_dir= defaultConfiguration.get(ConfigurationParameters.ReferenceDatabasePath);
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

		defaultConfiguration = loadVirusSpecificParameters(defaultConfiguration, reference_db);

		String outputPath = inputs.getString(CommandLineParameters.outputPrefix);
		File outputFile= new File(outputPath);
		if(outputFile.getParentFile().exists() &&outputFile.getParentFile().isDirectory()){

        }else{
		    throw new VigorException("Invalid -o parameter.Please provide valid output directory followed by prefix");
        }
		commandLineConfig.put(ConfigurationParameters.OutputPrefix,outputFile.getName());
        commandLineConfig.put(ConfigurationParameters.OutputDirectory,outputFile.getParentFile().getAbsolutePath());

		Integer min_gene_size = inputs.getInt(CommandLineParameters.minGeneSize);
		if (min_gene_size != null) {
			commandLineConfig.put(ConfigurationParameters.GeneMinimumSize, min_gene_size.toString());
		}

		String min_gene_coverage = inputs.getString(CommandLineParameters.minCoverage);
		if (min_gene_coverage != null ) {
			commandLineConfig.put(ConfigurationParameters.GeneMinimumCoverage, min_gene_coverage);
		}
		String frameshift_sensitivity = inputs.getString(CommandLineParameters.frameshiftSensitivity);
		if (frameshift_sensitivity != null ) {
			commandLineConfig.put(ConfigurationParameters.FrameShiftSensitivity, frameshift_sensitivity);
		}
		String candidate_selection = inputs.getString(CommandLineParameters.skipSelection);
		if (candidate_selection != null ) {
			commandLineConfig.put(ConfigurationParameters.CandidateSelection, candidate_selection);
		}
		Boolean use_locus_tags = inputs.getBoolean(CommandLineParameters.useLocusTags);
		if (use_locus_tags != null) {
			commandLineConfig.put(ConfigurationParameters.UseLocustags, use_locus_tags ? "1": "0");
		}
		Boolean ignore_reference_requirements = inputs.getBoolean(CommandLineParameters.ignoreReferenceRequirements);
		if (ignore_reference_requirements != null && ignore_reference_requirements) {
			commandLineConfig.put(ConfigurationParameters.CandidateMinimumSimilarity, "0");
			commandLineConfig.put(ConfigurationParameters.CandidateMinimumSubjectCoverage, "0");
			commandLineConfig.put(ConfigurationParameters.MaturePeptideMinimumCoverage, "0");
			commandLineConfig.put(ConfigurationParameters.MaturePeptideMinimumSimilarity, "0");
			commandLineConfig.put(ConfigurationParameters.MaturePeptideMinimumIdentity, "0");
			commandLineConfig.put(ConfigurationParameters.PseudoGeneMinimumIdentity, "0");
			commandLineConfig.put(ConfigurationParameters.PseudoGeneMinimumSimilarity, "0");
			commandLineConfig.put(ConfigurationParameters.PseudoGeneMinimumCoverage, "0");
		}

		String evalue = inputs.getString(CommandLineParameters.eValue);
		if (evalue != null) {
			commandLineConfig.put(ConfigurationParameters.CandidateEvalue, evalue);
		}
		Boolean jcvi_rules = inputs.getBoolean(CommandLineParameters.jcviRules);
		if (jcvi_rules != null) {
			commandLineConfig.put(ConfigurationParameters.JCVIRules, jcvi_rules ? "1": "0");
		}
		List<String> parameters = inputs.getList(CommandLineParameters.parameters);
		if (parameters != null) {
			final Pattern splitter = Pattern.compile("~~");
			Map<String, String> temp = parameters.stream()
												 .flatMap(p -> splitter.splitAsStream(p.trim()))
												 .map(s -> s.split("=", 2))
												 .collect(Collectors.toMap(a -> a[0], a -> a.length > 1 ? a[1] : ""));
			VigorConfiguration commandLineParametersConfig = LoadDefaultParameters.configurationFromMap("commandline-parameters",temp);

			// command line config overrides loaded configuration
			envConfiguration.setDefaults(defaultConfiguration);
			propConfiguration.setDefaults(envConfiguration);
			commandLineParametersConfig.setDefaults(propConfiguration);
			commandLineConfig.setDefaults(commandLineParametersConfig);
			defaultConfiguration = commandLineConfig;
		}
		form.setConfiguration(defaultConfiguration);
		form.setAlignmentEvidence(alignmentEvidence);
		return form;
	}

	/**
	 *
	 * @param vigorConfiguration:
	 *            Default vigor parameters from vigor.ini file
	 * @param reference_db
	 *            : Virus Specific reference_db
	 * @return configuration : Default vigor Parameters will be overridden
	 *         by virus specific parameters
	 */
	public VigorConfiguration loadVirusSpecificParameters(VigorConfiguration vigorConfiguration, String reference_db) throws VigorException {
	    String virusSpecificParametersPath = vigorConfiguration.get(ConfigurationParameters.VirusSpecificConfigurationPath);
		VigorConfiguration virusSpecificParameters = LoadDefaultParameters.loadVigorConfiguration(reference_db + " specific config",
				Thread.currentThread().getContextClassLoader().getResource(
				Paths.get(virusSpecificParametersPath , reference_db + ".ini").toString()));

		virusSpecificParameters.setDefaults(vigorConfiguration);
		return virusSpecificParameters;
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
