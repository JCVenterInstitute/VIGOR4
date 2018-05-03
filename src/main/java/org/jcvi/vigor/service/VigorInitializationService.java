package org.jcvi.vigor.service;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.apache.logging.log4j.core.tools.picocli.CommandLine;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.*;
import org.springframework.stereotype.Service;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
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

	public List<VigorConfiguration> getDefaultConfigurations() throws VigorException {
		List<VigorConfiguration> configurations = new ArrayList<>();
		VigorConfiguration defaultConfiguration = LoadDefaultParameters
				.loadVigorConfiguration("defaults",Thread.currentThread().getContextClassLoader().getResource(VigorUtils.getVigorParametersPath()));

		configurations.add(defaultConfiguration);

		VigorConfiguration propConfiguration = new VigorConfiguration("system-properties");
		String val;
		String key;
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			key = "vigor."+ param.configKey;
			val = System.getProperty(key);
			if (val != null) {
				if (param.hasFlag(ConfigurationParameters.Flags.VERSION_4)) {
					propConfiguration.put(param, val);
				} else {
					LOGGER.debug("Ignoring non VIGOR4 parameter {}={} set via system properties", param.configKey, val, key);
				}
			}
		}
		configurations.add(propConfiguration);


		VigorConfiguration envConfiguration = new VigorConfiguration("environment");
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			key = "VIGOR_" + param.configKey.toUpperCase();
			val = System.getenv(key);
			if (val != null) {
				if (param.hasFlag(ConfigurationParameters.Flags.VERSION_4)) {
					envConfiguration.put(param, val);
				} else {
					LOGGER.debug("ignoring non-VIGOR4 parameter {}={} set via environment variable {}", param.configKey, val, key);
				}
			}
		}
		configurations.add(envConfiguration);
		return configurations;
	}

	public VigorForm loadParameters(Namespace inputs) throws VigorException{

		List<VigorConfiguration> configurations = getDefaultConfigurations();
		configurations.addAll(getCommandLineConfiguration(inputs));

		String reference_db_dir = getConfigValue(ConfigurationParameters.ReferenceDatabasePath, configurations);

		String reference_db= inputs.getString(CommandLineParameters.referenceDB);
		if ("any".equals(reference_db)) {
			throw new VigorException("Auto-selecting reference database is not implemented");
		}else{
			File file = new File(reference_db);
			if(file.exists() && file.isFile() ){
				reference_db=file.getAbsolutePath();
				if (reference_db_dir == null || reference_db_dir.isEmpty()) {
					reference_db_dir = file.getParent();
				}
			}else if ( ! (reference_db_dir == null || reference_db_dir.isEmpty())){
				reference_db=Paths.get(reference_db_dir,reference_db).toString();
			}
		}

		if (reference_db_dir == null) {
			throw new VigorException("Reference database path is required");
		}

		LOGGER.debug("Reference_db is {}", reference_db);
		if ( reference_db == null || reference_db.isEmpty()) {
			throw new VigorException("reference database is required");
		}

		String virusSpecificConfig = getConfigValue(ConfigurationParameters.VirusSpecificConfiguration, configurations);
		String virusSpecificConfigPath = getConfigValue(ConfigurationParameters.VirusSpecificConfigurationPath, configurations);
		VigorConfiguration defaultConfiguration = configurations.get(0);
		defaultConfiguration = loadVirusSpecificParameters(defaultConfiguration, reference_db, virusSpecificConfigPath, virusSpecificConfig);
		configurations.set(0, defaultConfiguration);

		defaultConfiguration = mergeConfigurations(configurations);

		AlignmentEvidence alignmentEvidence = new AlignmentEvidence();
		VigorForm form = new VigorForm(defaultConfiguration);

		form.setConfiguration(defaultConfiguration);
		form.setAlignmentEvidence(alignmentEvidence);
		alignmentEvidence.setReference_db(reference_db);

		return form;
	}


	public VigorConfiguration mergeConfigurations(List<VigorConfiguration> configurations) {
		// now setup all the defaults

		VigorConfiguration previousConfig = null;

		for (VigorConfiguration config: configurations) {
			if (previousConfig == null) {
				previousConfig = config;
				continue;
			}
			config.setDefaults(previousConfig);
			previousConfig = config;
		}
		return previousConfig;
	}

	public List<VigorConfiguration> getCommandLineConfiguration(Namespace inputs) throws VigorException {

		List<VigorConfiguration> configurations = new ArrayList<>();

		if (inputs.get("config_file") != null) {
			VigorConfiguration configFileConfiguration = LoadDefaultParameters.loadVigorConfiguration("config-file",
					new File(inputs.getString("config_file")));
			configurations.add(configFileConfiguration);
		}

		VigorConfiguration commandLineConfig = new VigorConfiguration("commandline");

		String outputPath = inputs.getString(CommandLineParameters.outputPrefix);
		File outputFile= new File(outputPath);
		if( ! (outputFile.getParentFile().exists() &&outputFile.getParentFile().isDirectory()) ){
			throw new VigorException(String.format("Invalid output prefix %s. Please provide an existing output directory followed by a file prefix", outputPath));
		}
		commandLineConfig.put(ConfigurationParameters.OutputPrefix,outputFile.getName());
		commandLineConfig.put(ConfigurationParameters.OutputDirectory,outputFile.getParentFile().getAbsolutePath());

		commandLineConfig.put(ConfigurationParameters.OverwriteOutputFiles,
				inputs.getBoolean(CommandLineParameters.overwriteOutputFiles) ? "true": "false");

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

		String locus_tag = inputs.getString(CommandLineParameters.locusTag);
		if (locus_tag != null) {
			commandLineConfig.put(ConfigurationParameters.Locustag, locus_tag);
		}

		Boolean ignore_reference_requirements = inputs.getBoolean(CommandLineParameters.ignoreReferenceRequirements);
		if (ignore_reference_requirements != null && ignore_reference_requirements) {
			LOGGER.debug("Ignoring legacy parameter ignore reference requirement");
		}

		String evalue = inputs.getString(CommandLineParameters.eValue);
		if (evalue != null) {
			LOGGER.debug("Ignoring legacy parameter evalue");
		}

		Boolean jcvi_rules = inputs.getBoolean(CommandLineParameters.jcviRules);
		if (jcvi_rules != null) {
			LOGGER.debug("Ignoring legacy parameter evalue");
		}

		String reference_database_path = inputs.getString(CommandLineParameters.referenceDB_Path);
		if (reference_database_path != null) {
			commandLineConfig.put(ConfigurationParameters.ReferenceDatabasePath, reference_database_path);
		}

		String virus_specific_path = inputs.getString(CommandLineParameters.virusSpecificConfigPath);
		if (virus_specific_path != null) {
			commandLineConfig.put(ConfigurationParameters.VirusSpecificConfigurationPath, virus_specific_path);
		}

		virus_specific_path = inputs.getString(CommandLineParameters.virusSpecificConfig);
		if (virus_specific_path != null) {
			commandLineConfig.put(ConfigurationParameters.VirusSpecificConfiguration, virus_specific_path);
		}

		configurations.add(commandLineConfig);


		List<String> parameters = inputs.getList(CommandLineParameters.parameters);
		if (parameters != null) {
			final Pattern splitter = Pattern.compile("~~");
			Map<String, String> temp = parameters.stream()
												 .flatMap(p -> splitter.splitAsStream(p.trim()))
												 .map(s -> s.split("=", 2))
												 .collect(Collectors.toMap(a -> a[0], a -> a.length > 1 ? a[1] : ""));
			VigorConfiguration commandLineParametersConfig = LoadDefaultParameters.configurationFromMap("commandline-parameters", temp);
			configurations.add(commandLineParametersConfig);
		}

		return configurations;
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
	public VigorConfiguration loadVirusSpecificParameters(VigorConfiguration vigorConfiguration, String reference_db, String virusSpecificConfigPath, String virusSpecificConfig) throws VigorException {
		String configPath = virusSpecificConfig;
		if (configPath == null) {
			configPath = Paths.get(virusSpecificConfigPath, reference_db + ".ini").toString();
		}
		File configFile = new File(configPath);
		// config file was specified, but doesn't exist
		if (! configFile.exists()  && virusSpecificConfig != null) {
			throw new VigorException("Virus specific config file {} doesn't exist or is not readable");
		}
		// virus specific configuration files may not exist
		if (configFile.exists()) {
			VigorConfiguration virusSpecificParameters = LoadDefaultParameters.loadVigorConfiguration(configPath, new File(configPath));
			virusSpecificParameters.setDefaults(vigorConfiguration);
			return virusSpecificParameters;
		}
		return vigorConfiguration;
	}

	public void initiateReportFile(String outputDir, String outputPrefix ){
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        FileAppender fa = FileAppender.newBuilder().withName("mylogger").withAppend(false).withFileName(new File(outputDir, outputPrefix+".rpt").toString())
                .build();
        fa.start();
        lc.getConfiguration().addAppender(fa);
        lc.getLogger("org.jcvi.vigor").addAppender(lc.getConfiguration().getAppender(fa.getName()));
        lc.updateLoggers();
    }

    private String getConfigValue(ConfigurationParameters param, List<VigorConfiguration> configurations) {
		String value = null;
		for (int i = configurations.size() - 1; i >= 0; i--) {
			value = configurations.get(i).get(param);
			if (value != null) {
				break;
			}
		}
		return value;
	}
}
