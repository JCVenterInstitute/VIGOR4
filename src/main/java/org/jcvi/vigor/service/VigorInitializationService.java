package org.jcvi.vigor.service;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.jcvi.vigor.component.AlignmentEvidence;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.forms.VigorForm;
import org.jcvi.vigor.utils.*;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
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
            VigorConfiguration configuration =  form.getConfiguration();
            String outputDir = configuration.get(ConfigurationParameters.OutputDirectory);
            String outputPrefix = configuration.get(ConfigurationParameters.OutputPrefix);
            initiateReportFile(outputDir,outputPrefix, inputs.getInt(CommandLineParameters.verbose));
            form.getConfiguration().put(ConfigurationParameters.CircularGene, isCircular ? "1" : "0");
            form.getConfiguration().put(ConfigurationParameters.CompleteGene, isComplete ? "1" : "0");
            return form;

	}

	public List<VigorConfiguration> getDefaultConfigurations() throws VigorException {
		List<VigorConfiguration> configurations = new ArrayList<>();
		VigorConfiguration defaultConfiguration = LoadDefaultParameters
				.loadVigorConfiguration("defaults",
										Thread.currentThread().getContextClassLoader().getResource(VigorUtils.getDefaultConfigurationPath()));

		configurations.add(defaultConfiguration);

		VigorConfiguration propConfiguration = new VigorConfiguration("system-properties");
		String val;
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			val = System.getProperty(param.getSystemPropertyName());
			if (val != null) {
				if (param.hasFlag(ConfigurationParameters.Flags.VERSION_4)) {
					propConfiguration.put(param, val);
				} else {
					LOGGER.debug("Ignoring non VIGOR4 parameter {}={} set via system properties", param.configKey, val, param.getSystemPropertyName());
				}
			}
		}
		configurations.add(propConfiguration);


		VigorConfiguration envConfiguration = new VigorConfiguration("environment");
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			val = System.getenv(param.getEnvVarName());
			if (val != null) {
				if (param.hasFlag(ConfigurationParameters.Flags.VERSION_4)) {
					envConfiguration.put(param, val);
				} else {
					LOGGER.debug("ignoring non-VIGOR4 parameter {}={} set via environment variable {}", param.configKey, val, param.getEnvVarName());
				}
			}
		}
		configurations.add(envConfiguration);
		return configurations;
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

		List<VigorConfiguration> configurations = new ArrayList<>();
		configurations.addAll(getDefaultConfigurations());
		configurations.addAll(getCommandLineConfiguration(inputs));

		String reference_db_dir = getConfigValue(ConfigurationParameters.ReferenceDatabasePath, configurations);
		LOGGER.debug("Reference database path is {}", reference_db_dir);

		String reference_db= inputs.getString(CommandLineParameters.referenceDB);
		LOGGER.debug("reference database is {}", reference_db);

		if ("any".equals(reference_db)) {
			throw new VigorException("Auto-selecting reference database is not implemented");
		} else if (inputs.get(CommandLineParameters.genbankDB) != null) {
			throw new VigorException("Using genbank files as the reference database is not implemented");
		} else if (reference_db == null || reference_db.isEmpty()) {
			throw new VigorException("no reference database provided");
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
		File referenceDBFile = new File(reference_db);
		if (! referenceDBFile.exists() ) {
			throw new VigorException(String.format("Reference database file \"%s\" does not exist", reference_db));
		}

		if (! referenceDBFile.isFile()) {
			throw new VigorException(String.format("Reference database \"%s\" is not a file", reference_db));
		}

		if (! referenceDBFile.canRead()) {
			throw new VigorException(String.format("Reference database \"%s\" is not readable", reference_db));
		}

		String virusSpecificConfig = getConfigValue(ConfigurationParameters.VirusSpecificConfiguration, configurations);
		String virusSpecificConfigPath = getConfigValue(ConfigurationParameters.VirusSpecificConfigurationPath, configurations);
		if (virusSpecificConfigPath == null || virusSpecificConfigPath.isEmpty()) {
			virusSpecificConfigPath = reference_db_dir;
		}
		VigorConfiguration defaultConfiguration = configurations.get(0);
		String referenceDBName = Paths.get(reference_db).getFileName().toString();
		defaultConfiguration = loadVirusSpecificParameters(defaultConfiguration, referenceDBName, virusSpecificConfigPath, virusSpecificConfig);
		configurations.set(0, defaultConfiguration);

		defaultConfiguration = mergeConfigurations(configurations);

		String temporaryDirectory = defaultConfiguration.get(ConfigurationParameters.TemporaryDirectory);
		if ( temporaryDirectory == null || temporaryDirectory.isEmpty()) {
			throw new VigorException("temporary directory not set");
		}

		File tempDir = Paths.get(temporaryDirectory).toFile();
		if (tempDir.exists()) {
			if (! tempDir.isDirectory()) {
				throw new VigorException(String.format("temporary directory %s exists but is not a directory", temporaryDirectory));
			}
			if (! (tempDir.canRead() && tempDir.canWrite())) {
				throw new VigorException(String.format("temporary directory %s is not readable or not writable", temporaryDirectory));
			}
		} else {
			if (! tempDir.mkdirs()) {
				throw new VigorException(String.format("unable to create temporary directory %s", tempDir));
			}
		}
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

		String config_file = inputs.get(CommandLineParameters.configFile);
		if (config_file == null || config_file.isEmpty()) {
		    LOGGER.debug("checking environment variable VIGOR_CONFIG_FILE for config file");
		    config_file = System.getenv("VIGOR_CONFIG_FILE");
        }
		if (! (config_file == null || config_file.isEmpty()) ) {
			if (config_file.startsWith("~")) {
				try {
					config_file = VigorUtils.expandTilde(config_file);
				} catch (IOException e) {
					throw new VigorException(String.format("problem expanding ~ in config file %s", config_file), e);
				}
			}
			File config_path = new File(config_file).getAbsoluteFile();
			if (! config_path.exists()) {
				throw new VigorException(String.format("config file %s does not exist", config_path.toString()));
			} else if (! config_path.canRead()) {
				throw new VigorException(String.format("config file %s is not readable", config_path.toString()));
			}
		    LOGGER.debug("loading config file {}", config_file);
		    // use the file as configuration name so it's unambigious
			VigorConfiguration configFileConfiguration = LoadDefaultParameters.loadVigorConfiguration(config_file,config_path);
			configurations.add(configFileConfiguration);
		}

		VigorConfiguration commandLineConfig = new VigorConfiguration("commandline");

		String outputPath = inputs.getString(CommandLineParameters.outputPrefix);
		File outputFile= new File(outputPath);
		if (! (outputFile.getParentFile().exists() || outputFile.getParentFile().mkdirs()) ) {
			throw new VigorException(String.format("unable to create directory %s", outputFile.getParent()));
		}
		if( ! (outputFile.getParentFile().exists() && outputFile.getParentFile().isDirectory()) ){
			throw new VigorException(String.format("Invalid output prefix %s. Please provide a directory followed by a file prefix", outputPath));
		}
		commandLineConfig.put(ConfigurationParameters.OutputPrefix,outputFile.getName());
		commandLineConfig.put(ConfigurationParameters.OutputDirectory,outputFile.getParentFile().getAbsolutePath());

		// defaults to false on commandline, so set it in the commandline configuration if true, so we don't override
		// a value in the vigor.ini file
		if (inputs.getBoolean(CommandLineParameters.overwriteOutputFiles)) {
			commandLineConfig.put(ConfigurationParameters.OverwriteOutputFiles,
								  inputs.getBoolean(CommandLineParameters.overwriteOutputFiles) ? "true" : "false");
		}

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

		String reference_database = inputs.getString(CommandLineParameters.referenceDB);
		if (reference_database != null) {
			Path referenceDBPath = Paths.get(reference_database);
			if (referenceDBPath.isAbsolute()) {
				commandLineConfig.put(ConfigurationParameters.ReferenceDatabasePath, referenceDBPath.getParent().toString());
			}
		}

		String reference_database_path = inputs.getString(CommandLineParameters.referenceDB_Path);
		if (reference_database_path != null) {
			if (! (commandLineConfig.get(ConfigurationParameters.ReferenceDatabasePath) == null ||
					commandLineConfig.get(ConfigurationParameters.ReferenceDatabasePath) == reference_database_path)) {
				throw new VigorException(String.format("conflicting reference database paths db path %s reference db file %s",
						reference_database_path,
						reference_database));
			}
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
        commandLineConfig.put(ConfigurationParameters.Verbose,
                inputs.getInt(CommandLineParameters.verbose) > 0 ? "true": "false");

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
	public VigorConfiguration loadVirusSpecificParameters(VigorConfiguration vigorConfiguration,
														  String reference_db,
														  String virusSpecificConfigPath,
														  String virusSpecificConfig) throws VigorException {
		LOGGER.debug("checking virus specific config for reference db {}, virus specific config {}, virus specific config path {}", reference_db, virusSpecificConfig, virusSpecificConfigPath);
		String configPath = virusSpecificConfig;
		if (configPath == null && virusSpecificConfigPath != null) {
			LOGGER.debug("using default virus specific config path {}", virusSpecificConfigPath);

			configPath = Paths.get(virusSpecificConfigPath, reference_db + ".ini").toString();
		}
		if (configPath == null || configPath.isEmpty()) {
			LOGGER.warn("virus specific config path not specified");
			return vigorConfiguration;
		}

		File configFile = new File(configPath);
		LOGGER.debug("virus specific config file for reference_db {} is {}", reference_db, configPath);
		// config file was specified, but doesn't exist
		if ( (! configFile.exists())  && virusSpecificConfig != null) {
			throw new VigorException("Virus specific config file {} doesn't exist or is not readable");
		}
		// virus specific configuration files may not exist
		if (configFile.exists()) {
			VigorConfiguration virusSpecificParameters = LoadDefaultParameters.loadVigorConfiguration(configPath, configFile);
			virusSpecificParameters.setDefaults(vigorConfiguration);
			return virusSpecificParameters;
		}
		return vigorConfiguration;
	}

	public void initiateReportFile(String outputDir, String outputPrefix, int verbose){
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        FileAppender fa = FileAppender.newBuilder().withName("mylogger").withAppend(false).withFileName(new File(outputDir, outputPrefix+".rpt").toString())
                .build();
        fa.start();
        lc.getConfiguration().addAppender(fa);

		lc.getLogger("org.jcvi.vigor").addAppender(lc.getConfiguration().getAppender(fa.getName()));

		if (verbose > 0) {
			lc.getConfiguration().getLoggerConfig("org.jcvi.vigor").setLevel(verbose == 1 ? Level.DEBUG: Level.TRACE);
		}
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
