package org.jcvi.vigor.service;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.Filter;
import org.apache.logging.log4j.core.Layout;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.config.Configuration;
import org.apache.logging.log4j.core.filter.LevelRangeFilter;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.UserFacingException;
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

	public VigorConfiguration initializeVigor(Namespace inputs) throws VigorException {

            boolean isComplete = false;
            boolean isCircular = false;
            Boolean complete_gene = inputs.getBoolean(CommandLineParameters.completeGenome);
            if (complete_gene != null && complete_gene) {
                isComplete = true;
            }
            Boolean circular_gene = inputs.getBoolean(CommandLineParameters.circularGenome);
            if (circular_gene != null && circular_gene) {
                isComplete = true;
                isCircular = true;
            }
            VigorConfiguration configuration = loadParameters(inputs);
            configuration.putString(ConfigurationParameters.CircularGene, isCircular ? "1" : "0");
            configuration.putString(ConfigurationParameters.CompleteGene, isComplete ? "1" : "0");
            return configuration;

	}

	public List<VigorConfiguration> getDefaultConfigurations() throws VigorException {
		List<VigorConfiguration> configurations = new ArrayList<>();
		VigorConfiguration defaultConfiguration = LoadDefaultParameters
				.loadVigorConfiguration("defaults",
										Thread.currentThread().getContextClassLoader().getResource(VigorUtils.getDefaultConfigurationPath()),
										ConfigurationParameters.Flags.PROGRAM_CONFIG_SET);

		configurations.add(defaultConfiguration);

		VigorConfiguration propConfiguration = new VigorConfiguration("system-properties");
		String val;
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			val = System.getProperty(param.getSystemPropertyName());
			if (val != null) {
				if (param.hasFlag(ConfigurationParameters.Flags.VERSION_4)) {
					propConfiguration.putString(param, val);
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
					envConfiguration.putString(param, val);
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

	public VigorConfiguration loadParameters(Namespace inputs) throws VigorException{

		List<VigorConfiguration> configurations = new ArrayList<>();
		configurations.addAll(getDefaultConfigurations());
		configurations.addAll(getCommandLineConfiguration(inputs));

		String reference_db_dir = (String) getConfigValue(ConfigurationParameters.ReferenceDatabasePath, configurations);
		LOGGER.debug("Reference database path is {}", reference_db_dir);

		String reference_db= inputs.getString(CommandLineParameters.referenceDB);
		LOGGER.debug("reference database is {}", reference_db);

		if ("any".equals(reference_db)) {
			throw new UserFacingException("Auto-selecting reference database is not implemented");
		} else if (inputs.get(CommandLineParameters.genbankDB) != null) {
			throw new UserFacingException("Using genbank files as the reference database is not implemented");
		} else if (reference_db == null || reference_db.isEmpty()) {
			throw new UserFacingException("no reference database provided");
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
			throw new UserFacingException("Reference database path is required");
		}

		LOGGER.debug("Reference_db is {}", reference_db);
		if ( reference_db == null || reference_db.isEmpty()) {
			throw new UserFacingException("reference database is required");
		}

		try {
			VigorUtils.checkFilePath("Reference database file", reference_db,
									 VigorUtils.FileCheck.EXISTS,
									 VigorUtils.FileCheck.FILE,
									 VigorUtils.FileCheck.READ);
		} catch (VigorException e) {
			throw new UserFacingException(e.getMessage(), e);
		}

		configurations.get(configurations.size() -1).put(ConfigurationParameters.ReferenceDatabaseFile, reference_db);

		String virusSpecificConfig = (String) getConfigValue(ConfigurationParameters.VirusSpecificConfiguration, configurations);
		String virusSpecificConfigPath = (String) getConfigValue(ConfigurationParameters.VirusSpecificConfigurationPath, configurations);
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
			VigorUtils.checkFilePath("temporary directory", temporaryDirectory,
									 VigorUtils.FileCheck.DIRECTORY,
									 VigorUtils.FileCheck.READ,
									 VigorUtils.FileCheck.WRITE);
		} else {
			if (! tempDir.mkdirs()) {
				throw new VigorException(String.format("unable to create temporary directory %s", tempDir));
			}
		}


		return defaultConfiguration;
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

		String outputPath = inputs.getString(CommandLineParameters.outputPrefix);
		File outputFile= new File(outputPath).getAbsoluteFile();
		if (! (outputFile.getParentFile().exists() || outputFile.getParentFile().mkdirs()) ) {
			throw new VigorException(String.format("unable to create directory %s", outputFile.getParent()));
		}
		if( ! (outputFile.getParentFile().exists() && outputFile.getParentFile().isDirectory()) ){
			throw new VigorException(String.format("Invalid output prefix %s. Please provide a directory followed by a file prefix", outputPath));
		}
		initiateReportFile(outputFile.getParentFile().getAbsolutePath(), outputFile.getName(), inputs.getInt(CommandLineParameters.verbose));

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
			VigorUtils.checkFilePath("config file", config_file,
									 VigorUtils.FileCheck.EXISTS,
									 VigorUtils.FileCheck.READ,
									 VigorUtils.FileCheck.FILE);
		    LOGGER.debug("loading config file {}", config_file);
		    // use the file as configuration name so it's unambigious
			VigorConfiguration configFileConfiguration = LoadDefaultParameters.loadVigorConfiguration(config_file, new File(config_file),
																									  ConfigurationParameters.Flags.COMMANDLINE_SET,
																									  ConfigurationParameters.Flags.PROGRAM_CONFIG_SET);
			configurations.add(configFileConfiguration);
		}

		VigorConfiguration commandLineConfig = new VigorConfiguration("commandline");

		commandLineConfig.putString(ConfigurationParameters.OutputPrefix, outputFile.getName());
		commandLineConfig.putString(ConfigurationParameters.OutputDirectory, outputFile.getParentFile().getAbsolutePath());

		// defaults to false on commandline, so set it in the commandline configuration if true, so we don't override
		// a value in the vigor.ini file
		if (inputs.getBoolean(CommandLineParameters.overwriteOutputFiles)) {
			commandLineConfig.putString(ConfigurationParameters.OverwriteOutputFiles,
										inputs.getBoolean(CommandLineParameters.overwriteOutputFiles) ? "true" : "false");
		}

		String min_gene_coverage = inputs.getString(CommandLineParameters.minCoverage);
		if (min_gene_coverage != null ) {
			commandLineConfig.putString(ConfigurationParameters.GeneMinimumCoverage, min_gene_coverage);
		}
		String frameshift_sensitivity = inputs.getString(CommandLineParameters.frameshiftSensitivity);
		if (frameshift_sensitivity != null ) {
			commandLineConfig.putString(ConfigurationParameters.FrameShiftSensitivity, frameshift_sensitivity);
		}

		String locus_tag = inputs.getString(CommandLineParameters.locusTag);
		if (locus_tag != null) {
			commandLineConfig.putString(ConfigurationParameters.Locustag, locus_tag);
		}

		Boolean ignore_reference_requirements = inputs.getBoolean(CommandLineParameters.ignoreReferenceRequirements);
		if (ignore_reference_requirements != null && ignore_reference_requirements) {
			LOGGER.debug("Ignoring legacy parameter ignore reference requirement");
		}

		Boolean jcvi_rules = inputs.getBoolean(CommandLineParameters.jcviRules);
		if (jcvi_rules != null) {
			LOGGER.debug("Ignoring legacy parameter evalue");
		}

		String reference_database = inputs.getString(CommandLineParameters.referenceDB);
		if (reference_database != null) {
			Path referenceDBPath = Paths.get(reference_database);
			if (referenceDBPath.isAbsolute()) {
				commandLineConfig.putString(ConfigurationParameters.ReferenceDatabasePath, referenceDBPath.getParent().toString());
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
			commandLineConfig.putString(ConfigurationParameters.ReferenceDatabasePath, reference_database_path);
		}

		String virus_specific_path = inputs.getString(CommandLineParameters.virusSpecificConfigPath);
		if (virus_specific_path != null) {
			commandLineConfig.putString(ConfigurationParameters.VirusSpecificConfigurationPath, virus_specific_path);
		}

		virus_specific_path = inputs.getString(CommandLineParameters.virusSpecificConfig);
		if (virus_specific_path != null) {
			commandLineConfig.putString(ConfigurationParameters.VirusSpecificConfiguration, virus_specific_path);
		}
        commandLineConfig.putString(ConfigurationParameters.Verbose,
									inputs.getInt(CommandLineParameters.verbose) > 0 ? "true": "false");

		configurations.add(commandLineConfig);


		List<String> parameters = inputs.getList(CommandLineParameters.parameters);
		if (parameters != null) {
			final Pattern splitter = Pattern.compile("~~");
			Map<String, String> temp = parameters.stream()
												 .flatMap(p -> splitter.splitAsStream(p.trim()))
												 .map(s -> s.split("=", 2))
												 .collect(Collectors.toMap(a -> a[0], a -> a.length > 1 ? a[1] : ""));
			VigorConfiguration commandLineParametersConfig = LoadDefaultParameters.configurationFromMap("commandline-parameters", temp,
																										ConfigurationParameters.Flags.PROGRAM_CONFIG_SET,
																										ConfigurationParameters.Flags.COMMANDLINE_SET);
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
			throw new VigorException(String.format("Virus specific config file %s doesn't exist or is not readable", virusSpecificConfig));
		}
		// virus specific configuration files may not exist
		if (configFile.exists()) {
			VigorConfiguration virusSpecificParameters = LoadDefaultParameters.loadVigorConfiguration(configPath, configFile,
																									  ConfigurationParameters.Flags.VIRUS_SET,
																									  ConfigurationParameters.Flags.GENE_SET
																									  );
			LOGGER.info("loaded virus specific config from {}", configPath);
			if (virusSpecificParameters.hasSection(VigorConfiguration.METADATA_SECTION)) {
				String version = virusSpecificParameters.getOrDefault(VigorConfiguration.METADATA_SECTION, ConfigurationParameters.Version, "");
				String virusName = virusSpecificParameters.getOrDefault(VigorConfiguration.METADATA_SECTION, ConfigurationParameters.VirusName, "");
				if (! (version.isEmpty()  && virusName.isEmpty()) ) {
					LOGGER.info("Virus config for virus \"{}\" version \"{}\"",
								NullUtil.emptyOrElse(virusName, "not set"),
								NullUtil.emptyOrElse(version, "not set")
								);
				}

			}
			virusSpecificParameters.setDefaults(vigorConfiguration);
			return virusSpecificParameters;
		}
		return vigorConfiguration;
	}

	public void initiateReportFile(String outputDir, String outputPrefix, int verbose){
        LoggerContext lc = (LoggerContext) LogManager.getContext(false);
        Configuration config = lc.getConfiguration();

        FileAppender fa = FileAppender.newBuilder()
									  .withName("mylogger")
									  .withAppend(false)
									  .withFileName(new File(outputDir, outputPrefix+".rpt").toString())
									  .build();
        fa.start();
        config.addAppender(fa);

        Layout warningLayout = PatternLayout.newBuilder()
				.withConfiguration(config)
				.withPattern("%level %msg %exception{full}\n")
				.build();

		Filter warningFilter = LevelRangeFilter.createFilter(Level.FATAL, Level.WARN, Filter.Result.ACCEPT, Filter.Result.DENY);

		FileAppender warnings = FileAppender.newBuilder()
											.withName("__warnings")
											.withAppend(false)
											.withLayout(warningLayout)
											.withFileName(new File(outputDir, outputPrefix + ".warnings").toString())
											.build();
		warningFilter.start();
		warnings.addFilter(warningFilter);
		warnings.start();
		config.addAppender(warnings);

		lc.getLogger("org.jcvi.vigor").addAppender(warnings);
		lc.getLogger("org.jcvi.vigor").addAppender(fa);

		if (verbose > 0) {
			lc.getConfiguration().getLoggerConfig("org.jcvi.vigor").setLevel(verbose == 1 ? Level.DEBUG: Level.TRACE);
		}
        lc.updateLoggers();
    }

    private Object getConfigValue(ConfigurationParameters param, List<VigorConfiguration> configurations) {
		Object value = null;
		for (int i = configurations.size() - 1; i >= 0; i--) {
			value = configurations.get(i).get(param);
			if (value != null) {
				break;
			}
		}
		return value;
	}
}
