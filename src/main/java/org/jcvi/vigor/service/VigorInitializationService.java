package org.jcvi.vigor.service;

import net.sourceforge.argparse4j.inf.Namespace;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.exception.UserFacingException;
import org.jcvi.vigor.utils.*;
import org.springframework.stereotype.Service;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * This class is under service layer. The methods in this class has
 * functionality to initialize the vigor application
 */
@Service
public class VigorInitializationService {

	private static final Logger LOGGER = LogManager.getLogger(VigorInitializationService.class);
	private static final Function<String, EnumSet<ConfigurationParameters.Flags>> defaultConfigFlags = (section) -> {

		switch (section) {
			case VigorConfiguration.DEFAULT_GENE_SECTION:
				return EnumSet.of(ConfigurationParameters.Flags.GENE_SET);
			case VigorConfiguration.DEFAULT_VIRUS_SECTION:
				return EnumSet.of(ConfigurationParameters.Flags.VIRUS_SET);
			default:
				return EnumSet.of(ConfigurationParameters.Flags.PROGRAM_CONFIG_SET);
		}
	};

	public static final Function<String, EnumSet<ConfigurationParameters.Flags>> programConfigFlags = (section) -> {

		switch (section) {
			case VigorConfiguration.DEFAULT_GENE_SECTION:
				return EnumSet.of(ConfigurationParameters.Flags.GENE_SET);
			case VigorConfiguration.DEFAULT_VIRUS_SECTION:
				return EnumSet.of(ConfigurationParameters.Flags.VIRUS_SET);
			default:
				return EnumSet.of(ConfigurationParameters.Flags.PROGRAM_CONFIG_SET, ConfigurationParameters.Flags.COMMANDLINE_SET);
		}
	};

	public static final Function<String, EnumSet<ConfigurationParameters.Flags>> virusConfigFlags = (section) -> {

		switch (section) {
			case VigorConfiguration.METADATA_SECTION:
				return EnumSet.of(ConfigurationParameters.Flags.METADATA_SET);
			case VigorConfiguration.DEFAULT_GENE_SECTION:
				return EnumSet.of(ConfigurationParameters.Flags.GENE_SET);
			default:
				if (section != null && section.startsWith("gene:")) {
					return EnumSet.of(ConfigurationParameters.Flags.GENE_SET);
				}
				return EnumSet.of(ConfigurationParameters.Flags.VIRUS_SET);
		}
	};

	public class DatabaseInfo {
		public final File databaseFile;
		public final Optional<File> configFile;

		DatabaseInfo(File databaseFile, File configFile) {
			this.databaseFile = databaseFile;
			this.configFile = Optional.ofNullable(configFile);
		}
	}

	/**
	 * @param inputs:
	 *            User provided command line inputs Retrieve each Genomic
	 *            sequence from the input file and determine AlignmentEvidence
	 */

	public VigorConfiguration initializeVigor(Namespace inputs) throws VigorException {

            boolean isCircular = false;
            Boolean circular_gene = inputs.getBoolean(CommandLineParameters.circularGenome);
            if (circular_gene != null && circular_gene) {
                isCircular = true;
            }
            VigorConfiguration configuration = loadParameters(inputs);
            configuration.putString(ConfigurationParameters.CircularGene, isCircular ? "1" : "0");
            return configuration;

	}

	public List<VigorConfiguration> getDefaultConfigurations() throws VigorException {
		List<VigorConfiguration> configurations = new ArrayList<>();
		URL defaultConfigurationFile = Thread.currentThread().getContextClassLoader().getResource(VigorUtils.getDefaultConfigurationPath());
		Map<String,Map<String,String>> sectionMap = LoadDefaultParameters.configFileToSectionMap(defaultConfigurationFile);
		String configurationSource = "defaults";
		VigorConfiguration config = LoadDefaultParameters
				.configurationFromSectionMap(configurationSource, sectionMap, defaultConfigFlags);

		configurations.add(config);
		if (sectionMap.containsKey(VigorConfiguration.DEFAULT_VIRUS_SECTION)) {
			config = LoadDefaultParameters.configurationFromMap(String.format("%s:%s", configurationSource, VigorConfiguration.DEFAULT_VIRUS_SECTION),
																sectionMap.get(VigorConfiguration.DEFAULT_VIRUS_SECTION),
																section -> EnumSet.of(ConfigurationParameters.Flags.VIRUS_SET));
			configurations.add(config);
		}

		if (sectionMap.containsKey(VigorConfiguration.DEFAULT_GENE_SECTION)) {
			config = LoadDefaultParameters.configurationFromMap(String.format("%s:%s", configurationSource, VigorConfiguration.DEFAULT_GENE_SECTION),
																sectionMap.get(VigorConfiguration.DEFAULT_GENE_SECTION),
																section -> EnumSet.of(ConfigurationParameters.Flags.GENE_SET));
			configurations.add(config);
		}

		VigorConfiguration propConfiguration = new VigorConfiguration("system-properties");
		String val;
		for (ConfigurationParameters param: ConfigurationParameters.values()) {
			val = System.getProperty(param.getSystemPropertyName());
			if (val != null) {
				if (param.hasFlag(ConfigurationParameters.Flags.VERSION_4)) {
					propConfiguration.putString(param, val);
				} else {
					LOGGER.debug("Ignoring non-VIGOR4 parameter {}={} set via system properties", param.configKey, val, param.getSystemPropertyName());
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
					LOGGER.debug("Ignoring non-VIGOR4 parameter {}={} set via environment variable {}", param.configKey, val, param.getEnvVarName());
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
		LOGGER.debug("Reference database is {}", reference_db);

		if ("any".equals(reference_db)) {
			throw new UserFacingException("Auto-selecting reference database is not implemented");
		} else if (! NullUtil.isNullOrEmpty(reference_db)){
			File file = new File(reference_db);
			if(file.exists() && file.isFile() ){
				reference_db=file.getAbsolutePath();
				if (reference_db_dir == null || reference_db_dir.isEmpty()) {
					reference_db_dir = file.getParent();
				}
			} else if ( ! NullUtil.isNullOrEmpty(reference_db_dir)) {
				Path referenceDBFile = Paths.get(reference_db_dir, reference_db);
				if (referenceDBFile.toFile().exists()) {
					reference_db = Paths.get(reference_db_dir, reference_db).toString();
				} else {
					// maybe reference_db is an alias, eg "Influenza A"
					Optional<Path> aliasDatabase = findReferenceDBByAlias(reference_db_dir, reference_db);
					if (aliasDatabase.isPresent()) {
						reference_db = aliasDatabase.get().toString();
					}
				}
			}
			configurations.get(configurations.size() -1).put(ConfigurationParameters.ReferenceDatabaseFile, reference_db);
			LOGGER.debug("Reference_db is {}", reference_db);
			String virusSpecificConfig = (String) getConfigValue(ConfigurationParameters.VirusSpecificConfiguration, configurations);
			String virusSpecificConfigPath = (String) getConfigValue(ConfigurationParameters.VirusSpecificConfigurationPath, configurations);
			if (virusSpecificConfigPath == null || virusSpecificConfigPath.isEmpty()) {
				virusSpecificConfigPath = reference_db_dir;
			}
			String referenceDBName = Paths.get(reference_db).getFileName().toString();
			List<VigorConfiguration> virusConfigurations = loadVirusConfiguration(referenceDBName, virusSpecificConfigPath, virusSpecificConfig);
			configurations.addAll(1, virusConfigurations);
		}

		VigorConfiguration defaultConfiguration = mergeConfigurations(configurations);



		return defaultConfiguration;
	}

    public List<DatabaseInfo> getDatabaseInfo(String referenceDatabasePath) throws IOException, VigorException {
        LOGGER.debug("Looking for databases under {}", referenceDatabasePath);
        try {
			VigorUtils.checkFilePath("reference database path", referenceDatabasePath,
									 VigorUtils.FileCheck.EXISTS,
									 VigorUtils.FileCheck.DIRECTORY,
									 VigorUtils.FileCheck.READ);
		} catch (VigorException e) {
			throw new UserFacingException(e.getMessage(), e);
		}
        List<DatabaseInfo> databases = Files.list(Paths.get(referenceDatabasePath))
											.filter(f -> f.getFileName().toString().endsWith(("_db")))
											.map(f -> {
                 File configFile = f.resolveSibling(f.getFileName().toString() + ".ini").toFile();
                 configFile = configFile.exists() ? configFile : null;
                 return new DatabaseInfo(f.toFile(), configFile);
             }).collect(Collectors.toList());
        return databases;
    }


	// see if any config file has a defined alias
	private Optional<Path> findReferenceDBByAlias(String reference_db_dir, String alias) throws VigorException {
		LOGGER.debug("Checking under reference database path {} for database alias {}", reference_db_dir, alias);
		try {
			VigorConfiguration config;
			for (DatabaseInfo dbInfo: getDatabaseInfo(reference_db_dir)) {
				if (dbInfo.configFile.isPresent()) {
					config = mergeConfigurations(loadVirusConfiguration(dbInfo.configFile.get()));
					for (String configAlias: (List<String>) config.getOrDefault(VigorConfiguration.METADATA_SECTION, ConfigurationParameters.Alias, Collections.EMPTY_LIST)) {
						if (alias.equalsIgnoreCase(configAlias)) {
							return Optional.of(dbInfo.databaseFile.toPath());
						}
					}
				}
			}
		} catch (IOException e) {
			throw new VigorException(String.format("Problem reading configuration files from %s", reference_db_dir), e);
		}
		return Optional.empty();
	}


	public VigorConfiguration mergeConfigurations(List<VigorConfiguration> configurations) {
		// now setup all the defaults

		VigorConfiguration previousConfig = null;

		LOGGER.trace("Merging configurations: {}", () -> configurations.stream().map(VigorConfiguration::getSource).collect(Collectors.joining(" <- ")));
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
			VigorUtils.checkFilePath("config file", config_file,
									 VigorUtils.FileCheck.EXISTS,
									 VigorUtils.FileCheck.READ,
									 VigorUtils.FileCheck.FILE);
		    LOGGER.debug("loading config file {}", config_file);
		    Map<String,Map<String,String>> sectionMap = LoadDefaultParameters.configFileToSectionMap(new File(config_file));
		    // use the file as configuration name so it's unambigious
			VigorConfiguration config = LoadDefaultParameters.configurationFromSectionMap(config_file, sectionMap, programConfigFlags);
			configurations.add(config);
			if (sectionMap.containsKey(VigorConfiguration.DEFAULT_VIRUS_SECTION)) {
				config = LoadDefaultParameters.configurationFromMap(String.format("%s:%s",config_file, VigorConfiguration.DEFAULT_VIRUS_SECTION), sectionMap.get(VigorConfiguration.DEFAULT_VIRUS_SECTION),
																	section -> EnumSet.of(ConfigurationParameters.Flags.VIRUS_SET));
				configurations.add(config);
			}

			if (sectionMap.containsKey(VigorConfiguration.DEFAULT_GENE_SECTION)) {
				config = LoadDefaultParameters.configurationFromMap(String.format("%s:%s",config_file, VigorConfiguration.DEFAULT_GENE_SECTION),
																	sectionMap.get(VigorConfiguration.DEFAULT_GENE_SECTION),
																	section -> EnumSet.of(ConfigurationParameters.Flags.GENE_SET));
				configurations.add(config);
			}

		}

		VigorConfiguration commandLineConfig = new VigorConfiguration("commandline");

		String outputPath = inputs.getString(CommandLineParameters.outputPrefix);
		if (! NullUtil.isNullOrEmpty(outputPath)) {
			File outputFile = new File(outputPath).getAbsoluteFile();
			commandLineConfig.putString(ConfigurationParameters.OutputPrefix, outputFile.getName());
			commandLineConfig.putString(ConfigurationParameters.OutputDirectory, outputFile.getParentFile().getAbsolutePath());
		}

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
			VigorConfiguration commandLineParametersConfig = LoadDefaultParameters.configurationFromMap("commandline-parameters", temp, programConfigFlags);
			configurations.add(commandLineParametersConfig);
		}

		return configurations;
	}
	/**
	 *
	 * @param reference_db
	 *            : Virus Specific reference_db
	 * @return configuration : Default vigor Parameters will be overridden
	 *         by virus specific parameters
	 */
	public List<VigorConfiguration> loadVirusConfiguration(String reference_db,
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
			return Collections.EMPTY_LIST;
		}
		File configFile = new File(configPath);
		// config file was specified, but doesn't exist
		if (! (configFile.exists() || virusSpecificConfig == null)) {
			throw new VigorException(String.format("virus specific configfile doesn't exist %s", virusSpecificConfig));
		}
		return loadVirusConfiguration(configFile);

	}
	public List<VigorConfiguration> loadVirusConfiguration(File configFile) throws VigorException{
		List<VigorConfiguration> configurations = new ArrayList<>(2);
		LOGGER.debug("virus specific config file for reference_db {} is {}", configFile.getParent(), configFile.getName());
		// virus specific configuration files may not exist
		if (configFile.exists()) {
			Map<String,Map<String,String>> sectionMap = LoadDefaultParameters.configFileToSectionMap(configFile);
			VigorConfiguration virusSpecificParameters = LoadDefaultParameters.configurationFromSectionMap(configFile.getPath(), sectionMap, virusConfigFlags);
			configurations.add(virusSpecificParameters);
			LOGGER.info("loaded virus specific config from {}", configFile.getPath());
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
			Map<String,String> defaultGeneConfigMap = sectionMap.getOrDefault(VigorConfiguration.DEFAULT_GENE_SECTION, Collections.EMPTY_MAP);
			if (! defaultGeneConfigMap.isEmpty()) {
				VigorConfiguration defaultGeneConfig = LoadDefaultParameters.configurationFromMap(String.format("%s:%s", configFile.getPath(),
																												VigorConfiguration.DEFAULT_GENE_SECTION),
																								  defaultGeneConfigMap,
																								  section -> EnumSet.of(ConfigurationParameters.Flags.GENE_SET));
				configurations.add(defaultGeneConfig);
			}
		}
		return configurations;
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
