package org.jcvi.vigor.utils;

import java.io.File;
import java.net.URL;
import java.util.Comparator;
import java.util.Map;
import java.util.HashMap;
import java.util.stream.Collectors;

import org.apache.commons.configuration2.INIConfiguration;
import org.apache.commons.configuration2.builder.fluent.Configurations;
import org.apache.commons.configuration2.ex.ConfigurationException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.VigorInitializationService;

public class LoadDefaultParameters {
	private static final Logger LOGGER = LogManager.getLogger(VigorInitializationService.class);

	public static VigorConfiguration loadVigorConfiguration(String sourceName, File path) throws VigorException {
		try {
			Configurations configs = new Configurations();
			INIConfiguration iniConfig = configs.ini(path);
			return loadVigorConfiguration(sourceName, iniConfig);
		} catch (ConfigurationException e) {
			LOGGER.error(e.getMessage(),e);
			throw new VigorException(String.format("unable to load configuration file %s", path), e);
		}
	}

	public static VigorConfiguration loadVigorConfiguration(String sourceName, URL path) throws VigorException {
		try {
			Configurations configs = new Configurations();
			INIConfiguration iniConfig = configs.ini(path);
			return loadVigorConfiguration(sourceName, iniConfig);
		} catch (ConfigurationException e) {
			LOGGER.error(e.getMessage(),e);
			throw new VigorException(String.format("unable to load configuration file %s", path), e);
		}
	}

	public static VigorConfiguration loadVigorConfiguration(String sourceName, INIConfiguration iniConfig)  {
			final Map<String, String> parametersMap = new HashMap<>();
				iniConfig.getSections().stream().forEach(i -> iniConfig.getSection(i).getKeys()
																	   .forEachRemaining(n -> parametersMap.put(n, iniConfig.getSection(i).getString(n).replaceAll("\\s+",""))));
				return configurationFromMap(sourceName, parametersMap);
	}
	
	public static VigorConfiguration configurationFromMap(String sourceName, Map<String, String> configurationMap) {
		final VigorConfiguration configuration = new VigorConfiguration(sourceName);
		Map<String,String> configEntries = new HashMap<>();
		configEntries.putAll(configurationMap);
		for (ConfigurationParameters parameter: ConfigurationParameters.values()) {
			if (configEntries.containsKey(parameter.configKey)) {
				configuration.put(parameter, configEntries.remove(parameter.configKey));
			}
		}
		if (! configEntries.isEmpty()) {
			LOGGER.warn("In config {} unrecognized configuration entries:\n{}\n",
					sourceName,
					configEntries.entrySet().stream()
								 .sorted(Comparator.comparing(es -> es.getKey()))
								 .map(es -> String.format("\t%s: %s", es.getKey(), es.getValue()))
								 .collect(Collectors.joining("\n"))
			);
		}
		return configuration;
	}
}


