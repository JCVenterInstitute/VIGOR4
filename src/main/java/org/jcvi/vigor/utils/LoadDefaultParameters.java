package org.jcvi.vigor.utils;

import java.util.Map;
import java.util.HashMap;
import org.apache.commons.configuration2.INIConfiguration;
import org.apache.commons.configuration2.builder.fluent.Configurations;
import org.apache.commons.configuration2.ex.ConfigurationException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;
import org.jcvi.vigor.service.VigorInitializationService;

public class LoadDefaultParameters {
	private static final Logger LOGGER = LogManager.getLogger(VigorInitializationService.class);
	

	public static Map<String, String> loadVigorParameters(String path) {
		final Map<String, String> vigorParametersList = new HashMap<String, String>();
		try {
			Configurations configs = new Configurations();
			INIConfiguration iniConfig = configs.ini(Thread.currentThread().getContextClassLoader().getResource(path));

			iniConfig.getSections().stream().forEach(i -> iniConfig.getSection(i).getKeys()
					.forEachRemaining(n -> vigorParametersList.put(n, iniConfig.getSection(i).getString(n).replaceAll("\\s+",""))));

		} catch (ConfigurationException e) {
			VigorException.printExceptionMessage(e.getMessage());
			LOGGER.error(e.getMessage(),e);
		}
		return vigorParametersList;

	}
	
}


