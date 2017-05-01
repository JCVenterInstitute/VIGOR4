package com.vigor.utils;

import java.util.HashMap;
import org.apache.commons.configuration2.INIConfiguration;
import org.apache.commons.configuration2.builder.fluent.Configurations;
import org.apache.commons.configuration2.ex.ConfigurationException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import com.vigor.exception.VigorException;
import com.vigor.service.VigorInitializationService;

public class LoadDefaultParameters {
	private static final Logger LOGGER = LogManager.getLogger(VigorInitializationService.class);
	

	public static HashMap<String, String> loadVigorParameters() {
		final HashMap<String, String> vigorParametersList = new HashMap<String, String>();
		try {
			Configurations configs = new Configurations();
			INIConfiguration iniConfig = configs.ini(Thread.currentThread().getContextClassLoader().getResource(VigorUtils.getConfigIniPath()));

			iniConfig.getSections().stream().forEach(i -> iniConfig.getSection(i).getKeys()
					.forEachRemaining(n -> vigorParametersList.put(n, iniConfig.getSection(i).getString(n))));

			for (String key : vigorParametersList.keySet()) {
				System.out.println(key + " : " + vigorParametersList.get(key));
			}

		} catch (ConfigurationException e) {
			VigorException.printExceptionMessage(e.getMessage());
			LOGGER.error(e.getMessage(),e);
		}
		return vigorParametersList;

	}

}
