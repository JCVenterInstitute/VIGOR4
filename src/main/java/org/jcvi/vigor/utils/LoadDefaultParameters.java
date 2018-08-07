package org.jcvi.vigor.utils;

import java.io.File;
import java.net.URL;
import java.util.*;
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

    @SuppressWarnings("Duplicates")
    public static VigorConfiguration loadVigorConfiguration ( String sourceName, File path, ConfigurationParameters.Flags... flags ) throws VigorException {

        try {
            Configurations configs = new Configurations();
            INIConfiguration iniConfig = configs.ini(path);
            return loadVigorConfiguration(sourceName, iniConfig, flags);
        } catch (ConfigurationException e) {
            LOGGER.error(e.getMessage(), e);
            throw new VigorException(String.format("unable to load configuration file %s", path), e);
        }
    }

    @SuppressWarnings("Duplicates")
    public static VigorConfiguration loadVigorConfiguration (String sourceName, URL path, ConfigurationParameters.Flags ... flags)
            throws VigorException {

        try {
            Configurations configs = new Configurations();
            INIConfiguration iniConfig = configs.ini(path);
            return loadVigorConfiguration(sourceName, iniConfig, flags);
        } catch (ConfigurationException e) {
            LOGGER.error(e.getMessage(), e);
            throw new VigorException(String.format("unable to load configuration file %s", path), e);
        }
    }

    public static VigorConfiguration loadVigorConfiguration (String sourceName, INIConfiguration iniConfig, ConfigurationParameters.Flags ... flags) throws VigorException {

        final Map<String, Map<String,String>> parametersMap = new HashMap<>();
        Map<String,String> sectionMap;
        String key;
        String val;
        String configSection;
        for (String sectionName: iniConfig.getSections()) {
            configSection = sectionName == null || sectionName.isEmpty() ? VigorConfiguration.DEFAULT_SECTION : sectionName;
            sectionMap = parametersMap.computeIfAbsent(configSection, k -> new HashMap<>());
            Iterator<String> keyIter = iniConfig.getSection(sectionName).getKeys();
            while (keyIter.hasNext()) {
                key = keyIter.next();
                val = iniConfig.getSection(sectionName).getString(key);
                val = val.trim().replaceAll("\\s+"," ");
                sectionMap.put(key,val);
            }
        }
        return configurationFromSectionMap(sourceName, parametersMap, flags);
    }

    public static VigorConfiguration configurationFromMap(String sourceName, Map<String,String> configurationMap, ConfigurationParameters.Flags ... flags) throws VigorException {
       return ConfigurationUtils.configurationFromMap(sourceName, configurationMap, flags);
    }

    public static VigorConfiguration configurationFromSectionMap (String sourceName, Map<String, Map<String,String>> configurationMap, ConfigurationParameters.Flags ... flags)
            throws VigorException {
        return ConfigurationUtils.configurationFromSectionMap(sourceName, configurationMap, flags);
    }
}


