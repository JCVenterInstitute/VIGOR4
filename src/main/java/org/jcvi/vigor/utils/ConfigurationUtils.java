package org.jcvi.vigor.utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.exception.VigorException;

import java.lang.reflect.Parameter;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class ConfigurationUtils {

    private static final Logger LOGGER = LogManager.getLogger(ConfigurationUtils.class);

    public static VigorConfiguration configurationFromMap(String sourceName,
                                                          Map<String,String> values,
                                                          ConfigurationParameters.Flags ... flags) throws VigorException {
        return configurationFromMap(sourceName, p -> p.configKey, values, flags);
    }
    public static VigorConfiguration configurationFromMap(String sourceName,
                                                          Function<ConfigurationParameters, String> keyFunction,
                                                          Map<String,String> values,
                                                          ConfigurationParameters.Flags ... flags) throws VigorException {
        Map<String, Map<String, String>> sectionMap = new HashMap<>();
        sectionMap.put(VigorConfiguration.DEFAULT_SECTION, values);
        return configurationFromSectionMap(sourceName, keyFunction, sectionMap, flags);
    }

    public static VigorConfiguration configurationFromSectionMap (String sourceName ,Map<String, Map<String,String>> configurationMap, ConfigurationParameters.Flags ... flags) throws VigorException {
        return configurationFromSectionMap(sourceName, p -> p.configKey, configurationMap, flags);
    }

    public static VigorConfiguration configurationFromSectionMap (String sourceName, Function<ConfigurationParameters, String> keyFunction, Map<String, Map<String,String>> configurationMap, ConfigurationParameters.Flags ... flags)
            throws VigorException {

        EnumSet<ConfigurationParameters.Flags> flagSet;

        if (flags.length == 0) {
            flagSet = ConfigurationParameters.Flags.SET_FLAGS;
        } else {
            flagSet = EnumSet.copyOf(Arrays.asList(flags));
        }

        final VigorConfiguration configuration = new VigorConfiguration(sourceName);
        Map<String, Map<String,String>> configEntries = new HashMap<>();
        for (String section: configurationMap.keySet()) {
            configEntries.computeIfAbsent(section, k -> new HashMap<>()).putAll(configurationMap.get(section));
        }
        List<String> errors = new ArrayList<>();

        for (String section: configEntries.keySet()) {
            String sectionString = section == null || section.equals(VigorConfiguration.DEFAULT_SECTION) ? "" : String.format(" section \"%s\" ",section);
            for (ConfigurationParameters parameter : ConfigurationParameters.values()) {
                String key = keyFunction.apply(parameter);
                if (configEntries.get(section).containsKey(key)) {
                    String value = configEntries.get(section).remove(key);
                    if ( parameter.hasOneOrMoreFlags(ConfigurationParameters.Flags.IGNORE)) {
                        LOGGER.debug("Ignoring parameter \"{}\" with value \"{}\" in config source {}{}",
                                     key, value, sourceName, sectionString);
                    }
                    if (! parameter.hasOneOrMoreFlags(ConfigurationParameters.Flags.VERSION_4)) {
                        LOGGER.warn("Ignoring deprecated parameter \"{}\" with value \"{}\" in config source {}{}",
                                    key, value, sourceName, sectionString);
                        continue;
                    }
                    LOGGER.trace(() -> String.format("checking config %s%s key %s settable for flags %s: %s",
                                 sourceName,
                                 sectionString,
                                 key,
                                 flagSet.stream().map(e -> e.toString()).collect(Collectors.joining(",")),
                                 ! parameter.hasFlags(flagSet).isEmpty()));

                    if (parameter.hasFlags(flagSet).isEmpty()) {
                        errors.add(String.format("from source \"%s\"%s attempt to set config value %s with value %s; configuration flags for %s are %s",
                                                 sourceName,
                                                 sectionString,
                                                 key,
                                                 value,
                                                 key,
                                                 String.join(",", flagSet.stream().map(String::valueOf).collect(Collectors.joining(",")))));

                    }
                    try {
                        configuration.putString(section, parameter, value);
                    } catch (ConfigurationParameterFunctions.InvalidValue e) {
                        errors.add(String.format("from source \"%s\"%s got error: %s",
                                                 sourceName,
                                                 sectionString,
                                                 e.getMessage()));
                    }
                }
            }
        }
        StringBuilder unrecognized = new StringBuilder();
        for (String section: configEntries.keySet()) {
            if (configEntries.get(section).isEmpty()) {
                continue;
            }
            unrecognized.append("Unrecognized parameters in config ").append(sourceName);
            if (! (section == null || section.equals(VigorConfiguration.DEFAULT_SECTION) )) {
                unrecognized.append(" section ").append(section);
            }
            unrecognized.append("\n");
            configEntries.get(section).entrySet()
                         .stream()
                         .sorted(Comparator.comparing(es -> es.getKey()))
                         .map(es -> String.format("\t%s: %s", es.getKey(), es.getValue()))
                         .forEach(v -> unrecognized.append(v).append("\n"));

        }
        String unrecognizedWarning = unrecognized.toString();
        if (!unrecognizedWarning.isEmpty()) {
            LOGGER.warn(unrecognizedWarning);
        }

        if (! errors.isEmpty()) {
            throw new VigorException(errors.stream()
                                           .collect(Collectors.joining("\n", "Errors found in configuration:\n", ""))
            );
        }
        return configuration;
    }

    public static String getGeneSectionName(String geneSymbol) {
        return "gene:" + geneSymbol;
    }

}
