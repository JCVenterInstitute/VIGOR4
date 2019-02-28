package org.jcvi.vigor.utils;

import java.util.*;
import java.util.stream.Collectors;

public class VigorConfiguration {

    public static final String DEFAULT_SECTION = "__DEFAULT__" ;
    public static final String METADATA_SECTION = "Metadata";
    public static final String DEFAULT_VIRUS_SECTION="Default:virus";
    public static final String DEFAULT_GENE_SECTION="Default:gene";

    public static class ValueWithSource {

        public final Object value;
        public final String source;

        private ValueWithSource ( Object value, String source ) {

            this.value = value;
            this.source = Objects.requireNonNull(source);
        }

        public static ValueWithSource of ( Object value, String source ) {

            return new ValueWithSource(value, source);
        }

        @Override
        public String toString () {

            return String.format("\"%s\" [%s]", value, source);
        }
    }

    private VigorConfiguration defaults = null;
    private final Map<String,Map<ConfigurationParameters, Object>> values = new HashMap<>();
    private final String source;

    // New configuration with no defaults
    public VigorConfiguration ( String source ) {
        this.source = source;
    }

    // New configuration with existing defaults
    public VigorConfiguration ( String source, VigorConfiguration defaults ) {

        this(source);
        if (defaults != null) {
            this.defaults = new VigorConfiguration(defaults);
        }
    }

    // Copy constructor
    public VigorConfiguration ( VigorConfiguration defaults ) {

        this(defaults.source, defaults.defaults);
        for (String section: defaults.values.keySet()) {
            this.values.computeIfAbsent(section, k -> new HashMap<>())
                       .putAll(defaults.values.getOrDefault(section,Collections.EMPTY_MAP));
        }
    }

    public String getSource() {
        return this.source;
    }

    public <T> T get ( ConfigurationParameters parameter ) {
        return get(DEFAULT_SECTION, parameter);
    }

    public <T> T get (String section, ConfigurationParameters parameter) {

        if (this.values.getOrDefault(section,Collections.EMPTY_MAP).containsKey(parameter)) {
            return (T) this.values.get(section).get(parameter);
        }
        if (this.values.getOrDefault(DEFAULT_SECTION, Collections.EMPTY_MAP).containsKey(parameter)) {
            return (T) this.values.get(DEFAULT_SECTION).get(parameter);
        }
        if (this.defaults != null) {
            return this.defaults.get(section, parameter);
        }

        return null;
    }

    public void put(ConfigurationParameters param, Object value) {
        put(DEFAULT_SECTION, param, value);
    }

    public void putAll(String section, Map<ConfigurationParameters, Object> values) {
        this.values.computeIfAbsent(section, k -> new HashMap<>()).putAll(values);
    }

    public void putAll(Map<ConfigurationParameters, Object> values) {
        putAll(DEFAULT_SECTION, values);
    }

    // TODO atleast check type
    public void put(String section, ConfigurationParameters parameter, Object value) {
        this.values.computeIfAbsent(section, k-> new HashMap<>()).put(parameter, value);
    }

    public void putString(ConfigurationParameters parameter, String value ) {
        putString(DEFAULT_SECTION, parameter, value);
    }

    public void putString(String sectionName, ConfigurationParameters parameter, String value ) {
        put(sectionName, parameter, parameter.stringToValue(value));
    }

    public void putAllString(String sectionName, Map<ConfigurationParameters, String> other ) {
        for (Map.Entry<ConfigurationParameters,String> entry: other.entrySet()) {
            putString(sectionName, entry.getKey(), entry.getValue());
        }
    }

    public void putAllString(Map<ConfigurationParameters, String> other ) {
        putAllString(DEFAULT_SECTION, other);
    }

    public <T> T getOrDefault ( ConfigurationParameters parameter, Object defaultValue ) {
        return getOrDefault(DEFAULT_SECTION, parameter, defaultValue);
    }

    public <T> T getOrDefault (String sectionName, ConfigurationParameters parameter, Object defaultValue ) {

        if (containsKey(sectionName, parameter)) {
            return (T) get(sectionName, parameter);
        }
        return (T) defaultValue;
    }

    public void setDefaults ( VigorConfiguration defaults ) {

        this.defaults = new VigorConfiguration(defaults);
    }

    public Set<Map.Entry<ConfigurationParameters, Object>> entrySet () {
        return entrySet(DEFAULT_SECTION);
    }

    public Set<Map.Entry<ConfigurationParameters, Object>> entrySet (String section) {

        Set<ConfigurationParameters> seenParameters = new HashSet<>();
        Set<Map.Entry<ConfigurationParameters, Object>> entries = new HashSet<>();
        entries.addAll(values.getOrDefault(section, Collections.EMPTY_MAP).entrySet());
        entries.stream().map(es -> es.getKey()).forEach(seenParameters:: add);
        // add default section values
        ((Set< Map.Entry<ConfigurationParameters,Object>>) values.getOrDefault(DEFAULT_SECTION, Collections.EMPTY_MAP).entrySet())
                              .stream()
                              .filter(es -> !seenParameters.contains(es.getKey()))
                              .forEach(es -> {entries.add(es); seenParameters.add(es.getKey());});

        // add defaults
        if (defaults != null) {
            defaults.entrySet(section)
                    .stream()
                    .filter(es -> !seenParameters.contains(es.getKey()))
                    .forEach(entries:: add);
        }
        return entries;
    }

    public boolean containsKey ( ConfigurationParameters parameter ) {
        return containsKey(DEFAULT_SECTION, parameter);
    }

    public boolean containsKey (String section,  ConfigurationParameters parameter ) {

        return values.getOrDefault(section, Collections.EMPTY_MAP).containsKey(parameter) ||
                values.getOrDefault(DEFAULT_SECTION, Collections.EMPTY_MAP).containsKey(parameter) ||
                ( defaults != null && defaults.containsKey(section, parameter) );
    }

    public boolean hasSection(String section) {
        return values.containsKey(section) || (defaults != null && defaults.hasSection(section));
    }

    public Optional<ValueWithSource> getWithSource ( ConfigurationParameters parameter ) {
        return getWithSource(DEFAULT_SECTION, parameter);
    }
    // return value and configuration source for a given parameter
    public Optional<ValueWithSource> getWithSource ( String section, ConfigurationParameters parameter ) {

        if ( (!section.equals(DEFAULT_SECTION)) && values.getOrDefault(section, Collections.EMPTY_MAP).containsKey(parameter)) {
            return Optional.of(ValueWithSource.of(values.get(section).get(parameter), source + ":" + section));
        }
        if (values.getOrDefault(DEFAULT_SECTION, Collections.EMPTY_MAP).containsKey(parameter)) {
            return Optional.of(ValueWithSource.of(values.get(DEFAULT_SECTION).get(parameter), source));
        }
        if (defaults != null) {
            return defaults.getWithSource(section, parameter);
        }
        return Optional.empty();
    }


    public Set<ConfigurationParameters> keySet () {
        return keySet(DEFAULT_SECTION);
    }

    public Set<ConfigurationParameters> keySet (String section) {
        return entrySet(section).stream().map(es -> es.getKey()).collect(Collectors.toSet());
    }

    public Map<ConfigurationParameters, Object> getSectionConfig(String sectionName) {
        Map<ConfigurationParameters, Object> result = new HashMap<>();
        for (Map.Entry<ConfigurationParameters, Object> entry: entrySet(sectionName)) {
            result.put(entry.getKey(), entry.getValue());
        }
        return result;
    }


}
