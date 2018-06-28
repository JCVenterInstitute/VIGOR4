package org.jcvi.vigor.utils;

import java.util.*;
import java.util.stream.Collectors;

public class VigorConfiguration {

    public static class ValueWithSource {

        public final String value;
        public final String source;

        private ValueWithSource ( String value, String source ) {

            this.value = value;
            this.source = Objects.requireNonNull(source);
        }

        public static ValueWithSource of ( String value, String source ) {

            return new ValueWithSource(value, source);
        }

        @Override
        public String toString () {

            return String.format("\"%s\" [%s]", value, source);
        }
    }

    private VigorConfiguration defaults = null;
    private final Map<ConfigurationParameters, String> values = new HashMap<>();
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
        this.values.putAll(defaults.values);
    }

    public String get ( ConfigurationParameters parameter ) {

        if (this.values.containsKey(parameter)) {
            return this.values.get(parameter);
        }
        if (this.defaults != null) {
            return this.defaults.get(parameter);
        }
        return null;
    }

    public void put ( ConfigurationParameters parameter, String value ) {

        this.values.put(parameter, value);
    }

    public void putAll ( Map<ConfigurationParameters, String> other ) {

        this.values.putAll(other);
    }

    public String getOrDefault ( ConfigurationParameters parameter, String defaultValue ) {

        if (containsKey(parameter)) {
            return get(parameter);
        }
        return defaultValue;
    }

    public void setDefaults ( VigorConfiguration defaults ) {

        this.defaults = new VigorConfiguration(defaults);
    }

    // TODO
    public Set<Map.Entry<ConfigurationParameters, String>> entrySet () {

        Set<ConfigurationParameters> seenParameters = new HashSet<>();
        Set<Map.Entry<ConfigurationParameters, String>> entries = new HashSet<>();
        entries.addAll(values.entrySet());
        entries.stream().map(es -> es.getKey()).forEach(seenParameters:: add);
        if (defaults != null) {
            defaults.entrySet()
                    .stream()
                    .filter(es -> !seenParameters.contains(es.getKey()))
                    .forEach(entries:: add);
        }
        return entries;
    }

    public boolean containsKey ( ConfigurationParameters parameter ) {

        return values.containsKey(parameter) || ( defaults != null && defaults.containsKey(parameter) );
    }

    // return value and configuration source for a given parameter
    public Optional<ValueWithSource> getWithSource ( ConfigurationParameters parameter ) {

        if (values.containsKey(parameter)) {
            return Optional.of(ValueWithSource.of(values.get(parameter), source));
        }
        if (defaults != null) {
            return defaults.getWithSource(parameter);
        }
        return Optional.empty();
    }

    public Map<ConfigurationParameters, String> toMap () {

        return entrySet().stream()
                .collect(Collectors.toMap(es -> es.getKey(), es -> es.getValue()));
    }

    public Set<ConfigurationParameters> keySet () {

        return entrySet().stream().map(es -> es.getKey()).collect(Collectors.toSet());
    }
}
