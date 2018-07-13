package org.jcvi.vigor.utils;

import org.jcvi.vigor.testing.category.Fast;
import org.jcvi.vigor.testing.category.Isolated;
import org.junit.Test;
import org.junit.experimental.categories.Category;

import java.util.*;

import static org.hamcrest.CoreMatchers.*;
import static org.hamcrest.MatcherAssert.assertThat;

@Category({Fast.class, Isolated.class})
public class VigorConfigurationTest {

    /**
     * Test that default values are accessible and overrideable
     */
    @Test
    public void testDefaultValues () {

        VigorConfiguration defaults = new VigorConfiguration("defaults");

        String defaultOutputPrefix = "/usr/local/scratch/1234";
        defaults.put(ConfigurationParameters.OutputPrefix, defaultOutputPrefix);
        defaults.put(ConfigurationParameters.CircularGene, "0");
        defaults.put(ConfigurationParameters.MinimumMissingAASize, "19");

        VigorConfiguration overrides = new VigorConfiguration("overrides");
        overrides.put(ConfigurationParameters.CircularGene, "1");
        overrides.setDefaults(defaults);

        assertThat("Override value should mask default value",
                overrides.get(ConfigurationParameters.CircularGene), equalTo("1"));
        assertThat("Default values should be visible when no override is set",
                overrides.get(ConfigurationParameters.OutputPrefix), equalTo(defaultOutputPrefix));
        assertThat("When no key exists in either configuration null is returned",
                overrides.get(ConfigurationParameters.ExonMinimumSize), equalTo(null));

        String newDefaultOutputPrefix = "new/default/value";
        defaults.put(ConfigurationParameters.OutputPrefix, newDefaultOutputPrefix);
        assertThat("Changing default config after calling setDefaults has no affect on overrides",
                overrides.get(ConfigurationParameters.OutputPrefix), equalTo(defaultOutputPrefix));

        assertThat("override should report containing default keys",
                overrides.containsKey(ConfigurationParameters.OutputPrefix), equalTo(true));
        assertThat("override should report containing it's own keys",
                overrides.containsKey(ConfigurationParameters.CircularGene), equalTo(true));
        assertThat("override should report not containing keys not in it or default",
                overrides.containsKey(ConfigurationParameters.CompleteGene), equalTo(false));

        Set<ConfigurationParameters> keySet = overrides.keySet();
        assertThat("keySet should return unique set across override and defaults", keySet.size(), equalTo(3));
        assertThat("keySet should return unique set across override and default", keySet,
                allOf(hasItem(ConfigurationParameters.OutputPrefix),
                        hasItem(ConfigurationParameters.CircularGene),
                        hasItem(ConfigurationParameters.MinimumMissingAASize)
                        ));

        Optional<VigorConfiguration.ValueWithSource> valWithSource = overrides.getWithSource(ConfigurationParameters.ExonMinimumSize);
        assertThat("Keys not set in either override or default are not present", valWithSource.isPresent(), equalTo(false));
        valWithSource = overrides.getWithSource(ConfigurationParameters.OutputPrefix);
        assertThat("Set keys should be reported as present", valWithSource.isPresent(), equalTo(true));
        assertThat("Source should be reported correctly", valWithSource.get().source, equalTo("defaults"));

        String overrideOutputPrefix = "not/here/or/there";
        overrides.put(ConfigurationParameters.OutputPrefix, overrideOutputPrefix );
        assertThat("override.put should override default",
                overrides.get(ConfigurationParameters.OutputPrefix), equalTo(overrideOutputPrefix));
        assertThat("override.put should not affect default",
                defaults.get(ConfigurationParameters.OutputPrefix), equalTo(newDefaultOutputPrefix));


        Map<ConfigurationParameters, String> newValues = new HashMap<>();
        newValues.put(ConfigurationParameters.CompleteGene,"0");
        overrides.putAll(newValues);

        assertThat("putall only affects override",
                overrides.get(ConfigurationParameters.CompleteGene), equalTo("0"));
        assertThat("putall does not affect not default",
                defaults.get(ConfigurationParameters.CompleteGene), nullValue());


        Set<Map.Entry<ConfigurationParameters,String>> entrySet = overrides.entrySet();
        assertThat("entrySet contains unique key,value pairs from override and default",
                entrySet.size(), equalTo(4));
        assertThat("entrySet contains unique key,value pairs from override and default", entrySet,
                allOf(
                        hasItem(equalTo(entryOf(ConfigurationParameters.OutputPrefix, overrideOutputPrefix))),
                        hasItem(equalTo(entryOf(ConfigurationParameters.CompleteGene, "0"))),
                        hasItem(equalTo(entryOf(ConfigurationParameters.CircularGene, "1"))),
                        hasItem(equalTo(entryOf(ConfigurationParameters.MinimumMissingAASize, "19")))
                ));
    }

    private Map.Entry<ConfigurationParameters, String> entryOf ( ConfigurationParameters parameter, String value ) {

        return new AbstractMap.SimpleEntry(parameter, value);
    }
}
