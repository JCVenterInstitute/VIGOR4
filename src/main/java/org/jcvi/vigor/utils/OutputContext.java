package org.jcvi.vigor.utils;

import java.util.*;

public class OutputContext {

    public enum Key {
        GENOME,
        GENE,
        PEP
    }

    private final Map<Key, String> context = new EnumMap(Key.class);

    public OutputContext () {
    }

    public OutputContext addContext(Key key, String value) {
        context.put(key, value);
        return this;
    }

    public Optional<String> getContext(Key key) {
        return Optional.ofNullable(context.get(key));
    }

    public Optional<String> removeContext(Key key) {
        return Optional.ofNullable(context.remove(key));
    }

    public Set<Key> keySet() {
        return context.keySet();
    }


}
