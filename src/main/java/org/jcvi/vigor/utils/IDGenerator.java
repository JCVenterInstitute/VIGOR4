package org.jcvi.vigor.utils;

import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;

public class IDGenerator implements Iterable<String>, Iterator<String> {

    private final String seed;
    private AtomicInteger counter = new AtomicInteger(0);

    public IDGenerator ( String seed ) {

        this.seed = seed;
    }

    @Override
    public Iterator<String> iterator () {

        return this;
    }

    @Override
    public boolean hasNext () {

        return true;
    }

    @Override
    public String next () {

        return seed + "." + counter.incrementAndGet();
    }

    public static IDGenerator of ( String seed ) {

        return new IDGenerator(seed);
    }
}
