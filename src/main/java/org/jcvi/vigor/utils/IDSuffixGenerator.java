package org.jcvi.vigor.utils;

import java.util.Arrays;
import java.util.Iterator;

public class IDSuffixGenerator implements Iterable<String>, Iterator<String> {

    private static final int start = (int) 'a';
    private static final int end = (int) 'z';

    private int[] current;
    public IDSuffixGenerator() {
        current = new int[] {};
    }

    @Override
    public Iterator<String> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return true;
    }

    @Override
    public String next() {
        increment();
        return currentValue();
    }

    private void increment() {
        int val;

        boolean rollover = true;
        for (int i=current.length -1 ; i>= 0; i--) {
            val = current[i];
            if (rollover) {
                val += 1;
                current[i] = val;
            }
            rollover = false;
            if (val > end ) {
                val = start;
                current[i] =  val;
                rollover = true;
            }
        }
        if (rollover) {
            current = Arrays.copyOf(current, current.length + 1);
            current[current.length -1 ] = start;
        }
    }

    private String currentValue() {
        StringBuffer value = new StringBuffer();
        for (int val: current) {
            value.append( (char) val);
        }
        return value.toString();
    }
}
