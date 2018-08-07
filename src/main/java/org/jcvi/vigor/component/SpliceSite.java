package org.jcvi.vigor.component;

import com.google.common.collect.ImmutableList;

import java.util.List;

public class SpliceSite {

    public static final List<SpliceSite> DEFAULT_SPLICE_SITES = ImmutableList.of(new SpliceSite("GT", "AG"));

    public final String donor;
    public final String acceptor;

    public SpliceSite(String donor, String acceptor ) {

        this.donor = donor;
        this.acceptor = acceptor;
    }
}
