package org.jcvi.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import java.util.Collections;
import java.util.List;

/**
 * Created by snettem on 5/8/2017.
 */
@Component
@Scope("prototype")
@Data
public class Splicing {

    private final boolean isSpliced;
    private final List<SpliceSite> nonCanonical_spliceSites;
    private final String spliceform;

    public class SpliceSite {

        public final String donor;
        public final String acceptor;

        public SpliceSite ( String donor, String acceptor ) {

            this.donor = donor;
            this.acceptor = acceptor;
        }
    }

    public Splicing ( boolean splicing, List<SpliceSite> spliceSites, String spliceform ) {

        this.isSpliced = splicing;
        this.nonCanonical_spliceSites = spliceSites;
        this.spliceform = spliceform;
    }

    public static Splicing NO_SPLICING = new Splicing(false, Collections.EMPTY_LIST, "");
}

