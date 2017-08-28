package org.jcvi.vigor.component;

import lombok.Data;

import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;
import java.util.Map;


/**
 * Created by snettem on 5/8/2017.
 */
@Component
@Scope("prototype")
@Data
public class Splicing {

    private boolean isSpliced=false;
    private Map<String,String> nonCanonical_spliceSites; //Donor+Acceptor
    private String spliceform;


}
