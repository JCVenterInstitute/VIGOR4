package com.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import java.util.HashMap;
import java.util.List;

/**
 * Created by snettem on 5/9/2017.
 */
@Component
@Scope("prototype")
@Data
public class StructuralSpecifications {

    private String shared_cds;
    private String is_optional;
    private String is_required;
    private List<String> excludes_gene;
    private HashMap<String,Integer> tiny_exon3;
    private HashMap<String,Integer> tiny_exon5;
    private long minFunctionalLength;

}
