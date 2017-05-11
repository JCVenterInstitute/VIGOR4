package com.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

/**
 * Created by snettem on 5/8/2017.
 */
@Component
@Scope("prototype")
@Data
public class Ribosomal_Slippage {

    private boolean is_ribosomal_slippage=false;
    private String slippage_motif;
    private int slippage_offset;
    private int slippage_frameshift;

}
