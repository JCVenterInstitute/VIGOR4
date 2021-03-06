package org.jcvi.vigor.component;

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

    private final boolean has_ribosomal_slippage;
    private final String slippage_motif;
    private final int slippage_offset;
    private final int slippage_frameshift;

    public Ribosomal_Slippage ( boolean has_slippage, String motif, int offset, int frameshift ) {

        this.has_ribosomal_slippage = has_slippage;
        this.slippage_motif = motif;
        this.slippage_offset = offset;
        this.slippage_frameshift = frameshift;
    }

    public static Ribosomal_Slippage parseFromString(String slippageString) throws IllegalArgumentException {
        String[] temp = slippageString.split("/");
        if (temp.length == 3) {
            return new Ribosomal_Slippage(true, temp[2], Integer.parseInt(temp[0]), Integer.parseInt(temp[1]));
        }
        throw new IllegalArgumentException(
                String.format("Invalid ribsomal slippage format\"%s\". Format is OFFSET/FRAMESHIFT/MOTIF_REGEX", slippageString)
        );
    }

    public final static Ribosomal_Slippage NO_SLIPPAGE = new Ribosomal_Slippage(false, "", 0, 0);
}
