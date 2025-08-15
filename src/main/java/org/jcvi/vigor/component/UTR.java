package org.jcvi.vigor.component;

import lombok.Data;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

@Component
@Scope("prototype")
@Data
public class UTR {

    public final String prime;
    public final String utrStrand;
    public final int startCoordinate;
    public final int endCoordinate;

    public UTR(String prime) {

        this.prime = prime;
        this.utrStrand = "";
        this.startCoordinate = 0;
        this.endCoordinate = 0;
    }

    public UTR(String prime, String utrStrand, int startCoordinate, int endCoordinate) {

        this.prime = prime;
        this.utrStrand = utrStrand;
        this.startCoordinate = startCoordinate;
        this.endCoordinate = endCoordinate;
    }

    public static UTR parseFromString(String utrString) throws IllegalArgumentException {
        if (utrString == null || utrString.isEmpty()) {
            throw new IllegalArgumentException(String.format("Invalid UTR format\"%s\". Format is BP/Sequence", utrString));
        }

        String[] temp = utrString.split("/");
        String prime = temp[0];
        int start = Integer.parseInt(temp[1]);
        int end = Integer.parseInt(temp[2]);
        String utr = temp[3];

        return new UTR(prime, utr, start, end);
    }

    public static final UTR NO_UTR3 = new UTR("3'");
    public static final UTR NO_UTR5 = new UTR("5'");
}
