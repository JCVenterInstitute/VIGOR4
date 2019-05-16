package org.jcvi.vigor.utils;


import java.util.HashMap;
import java.util.Map;
import java.util.function.Supplier;

public class OutputWriters {

    private OutputWriters() {
    }

    public static Map<String, Supplier<IOutputWriter>> Writers = new HashMap<>();
    static {
        Writers.put("TBL", TBLWriter::new);
        Writers.put("PEP", PEPWriter::new);
        Writers.put("CDS", CDSWriter::new);
        Writers.put("ALN", AlignmentWriter::new);
        Writers.put("SUM", SUMWriter::new);
        Writers.put("GFF3", GFF3Writer::new);
    }
}
