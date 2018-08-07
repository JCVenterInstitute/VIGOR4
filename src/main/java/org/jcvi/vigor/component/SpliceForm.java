package org.jcvi.vigor.component;

import org.jcvi.vigor.exception.VigorRuntimeException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class SpliceForm {
    public enum SpliceType {EXON, INTRON}

    public final SpliceType type;
    public final long length;

    public SpliceForm(SpliceType type, long length) {
        this.type = type;
        this.length = length;
    }

    private final static Pattern spliceFormPattern = Pattern.compile("(?<spec>[ie]-?\\d+)");

    public static SpliceForm fromString(String spliceForm) {
        spliceForm = spliceForm.toLowerCase();
        return new SpliceForm(spliceForm.charAt(0) == 'e' ? SpliceType.EXON : SpliceType.INTRON,
                              Long.parseLong(spliceForm.substring(1)));
    }

    public String toString() {
        return (type == SpliceType.EXON ? "e" : "i") + Long.toString(this.length);
    }

    public static String spliceFormsToString(List<SpliceForm> spliceForms) {
        return String.join("", spliceForms.stream().map(SpliceForm::toString).collect(Collectors.toList()));
    }

    public static List<SpliceForm> parseFromString(String spliceForm) {
        if (spliceForm == null || spliceForm.isEmpty()) {
            return Collections.emptyList();
        }
        Matcher matcher = spliceFormPattern.matcher(spliceForm);
        List<SpliceForm> spliceForms = new ArrayList<>();
        int end = 0;
        while (matcher.find()) {
            if (matcher.start() != end) {
                throw new VigorRuntimeException(String.format("malformed splice_form: %s problem substring %s",
                                                              spliceForm,
                                                              spliceForm.substring(end, matcher.start())));
            }
            String spec = matcher.group("spec");
            end = matcher.end();
            spliceForms.add(fromString(spec));
        }
        return spliceForms;
    }
}
