package org.jcvi.vigor.utils;

import java.util.Objects;

public class NullUtil {

    public static <T> T nullOrElse(T obj, T defaultVal) {
        return obj != null ? obj : Objects.requireNonNull(defaultVal, "default value may not be null");
    }

    public static  boolean isNullOrEmpty(String str) {
        return str == null || str.trim().isEmpty();
    }

    public static String emptyOrElse(String str, String defaultVal) {
        return ! isNullOrEmpty(str) ? str : Objects.requireNonNull(defaultVal, "default value may not be null");
    }
}
