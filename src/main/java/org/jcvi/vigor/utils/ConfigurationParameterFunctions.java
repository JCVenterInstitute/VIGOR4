package org.jcvi.vigor.utils;

import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.vigor.component.SpliceForm;
import org.jcvi.vigor.component.StopTranslationException;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

public class ConfigurationParameterFunctions {

    public static class ValueFunction {
        final Class<?> valueClass;
        final Function<String, ?> valueFunction;

        public ValueFunction(Class<?> valueClass, Function<String,?> valueFunction ) {
            this.valueClass = valueClass;
            this.valueFunction = valueFunction;
        }


    }

    public static ValueFunction of(Class<?> valueClass, Function<String,?> valueFunction) {
        return new ValueFunction(valueClass, valueFunction);
    }

    public static class InvalidValue extends RuntimeException {
        public InvalidValue(String message) {
            super(message);
        }
    }

    public static ValueFunction toListOfStrings = of(List.class, s -> Arrays.stream(s.split(","))
                                                                            .map(i -> i.trim())
                                                                            .collect(Collectors.toList()));
    public static ValueFunction toSpliceForms = of (List.class, SpliceForm::parseFromString);
    public static ValueFunction toStopException = of(StopTranslationException.class,
                                                     s -> {
                                                         String[] temp = s.split("/");
                                                         return new StopTranslationException(true,
                                                                                             AminoAcid.parse(temp[ 1 ]),
                                                                                             temp[ 2 ],
                                                                                             Integer.parseInt(temp[ 0 ]));

                                                     });

    public static ValueFunction toPercent = toBoundedDouble(0,100);

    private static final Set<String> booleanYes = new HashSet<>();
    private static final Set<String> booleanNo = new HashSet<>();

    static {
        booleanYes.addAll(Arrays.asList("true","t","yes","1"));
        booleanNo.addAll(Arrays.asList("false","f","no","0"));
    }

    public static ValueFunction toBoolean = of(Boolean.class, s -> {
        String test = s.toLowerCase();
        if (booleanYes.contains(test)){
            return Boolean.TRUE;
        } else if (booleanNo.contains(test)) {
            return Boolean.FALSE;
        }
        throw new InvalidValue(String.format("acceptable values are %s or %s",
                                             String.join(",", booleanYes),
                                             String.join(",", booleanNo)));
    });

    public static ValueFunction toInteger = of(Integer.class, Integer::parseInt);
    public static ValueFunction toDouble = of(Double.class, Double::parseDouble);

    public static ValueFunction toBoundedInteger(int min, int max) {
        return of(Integer.class,
                  s -> {
                      Integer i = Integer.parseInt(s);
                      if (i < min || i > max) {
                          throw new InvalidValue(String.format("Value %s must be between %s and %s inclusive", s, min, max));
                      }
                      return i;
                  });
    }

    public static ValueFunction toPositiveInteger() {
        return toBoundedInteger(0, Integer.MAX_VALUE);
    }

    public static ValueFunction toLong = of(Long.class, s -> Long.parseLong(s));
    public static ValueFunction toBoundedLong(long min, long max) {
        return of(Long.class,
                  s -> {
                      Long i = Long.parseLong(s);
                      if (i < min || i > max) {
                          throw new InvalidValue(String.format("Value %s must be between %s and %s inclusive", s, min, max));
                      }
                      return i;
                  });
    }

    public static ValueFunction toBoundedDouble (double min, double max) {
        return of(Double.class,
                  s -> {
                      Double i = Double.parseDouble(s);
                      if (i < min || i > max) {
                          throw new InvalidValue(String.format("Value %s must be between %s and %s inclusive", s, min, max));
                      }
                      return i;
                  });
    }

    public static ValueFunction isMemberOfSet(Object ... members) {
        return isMemberOfSet(String.class, s -> s, members);
    }

    public static ValueFunction isMemberOfSet(Class valueClass, Function<String,Object> valueFunction, Object ... members) {
        final Set<Object> validMembers = new HashSet<>();
        validMembers.addAll(Arrays.asList(members));
        return of(valueClass,
                  s -> {
                      Object v = valueFunction.apply(s);
                      if (! validMembers.contains(s)) {
                          throw new InvalidValue(String.format("value \"%s\" must be one of %s", s,
                                                               Arrays.stream(members).map(String::valueOf).collect(Collectors.joining(",")))
                          );
                      }
                      return v;
                  });
    }
}
