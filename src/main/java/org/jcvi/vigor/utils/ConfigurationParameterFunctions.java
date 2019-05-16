package org.jcvi.vigor.utils;

import com.google.common.collect.Sets;
import org.jcvi.jillion.core.residue.aa.AminoAcid;
import org.jcvi.vigor.component.SpliceForm;
import org.jcvi.vigor.component.StopTranslationException;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

public class ConfigurationParameterFunctions {

    public static class ValueFunction {
        final Class<?> valueClass;
        final Function<String, ?> valueFunction;
        final Function<Object, String> stringFunction;

        public ValueFunction(Class<?> valueClass, Function<String,?> valueFunction, Function<Object, String> stringFunction) {
            this.valueClass = valueClass;
            this.valueFunction = valueFunction;
            this.stringFunction = stringFunction;
        }

        public ValueFunction(Class<?> valueClass, Function<String,?> valueFunction ) {
            this(valueClass, valueFunction, String::valueOf);
        }


    }

    public static ValueFunction of(Class<?> valueClass, Function<String,?> valueFunction) {
        return new ValueFunction(valueClass, valueFunction);
    }

    public static ValueFunction of(Class<?> valueClass, Function<String,?> valueFunction, Function<Object, String> stringFunction) {
        return new ValueFunction(valueClass, valueFunction, stringFunction);
    }

    public static class InvalidValue extends RuntimeException {
        public InvalidValue(String message) {
            super(message);
        }
    }

    public static Function<String,List<String>> stringToList = string -> Arrays.stream(string.split(","))
                                                                               .map(String::trim)
                                                                               .filter(s -> ! s.isEmpty())
                                                                               .collect(Collectors.toList());
    public static ValueFunction toListOfStrings = of(List.class, stringToList,
                                                     v -> String.join(",", (List) v));
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
    public static ValueFunction toTinyExonMap = of(Map.class, (s) ->
        Arrays.stream(s.split(","))
              .map(i -> i.split(":", 2))
              .collect(Collectors.toMap(t -> t[0].trim(),
                                        t -> t.length > 1 ? Integer.parseInt(t[1].trim()): 0))
    );
    private static final Set<String> booleanYes = new HashSet<>();
    private static final Set<String> booleanNo = new HashSet<>();

    static {
        booleanYes.addAll(Arrays.asList("true","t","yes","y","1"));
        booleanNo.addAll(Arrays.asList("false","f","no","n","0"));
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

    public static ValueFunction isPresent = of(Boolean.class, (s) -> {
        if (! (s == null || s.isEmpty())) {
            throw new InvalidValue("No value allowed");
        }
        return Boolean.TRUE;
    });

    public static ValueFunction isPresentOrBoolean = of(Boolean.class, (s) -> {
        if ( s == null || s.isEmpty()) {
            return Boolean.TRUE;
        }
        return toBoolean.valueFunction.apply(s);
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

    public static ValueFunction toPositiveInteger = toBoundedInteger(0, Integer.MAX_VALUE);


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

    public static ValueFunction areMembersOfSet(String ... members) {

        Set<String> validMembers = Arrays.stream(members).map(String::toUpperCase).collect(Collectors.toSet());
        return of(Set.class,
                  s -> {
                      Set<String> value = stringToList.apply(s).stream().map(String::toUpperCase).collect(Collectors.toSet());
                      Set<String> unknownFormats = Sets.difference(value, validMembers);
                      if (!unknownFormats.isEmpty()) {
                          throw new InvalidValue(String.format("Unexpected value(s): %s. Value(s) must be one of: %s",
                                                               unknownFormats.stream()
                                                                             .sorted(String.CASE_INSENSITIVE_ORDER)
                                                                             .collect(Collectors.joining(",")),
                                                               validMembers.stream()
                                                                           .sorted(String.CASE_INSENSITIVE_ORDER)
                                                                           .collect(Collectors.joining(",")))
                          );
                      }
                      return value;
                  },
                  set -> ((Set<String>)set).stream().sorted(String.CASE_INSENSITIVE_ORDER).collect(Collectors.joining(","))
        );
    }
}
