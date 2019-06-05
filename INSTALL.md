# Building and Installing VIGOR4

## Build Dependencies

### Maven

VIGOR4 uses Maven to build and package the program. Version 3.5 or
later is recommended.

## Runtime Dependencies

### A Unix environment

Although VIGOR4 may work on other operating systems, it has only been
tested in a linux environment.

### Java 8 or above

VIGOR4 uses features, such as lambda expressions and the `Stream` API,
that are only available in Java 8 or above.

### Exonerate

By default, VIGOR4 uses exonerate to generate its initial
alignments. Exonerate is licensed under the GPL and is available from
the European Bioinformatics Institute at
[Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate).
VIGOR4 has been tested with the latest released exonerate version
2.2.0 and with the development version 2.4.0. 

Note: There is a bug in exonerate that is triggered by VIGOR when
exonerate is compiled with assertions enabled. Add --disable-assert to
configure options when building from source.

### Vigor Viral Database

VIGOR4 requires the VIGOR viral database. It is distributed separately and is available at
[VIGOR DB](https://github.com/JCVenterInstitute/VIGOR_DB)

## Building

From the root folder where the `pom.xml` file is located, type the
following on the commandline

```
%mvn clean package -DskipTests
```

This will create the zip file referenced in the Installing and
Configuring section

Running the tests during the build requires the path to the viral
database be set. This can be done by setting the
vigor4.reference_database_path system property on the commandlline

```
%mvn clean package -Dvigor4.reference_database_path=PATH_TO_DATABASE
```


## Installing and Configuring

Vigor is distributed as a zip file. Unzipping will create a directory
structure with a root directory vigor-[VERSION], for example:

```
    vigor-4.0.0/
    vigor-4.0.0/bin
    vigor-4.0.0/config
    vigor-4.0.0/lib
```

The default location for the configuration file is config/vigor.ini
and a skeleton configuration file is distributed with vigor4. The
location to the configuration file may be passed on the command line
using the `--config-file` option or set via the `VIGOR_CONFIG_FILE`
environment variable.

The minimum configuration requires setting the location of the
exonerate binary, the path to the vigor viral database and the
directory under which VIGOR4 can create temporary files and
directories. An example might look like

```
    reference_database_path=/data/VIGOR-DB/Reference_DBs/
    exonerate_path=/usr/local/exonerate-2.2.0/bin/exonerate
    temporary_directory=/tmp/vigor-temp
```

For a full listing of the configuration parameters with a description,
type the following on the commandline:

```
% VIGOR-4.0.0/bin/vigor4 --list-config-parameters
```

Configuration parameters may be set in a number of ways:

- as environment variables
- as java system properties passed via the JAVA_OPTS enviroment variable
- on the commandline by using the -P or --parameter options
- in the configuration file

Where conflicting configuration parameters are set, values are resolved in the following order:

- Parameters set via the command line
- Parameters set via environment variables
- Parameters set via system properties
- Parameters set in the viral specific config file, if any
- Parameters set in the general config file

VIGOR4 will print a warning for any unrecognized configuration parameters.

