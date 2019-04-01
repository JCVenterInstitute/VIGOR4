package org.jcvi.vigor.service;

import com.google.common.base.Joiner;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.jcvi.vigor.utils.NullUtil;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

@Service
public class VigorInputValidationService {

	private static final Logger LOGGER = LogManager.getLogger(VigorInputValidationService.class);

	public Namespace processInput(String... inputParams)  {
		ArgumentParser parser = getArgumentParser();
		return parser.parseArgsOrFail(inputParams);
	}

	private ArgumentParser getArgumentParser() {

		ArgumentParser parser = ArgumentParsers.newFor("vigor4")
											   .build()
											   .usage("${prog} -i inputfasta -o outputprefix [ -d refdb ]")
				.epilog(String.join("\n","Outputs:",
						"  outputprefix.rpt   - summary of program results",
						"  outputprefix.cds   - fasta file of predicted CDSs",
						"  outputprefix.pep   - fasta file of predicted proteins",
						"  outputprefix.tbl   - predicted features in GenBank tbl format",
						"  outputprefix.aln   - alignment of predicted protein to reference, and reference protein to genome"
						));

		MutuallyExclusiveGroup referenceGroup = parser.addMutuallyExclusiveGroup("reference database");
		MutuallyExclusiveGroup locusGroup = parser.addMutuallyExclusiveGroup("locus tag usage");

		ArgumentGroup unImplementedGroup = parser.addArgumentGroup("Unimplemented/Ignored for backward compatibilty");

		parser.addArgument("-i","--input-fasta")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.inputFile)
			  .metavar("<input fasta>")
			  .help("path to fasta file of genomic sequences to be annotated.");

		parser.addArgument("-o","--output-prefix")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.outputPrefix)
			  .metavar("<output prefix>")
			  .help("prefix for outputfile files, e.g. if the ouput prefix is /mydir/anno VIGOR will create output files /mydir/anno.tbl, /mydir/anno.stats, etc. An output prefix without a directory element will create the output files in the current working directory.");

		referenceGroup.addArgument("-d", "--reference-database")
					  .dest(CommandLineParameters.referenceDB)
					  .action(Arguments.store())
					  .metavar("<ref db>")
					  .help("specify the reference database to be used");

		parser.addArgument("-c", "--min-coverage")
			  .dest(CommandLineParameters.minCoverage)
			  .action(Arguments.store())
			  .help("minimum coverage of reference product (0-100) required to report a gene, by default coverage is ignored");

		unImplementedGroup.addArgument("-0", "--circular")
				   .dest(CommandLineParameters.circularGenome)
				   .help("complete circular genome (allows gene to span origin). This feature is currently unimplemented")
				   .action(Arguments.storeTrue());
		unImplementedGroup.addArgument("-f", "--frameshift-sensitivity")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.frameshiftSensitivity)
			  .choices("0","1","2")
			  .setDefault("1")
			  .help("frameshift sensitivity, 0=ignore frameshifts, 1=normal (default), 2=sensitive. ");

		// use storeConst rather than storeTrue to avoid automatically setting a default value
		locusGroup.addArgument("-l", "--no-locus-tags")
				  .dest(CommandLineParameters.locusTag)
				  .action(Arguments.storeConst())
				  .setConst("")
				  .help("do NOT use locus_tags in TBL file output (incompatible with -L)");
		locusGroup.addArgument("-L", "--locus-tags")
				  .dest(CommandLineParameters.locusTag)
				  .action(Arguments.store())
				  .nargs("?")
				  // default is used when the option is not present
				  .setDefault("vigor_")
				  // const is used when the option is present, but no argument is passed
				  .setConst("vigor_")
				  .metavar("<locus_tag_prefix>")
				  .help("USE locus_tags in TBL file output (incompatible with -l). If no prefix is provided, the prefix \"vigor_\" will be used.");

		parser.addArgument("-P","--parameter")
			  .action(Arguments.append())
			  .dest(CommandLineParameters.parameters)
			  .metavar("<parameter=value~~...~~parameter=value>")
			  .help("~~ separated list of VIGOR parameters to override default values. Use --list-config-parameters to see settable parameters.");
		unImplementedGroup.addArgument("-m","--ignore-reference-requirements")
			  .action(Arguments.storeTrue())
			  .dest(CommandLineParameters.ignoreReferenceRequirements)
			  .help("ignore reference match requirements (coverage/identity/similarity), sometimes useful when running VIGOR to evaluate raw contigs and rough draft sequences");
		parser.addArgument("-v","--verbose")
			  .action(Arguments.count())
			  .dest(CommandLineParameters.verbose)
			  .help("verbose logging (default=terse)");

		unImplementedGroup.addArgument("-x", "--ignore-refID")
			  .action(Arguments.append())
			  .dest(CommandLineParameters.ignoreRefID)
			  .metavar("<ref_id,...,ref_id>")
			  .help("comma separated list of reference sequence IDs to ignore (useful when debugging a reference database). Not currently implemented");

		parser.addArgument("--list-config-parameters")
			  .action(new ListConfigurations())
			  .choices("all","current")
			  .setDefault("current")
			  .setConst("current")
			  .nargs("?")
			  .help("list available configuration parameters and exit. By default only lists description, use the verbose option before this option to list more information");

		parser.addArgument("--list-databases")
			  .action(Arguments.storeTrue())
			  .dest(CommandLineParameters.listDatabases)
			  .help(String.format("list the names and other information about the found vigor compatible databases. Requires reference database path to be set either by passing the --reference-database-path command line parameter or setting %s in the configuration file", ConfigurationParameters.ReferenceDatabasePath.configKey));

		parser.addArgument("--version")
			  .action(new PrintVersion())
			  .help("print version information");

		parser.addArgument("--config-file")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.configFile)
			  .help("config file to use");

		parser.addArgument("--reference-database-path")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.referenceDB_Path)
			  .help("reference database path");

		parser.addArgument("--virus-config")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.virusSpecificConfig)
			  .help("Path to virus specific configuration");

		parser.addArgument("--virus-config-path")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.virusSpecificConfigPath)
			  .help("Path to directory containing virus specific config files.");

        parser.addArgument("--overwrite-output")
              .action(Arguments.storeTrue())
              .dest(CommandLineParameters.overwriteOutputFiles)
              .help("overwrite existing output files if they exist");

        parser.addArgument("--temporary-directory")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.temporaryDirectory)
			  .help("Root directory to use for temporary directories");

		return parser;
	}

	private static abstract class CustomNoArgumentAction implements ArgumentAction {

		@Override
		public void onAttach(Argument argument) {
		}

		@Override
		public boolean consumeArgument() {
			return false;
		}

	}

	private static class PrintVersion extends CustomNoArgumentAction {

		@Override
		public void run(ArgumentParser argumentParser, Argument argument, Map<String, Object> map, String s, Object o) throws ArgumentParserException {
			try {
				Properties buildProperties = new Properties();
				buildProperties.load(this.getClass().getResourceAsStream("/build.properties"));
				String branch = NullUtil.emptyOrElse(buildProperties.getProperty("branch"), "master");
				String host = NullUtil.nullOrElse(buildProperties.getProperty("build.host"),"").trim();
				String buildTimeString = NullUtil.nullOrElse(buildProperties.getProperty("build.time"),"").trim();

				host = host.contains("localhost") ? "" : host;
				if (! buildTimeString.isEmpty()) {
					LocalDateTime buildDateTime = LocalDateTime.of(Integer.parseInt(buildTimeString.substring(0,4)),
																   Integer.parseInt(buildTimeString.substring(4,6)),
																   Integer.parseInt(buildTimeString.substring(6,8)),
																   Integer.parseInt(buildTimeString.substring(9,11)),
																   Integer.parseInt(buildTimeString.substring(11,13)),
																   Integer.parseInt(buildTimeString.substring(13,15))
					);
					buildTimeString = buildDateTime.format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss", Locale.getDefault()));
				}
				String buildInfo = "";
				if (! (host.isEmpty() && buildTimeString.isEmpty()) ) {
					buildInfo = String.format(" ( built %s%s )",
											  ! host.isEmpty() ? String.format(" on host %s", host) : "",
											  ! buildTimeString.isEmpty() ? String.format("at %s", buildTimeString): "");
				}
				System.out.println(String.format("%s%s%s",
												 buildProperties.getProperty("vigor.version"),
												 "master".equals(branch) ? "" : String.format(" (branch %s)", branch),
												 buildInfo)
				);
				System.exit(0);
			} catch (IOException e ) {
				throw new ArgumentParserException(argumentParser);
			}
		}
	}
	private static class ListConfigurations implements ArgumentAction {
		@Override
		public void onAttach(Argument argument) {
		}

		@Override
		public boolean consumeArgument() {
			return true;
		}

		@Override
		public void run(ArgumentParser argumentParser, Argument argument, Map<String, Object> map, String s, Object o) throws ArgumentParserException {
			Function<ConfigurationParameters.Flags, String> flagToString = (flag) -> {
				switch (flag) {
					case VIRUS_SET:
						return "virus config";
					case GENE_SET:
						return "gene config/variant";
					case PROGRAM_CONFIG_SET:
						return "program config";
					case COMMANDLINE_SET:
						return "commandline";
					default:
						return null;
				}
			};

			String selectionString = String.valueOf(o);
			EnumSet<ConfigurationParameters.Flags> settableFlags = EnumSet.of(ConfigurationParameters.Flags.GENE_SET,
																			  ConfigurationParameters.Flags.VIRUS_SET,
																			  ConfigurationParameters.Flags.PROGRAM_CONFIG_SET,
																			  ConfigurationParameters.Flags.COMMANDLINE_SET);
			Map<Boolean,List<ConfigurationParameters>> partitionedParameters =
					Arrays.stream(ConfigurationParameters.values())
						  .sorted(Comparator.comparing(p-> p.configKey, String.CASE_INSENSITIVE_ORDER))
						  .collect(Collectors.partitioningBy( p ->	p.hasFlag(ConfigurationParameters.Flags.VERSION_4) &&
								  ! p.hasOneOrMoreFlags(ConfigurationParameters.Flags.IGNORE,
														ConfigurationParameters.Flags.METADATA)));
			boolean[] parametersToList = selectionString.equals("all") ? new boolean[] {true, false} : new boolean[] {true};
			boolean verbose = (int) map.getOrDefault(CommandLineParameters.verbose,0) > 0;
			for (boolean currentParameters: parametersToList) {

				if (currentParameters) {
					System.out.println("\nCONFIGURATION PARAMETERS\n");
				} else {
					System.out.println("\n\nDEPRECATED/IGNORED PARAMETERS\n");
				}

				for (ConfigurationParameters param : partitionedParameters.get(currentParameters)) {
					boolean printedSomething = false;
					System.out.println(param.configKey);
					System.out.println();
					if (currentParameters && verbose) {
						printedSomething = true;
						System.out.println(String.format("\t%-30s VIGOR_%s", "Environment variable:", param.configKey.toUpperCase()));
						System.out.println(String.format("\t%-30s vigor.%s", "System.property:", param.configKey));
						System.out.println(String.format("\t%-30s %s", "Settable levels", Joiner.on(",")
																								.skipNulls()
																								.join(param.hasFlags(settableFlags).stream()
																										   .sorted()
																										   .map(flagToString)
																										   .collect(Collectors.toList()))));
						if (param.hasFlag(ConfigurationParameters.Flags.UNIMPLEMENTED)) {
							System.out.println(String.format("\t%-30s This parameters is currently unimplemented","Note:"));
						}

					}
					if (!(param.description == null || param.description.isEmpty())) {
						if (printedSomething || verbose) {
							System.out.println(String.format("\t%-30s %s", "Description:", param.description));
						} else {
							System.out.print("\t");
							System.out.println(param.description);
						}
						printedSomething = true;

					}
					if (printedSomething) {
						System.out.println();
					}
				}
			}
			System.exit(0);
		}
	}
}
