package org.jcvi.vigor.service;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.tools.picocli.CommandLine;
import org.jcvi.vigor.utils.ConfigurationParameters;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Map;
import java.util.stream.Collectors;

@Service
public class VigorInputValidationService {

	private static final Logger LOGGER = LogManager.getLogger(VigorInputValidationService.class);

	@Autowired
	private VigorInitializationService vigorInitializationServiceObj;

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
						"  outputprefix.stats - run statistics (per genome sequence) in tab-delimited format",
						"  outputprefix.cds   - fasta file of predicted CDSs",
						"  outputprefix.pep   - fasta file of predicted proteins",
						"  outputprefix.tbl   - predicted features in GenBank tbl format",
						"  outputprefix.aln   - alignment of predicted protein to reference, and reference protein to genome",
						"  outputprefix.fs    - subset of aln report for those genes with potential sequencing issues",
						"  outputprefix.at    - potential sequencing issues in tab-delimited format"));

		MutuallyExclusiveGroup inputGroup = parser.addMutuallyExclusiveGroup().required(true);
		inputGroup.addArgument("-i","--input-fasta")
				  .action(Arguments.store())
				  .dest(CommandLineParameters.inputFile)
				  .metavar("<input fasta>")
				  .help("path to fasta file of genomic sequences to be annotated, (-I is a synonym for this option)");

		inputGroup.addArgument("-I")
				  .action(Arguments.store())
				  .dest(CommandLineParameters.inputFile)
				  .metavar("<input fasta>")
				  .help("synonym for -i/--input-fasta)");

		MutuallyExclusiveGroup outputGroup = parser.addMutuallyExclusiveGroup().required(true);
		outputGroup.addArgument("-o","--output-prefix")
				   .action(Arguments.store())
				   .dest(CommandLineParameters.outputPrefix)
				   .metavar("<output prefix>")
				   .help("prefix for outputfile files, e.g. if the ouput prefix is /mydir/anno VIGOR will create output files /mydir/anno.tbl, /mydir/anno.stats, etc., (-O is a synonym for this option)");

		outputGroup.addArgument("-O")
				   .action(Arguments.store())
				   .dest(CommandLineParameters.outputPrefix)
				   .metavar("<output prefix>")
				   .help("synonym for -o/--output-prefix");


		MutuallyExclusiveGroup referenceGroup = parser.addMutuallyExclusiveGroup("reference database");


		referenceGroup.addArgument("-a", "--autoselect-reference")
					  .dest(CommandLineParameters.referenceDB)
					  .action(Arguments.storeConst())
					  .setConst("any")
					  .help("auto-select the reference database, equivalent to '-d any ', default behavior unless overridden by -d or -G, (-A is a synonym for this option)");
		referenceGroup.addArgument("-A")
					  .dest(CommandLineParameters.referenceDB)
					  .action(Arguments.storeConst())
					  .setConst("any")
					  .help("synonym for -a/--autoselect-reference");

		referenceGroup.addArgument("-d", "--reference-database")
					  .dest(CommandLineParameters.referenceDB)
					  .action(Arguments.store())
					  .metavar("<ref db>")
					  .help("specify the reference database to be used, (-D is a synonym for this option)");
		referenceGroup.addArgument("-D")
					  .dest(CommandLineParameters.referenceDB)
					  .action(Arguments.store())
					  .metavar("<ref db>")
					  .help("synonym for -d/--reference-database");
		referenceGroup.addArgument("-G", "--genbank-reference")
					  .metavar("<genback file>")
					  .dest(CommandLineParameters.genbankDB)
					  .action(Arguments.store())
					  .help("use a genbank file as the reference database, caution: VIGOR genbank parsing is fairly rudimentary and many genbank files are unparseable.  Partial genes will be ignored. Note: genbank files do not record enough information to handle RNA editing");

		// TODO what are acceptable values for this
		parser.addArgument("-e", "--evalue")
			  .action(Arguments.store())
			  .type(Integer.class)
			  .dest(CommandLineParameters.eValue)
			  .help("<evalue>, override the default evalue used to identify potential genes, the default is usually 1E-5, but varies by reference database");

		// TODO add validation
		parser.addArgument("-c", "--min-coverage")
			  .dest(CommandLineParameters.minCoverage)
			  .action(Arguments.store())
			  .help("minimum coverage of reference product (0-100) required to report a gene, by default coverage is ignored");

		parser.addArgument("-C", "--complete")
				   .help("complete (linear) genome (do not treat edges as gaps)")
				   .dest(CommandLineParameters.completeGene)
				   .action(Arguments.storeTrue());
		parser.addArgument("-0", "--circular")
				   .dest(CommandLineParameters.circularGene)
				   .help("complete circular genome (allows gene to span origin)")
				   .action(Arguments.storeTrue());
		parser.addArgument("-f", "--frameshift-sensitivity")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.frameshiftSensitivity)
			  .choices("0","1","2")
			  .setDefault("1")
			  .help("frameshift sensitivity, 0=ignore frameshifts, 1=normal (default), 2=sensitive");

		parser.addArgument("-K", "--skip-candidate-selection")
			  .choices("0", "1")
			  .setDefault("1")
			  .dest(CommandLineParameters.skipSelection)
			  .metavar("<value>")
			  .help("value=0 skip candidate selection (default=1)");

		MutuallyExclusiveGroup locusGroup = parser.addMutuallyExclusiveGroup("locus tag usage");
		// use storeConst rather than storeTrue to avoid automatically setting a default value
		locusGroup.addArgument("-l", "--no-locus-tags")
				  .dest(CommandLineParameters.useLocusTags)
				  .action(Arguments.storeConst())
				  .setConst(true)
				  .help("do NOT use locus_tags in TBL file output (incompatible with -L)");
		locusGroup.addArgument("-L", "--locus-tags")
				  .dest(CommandLineParameters.useLocusTags)
				  .action(Arguments.storeConst())
				  .setConst(false)
				  .help("USE locus_tags in TBL file output (incompatible with -l)");

		parser.addArgument("-P","--parameter")
			  .action(Arguments.append())
			  .dest(CommandLineParameters.parameters)
			  .metavar("<parameter=value~~...~~parameter=value>")
			  .help("~~ separated list of VIGOR parameters to override default values");
		parser.addArgument("-j", "--jcvi-rules-off")
			  .action(Arguments.storeFalse())
			  .dest(CommandLineParameters.jcviRules)
			  .setDefault(true)
			  .help("turn off JCVI rules, JCVI rules treat gaps and ambiguity codes conservatively, use this option to relax these constraints and produce a more speculative annotation");
		parser.addArgument("-m","--ignore-reference-requirements")
			  .action(Arguments.storeTrue())
			  .dest(CommandLineParameters.ignoreReferenceRequirements)
			  .help("ignore reference match requirements (coverage/identity/similarity), sometimes useful when running VIGOR to evaluate raw contigs and rough draft sequences");
		parser.addArgument("-s","--min-gene-size")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.minGeneSize)
			  .type(Integer.class)
			  .metavar("<gene size>")
			  .help("minimum size (aa) of product required to report a gene, by default size is ignored");
		parser.addArgument("-v","--verbose")
			  .action(Arguments.storeTrue())
			  .dest(CommandLineParameters.verbose)
			  .help("verbose logging (default=terse)");

		parser.addArgument("-x", "--ignore-refID")
			  .action(Arguments.append())
			  .dest(CommandLineParameters.ignoreRefID)
			  .metavar("<ref_id,...,ref_id>")
			  .help("comma separated list of reference sequence IDs to ignore (useful when debugging a reference database)");

		parser.addArgument("--list-config-parameters")
			  .action(new ListConfigurations())
			  .dest(CommandLineParameters.listConfigParameters)
			  .help("list available configuration parameters and exit");

		parser.addArgument("--config-file")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.configFile)
			  .help("config file to use");

		parser.addArgument("--reference-database-path")
			  .action(Arguments.store())
			  .dest(CommandLineParameters.referenceDB_Path)
			  .help("reference database path");

		return parser;
	}

	private static class ListConfigurations implements ArgumentAction {

		@Override
		public void run(ArgumentParser argumentParser, Argument argument, Map<String, Object> map, String s, Object o) throws ArgumentParserException {
			System.out.println("Configuration Parameters\n");
			for (ConfigurationParameters param: Arrays.stream(ConfigurationParameters.values())
													  .sorted(Comparator.comparing(p-> p.configKey, String.CASE_INSENSITIVE_ORDER))
													  .collect(Collectors.toList())) {
				if (! param.hasFlag(ConfigurationParameters.Flags.VERSION_4)) {
					continue;
				}
				System.out.println(param.configKey);
				System.out.println();
				System.out.println(String.format("\t%-30s VIGOR_%s","Environment variable:",param.configKey.toUpperCase()));
				System.out.println(String.format("\t%-30s vigor.%s","System.property:", param.configKey));
				if (! (param.description == null || param.description.isEmpty()) ) {
					System.out.println(String.format("\t%-30s %s", "Description:", param.description));
				}
				System.out.println();
			}
			System.exit(0);
		}

		@Override
		public void onAttach(Argument argument) {

		}

		@Override
		public boolean consumeArgument() {
			return false;
		}
	}
}
