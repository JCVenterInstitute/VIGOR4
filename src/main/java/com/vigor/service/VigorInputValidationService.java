package com.vigor.service;


import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.UnrecognizedOptionException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;


import com.vigor.exception.VigorException;
import com.vigor.forms.VigorForm;
import com.vigor.utils.VigorUtils;

@Service
public class VigorInputValidationService {

	private static final Logger LOGGER = LogManager.getLogger(VigorInputValidationService.class);
	

	@Autowired
	private VigorInitializationService vigorInitializationServiceObj;

	public void processInput(String... inputParams) {

		CommandLine inputs = validateInput(inputParams);
		try{
		if(!(inputs == null)){
		vigorInitializationServiceObj.initializeVigor(inputs);
		}
		}
		catch(Exception e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}

	}

	// Retrieve all the input parameters and validate
	public CommandLine validateInput(String[] inputParams) {
		CommandLineParser parser = new DefaultParser();
		CommandLine inputs = null;
		try {
			Options options = createOptions();
			inputs = parser.parse(options, inputParams);
			if(inputs.hasOption('i')){
				if (VigorUtils.isNullOrEmpty(inputs.getOptionValue('i'))) {
					VigorException.printExceptionMessage("File Does Not Exist");
				}
			}
			else
			{
				System.out.println("Input file is missing! Please provide input file.");
			}
			if(!(inputs.hasOption('o')))
			{
				System.out.println("Prefix for output file is missing!");
			}
			if (inputs.hasOption('F')) {
				int f = Integer.parseInt(inputs.getOptionValue('F'));
				if (!(f >= 0 && f < 3)) {
					VigorException.printExceptionMessage("Invalid value for Frameshift Sensitivity");
				}
			}
			if (inputs.hasOption('K')) {
				int k = Integer.parseInt(inputs.getOptionValue('K'));
				if (!(k >= 0 && k <= 1)) {
					VigorException.printExceptionMessage("Invalid value for Candidate Selection");
				}
			}

			if (inputs.hasOption('e')) {
				int e = Integer.parseInt(inputs.getOptionValue('e'));
				if (!(e > 0)) {
					VigorException.printExceptionMessage("E-value must be > 0 (-e " + e + ")");
				}
			}

		} catch (UnrecognizedOptionException e) {
			LOGGER.error(e.getMessage(), e);
			VigorException.printExceptionMessage(e.getMessage());

		} catch (MissingArgumentException e) {
			LOGGER.error(e.getMessage(), e);
			VigorException.printExceptionMessage(e.getMessage());

		} catch (ParseException e) {
			LOGGER.error(e.getMessage(), e);
			VigorException.printExceptionMessage(e.getMessage());
		}

		return inputs;

	}

	// CommandLine options and database related information to display to the user
	public void printHelp() {
		Options options = createOptions();

		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp("CommandLineOptions", options);
	}

	// Set the command line parameters into Options object
	public Options createOptions() {
		Options options = new Options();
		// Option h ;

		options.addOption("a", true,
				"auto-select the reference database, equivalent to '-d any ', default behavior unless overridden by -d or -G, (-A is a synonym for this option)");
		options.addOption("d", true,
				"<ref db>, specify the reference database to be used, (-D is a synonym for this option)");
		options.addOption("e", true,
				"<evalue>, override the default evalue used to identify potential genes, the default is usually 1E-5, but varies by reference database");
		options.addOption("c", true,
				"minimum coverage of reference product (0-100) required to report a gene, by default coverage is ignored");
		options.addOption("C", true,
				"complete (linear) genome (do not treat edges as gaps)");
		options.addOption("0", true, "complete circular genome (allows gene to span origin)");
		options.addOption("f", true,
				"<0, 1, or 2>, frameshift sensitivity, 0=ignore frameshifts, 1=normal (default), 2=sensitive");
		options.addOption("G", true,
				"<genbank file>, use a genbank file as the reference database, caution: VIGOR genbank parsing is fairly rudimentary and many genbank files are unparseable.  Partial genes will be ignored. Note: genbank files do not record enough information to handle RNA editing");
		options.addOption("i", "I", true,
				"<input fasta>, path to fasta file of genomic sequences to be annotated, (-I is a synonym for this option)");
		options.addOption("K", true, "<value>, value=0 skip candidate selection (default=1)");
		options.addOption("l", true, "do NOT use locus_tags in TBL file output (incompatible with -L)");
		options.addOption("L", true, "USE locus_tags in TBL file output (incompatible with -l)");
		options.addOption("o", "O", true,
				"<output prefix>, prefix for outputfile files, e.g. if the ouput prefix is /mydir/anno VIGOR will create output files /mydir/anno.tbl, /mydir/anno.stats, etc., (-O is a synonym for this option)");
		options.addOption("P", true,
				"<parameter=value~~...~~paramaeter=value>, override default values of VIGOR parameters");
		options.addOption("j", true,
				"turn off JCVI rules, JCVI rules treat gaps and ambiguity codes conservatively, use this option to relax these constraints and produce a more speculative annotation");
		options.addOption("m", true,
				"ignore reference match requirements (coverage/identity/similarity), sometimes useful when running VIGOR to evaluate raw contigs and rough draft sequences");
		options.addOption("s", true,
				"<gene size> minimum size (aa) of product required to report a gene, by default size is ignored");
		options.addOption("v", true, "verbose logging (default=terse)");
		options.addOption("x", true,
				"<ref_id,...,ref_id> comma separated list of reference sequence IDs to ignore (useful when debugging a reference database)");
		options.addOption("h", false, "For Help");

		return options;

	}

}
