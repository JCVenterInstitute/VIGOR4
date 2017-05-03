package com.vigor.service;


import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;
import com.vigor.exception.VigorException;
import com.vigor.utils.VigorUtils;

@Service
public class VigorInputValidationService {

    private static final Logger LOGGER = LogManager.getLogger ( VigorInputValidationService.class );


    @Autowired
    private VigorInitializationService vigorInitializationServiceObj;

    public void processInput(String... inputParams) {
        CommandLine inputs = null;
        if (inputParams.length >= 1) {
            inputs = validateInput ( inputParams );
        } else {
            System.out.println ( "Missing input parameters" );
            printHelp ();
        }
        try {
            if (!(inputs == null)) {
                vigorInitializationServiceObj.initializeVigor ( inputs );
            }
        } catch (Exception e) {
            LOGGER.error ( e.getMessage (), e );
        }

    }

    // Retrieve all the input parameters and validate
    public CommandLine validateInput(String[] inputParams) {
        CommandLineParser parser = new DefaultParser ();
        CommandLine inputs = null;

        try {
            Options options = createOptions(false);

            OptionGroup synonymOptionGroup1 = new OptionGroup ();
            synonymOptionGroup1.setRequired ( true );
            synonymOptionGroup1.addOption ( options.getOption ( "I" ) );
            synonymOptionGroup1.addOption ( options.getOption ( "i" ) );
            OptionGroup synonymOptionGroup2 = new OptionGroup();
            synonymOptionGroup2.setRequired ( true );
            synonymOptionGroup2.addOption ( options.getOption ( "o" ) );
            synonymOptionGroup2.addOption ( options.getOption ( "O" ) );
            OptionGroup synonymOptionGroup3 = new OptionGroup();
            synonymOptionGroup3.setRequired ( false );
            synonymOptionGroup3.addOption ( options.getOption ( "d" ) );
            synonymOptionGroup3.addOption ( options.getOption ( "D" ) );
            options.addOptionGroup ( synonymOptionGroup1 );
            options.addOptionGroup ( synonymOptionGroup2 );
            options.addOptionGroup ( synonymOptionGroup3 );
            OptionGroup conflictOptionGroup1 = new OptionGroup ();
            conflictOptionGroup1.addOption ( options.getOption ( "l" ) );
            conflictOptionGroup1.addOption ( options.getOption ( "L" ) );
            options.addOptionGroup ( conflictOptionGroup1 );

            inputs = parser.parse ( options, inputParams );


            if (inputs.hasOption ( 'F' )) {
                int f = Integer.parseInt ( inputs.getOptionValue ( 'F' ) );
                if (!(f >= 0 && f < 3)) {
                    VigorException.printExceptionMessage ( "Invalid value for Frameshift Sensitivity" );
                    System.exit ( 0 );
                }
            }
            if (inputs.hasOption ( 'K' )) {
                int k = Integer.parseInt ( inputs.getOptionValue ( 'K' ) );
                if (!(k >= 0 && k <= 1)) {
                    VigorException.printExceptionMessage ( "Invalid value for Candidate Selection" );
                    System.exit ( 0 );
                }
            }

            if (inputs.hasOption ( 'e' )) {
                int e = Integer.parseInt ( inputs.getOptionValue ( 'e' ) );
                if (!(e > 0)) {
                    VigorException.printExceptionMessage ( "E-value must be > 0 (-e " + e + ")" );
                    System.exit ( 0 );
                }
            }

        } catch(AlreadySelectedException e){
            LOGGER.error(e.getMessage (),e);
            VigorException.printExceptionMessage ( e.getMessage ()+". Please select either of the two options");

        }
        catch(MissingOptionException e)
        {
            LOGGER.error(e.getMessage (),e);
            VigorException.printExceptionMessage ( e.getMessage () );

        }

        catch (UnrecognizedOptionException e) {
            LOGGER.error ( e.getMessage (), e );
            VigorException.printExceptionMessage ( e.getMessage () );

        } catch (MissingArgumentException e) {
            LOGGER.error ( e.getMessage (), e );
            VigorException.printExceptionMessage ( e.getMessage () );

        } catch (ParseException e) {
            LOGGER.error ( e.getMessage (), e );
            VigorException.printExceptionMessage ( e.getMessage () );
        }

        return inputs;

    }

    // CommandLine options and database related information to display to the user
    public void printHelp() {
        Options options = createOptions ( true );

        HelpFormatter formatter = new HelpFormatter ();
        formatter.printHelp ( "CommandLineOptions", options );

    }

    // Set the command line parameters into Options object
    public Options createOptions(boolean helpOptions) {
        Options helperOptions = new Options ();
        Options allOptions = new Options();
        // Option h ;
        Option option1 = Option.builder ( "a" ).hasArg ().desc ( "auto-select the reference database, equivalent to '-d any ', default behavior unless overridden by -d or -G, (-A is a synonym for this option)" ).build ();
        Option option2 = Option.builder ( "d" ).hasArg ().desc (
                "<ref db>, specify the reference database to be used, (-D is a synonym for this option)" ).build ();
        Option option2_1= Option.builder ( "D" ).hasArg ().desc (
                "<ref db>, specify the reference database to be used, (-d is a synonym for this option)" ).build ();
        Option option3 = Option.builder ( "e" ).hasArg ().desc (
                "<evalue>, override the default evalue used to identify potential genes, the default is usually 1E-5, but varies by reference database" ).build ();
        Option option4 = Option.builder ( "c" ).hasArg ().desc (
                "minimum coverage of reference product (0-100) required to report a gene, by default coverage is ignored" ).build ();
        Option option5 = Option.builder ( "C" ).hasArg ().desc (
                "complete (linear) genome (do not treat edges as gaps)" ).build ();
        Option option6 = Option.builder ( "0" ).hasArg ().desc ( "complete circular genome (allows gene to span origin)" ).build ();
        Option option7 = Option.builder ( "f" ).hasArg ().desc (
                "<0, 1, or 2>, frameshift sensitivity, 0=ignore frameshifts, 1=normal (default), 2=sensitive" ).build ();
        Option option8 = Option.builder ( "G" ).hasArg ().desc (
                "<genbank file>, use a genbank file as the reference database, caution: VIGOR genbank parsing is fairly rudimentary and many genbank files are unparseable.  Partial genes will be ignored. Note: genbank files do not record enough information to handle RNA editing" ).build ();
        Option option9 = Option.builder ( "i" ).hasArg ().desc (
                "<input fasta>, path to fasta file of genomic sequences to be annotated, (-I is a synonym for this option)" ).build ();
        Option option9_1 = Option.builder ( "I" ).hasArg ().desc (
                "<input fasta>, path to fasta file of genomic sequences to be annotated, (-i is a synonym for this option)" ).build ();
        Option option10 = Option.builder ( "K" ).hasArg ().desc ( "<value>, value=0 skip candidate selection (default=1)" ).build ();
        Option option11 = Option.builder ( "l" ).hasArg ().desc ( "do NOT use locus_tags in TBL file output (incompatible with -L)" ).build ();
        Option option12 = Option.builder ( "L" ).hasArg ().desc ( "USE locus_tags in TBL file output (incompatible with -l)" ).build ();
        Option option13 = Option.builder ( "o" ).hasArg ().desc (
                "<output prefix>, prefix for outputfile files, e.g. if the ouput prefix is /mydir/anno VIGOR will create output files /mydir/anno.tbl, /mydir/anno.stats, etc., (-O is a synonym for this option)" ).build ();
        Option option13_1 = Option.builder ( "O" ).hasArg ().desc (
                "<output prefix>, prefix for outputfile files, e.g. if the ouput prefix is /mydir/anno VIGOR will create output files /mydir/anno.tbl, /mydir/anno.stats, etc., (-o is a synonym for this option)" ).build ();
        Option option14 = Option.builder ( "P" ).hasArg ().desc (
                "<parameter=value~~...~~paramaeter=value>, override default values of VIGOR parameters" ).build ();
        Option option15 = Option.builder ( "j" ).hasArg ().desc (
                "turn off JCVI rules, JCVI rules treat gaps and ambiguity codes conservatively, use this option to relax these constraints and produce a more speculative annotation" ).build ();
        Option option16 = Option.builder ( "m" ).hasArg ().desc (
                "ignore reference match requirements (coverage/identity/similarity), sometimes useful when running VIGOR to evaluate raw contigs and rough draft sequences" ).build ();
        Option option17 = Option.builder ( "s" ).hasArg ().desc (
                "<gene size> minimum size (aa) of product required to report a gene, by default size is ignored" ).build ();
        Option option18 = Option.builder ( "v" ).hasArg ().desc ( "verbose logging (default=terse)" ).build ();
        Option option19 = Option.builder ( "x" ).hasArg ().desc (
                "<ref_id,...,ref_id> comma separated list of reference sequence IDs to ignore (useful when debugging a reference database)" ).build ();
        Option option20 = Option.builder ( "h" ).hasArg ().desc ( "For Help" ).build ();
        helperOptions.addOption ( option1 ).addOption ( option2 ).addOption ( option3 ).addOption ( option4 ).addOption ( option5 ).addOption ( option6 ).
                addOption ( option7 ).addOption ( option8 ).addOption ( option9 ).addOption ( option10 ).addOption ( option11 ).addOption ( option12 )
                     .addOption ( option13 ).addOption ( option14 ).addOption ( option15 ).addOption ( option16 ).addOption ( option17 ).addOption ( option18 ).
                             addOption ( option19 ).addOption ( option20 );
        allOptions.addOption ( option1 ).addOption ( option2 ).addOption ( option3 ).addOption ( option4 ).addOption ( option5 ).addOption ( option6 ).
                addOption ( option7 ).addOption ( option8 ).addOption ( option9 ).addOption ( option10 ).addOption ( option11 ).addOption ( option12 )
                     .addOption ( option13 ).addOption ( option14 ).addOption ( option15 ).addOption ( option16 ).addOption ( option17 ).addOption ( option18 ).
                             addOption ( option19 ).addOption ( option20 ).addOption ( option2_1 ).addOption ( option9_1 ).addOption ( option13_1 );
       if(helpOptions) {
           return helperOptions;}

       else {return allOptions;}

    }

}
