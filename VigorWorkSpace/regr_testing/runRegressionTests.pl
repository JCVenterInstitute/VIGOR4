#!/usr/local/bin/perl

# File: runRegressionTests.pl
# Author: Paolo Amedeo
# Created: November,24 2014
#
# $Author: $
# $Date:  $
# $Revision:  $
# $HeadURL: $

#######################################################################################
#
# Copyright (c) 2009 - 2015 J. Craig Venter Institute.
#   This file is part of JCVI VIGOR
# 
#   JCVI VIGOR is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#   
#   JCVI VIGOR is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with JCVI VIGOR.  If not, see <http://www.gnu.org/licenses/>.
# 
# Contributors:
#     Shiliang Wang - Initial idea and implementation.
#     Jeff Hoover - Redesigning, refactoring, and expanding the scope.
#     Susmita Shrivastava and Neha Gupta - Creation and curation of sequence databases.
#     Paolo Amedeo and Danny Katzel - Maintenance and further improvements.
#
#######################################################################################

# runRegressionTests.pl is a a wrapper script that runs all the available regression tests for Vigor

=head1 NAME
    
    runRegressionTests.pl
    
=head1 USAGE

    runRegressionTests.pl [-]-vigor_test_dir <vigor_test> [-]-o[utput_file] <report>  [options]

=head1 REQUIRED ARGUMENTS

=over

=item [-]-vigor_test_dir <vigor_test>

Path to the new version of Vigor to be tested

=for Euclid:
    vigor_test.type: string

=item [-]-o[utput_file] <report> 

Name for the report file summarizing all the tests results.

=for Euclid:
    report.type: writeable

=back

=head1 OPTIONS

=over

=item [-]-vigor_prod_dir <path_to_production>

Path to the current production VIGOR.pl (default: /usr/local/devel/VIRIFX/software/VIGOR3/prod3)

=for Euclid:
    path_to_production.type:     string
    path_to_production.default: '/usr/local/devel/VIRIFX/software/VIGOR3/prod3'

=item [-]-test_area <path_to_output>

Area where the output files are written (default: vigorscratch/RegrTest)

=for Euclid:
    path_to_output.type:     string
    path_to_output.default: "$::vigor_scratch/RegrTest"

=item [-]-debug [<log_level>]

Logging threshold for the error log file.
Valid values: TRACE, DEBUG, INFO, WARN, ERROR, and FATAL
Not specifying any value, it will default to 'DEBUG'

=for Euclid:
    log_level.type:        string
    log_level.default:     'NOT_SET'
    log_level.opt_default: 'DEBUG'
    
=item [-]-log[_file] <error_log_file>

Local logging file, used to monitor the process.

=for Euclid:
    error_log_file.type:   writeable

=item --help

    Prints this documentation and quit
    
=back

=head1 DESCRIPTION

This is just a working template for creating scripts. If you don't replace the documentation with your real documentation, you'll entartain the users seeking help.

=cut

BEGIN {
    use Cwd (qw(abs_path getcwd));
    use FindBin;
    use ProcessingObjects::SafeIO;
    $::vigor_scratch = "$FindBin::Bin/../prod3/vigorscratch";
    
    unless (-e $::vigor_scratch) {
        mk_tree_safe($::vigor_scratch, 0777) || die "Impossible to create vigorscratch directory ($::vigor_scratch).";
    }
    $::cmd = join(' ', $0, @ARGV);
    $::working_dir = getcwd();
}

use strict;
use warnings;
use lib ("$FindBin::Bin/../prod3");
use Getopt::Euclid 0.2.4 (qw(:vars));
## Commonly used modules (remove whatever doesn't apply):
#use Data::Dumper;
use File::Basename;
use File::Path;
use Log::Log4perl (qw(get_logger));
use L4pTools;

## Constants declaration
#
use constant VIGOR_EXECUTABLE  => 'VIGOR3.pl';

our ($ARGV_vigor_test_dir, $ARGV_output_file, $ARGV_vigor_prod_dir, $ARGV_test_area, $ARGV_debug, $ARGV_log_file);

my $jlt = L4pTools->init(ARGV_log_file => $ARGV_log_file, ARGV_debug => $ARGV_debug);
my $logger = get_logger(basename($0));

$logger->info("Command line: $::cmd\nInitial working directory: $::working_dir\nDebug level: \"$ARGV_debug\"");

$ARGV_vigor_test_dir =~ s#/$##;
$ARGV_vigor_prod_dir =~ s#/$##;
my ($exec_prod, $exec_test) = ($ARGV_vigor_prod_dir .'/'. VIGOR_EXECUTABLE, $ARGV_vigor_test_dir .'/'. VIGOR_EXECUTABLE);
my ($same, $diff, $error) = (0) x 3;

unless (-d $ARGV_vigor_test_dir && -x $exec_test) {
    $logger->logdie("Unable to find the Vigor test directory \"$ARGV_vigor_test_dir\" or \"$exec_test\" is not executable.");
}
unless (-d $ARGV_vigor_prod_dir && -x $exec_prod) {
    $logger->logdie("Unable to find the Vigor production directory \"$ARGV_vigor_prod_dir\" or \"$exec_prod\" is not executable.");
}

unless (-d $ARGV_test_area) {
    mk_tree_safe($ARGV_test_area, 0777) || $logger->logdie("Impossible to create the test directory (\"$ARGV_test_area\").");
}

open(my $report, ">$ARGV_output_file") || $logger->logdie("Impossible to open the output report file \"$ARGV_output_file\".");

## Hash containing all the parameters for all the tests

my %vigor_params = (
    flu   => {i => "$FindBin::Bin/flu.fasta"},
    gcv   => {i => "$FindBin::Bin/gcv.fasta"},
    hadv  => {i => "$FindBin::Bin/hadv.fasta"},
    hhv4  => {i => "$FindBin::Bin/hhv4.fasta",
              0 => ''},
    hpiv1 => {i => "$FindBin::Bin/hpiv1.fasta"},
    hpiv3 => {i => "$FindBin::Bin/hpiv3.fasta"},
    hrv   => {i => "$FindBin::Bin/hrv.fasta"},
    mmp   => {i => "$FindBin::Bin/mmp.fasta"},
    mpv   => {i => "$FindBin::Bin/mpv.fasta"},
    msl   => {i => "$FindBin::Bin/msl.fasta"},
    norv  => {i => "$FindBin::Bin/norv.fasta"},
    rbl   => {i => "$FindBin::Bin/rbl.fasta"},
    rsv   => {i => "$FindBin::Bin/rsv.fasta"},
    rtv   => {i => "$FindBin::Bin/rtv.fasta"},
    var   => {i => "$FindBin::Bin/var.fasta"},
    veev  => {i => "$FindBin::Bin/veev.fasta"},
    yfv   => {i => "$FindBin::Bin/yfv.fasta"});

## Now performongthe tests

while (my ($db, $data) = each(%vigor_params)) {
    my $msg = "Now processing database \"$db\".";
    $logger->info($msg);
    print "\n$msg\n";
    print {$report} "$msg\n";
    my ($prod_dir, $test_dir) = ("$ARGV_test_area/$db/PRODUCTION", "$ARGV_test_area/$db/TEST");
    
    unless (-d "$ARGV_test_area/$db" && -d $prod_dir) {
        unless (mk_tree_safe($prod_dir, 0777)) {
            $msg = "Impossible to create the test directory \"$prod_dir\". - Skipping this test";
            $logger->error($msg);
            print {$report} "$msg\n";
            ++$error;
            &endReport();
            next;
        }
    }
    unless (-d $test_dir) {
        unless (mk_tree_safe($test_dir, 0777)) {
            $msg = "Impossible to create the test directory \"$test_dir\". - Skipping this test";
            $logger->error($msg);
            print {$report} "$msg\n";
            ++$error;
            &endReport();
            next;
        }
    }
    ## Running the production version of Vigor
    
    my $cmd = "$exec_prod -D $db -O $prod_dir/$db";
    
    while (my ($opt, $val) = each(%{$data})) {
        $cmd .= " -$opt $val";
    }
    $cmd .= " > $prod_dir/$db.log";
    $logger->info("CMD: \"$cmd\"");
    
    if (system($cmd)) {
        $msg = "Problems with running the production vesrion of Vigor on database \"$db\" - Skipping the remainder of the test on this database.";
        $logger->error($msg);
        print {$report} "$msg\n";
        ++$error;
        &endReport();
        next;
    }
    ## Running the production version of Vigor
    
    $cmd = "$exec_test -D $db -O $test_dir/$db";
    
    while (my ($opt, $val) = each(%{$data})) {
        $cmd .= " -$opt $val";
    }
    $cmd .= " > $test_dir/$db.log";
    $logger->info("CMD: \"$cmd\"");
    
    if (system($cmd)) {
        $msg = "Problems with running the test vesrion of Vigor on database \"$db\" - Skipping the evaluation of the results of this test.";
        $logger->error($msg);
        print {$report} "$msg\n";
        ++$error;
        &endReport();
        next;
    }
    ## Removing the parameter listing portion of the .rpt files (which is always different in at least date, executable, and temp space)
    
    foreach my $file ("$prod_dir/$db.rpt",  "$test_dir/$db.rpt") {
        my $shorten = $file . '_short';
        open(my $rpt, $file) || $logger->logdie("Impossible to open the file \"$file\" for reading.");
        open(my $short, ">$shorten") || $logger->logdie("Impossible to open the file \"$shorten\" for writing.");
        my $print = 0;
        
        while (<$rpt>) {
            if (/^Sequence/) {
                ++$print;
            }
            if ($print) {
                print {$short} $_;
            }
        }
        close($rpt);
        close($short);
    }
    my $diff_output = `diff $prod_dir/$db.rpt_short $test_dir/$db.rpt_short`;
    
    if ($diff_output =~ /\S/) {
        $msg = "The two versions produce different output files.";
        print "$msg\n\n";
        print {$report} "$msg\n\n", '-'x80, $diff_output, '-'x80, "\n\n";
        ++$diff;
    }
    else {
        $msg = "The two versions produce the same output.";
        print "$msg\n\n";
        print {$report} "$msg\n\n";
        ++$same;
    }
    &endReport();
}
print "\n\nDone.\n\n$same identical results.\n$diff different results.\n$error errors.\n\n";

sub endReport {
    print {$report} '=' x 80, "\n\n";
}
