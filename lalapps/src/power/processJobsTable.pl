#!/usr/bin/perl -w
#-----------------------------------------------------------------------------------
# processJobsTable.pl 
#-----------------------------------------------------------------------------------
# -Processes records in the jobs table created by script 
#   createJobsTable.pl
# - If status is P - pending, it submits job to condor and changes 
#   the status to R - running
# - If status is R - running, it checks to see if the output XML file 
# is bigger than 3 KB. If it is, the status
# is changed to C - completed. If it exists but is smaller than 3 KB, the status
# is changed to U - user review required.  NF is returned if the output file is missing.
#
# Revision History
# July, 2003 - Dennis Mackin <dsmackin@stupidlinux.com>  - Created the original version of the script
# $Id$
#-----------------------------------------------------------------------------------
use strict;

#Add the program path and the users entire path to @INC (where perl searches for modules)
my $programPath = $0;
# now clean the filename off the path
$programPath =~ s/(\w+)\.pl$//;

use lib split(":", $ENV{'PATH'}), qw($programPath);

#load the power tools module which contains functions shared by the scripts
use power_tools;

my %STATUS = (
					P => "Pending",
					R => "Running",
					C => "Complete",
					U => "User Review Required",
					E => "Error",
					BC => "Bad Cache File",
					NF => "Output file not found");
					
# Check to make sure date of table is included as arg
my $USAGE = "\nusage: processJobsTable.pl  PARAMETERS_FILE YYMMDD-RUN_NUM \n\n";
if(! $ARGV[1])
{
	die $USAGE;
}

#get the date and runNum
# NEED TO ADD ERROR CHECKING HERE !!!!
my ($DATE, $RUN_NUM) = split '-', $ARGV[1];	
my $RUN_ID = $ARGV[1];

my $parametersFile = $ARGV[0];

my $params = f_parseParametersFile($parametersFile);

#-----------------------------------------------------------------------------------
#  MAIN
#-----------------------------------------------------------------------------------

my $startTime = time();

#set the path for all program output files
my $runPath = ${$params}{OUTPUT_PATH}  . "/$RUN_ID";

#build the parameter dependent globals
my $JOBS_TABLE = "$runPath/${$params}{'JOBS_TABLE'}";
my $CONDOR_SUBMIT_FILE = "$runPath/power_jobs.sub";

open LOG, ">>$runPath/" . ${$params}{'LOG'};

my $t = localtime();
print LOG "Starting $0 at: $t\n"; # $0 stores program name
print "Starting $0 at: $t\n";

#set up the condor submit file by writing executable line 
# and other statements relevant to all the jobs to be submitted
f_writeCondorSubmitFileHeaders( );

# now create the submit script using quality data file and the 
# the playground seconds hash array. The return code is one
# if jobs need to be submitted to condor. 0 otherwise
my $returnCode = f_processJobsTable ($JOBS_TABLE);

if ($returnCode){					
	f_submitJobs($CONDOR_SUBMIT_FILE);
}

# update the notebook page that give stats for this run
f_updateNotebookPage( 
					$params,
					$runPath,
					$RUN_ID);

$t = localtime();
print "$0 finished at $t\n";
print LOG "$0 finished at $t\n";
close LOG;		
									

#-----------------------------------------------------------------------------------
# f_processJobsTable
#-----------------------------------------------------------------------------------
#  - Reads through the records in the jobs table. if a record
#    is P then the job is added to the condor submit script
#   and the value of $submitCondor is set to 1 (true) to tell the
#   calling part of the code that there are jobs to submit to condor
#    
#-----------------------------------------------------------------------------------
# Returns $submitCondor 
#-----------------------------------------------------------------------------------
sub f_processJobsTable {

	my ($jobsTableFile) = @_;
	open JOBS_TABLE, $jobsTableFile
			or die "In f_processJobsTable: Couldn't open $jobsTableFile." ;	
				
	my $tmpTableFile = "$jobsTableFile.tmp";
	open TMP_TABLE, ">$tmpTableFile"
			or die "In f_processJobsTable: Couldn't open $tmpTableFile." ;	
	
	my $submitCondor = 0;
	
	while(<JOBS_TABLE>){
		chomp;
			
		#read in fields by splitting line on spaces
		my @fields = split "\t";
		
		if (scalar @fields != ${$params}{'JOBS_TABLE_FIELDS'}) {
			die "Expected " . ${$params}{'JOBS_TABLE_FIELDS'} . " fields in $JOBS_TABLE. Found " . scalar(@fields)  . ".\n";
		}
		
		#read columns into variables
		my ($statusCode, $startSec, $stopSec, $framecache, $outfile)  = ($fields[0],$fields[2],$fields[3],$fields[4],$fields[5]);

		if($statusCode eq "C" or $statusCode eq "E" or $statusCode eq "BC"){
			#nothing to do now; may add code later
		}elsif ($statusCode eq "P" ){
			f_writeJobToCondorSubmitFile($startSec,  $stopSec, $framecache, $outfile);
			$submitCondor = 1;
			$statusCode = "R";
		} elsif ($statusCode eq "R") { #Check output file for completion
			$statusCode = f_checkForProgramCompletion($outfile);
		} elsif ($statusCode eq "NF"){ #If there's no output, try re-running it
		
			#rebuild the cache file with more data. The extra data 
			# should get rid of the "ABORT: End of Frame Data" errors
			if(-f $framecache){unlink($framecache);}
			#if the rebuild failes mark it bad cache - Can't create the cachefile
			unless(f_buildCacheFile($startSec,  $stopSec + 5000,$framecache, ${$params}{'INSTRUMENT'})){
				$statusCode = "BC";
				next;
			}
			
			#Now add it to the .sub so it will run again
			f_writeJobToCondorSubmitFile($startSec,  $stopSec, $framecache, $outfile);
			$submitCondor = 1;
			$statusCode = "R";
		} else {
			print "Unknown status code: $statusCode\n";
			print LOG "Unknown status code: $statusCode\n";
		}
		my $record = "$statusCode\t$STATUS{$statusCode}\t$startSec\t$stopSec\t$framecache\t$outfile\n";
		print TMP_TABLE $record;
	}
	close TMP_TABLE;
	close JOBS_TABLE;
	
	system "cp $jobsTableFile $jobsTableFile.bak";
	rename($tmpTableFile, $jobsTableFile);
	return $submitCondor;
}

#-----------------------------------------------------------------------------------
#   f_writeCondorSubmitFileHeaders
#-----------------------------------------------------------------------------------
#  - Prints the general condor statements that apply to all
#     the jobs that are being submitted. e.g. Getenv which
#     tells condor to use the environment variables of the user
#     who submitted the job
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub f_writeCondorSubmitFileHeaders{

	#check for parameter CONDOR for backwards compatibility
	my $universe = "Vanilla";
	if(${$params}{'CONDOR_UNIVERSE'}) { 
		$universe = ${$params}{'CONDOR_UNIVERSE'};
	}
	
	my $stmts = << "STATEMENTS";
Getenv = True
Universe = $universe
Executable = ${$params}{'EXECUTABLE'}
STATEMENTS

	my $bak = $CONDOR_SUBMIT_FILE;
	my $i = 1;
	#make a backup of existing submit file
	if(-f $CONDOR_SUBMIT_FILE){
		#backups of submit are numbered, find next number
		my $bak = $CONDOR_SUBMIT_FILE . ".bak$i";
		while(-f $bak){
			$i++;
			$bak = $CONDOR_SUBMIT_FILE . ".bak$i";
		}
		#move exisiting submit file to backup
		rename $CONDOR_SUBMIT_FILE, $bak;
	}
	
	open ("CONDOR_SUB", ">$CONDOR_SUBMIT_FILE") or die "Couldn't open $CONDOR_SUBMIT_FILE.";
	print CONDOR_SUB $stmts;
	
	close CONDOR_SUB;
}


#-----------------------------------------------------------------------------------
#   f_checkForProgramCompletion
#-----------------------------------------------------------------------------------
#  - The best way to check to see if the lalapps_power has
#    completed is to check the XML output file. If it is complete 
#    and has records, then it is assumed the program has 
#    completed. If it is complete without records then the user
#    needs to review it to see if there is an error; a U is returned.
#    If the file does not exist, then it is assumed to still be 
#    running and an R is returned.
#-----------------------------------------------------------------------------------
#  Returns status code
#-----------------------------------------------------------------------------------
sub f_checkForProgramCompletion{
	my $outfile = shift;
	#if the file exists, try to parse it 
	if( -f $outfile){
		# Parsing XML file is better test but it
		#  is too slow. I wrapped this code in a sub
		# so it won't execute and then switched to a
		# faster -s test. dsmackin 07/30/03
		sub no_longer_used { 
			#use XML::DOM;
			#use XML::DOM::NodeList;
			
			#my $parser = new XML::DOM::Parser;
			#my $doc = $parser->parsefile ($outfile);
			
			#my $data = $doc->getElementsByTagName ("Stream")->item(0)->getFirstChild->getNodeValue;
			
			#clean the whitespace characters off the end of the string
			#$data =~ s/^\s*//;
			#$data =~ s/\s*$//;
			
			#put the rows into an array. Use the array length as the event count
			#my @events = split "\n", $data;
			#my $numEvents = scalar(@events);
			
			#print "$numEvents found for $outfile.\n";
			#$doc->dispose();
	
			#if($numEvents > 0 ) { #return status complete
			#	return "C";
			#} else { #Need the user to look into this job. There's probably an error.
			#	return "U";
			#}
		}#READ NOTE ABOVE !
		
		if( -s $outfile < 3000){ # then file seems too small; might be an error
			return "U";
		}else{
			return "C";
		}
	}
	# No output. Return NF for File not found
	return "NF";
}

#-----------------------------------------------------------------------------------
#   f_writeJobToCondorSubmitFile
#-----------------------------------------------------------------------------------
#  - takes $seed, $epoch, and $nset
#  - updates the dynamic parameters
#  - builds the condor command
#  - writes the command to the submit file 
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub f_writeJobToCondorSubmitFile {

	my ($startSec, $stopSec, $framecache, $outfile) = @_;
	my $duration = $stopSec - $startSec;
	my $npts = ${$params}{'SRATE'} * ${$params}{'T'};
	my $epoch = $startSec;	
	#number of segments = (2*number of seconds) - 1
	my $nseg = 2*($stopSec - $startSec) - 1;
	my $olap = $npts/2;
	my $numpts = $npts + ($npts - $olap)*( $nseg-1)+ 2*$olap;

	#replace spaces in comment w/ underscores
	${$params}{'COMMENT'} =~ s/ /_/g;

	# BUILD THE ARGS
	my $args = << "POWER_ARGS";
		--npts $npts
		--nseg $nseg
		--olap $olap
		--olapfctr ${$params}{'OLAPFCTR'}
		--minfbin ${$params}{'MINFBIN'}
		--mintbin ${$params}{'MINTBIN'}
		--flow ${$params}{'FLOW'}
		--delf ${$params}{'DELF'}
		--lngth ${$params}{'LNGTH'}
		--nsigma ${$params}{'NSIGMA'}
		--alphdef ${$params}{'ALPHADEF'}
		--segdcle ${$params}{'SEGDCLE'}
		--threshold ${$params}{'THRESHOLD'}
		--etomstr ${$params}{'ETOMSTR'}
		--channel ${$params}{'CHANNEL'}
		--framecache $framecache		
		--simtype ${$params}{'SIMTYPE'}
		--srate ${$params}{'SRATE'}
		--spectype ${$params}{'SPECTYPE'}
		--window ${$params}{'WINDOW'}
		--start_time $startSec
      --start_time_ns 0 
		--numpts $numpts
		--comment ${$params}{'COMMENT'}
		--dbglevel ${$params}{'DBGLEVEL'}
		--verbose
POWER_ARGS

	#clean the args string so that it is only one line
	$args =~ s/(\n|\t)/ /g;

	#create the condor command
	my $condorCmd =<< "ARGUMENTS";
		Arguments  = $args
		Output     = $runPath/out/$epoch-$duration.out
		Error      = $runPath/err/$epoch-$duration.err
		Log        = $runPath/log/$epoch-$duration.log
		Queue

ARGUMENTS

	open ("CONDOR_SUB", ">>$CONDOR_SUBMIT_FILE") or die "Couldn't open $CONDOR_SUBMIT_FILE.";
	print CONDOR_SUB $condorCmd;
	close CONDOR_SUB;
	
	return;
}

#-----------------------------------------------------------------------------------
#   f_submitJobs
#-----------------------------------------------------------------------------------
#  - makes a system call to submit the jobs to condor 
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub f_submitJobs {
	my $jobsSubmitScript = shift;
	
	#change directory to outpath before submitting
	chdir ${$params}{'OUTPUT_PATH'} . "/$RUN_ID/xml/" 
		or die "Couldn't change to dir ${$params}{'OUTPUT_PATH'}/$RUN_NUM/xml/ in f_submitJobs.\n" ;
	
	#submit the jobs
	system("condor_submit $jobsSubmitScript");
	
	#call rescedule to minimize delay before jobs are started
	system("/opt/condor/sbin/condor_reschedule");
}