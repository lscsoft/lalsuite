#!/usr/bin/perl -w
use strict;
#-----------------------------------------------------------------------------------
# processJobsTable.pl 
#-----------------------------------------------------------------------------------
# -Processes records in the jobs table created by script 
#   createJobsTable.pl
# - If status is P - pending, it submits job to condor and changes 
#   the status to R - running
# - If status is R - running, it checks to see if the output XML file 
# is complete with records. If it is complete with records, the status
# is changed to C - completed. If it is complete w/o records, the status
# is changed to U - user review required.  E is used for error conditions.
#
# Revision History
# July, 2003 - Dennis Mackin <dsmackin@stupidlinux.com>  - Created the original version of the script
# $Id$
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#  GLOBAL VARIABLES
#-----------------------------------------------------------------------------------
# Check to make sure date of table is included as arg
my $USAGE = "\nusage: countJobsTable.pl  YYMMDD \n\n";
if(! $ARGV[0])
{
	die $USAGE;
}
#my $DATE = f_getDateYYMMDD();
my $DATE =$ARGV[0];	

my %STATUS = (
					P => "Pending",
					R => "Running",
					C => "Complete",
					U => "User Review Required",
					E => "Error");

#my $TIME_CHUNK_SIZE = 64;
#my $TIME_CHUNK_MIN_SIZE = 16;

my $DATA_SET_NAME = "S2H1v02";
#my $DATA_QUALITY_FILE = "/home/dsmackin/lal/S2H1v02_segs-TEST.txt";
#my $PLAYGROUND_FILE  = "/home/dsmackin/lal/s2-playground-TEST.txt";

#path for CACHE_FILES
my $CACHE_PATH = "/home/dsmackin/tmp/cache";

#OUTPUT FILES
my $OUTPUT_PATH = "/home/dsmackin/lal/tests";
my $OUTPUT_FILE_ROOT  =  "search$DATE-EPOCH";
my $CONDOR_SUBMIT_FILE = "Search-$DATA_SET_NAME-$DATE.sub";
my $JOBS_TABLE = "/scratch/power/power_jobs_$DATE.tbl";

my $JOBS_TABLE_FIELDS = 6;

my $INSTRUMENT = "H";	
my $TYPE = "RDS_R_L1";

#-----------------------------------------------------------------------------------
# LALApps_power  PARAMETERS
#-----------------------------------------------------------------------------------

# $SEED is used for the random number generator used to create random noise.
#my $SEED = 4;

# number of points to analyze per second
# do not change this value
my $NPTS = 16384;

# Number of segments to analyze
my $NSEG = 16;

# do not change, must be $NPTS/2
my $OLAP =  $NPTS/2;

# leave these 3 values as they are
my $OLAPFCTR = 3;
my $MINFBIN = 2;
my $MINTBIN = 2;

#the lower end of the frequency spectrum to analyze; set to 60 for S2 data
my $FLOW = 60.0;

#do not change
my $DELF = 1.0;

#Upper limit on frequency.
# This value must be a power of 2
my $LNGTH = 512;

# Standard diviation; can change this value but make bigger not smaller
my $NSIGMA = 2.0;

#can change but not very relevant at this stage
my $ALPHADEF = 0.5;

# Segment duty cycle. Data should take $SEGDCLE segments at a go
my $SEGDCLE = 64;

# If the probability of an event occuring from background noise
# is less than the threshold, then an event is considered an Event
# trigger generator (ETG). This number should be between 1.0e-09 and
# 1.0e-05
my $THRESHOLD = "10.0e-15";

# Controls how many events get written out per segment
my $ETOMSTR = 10;

# The source of the data
my $CHANNEL = "H1:LSC-AS_Q";
my $FRAMECACHE = "/scratch/LALcaches/S2-RDS-LHO.cache	";

# does not currently do anything but keep it set to 0
my $SIMTYPE = 0;

# can be useMedian or useMean
my $SPECTYPE = "useMedian";

# Must be an integer. 2 works well.
my $WINDOW = 2;

# start time for data analysis
my $EPOCH = 729273613;

#numpts = (number of points in first segment) + (number of points in each offset)(number of offsets)
# Do not change this value
# NUMPTS is dynamically set when needed
##my $NUMPTS = $NPTS + ($NPTS - $OLAP)*($NSEG-1);

# If you want to create a file with the print spectrum, uncomment the line with "--printSpectrum"
#	If you don't want the print spectrum output file, uncomment the line without "--printSpectrum"
my $PRINT_SPECTRUM = "";
#my $PRINT_SPECTRUM = "--printSpectrum";

# Name of the XML output file
#my $OUTFILE = "REALSEARCH-$date-SEED_$SEED-GPS_$EPOCH-LENGTH_$LNGTH.xml";

#-----------------------------------------------------------------------------------
#  MAIN
#-----------------------------------------------------------------------------------

#set up the condor submit file by writing executable line 
# and other statements relevant to all the jobs to be submitted
f_writeCondorSubmitFileHeaders( );

#playgroundSeconds is a reference to a hash array that contains all the playground
# seconds.
#my $playgroundSeconds = lf_getPlaygroundSeconds($PLAYGROUND_FILE);
#my $countPlaygroundSecs =0;


# now create the submit script using quality data file and the 
# the playground seconds hash array
f_processJobsTable ($JOBS_TABLE);
									
f_submitJobs($CONDOR_SUBMIT_FILE);							
									

#-----------------------------------------------------------------------------------
# f_processJobsTable
#-----------------------------------------------------------------------------------
#  
#-----------------------------------------------------------------------------------
# Returns
#-----------------------------------------------------------------------------------
sub f_processJobsTable {

	my ($jobsTableFile) = @_;
	open JOBS_TABLE, $jobsTableFile
			or die "In f_processJobsTable: Couldn't open $jobsTableFile." ;	
				
	my $tmpTableFile = "jobsTable.tmp";
	open TMP_TABLE, ">$jobsTableFile"
			or die "In f_processJobsTable: Couldn't open $tmpTableFile." ;	
	
	while(<JOBS_TABLE>){
		chomp;
			
		#read in fields by splitting line on spaces
		my @fields = split "\t";
		
		if (scalar @fields != $JOBS_TABLE_FIELDS) {
			die "Expected $JOBS_TABLE_FIELDS fields in $JOBS_TABLE. Found " . scalar(@fields) .".\n";
		}
		
		#read columns into variables
		my ($statusCode, $startSec, $stopSec, $framecache, $outfile)  = ($fields[0],$fields[2],$fields[3],$fields[4],$fields[5]);

		if ($statusCode eq "P"){
			f_writeJobToCondorSubmitFile($startSec,  $stopSec, $framecache, $outfile);
			$statusCode = "R";
		} elsif ($statusCode eq "R") { #Check output file for completion
			$statusCode = f_checkForProgramCompletion($outfile);
		} else {
			die "Unknown status code: $statusCode\n";
		}
		
		 my $record = "$statusCode\t$STATUS{$statusCode}\t$startSec\t$stopSec\t$framecache\t$outfile\n";
		print TMP_TABLE $record;
	}
	close TMP_TABLE;
	close JOBS_TABLE;
	unlink "$jobsTableFile";
	rename $tmpTableFile, $jobsTableFile;
	return ;
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
	my $stmts = << "STATEMENTS";
			Executable = /home/dsmackin/bin/lalapps_power
			Getenv = True
STATEMENTS

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
		use XML::DOM;
		use XML::DOM::NodeList;
		
		my $parser = new XML::DOM::Parser;
		my $doc = $parser->parsefile ($outfile);
		
		my $data = $doc->getElementsByTagName ("Stream")->item(0)->getFirstChild->getNodeValue;
		
		#clean the whitespace characters off the end of the string
		$data =~ s/^\s*//;
		$data =~ s/\s*$//;
		
		#put the rows into an array. Use the array length as the event count
		my @events = split "\n", $data;
		my $numEvents = scalar(@events);
		
		print "$numEvents found for $outfile.\n";
		$doc->dispose();
		
		if($numEvents > 0 ) { #return status complete
			return "C";
		} else { #Need the user to look into this job. There's probably an error.
			return "U";
		}
	}

	# The jobs not complete yet. Return R because it's still running
	return "R";
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
	
	#number of segments = (2*number of seconds) - 1
	my $epoch = $startSec;
	my $nseg = 2*($stopSec - $startSec) - 1;
	my $numpts = $NPTS + ($NPTS - $OLAP)*( $nseg-1);
	my $condorCmd;

	# BUILD THE ARGS
	my $args = << "POWER_ARGS";
		--npts $NPTS
		--nseg $nseg
		--olap $OLAP
		--olapfctr $OLAPFCTR
		--minfbin $MINFBIN
		--mintbin $MINTBIN
		--flow $FLOW
		--delf $DELF
		--lngth $LNGTH
		--nsigma $NSIGMA
		--alphdef $ALPHADEF
		--segdcle $SEGDCLE
		--threshold $THRESHOLD
		--etomstr $ETOMSTR
		--channel $CHANNEL
		--framecache $framecache		
		--simtype $SIMTYPE
		--spectype $SPECTYPE
		--window $WINDOW
		--epoch $epoch 0
		--numpts $numpts
		$PRINT_SPECTRUM
		--outfile  $outfile
POWER_ARGS

	#if cmd is to include sine wave, add it to the args
	#if($INCLUDE_SINE_WAVE){ $args .= " --sine $FREQ $DELTA $AMPLITUDE $WIDTH";}

	#clean the args string so that it is only one line
	$args =~ s/(\n|\t)/ /g;

	#create the condor command
	$condorCmd =<< "ARGUMENTS";
		Arguments      = $args
		Output     = $OUTPUT_PATH/out/$OUTPUT_FILE_ROOT$epoch.out
		Error      = $OUTPUT_PATH/err/$OUTPUT_FILE_ROOT$epoch.err
		Log        = $OUTPUT_PATH/log/$OUTPUT_FILE_ROOT$epoch.log
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
	
	#submit the jobs
	system("condor_submit $jobsSubmitScript");
	
	#call rescedule to minimize delay before jobs are started
	system("/opt/condor/sbin/condor_reschedule");
}

#-----------------------------------------------------------------------------------
#  f_getDateYYMMDD
#-----------------------------------------------------------------------------------
#  Returns date as YYMMDD
#-----------------------------------------------------------------------------------
sub f_getDateYYMMDD{
    my %month = (
                Jan => '01',
                Feb => '02',
                Mar => '03',
                Apr => '04',
                May => '05',
                Jun => '06',
                Jul => '07',
                Aug => '08',
                Sep => '09',
                Oct => '10',
                Nov => '11',
                Dec => '12'
                );
    my @timeparts = split(" ", localtime());

    #add in leading 0 for day if needed
    if ($timeparts[2]<10) {$timeparts[2] = "0" . $timeparts[2];};

	 #change YYYY to YY
	 $timeparts[4] =~ s/^20//;

    return "$timeparts[4]$month{$timeparts[1]}$timeparts[2]";
}