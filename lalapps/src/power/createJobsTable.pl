#!/usr/bin/perl -w
use strict;
#-----------------------------------------------------------------------------------
# createJobsTable.pl 
#-----------------------------------------------------------------------------------
#  -Takes a data quality file and playground data file as inputs.
#  -Creates a list of chunks of time chunks that have BAD_CODE==0
#   and are part of the playground data.
# 
# - This script works in conjunction with runJobsInJobTable.pl which
#   submits the jobs from the table to condor and checks on their progress
#   
# *****WARNING *************
#  This script does not properly overlap chunks of time. When a chunk of time starts next to another 
# chunk of time, they should overlap by one second.
# *****END OF WARNING *****
#
# Revision History
# July, 2003 - Dennis Mackin <dsmackin@stupidlinux.com>  - Created the original version of the script
# $Id$
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#  GLOBAL VARIABLES
#-----------------------------------------------------------------------------------

# use f_getDateYYMMDD to get current date
my $DATE = f_getDateYYMMDD();

my $TIME_CHUNK_SIZE = 600;
my $TIME_CHUNK_MIN_SIZE = 16;

my $DATA_SET_NAME = "S2H1v02";
my $DATA_QUALITY_FILE = "/scratch/power/input/$DATA_SET_NAME" . "_segs.txt";

#To allow data outside of playground, set USE_PLAYGROUND to 0
my $USE_PLAYGROUND =1;
my $PLAYGROUND_FILE  = "/scratch/power/input/s2-playground.txt";

#path for CACHE_FILES
my $CACHE_PATH = "/scratch/power/cache";

#OUTPUT FILES
my $OUTPUT_PATH = "/scratch/power/tests";
my $OUTPUT_FILE_ROOT  =  "search$DATE-EPOCH";
my $JOBS_TABLE = "power_jobs_$DATE.tbl";
my $LOG = "create.log";

my $INSTRUMENT = "H";	
my $TYPE = "RDS_R_L1";


#-----------------------------------------------------------------------------------
#  MAIN
#-----------------------------------------------------------------------------------


my $startTime = time();
#get the next runNumber as string NNN
my $runNum = f_getRunNumber($OUTPUT_PATH,$DATE);

#call function to setup output dirs
 f_setupOutputDirs($OUTPUT_PATH,$DATE, $runNum);
 
 #set the path for all program output files
 my $runPath = "$OUTPUT_PATH/$DATE-$runNum";

open LOG, ">$runPath/$LOG";
my $t = localtime();
print LOG "Start time: $t\n";
 
#if the jobs table exists, move it to a backup file
if (-f "$runPath/$JOBS_TABLE"){ rename "$runPath/$JOBS_TABLE", "$runPath/$JOBS_TABLE.bak";}

#playgroundSeconds is a reference to a hash array that contains all the playground
# seconds.
my $playgroundSeconds = lf_getPlaygroundSeconds($PLAYGROUND_FILE);

# For debug purposes, count and display the number of playground seconds
my $countPlaygroundSecs =0;
foreach (keys %{$playgroundSeconds}){
	$countPlaygroundSecs++;
}
print "The Playground is $countPlaygroundSecs seconds.\n\n";

# now create the submit script using quality data file and the 
# the playground seconds hash array
f_createJobsForGoodData (
									$runPath,
									$DATA_QUALITY_FILE, 
									$TIME_CHUNK_MIN_SIZE, 
									$TIME_CHUNK_SIZE,
									$playgroundSeconds);

my $totalTime = time() - $startTime;
print "Completed in 	$totalTime seconds.\n";
print LOG "Completed in $totalTime seconds\n.";
print "The jobs table for run $DATE-$runNum has been created. \n\n".
			"To submit the jobs to condor, run:\n\n processJobsTable.pl $DATE $runNum\n\n";

#-----------------------------------------------------------------------------------
#  f_createJobsForGoodData
#-----------------------------------------------------------------------------------
#  - Reads the data quality file
# - throws away bad chunks
#  - check each good chunk to see if it contains playground data
#-----------------------------------------------------------------------------------
#  lists chunks of good data from the  data quality file that are in the playground
#-----------------------------------------------------------------------------------
sub f_createJobsForGoodData {

	my ($runPath, $dataQualityFile, $minChunkSize, $maxChunkSize, $playgroundSecondsHash) = @_;
	my ($chunkStart, $chunkStop);
	
	open DATA_QUALITY_FILE, $dataQualityFile
			or die "In f_createJobsForGoodData: Couldn't open $dataQualityFile." ;	
	
	while(<DATA_QUALITY_FILE>){
		chomp;
			
		#read in fields by splitting line on spaces
		my @fields = split " ";
		
		#read columns into variables
		my ($startSecond, $stopSecond, $badCode)  = ($fields[1],$fields[2],$fields[4]);

		#if line is a bad chunk, skip it
		if ($badCode) { next;}
		
		#set the beging of the first chunk to analyze to the start second of the chunk from the
		# data quality file
		$chunkStart = $startSecond;
		
		#advance chunkStart to the first second in playground
		while( $USE_PLAYGROUND and not exists(${$playgroundSecondsHash}{$chunkStart}) and $chunkStart < $stopSecond){ 
			$chunkStart++;
		}
		
		#Loop through data quality file chunk
		for($chunkStop=$chunkStart; $chunkStop < $stopSecond; $chunkStop++){
			#print "chunksize = ", $chunkStop - $chunkStart , "\n";
			if(($USE_PLAYGROUND and not exists ${$playgroundSecondsHash}{$chunkStop})  or  ($chunkStop - $chunkStart)  >= $maxChunkSize ){
				
				f_processChunk($runPath,$chunkStart, $chunkStop);
				
				#fast forward to the next good second 
				while(($USE_PLAYGROUND and not exists ${$playgroundSecondsHash}{$chunkStop})and $chunkStop < $stopSecond){ 
					$chunkStop++;
				}
				
				$chunkStart = $chunkStop;
			}
		}		
		
		#process last chunk in data segment
		#print "Processing last chunk.\n";
		f_processChunk($runPath,$chunkStart, $chunkStop);
	}
	
	close DATA_QUALITY_FILE;
	return ;
}

#-----------------------------------------------------------------------------------
#   f_processChunk
#-----------------------------------------------------------------------------------
#  - Performs processing on a chunk once it has been identified.
# 
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub  f_processChunk {
	my ($runPath,$startSec, $stopSec) = @_;
	
	my $duration = $stopSec - $startSec;
	#print "Processing chunk $startSec, $stopSec.\n";
	#make sure chunk is larger than minimum size
	if($stopSec - $startSec > $TIME_CHUNK_MIN_SIZE){
	
		#build output file string
		my $outfile = "$runPath/xml/$OUTPUT_FILE_ROOT$startSec-$stopSec.xml";
		my $framecache = "$CACHE_PATH/$INSTRUMENT-$startSec-$duration";
		
		f_findData($startSec,  $stopSec,$framecache);
		f_writeJobsTable($runPath,$startSec,  $stopSec, $framecache, $outfile, $JOBS_TABLE);
	}
}

#-----------------------------------------------------------------------------------
#   f_findData
#-----------------------------------------------------------------------------------
#  - Takes the stop and start
# 
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub f_findData {
	my ($startSec, $stopSec, $framecache) = @_;
	
	#only call LALdataFind if the cache file doesn't currently exist.
	if (! -f $framecache){
		my $cmd =  "LALdataFind --lal-cache --instrument $INSTRUMENT --type RDS_R_L1 " .
		 		" --start $startSec --end $stopSec > $framecache";
		print "$cmd\n";
		system $cmd;
	}
}

#-----------------------------------------------------------------------------------
#   f_writeJobsTable
#-----------------------------------------------------------------------------------
#  - Creates a file based table that contains the following fields
#			--  run status code (use P for pending)
#			--  run status description
#			--  start second 
#			--  finish second 
#			--  cachfile
#			--  outputfile XML file
# 
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub  f_writeJobsTable {
	my ($runPath,$startSec, $stopSec, $cachefile, $xmlFile, $tableFile) = @_;
	
	open TABLE, ">>$runPath/$tableFile" or die "Couldn't open $tableFile.";
	#print  "P\tpending\tstart=$startSec\tstop=$stopSec\tcache=$cachefile\t$xmlFile\n";
	print TABLE  "P\tpending\t$startSec\t$stopSec\t$cachefile\t$xmlFile\n";
	close TABLE;
}

#-----------------------------------------------------------------------------------
#   lf_getPlaygroundSeconds 
#-----------------------------------------------------------------------------------
#  - takes path to playground file as a parameter.
#  - reads file line by line
#  - builds hash array with all of the playground seconds as keys
#-----------------------------------------------------------------------------------
#  RETURNS : the address of the playground seconds hash array 
#-----------------------------------------------------------------------------------
sub lf_getPlaygroundSeconds{

	my $playgroundFilePath = shift ;
	
	open PLAYGROUND_FILE, $playgroundFilePath or die "In  lf_getPlaygroundSeconds: Couldn't open $playgroundFilePath.";
	my %playgroundSeconds = ();
	
	# $sumPlaygroundSecs, the sum of the number of bad seconds, is used only for testing
	my $sumPlaygroundSecs = 0;
	
	while(<PLAYGROUND_FILE>){

		# if the first field is not a number, skip the line, header row
		if ( ! /^\d/ ){ next;} 
		
		#remove the newline character
		chomp;
		
		my ($startSecond, $stopSecond)  = split " ";
			
		#Create a hash array with all of the playground seconds
		for(my $t = $startSecond; $t <$stopSecond; $t++){
			$playgroundSeconds{$t}=1;
		}
		
		$sumPlaygroundSecs += $stopSecond - $startSecond;
	}
	
	close PLAYGROUND_FILE;
	return \%playgroundSeconds;
}

#-----------------------------------------------------------------------------------
#   f_setupOutputDirs()
#-----------------------------------------------------------------------------------
#  - checks to see if output dirs exist. if not, creates them
#-----------------------------------------------------------------------------------
#  Returns 
#-----------------------------------------------------------------------------------
sub f_setupOutputDirs {
	my ($path, $date, $runNum) = @_;
	
	$path .= "/$date-$runNum";
	
	my $xmldir = "$path/xml/";
	my $logdir = "$path/log/";
	my $errdir = "$path/err/";
	my $outdir = "$path/out/";
	
	if(! -d $path) { mkdir $path or die "Couldn't create $path.\n";}
	if(! -d $xmldir) { mkdir $xmldir or die "Couldn't create $xmldir.\n";}
	if(! -d $logdir) { mkdir $logdir or die "Couldn't create $logdir\n";}
	if(! -d $errdir) { mkdir $errdir or die "Couldn't create $errdir\n";}
	if(! -d $outdir) { mkdir $outdir or die "Couldn't create $outdir\n";}			
	
}
#-----------------------------------------------------------------------------------
#  f_getRunNumber
#-----------------------------------------------------------------------------------
#  - Looks in the run output directory to determine the next run
#    number for the current day
#-----------------------------------------------------------------------------------
#  Returns the next run number as NNN
#-----------------------------------------------------------------------------------
sub f_getRunNumber {
	my ($path, $date) = @_;
	opendir PATH, $path	;
	
	#get a bottom to top list of the files that end with DATE-XXX
	my $searchString = $date . '-\d{3}';
	my @dirs = reverse sort grep /$\$searchString/, readdir PATH;
	closedir PATH;
	
	if(scalar @dirs > 0){
		my @parts  = reverse split "-", $dirs[0];
		my $maxNum = $parts[0];
		my $nextNum = $maxNum + 1;
		
		#make sure nextNum is like NNN, not NN or N
		while(length($nextNum) < 3){ $nextNum = "0" . $nextNum;}
		
		return $nextNum;
	}

	#by default return 001
	return "001";
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
