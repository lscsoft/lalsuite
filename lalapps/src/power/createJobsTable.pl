#!/usr/bin/perl -w
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
use strict;

#Add the user's entire path and the program pathto @INC (where perl searches for modules)
my $program_path = $0;
# now clean the filename off the path
$program_path =~ s/\/(\w+)\.pl$//;
use lib split(":", $ENV{'PATH'}), qw($program_path);

#load the power tools module which contains functions shared by the scripts
use power_tools;

#-----------------------------------------------------------------------------------
#  GLOBAL VARIABLES
#-----------------------------------------------------------------------------------

# use f_getDateYYMMDD to get current date
my $DATE = f_getDateYYMMDD();

# Check to make sure  parms package is included as arg
my $USAGE = "\nusage: createJobsTable.pl  PARAMETERS_FILE \n\n";
if(! $ARGV[0])
{
	die $USAGE;
}

my $parametersPath = $ARGV[0];
my @parts = reverse split('/', $parametersPath);
my $parametersFile = $parts[0];

my $params = f_parseParametersFile($parametersPath);

#foreach(sort keys %{$params}){
#	print $_, "=" , ${$params}{$_}, "\n";
#}

my $startTime = time();
#get the next runNumber as string NNN
my $runNum = f_getRunNumber(${$params}{'OUTPUT_PATH'},$DATE);

#call function to setup output dirs
 f_setupOutputDirs(${$params}{'OUTPUT_PATH'}, $DATE, $runNum);
 
 #set the path for all program output files
my $outputPath = ${$params}{OUTPUT_PATH};
my $runID = "$DATE-$runNum";
my $runPath = $outputPath  . "/$runID";
 
#build the parameter dependent globals
my $JOBS_TABLE = "$runPath/${$params}{'JOBS_TABLE'}";

my %STATUS = (
					P => "Pending",
					R => "Running",
					C => "Complete",
					U => "User Review Required",
					E => "Error",
					BC => "Bad Cache File",
					NF => "Output file not found");

#-----------------------------------------------------------------------------------
#  MAIN
#-----------------------------------------------------------------------------------

open LOG, ">>$runPath/${$params}{'LOG'}" or die "Couldn't open log $runPath/${$params}{'LOG'}\n";
my $t = localtime();
print LOG "Starting $0 at: $t\n"; # $0 stores program name
print "Starting $0 at: $t\n";

#copy the parameters file to the output directory
print "cp $parametersPath $runPath\n";
system ("cp $parametersPath $runPath"); 

#Get the description for this run
print "\n\nEnter a description for this run: ";
my $description;
while(<STDIN>){
	chomp;
	$description = $_;
	last;
}
close STDIN;

# write the description to the runs table
open RUN_TABLE, ">>$outputPath/${$params}{'RUN_TABLE'}" or die "Couldn't open " . ${$params}{'RUN_TABLE'} . ".\n";
print RUN_TABLE "$runID\t$ENV{'USER'}\t$description\n"; 
close RUN_TABLE;

#if the jobs table exists, move it to a backup file
if (-f "$runPath/$JOBS_TABLE"){ rename "$runPath/$JOBS_TABLE", "$runPath/$JOBS_TABLE.bak";}

#playgroundSeconds is a reference to a hash array that contains all the playground
# seconds.
my $playgroundSeconds;
if(${$params}{'USE_PLAYGROUND'}){
	$playgroundSeconds = lf_getPlaygroundSeconds(${$params}{'PLAYGROUND_FILE'});
	
	# For debug purposes, count and display the number of playground seconds
	my $countPlaygroundSecs =0;
	foreach (keys %{$playgroundSeconds}){
		$countPlaygroundSecs++;
	}
	print "The Playground is $countPlaygroundSecs seconds.\n\n";
}

# now create the submit script using quality data file and the 
# the playground seconds hash array
lf_createJobsForGoodData (
									$runPath,
									${$params}{'DATA_QUALITY_FILE'}, 
									scalar(${$params}{'TIME_CHUNK_MIN_SIZE'}), 
									scalar(${$params}{'TIME_CHUNK_MAX_SIZE'}),
									$playgroundSeconds);

# update the html page that lists the job runs
f_updateNotebook( $params);
					
# update the notebook page that give stats for this run
f_updateNotebookPage( 
					$params,
					$runPath,
					$runID);

my $totalTime = time() - $startTime;
print "Completed in $totalTime seconds.\n";
print LOG "Completed in $totalTime seconds\n.";
print "The jobs table for run $DATE-$runNum has been created. \n\n".
			"To submit the jobs to condor, run:\n\nprocessJobsTable.pl $runPath/$parametersFile $DATE-$runNum\n\n";

#-----------------------------------------------------------------------------------
#  lf_createJobsForGoodData
#-----------------------------------------------------------------------------------
#  - Reads the data quality file
# - throws away bad chunks
#  - check each good chunk to see if it contains playground data
#-----------------------------------------------------------------------------------
#  lists chunks of good data from the  data quality file that are in the playground
#-----------------------------------------------------------------------------------
sub lf_createJobsForGoodData {

	my ($runPath, $dataQualityFile, $minChunkSize, $maxChunkSize, $playgroundSecondsHash) = @_;
	my ($chunkStart, $chunkStop);
	
	open DATA_QUALITY_FILE, $dataQualityFile
			or die "In lf_createJobsForGoodData: Couldn't open $dataQualityFile." ;	
	
	while(<DATA_QUALITY_FILE>){
		chomp;
		
		#if first character is # is then it is a comment, skip line
		if (/^\#/){next;}
			
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
		while( scalar(${$params}{'USE_PLAYGROUND'}) and not exists(${$playgroundSecondsHash}{$chunkStart}) and $chunkStart < $stopSecond){ 
			$chunkStart++;
		}
		
		#Loop through data quality file chunk
		for($chunkStop=$chunkStart; $chunkStop < $stopSecond; $chunkStop++){
			if((${$params}{'USE_PLAYGROUND'} and not exists ${$playgroundSecondsHash}{$chunkStop})  or  ($chunkStop - $chunkStart)  >= $maxChunkSize ){
				
				f_processChunk($runPath,$chunkStart, $chunkStop);
				
				#fast forward to the next good second 
				while((${$params}{'USE_PLAYGROUND'} and not exists ${$playgroundSecondsHash}{$chunkStop})and $chunkStop < $stopSecond){ 
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
	#print "Processing chunk $runPath, $startSec, $stopSec.\n";
	#make sure chunk is larger than minimum size
	if($stopSec - $startSec > ${$params}{'TIME_CHUNK_MIN_SIZE'}){
	
		#build output file string
		my $outfile = "$runPath/xml/$startSec-$duration.xml";
		my $framecache = ${$params}{'CACHE_PATH'} . "/" . ${$params}{'INSTRUMENT'} . "-$startSec-$duration";
		
		my $statusCode = "P"; 
		
		#build the cache file. if it is 0 length, then change status to BC - Bad Caceh File
		if (f_buildCacheFile($startSec,  $stopSec,$framecache, ${$params}{'INSTRUMENT'}) == 0 ){$statusCode = "BC";}
		f_writeJobsTable($runPath, $statusCode, $startSec,  $stopSec, $framecache, $outfile, $JOBS_TABLE);
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
	my ($runPath,$statusCode, $startSec, $stopSec, $cachefile, $xmlFile, $tableFile) = @_;
	
	open TABLE, ">>$tableFile" or die "Couldn't open $tableFile.";
	#print   "P\tpending\t$startSec\t$stopSec\t$cachefile\t$xmlFile\n";
	print TABLE  "$statusCode\t$STATUS{$statusCode}\t$startSec\t$stopSec\t$cachefile\t$xmlFile\n";
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

