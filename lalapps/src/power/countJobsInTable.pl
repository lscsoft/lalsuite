#!/usr/bin/perl  -w
#-----------------------------------------------------------------------------------
# countJobsTable.pl 
#-----------------------------------------------------------------------------------
#  -  utility to print the number of jobs in the jobs table
#
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

# Check to make sure date of table is included as arg
my $USAGE = "\nusage: countJobsTable.pl PARAMETERS_FILE YYMMDD-RUN_NUM  \n\n";
if(! $ARGV[1])
{
	die $USAGE;
}

#get the date and runNum
# NEED TO ADD ERROR CHECKING HERE !!!!
my ($DATE, $RUN_NUM) = split '-', $ARGV[1];	

my $parametersPath = $ARGV[0];
my $params = f_parseParametersFile($parametersPath);

#set the path for all program output files
my $runPath = ${$params}{OUTPUT_PATH}  . "/$DATE-$RUN_NUM";
#build the parameter dependent globals
my $JOBS_TABLE = "$runPath/power_jobs.tbl";

print "\nCounting jobs in table $JOBS_TABLE.\n\n";	
open TABLE, "$JOBS_TABLE" or die "Couldn't open $JOBS_TABLE.\n";

my $i=0;
my %statusCount = ( );
my $amountOfData = 0;
while(<TABLE>){
	chomp;
	my ($statusCode, $statusDescription, $startSec, $endSec, $cache, $xml) = split "\t";
	$i++;
	if ($statusCount{$statusDescription}){
		$statusCount{$statusDescription}++;
	}else{
		$statusCount{$statusDescription} = 1;
	}
	
	#get the length of data 
	#my @parts = reverse split "-", $fields[4];
	#$amountOfData += scalar($parts[0]);
	$amountOfData += $endSec - $startSec;
}

foreach(sort keys %statusCount){
	print "$_ =", $statusCount{$_},"\n";
}
print "\nTOTAL = $i \n";
print "Total Seconds Analyzed: $amountOfData\n\n";