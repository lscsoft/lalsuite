#!/usr/bin/perl -w
#-----------------------------------------------------------------------------------
# countJobsTable.pl 
#-----------------------------------------------------------------------------------
#  -  utility to print the number of jobs in the jobs table
#
# July, 2003 - Dennis Mackin <dsmackin@stupidlinux.com>  - Created the original version of the script
# $Id$
#-----------------------------------------------------------------------------------
use strict;

# Check to make sure date of table is included as arg
my $USAGE = "\nusage: countJobsTable.pl  YYMMDD \n\n";
if(! $ARGV[0])
{
	die $USAGE;
}

#my $DATE = f_getDateYYMMDD();
my $DATE =$ARGV[0];	
my $JOBS_TABLE = "/home/dsmackin/lal/power_jobs_$DATE.tbl";
open TABLE, $JOBS_TABLE;

my $i=0;
while(<TABLE>){
	$i++;
}

print "$i jobs in table.\n";

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