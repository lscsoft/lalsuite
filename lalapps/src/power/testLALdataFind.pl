#!/usr/bin/perl -w

use strict;
#-----------------------------------------------------------------------------------
#
# testLALdataFind.pl - Program tests the speed of LALdataFind for different lengths of
# data
#
#
# Revision History
# June, 2003 - Dennis Mackin <dsmackin@stupidlinux.com>  - Created the original version of the script
# $Id$
#-----------------------------------------------------------------------------------

my $DATE = f_getDateYYMMDD();
my $LOG_FILE = "/scratch/power/LALdataFind-TimeTrials-$DATE.log";
my $DATA_FILE = "/scratch/power/LALdataFind-TimeTrials-$DATE.dat";

my $NSEG = 127;

open LOG, ">$LOG_FILE" or die "Couldn't open $LOG_FILE.";
open DATA, ">$DATA_FILE" or die "Couldn't open $DATA_FILE.";

#my @lengths = qw(10 100 1000 10000 100000 1000000);
my @lengths = qw(500 1000 5000 10000 50000 100000);


my $startEpoch = 729273613;

foreach(@lengths){
	my $lt = localtime();
	print LOG ("Starting test of length $_ at $lt.\n");
	print "Starting test of length $_ at at $lt.\n\n";
	
	my $cacheFile = "/scratch/power/segment-" . $startEpoch . "-" . ($startEpoch + $_);
	my $outputFile = "/scratch/power/search_$startEpoch-" . ($startEpoch + $_) . ".xml";

	my $cmd = "LALdataFind --lal-cache --instrument H --type RDS_R_L1  " .
				" --start " . $startEpoch .
				" --end " . ($startEpoch + $_) . " >$cacheFile";
	print $cmd, "\n\n";
	my $startTime = time();
	system $cmd;
	my $totalTime = time() - $startTime;
	$lt = localtime();
	
	print LOG ("Completed test of length $_ at $lt.\n\n");
	print "Completed test of length $_ at $lt.\n\n";
	print "Length $_ took ", $totalTime, " seconds.\n";
	print LOG "Length $_ took $totalTime seconds.\n\n";
	print DATA "$_ $totalTime";

	#TIME RUNNING LALAPPS
	$startTime = time();
	$cmd = "lalapps_power --npts 16384   --nseg $NSEG " .
				"--olap 8192   --olapfctr 3   --minfbin 2 " .
				"--mintbin 2   --flow 60   --delf 1 " .
				"--lngth 512   --nsigma 2   --alphdef 0.5 " .
				"--segdcle 64   --threshold 10.0e-15 " .
				"--etomstr 10   --channel H1:LSC-AS_Q " .
				"--framecache $cacheFile " .
				"--simtype 0   --spectype useMedian   --window 2 " .
				"--epoch $startEpoch 0 " .
				"--numpts 1048576 " .
				"--outfile $outputFile";
	print $cmd, "\n\n";
	system $cmd;
	$totalTime = time() - $startTime;
	
	print DATA " $totalTime\n";
}

close LOG;


#PLOT THE RESULTS
my $plotFile = "/scratch/power/plot.tmp";
open PLOT_FILE, ">$plotFile";

#make the xtics labels from the lengths values
my $xtics;
my $divider = "";
foreach (@lengths){ $xtics .= $divider . $_; $divider = ",";}

print PLOT_FILE << "PLOT_COMMANDS";
	set term png small color
	set xlabel "Number of Seconds of Data"
	set ylabel "Time to Run LALdataFind"
	set data style line
	set size 1.5,1.5
	set logscale x
	set timestamp
	set output "/home/dsmackin/public_html/plots/time_to_run_laldatafind.png"
	set title "Cache File Time Lenght Optimization (first second in test $startEpoch; tested on Hydra)"
	set xlabel "Number of Seconds in Cache File"
	set xtics ($xtics)
	set ylabel "Time to Run"
	set data style lp
	plot "$DATA_FILE" using 1:2 t "LALdataFind", "$DATA_FILE" using 1:3 t "lalapps_power ($NSEG segments)"

PLOT_COMMANDS

system "gnuplot $plotFile";
#unlink ($plotFile);


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
