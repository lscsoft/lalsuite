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
my $PLOT_FILE = "/home/dsmackin/public_html/plots/powerTrial_$DATE.png";


open LOG, ">$LOG_FILE" or die "Couldn't open $LOG_FILE.";
open DATA, ">$DATA_FILE" or die "Couldn't open $DATA_FILE.";

#my @lengths = qw(10 100 1000 10000 100000 1000000);
#my @lengths = qw(1000 5000 10000 15000 100000);
my @lengths = qw(1000 1200 1400 1600 1800 2000);


my $startEpoch = 729439471;

foreach(@lengths){
	my $lt = localtime();
	print LOG ("Starting test of length $_ at $lt.\n");
	print "Starting test of length $_ at at $lt.\n\n";
	
	my $cacheFile = "/scratch/power/segment-" . $startEpoch . "-$_";
	my $outputFile = "/scratch/power/search_$startEpoch-$_.xml";

	my $cmd;
	my $startTime = time();
	my $totalTime = time() - $startTime;
	if(! -f $cacheFile){
		$cmd = "LALdataFind --lal-cache --instrument H --type RDS_R_L1  " .
					" --start " . $startEpoch .
					" --end " . ($startEpoch + $_) . " >$cacheFile";
		print $cmd, "\n\n";

		system $cmd;
		
		$totalTime = time() - $startTime;
		$lt = localtime();
		
		print LOG ("Completed test of length $_ at $lt.\n\n");
		print "Completed test of length $_ at $lt.\n\n";
		print "Length $_ took ", $totalTime, " seconds.\n";
		print LOG "Length $_ took $totalTime seconds.\n\n";
	}
	
	print DATA "$_ $totalTime";
	
	my $nseg = 2*($_) - 1;
	my $numpts = 16384 + (8192)*( $nseg-1)+ 2*(8192);

	#TIME RUNNING LALAPPS
	$startTime = time();
	$cmd = "lalapps_power --npts 16384   --nseg $nseg " .
				"--olap 8192   --olapfctr 3   --minfbin 2 " .
				"--mintbin 2   --flow 60   --delf 1 " .
				"--lngth 1024   --nsigma 2   --alphdef 0.5 " .
				"--segdcle 64   --threshold 10.0e-10 " .
				"--etomstr 10   --channel H1:LSC-AS_Q " .
				"--framecache $cacheFile " .
				"--simtype 0   --spectype useMedian   --window 2 " .
				"--epoch $startEpoch 0 " .
				"--numpts $numpts " .
				"--outfile $outputFile";
	print $cmd, "\n\n";
	system $cmd;
	$totalTime = time() - $startTime;
	
	print DATA " $totalTime\n";
}

close LOG;


#PLOT THE RESULTS
my $plotCmd = "/scratch/power/plot.tmp";
open PLOT_FILE, ">$plotCmd";

#make the xtics labels from the lengths values
my $xtics;
my $divider = "";
foreach (@lengths){ $xtics .= $divider . $_; $divider = ",";}

print PLOT_FILE << "PLOT_COMMANDS";
	set term png small color
	set xlabel "Number of Seconds of Data"
	set ylabel "Time in Seconds for Program to Run"
	set data style line
	set size 1.5,1.5
	set logscale x
	set logscale y
	set timestamp
	set output "$PLOT_FILE"
	set title "lalapps_power Optimization (first second in test $startEpoch; tested on Hydra)"
	set xlabel "Number of Seconds of Data Analyzed"
	set xtics ($xtics)
	set ylabel "Time to Run"
	set data style lp
	plot "$DATA_FILE" using 1:3 t "lalapps_power"

PLOT_COMMANDS

system "gnuplot $plotCmd";
#unlink ($plotCmd);


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
