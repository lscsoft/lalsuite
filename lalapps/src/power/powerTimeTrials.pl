#!/usr/bin/perl -w
use strict;
#-----------------------------------------------------------------------------------
#
# powerTimeTrials.pl - Program tests the speed of lalapps_power
# for multiple lengths of data and then plots the result.
#
# Revision History
# $Id$
#-----------------------------------------------------------------------------------

my $DATE = f_getDateYYMMDD();
my $LOG_FILE = "/scratch/power/Power-TimeTrials-$DATE.log";
my $DATA_FILE = "/scratch/power/Power-TimeTrials-$DATE.dat";
my $PLOT_FILE = "/home/dsmackin/public_html/plots/power_time_trials-$DATE.png";

open LOG, ">$LOG_FILE" or die "Couldn't open $LOG_FILE.";
open DATA, ">$DATA_FILE" or die "Couldn't open $DATA_FILE.";

#my @lengths = qw(10 100 1000 10000 100000 1000000);
my @lengths = qw(512 1024 2048 4096);

my $cacheFileLength = 10000;

my $startEpoch = 729283613;
my $cacheFile = "/scratch/power/segment-" . $startEpoch . "-" . ($startEpoch + $cacheFileLength );

foreach(@lengths){
	my $lt = localtime();
	print LOG ("Starting test of length $_ at $lt.\n");
	print "Starting test of length $_ at at $lt.\n\n";
	
	my $nseg = 2*$_ -1;
	my $numpts = 16384 + (16384/2)*($nseg-1) + 2*8192;
	my $outputFile = "/scratch/power/search_$startEpoch-" . ($startEpoch + $cacheFileLength) ."-nseg_$nseg.xml";

	#TIME RUNNING LALAPPS
	my $startTime = time();
	my $cmd = "lalapps_power --npts 16384   --nseg $nseg " .
				"--olap 8192   --olapfctr 3   --minfbin 2 " .
				"--mintbin 2   --flow 60   --delf 1 " .
				"--lngth 512   --nsigma 2   --alphdef 0.5 " .
				"--segdcle 64   --threshold 10.0e-15 " .
				"--etomstr 10   --channel H1:LSC-AS_Q " .
				"--framecache $cacheFile " .
				"--simtype 0   --spectype useMedian   --window 2 " .
				"--epoch $startEpoch 0 " .
				"--numpts  $numpts " .
				"--outfile $outputFile";
	print $cmd, "\n\n";
	print LOG $cmd, "\n\n";
	system $cmd;
	my $totalTime = time() - $startTime;
	
	print DATA "$_ $totalTime\n";
}

close LOG;
close DATA;


#PLOT THE RESULTS
my $plotFile = "/scratch/power/power_plot.tmp";
open PLOT_FILE, ">$plotFile";

#make the xtics labels from the lengths values
my $xtics;
my $divider = "";
foreach (@lengths){ $xtics .= $divider . $_; $divider = ",";}

print PLOT_FILE << "PLOT_COMMANDS";
	set term png small color
	set xlabel "Number of Seconds of Data"
	set ylabel "Time to Run lalapps_power"
	set data style line
	set size 1.5,1.5
	set logscale x
	set logscale y
	set timestamp
	set output "$PLOT_FILE"
	set title "lalapps_power Length Optimization (first second in test $startEpoch; tested on Hydra)"
	set xtics ($xtics)
	set data style lp
	plot "$DATA_FILE" using 1:2 t "lalapps_power (2*length -1  segments)"

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
