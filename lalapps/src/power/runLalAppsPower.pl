#!/usr/bin/perl -w
use strict;
#-----------------------------------------------------------------------------------
#
# runLalAppsPower.pl - This script was created to make it easier to test and run lalApps_power
#	with various combination of the parameter values
#
# Revision History
# June, 2003 - Dennis Mackin <dsmackin@stupidlinux.com>  - Created the original version of the script
# $Id$
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#  PARAMETERS
#-----------------------------------------------------------------------------------

# use f_getDateYYMMDD to get current date
my $date = f_getDateYYMMDD( );

# $SEED is used for the random number generator used to create random noise.
my $SEED = 4;

# number of points to analyze per second
# do not change this value
my $NPTS = 16384;

# Number of segments to analyze
my $NSEG = 128;

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
my $SEGDCLE = 16;

# If the probability of an event occuring from background noise
# is less than the threshold, then an event is considered an Event
# trigger generator (ETG). This number should be between 1.0e-09 and
# 1.0e-05
my $THRESHOLD = "10.0e-15";

# Controls how many events get written out per segment
my $ETOMSTR = 10;

# The source of the data
my $CHANNEL = "H1:LSC-AS_Q";
my $FRAMECACHE = "/home/saikat/powerscript/frcache_plgrnd_1";

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
my $NUMPTS = $NPTS + ($NPTS - $OLAP)*($NSEG-1) + 2*$OLAP;

# If you want to create a file with the print spectrum, uncomment the line with "--printSpectrum"
#	If you don't want the print spectrum output file, uncomment the line without "--printSpectrum"
my $PRINT_SPECTRUM = "";
#my $PRINT_SPECTRUM = "--printSpectrum";

# Do not change this value. Changing the seed will automatically change this.
my $NOISE = "1 $SEED";

# Name of the XML output file
my $OUTFILE = "REALSEARCH-$date-SEED_$SEED-GPS_$EPOCH-LENGTH_$LNGTH-THRESHOLD_$THRESHOLD.xml";

#==============================
#GAUSSIAN SINE WAVE PARAMETERS
#==============================
# To inject a sine wave, set $INCLUDE_SINE_WAVE = 1. Otherwise, set it to 0
my $INCLUDE_SINE_WAVE = 0;

#Set the frequency of the sine
my $FREQ = 450;
my $DELTA = 16;

#Set the amplitude of the sine wave
my $AMPLITUDE = 50;

#Set the width of the gaussian sine wave
my $WIDTH = .5;


#==============================
# main program
#==============================
# RUN THE POWER CODE
my $cmd = << "RUN_POWER";
	lalapps_power
		--npts $NPTS
		--nseg $NSEG
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
		--simtype $SIMTYPE
		--spectype $SPECTYPE
		--window $WINDOW
		--epoch $EPOCH 0
		--numpts $NUMPTS
		$PRINT_SPECTRUM
		--framecache $FRAMECACHE
		--outfile $OUTFILE
RUN_POWER

if($INCLUDE_SINE_WAVE){ $cmd .= "--sine $FREQ $DELTA $AMPLITUDE $WIDTH";}

$cmd =~ s/(\n|\t)/ /g;
print "\n\n\n",$cmd, "\n\n\n";
system $cmd and die ("Error running lalapps_power.");


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

