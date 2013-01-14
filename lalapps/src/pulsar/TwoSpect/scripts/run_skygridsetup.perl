#!/usr/bin/perl

use strict;
use warnings;

my $band = $ARGV[0];
my $fstart = $ARGV[1];
my $fstep = $ARGV[2];
my $numberBands = $ARGV[3];
my $ifo = $ARGV[4];
my $outputPathBase = $ARGV[5];
my $ephemdir = $ARGV[6];
my $ephemyear = "08-11-DE405";
my $Tcoh = 1800;
my $sftoverlap = 900;
my $sky = "allSky";

if ($ifo eq "LHO") {
   $ifo = "H1";
} elsif ($ifo eq "LLO") {
   $ifo = "L1";
}

system("mkdir $outputPathBase/$band\_$ifo");
die "mkdir failed: $?" if $?;

for(my $ii=0; $ii<$numberBands; $ii++) {
   my $fmin = sprintf("%.3f", $fstart+$ii*$fstep);
   my $outfile = "$outputPathBase/$band\_$ifo/skygrid-${fmin}\-${fstep}HzBand.dat";
   
   system("../helperprograms/skygridsetup --fmin=$fmin --fspan=$fstep --Tcoh=$Tcoh --IFO=$ifo --ephemDir=$ephemdir --ephemYear=$ephemyear --skyRegion=$sky --SFToverlap=$sftoverlap --outfilename=$outfile");
   die "skygridsetup failed: $?" if $?;
}



