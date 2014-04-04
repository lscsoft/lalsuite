#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $band = '';
my $fstart = 0.0;
my $fstep = 0.0;
my $numberBands = 0;
my $ifo = '';
my $outputPathBase = '';
my $ephemdir = "/home/egoetz/TwoSpect/S6";
my $ephemyear = "08-11-DE405";
my $Tcoh = 1800;
my $sftoverlap = 900;
my $sky = "allSky";

GetOptions('band=s' => \$band, 'fstart=f' => \$fstart, 'fstep=f' => \$fstep, 'numberBands=i' => \$numberBands, 'ifo=s' => \$ifo, 'outputPathBase=s' => \$outputPathBase, 'Tcoh:f' => \$Tcoh, 'sftOverlap:f' => \$sftoverlap, 'skyRegion:s' => \$sky);

if ($ifo eq "LHO") {
   $ifo = "H1";
} elsif ($ifo eq "LLO") {
   $ifo = "L1";
} elsif ($ifo eq "Virgo") {
   $ifo = "V1";
}

system("mkdir $outputPathBase/$band\_$ifo");
die "mkdir failed: $?" if $?;

for(my $ii=0; $ii<$numberBands; $ii++) {
   my $fmin = sprintf("%.3f", $fstart+$ii*$fstep);
   my $outfile = "$outputPathBase/$band\_$ifo/skygrid-${fmin}\-${fstep}HzBand.dat";
   
   system("../skygridsetup --fmin=$fmin --fspan=$fstep --Tcoh=$Tcoh --IFO=$ifo --ephemDir=$ephemdir --ephemYear=$ephemyear --skyRegion=$sky --SFToverlap=$sftoverlap --outfilename=$outfile");
   die "skygridsetup failed: $?" if $?;
}



