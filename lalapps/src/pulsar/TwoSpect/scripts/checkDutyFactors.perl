#!/usr/bin/perl

use strict;
use warnings;
use POSIX;

my $fmin = 50;
my $fend = 50.25;
my $fstep = 0.25;
my $ifo = "LHO";

while ($fmin < $fend) {
   my $lowestFneeded = $fmin - 0.1 - 0.1 - 4e-3;
   my $lowestFinteger = int(floor($lowestFneeded - fmod($lowestFneeded, 2.0) + 0.5)) + 1;
   if ($lowestFinteger>$lowestFneeded) {
      $lowestFinteger -= 2;
   }
   my $highestFinteger = $lowestFinteger + 3;

   my $sftFile = "";
   my $sftType = "";
   if ($lowestFinteger<173 && ($ifo eq "LHO" || $ifo eq "LLO")) {
      $sftType = "vladimir";
      $sftFile = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz/*.sft";
   } else {
      $sftType = "standard";
      $sftFile = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz/*.sft";
   }

   my $t0 = 0;
   my $ifokey = '';
   if ($ifo eq "LHO") {
      $t0 = 931081500;
      $ifokey = "H1";
   } elsif ($ifo eq "LLO") {
      $t0 = 931113900;
      $ifokey = "L1";
   } else {
      $t0 = 931131900;
      $ifokey = "V1";
   }

   my $tobs = 40551300;
   my $tsft = 1800;
   my $sftOverlap = 900;

   system("/home/egoetz/opt/lscsoft-master/bin/lalapps_TwoSpect --Tobs=$tobs --Tcoh=$tsft --SFToverlap=$sftOverlap --t0=$t0 --fmin=$fmin --fspan=$fstep --IFO=$ifokey --avesqrtSh=1e-22 --blksize=101 --sftType=$sftType --sftFile=$sftFile --Pmin=7200 --Pmax=8110260 --dfmin=0.0 --dfmax=0.1 --skyRegion=\"(0,0)\" --ihsfar=1.0 --ihsfomfar=1.0 --tmplfar=1.0 --markBadSFTs --lineDetection=1.5 --fastchisqinv --useSSE --templateTest --templateTestF=$fmin --templateTestP=8110260 --templateTestDf=0.001");
   die "system lalapps_TwoSpect failed: $?" if $?;

   open(LOGFILE, "/atlas/user/atlas3/egoetz/lalsuite-master/lalapps/src/pulsar/TwoSpect/scripts/output/logfile.txt");
   open(OUTPUT, ">>/home/egoetz/TwoSpect/S6/$ifo.dutyfactor.dat");
   while (my $line=<LOGFILE>) {
      chomp($line);
      if ($line =~ /^Loading in SFTs... Duty factor = (\d.\d+)/) {
         print OUTPUT "$fmin $1 ";
      } elsif ($line =~ /^Assessing background... Duty factor of usable SFTs = (\d.\d+)/) {
         print OUTPUT "$1\n";
      }
   }
   close(LOGFILE);
   close(OUTPUT);

   $fmin = $fmin + $fstep;
}

