#!/usr/bin/perl

use strict;
use warnings;

my $band = $ARGV[0];
my $fstart = $ARGV[1];
my $fstep = $ARGV[2];
my $numberBands = $ARGV[3];
my $ifo = $ARGV[4];
my $analysisdate = $ARGV[5];
$outdirectory0 = "/home/egoetz/TwoSpect/S6/$band\_$ifo/output/";
my $Tcoh = 1800;
my $Tobs = 40551300;
my $Pmin = 7200.0;
my $Pmax = 8110260.0;
my $dfmin = 0.0002;
my $dfmax = 0.1;
my $fspan = 0.25;
my $ihsfar = 1.0e-14;
my $ihsfomfar = 1.0;
my $tmplfar = 1.0e-18;
my $t0 = 0.0;
my $ifokey = "";
if ($ifo == "LHO" || $ifo == "H1") {
   $t0 = 931081500;  #H1 start
   $ifo = "LHO"
   $ifokey = "H1";
}
elsif ($ifo == "LLO") {
   $t0 = 931113900;  #L1 start
   $ifo == "LLO"
   $ifokey = "L1";
}
my $blksize = 101;
my $minTemplateLength = 1;
my $maxTemplateLength = 500;
my $avesqrtSh = 1.0e-22;
my $sftOverlap = 900;
my $sfttype = "vladimir";
my $FFTplanFlag = 3;
my $ephemDir = "/home/egoetz/TwoSpect/S6";
my $ephemYear = "08-11-DE405";
my $directorynumber = 0;

system("mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate");
die "mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate failed: $?" if $?;
system("mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/in");
die "mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/in failed: $?" if $?;
system("mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/out");
die "mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/out failed: $?" if $?;
system("mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/err");
die "mkdir /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/err failed: $?" if $?;

my $maxnumperdag = 200;

open(DAG,">>/home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/dag") or die "Cannot write to /home/egoetz/TwoSpect/S6/$band_$ifo/$analysisdate/dag $!";
for($ii=0; $ii<$numberBands; $ii++) {
   open(INPUT, "skygrid-".sprintf("%.3f", $fstart+$ii*$fstep)."-${fstep}HzBand.dat") or die "Cannot open skygrid-".sprintf("%.3f", $fstart+$ii*$fstep)."-${fstep}HzBand.dat $!";
   my $numberindag = 0;
   my $firstpointset = 0;
   my $startbin = int(($fstart+$ii*$fstep)*$Tcoh + 0.5);
   
   while($line=<INPUT>) {
      if ($line =~ /^(\d+.\d+)\s(-?\d+.\d+)/) {
         if ($firstpointset==0) {
            my $skyRegionFile = "/home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/in/skygrid.$directorynumber";
            
            open(SKYFILE,">$skyRegionFile") or die "Cannot write to $skyRegionFile $!";
            print SKYFILE "$1 $2\n";
            $numberindag++;
            $firstpointset = 1;
            
         } elsif ($numberindag==$maxnumperdag-1) {
            print SKYFILE "$1 $2\n";
            $numberindag++;
            close(SKYFILE);
            
            my $lowestFneeded = ($fstart+$ii*$fstep) - 0.1 - 0.1 - 4e-3;
            my $lowestFinteger = int(floor($lowestFneeded - ($lowestFneeded % 2) + 0.5)) + 1;
            if ($lowestFinteger>$lowestFneeded) $lowestFinteger -= 2;
            my $highestFinteger = $lowestFinteger + 3;
            my $sftDir = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz";
                        
            my $outdirectory = $outdirectory0.$directorynumber;
            open(OUT,">/home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/in/$directorynumber") or die "Cannot write to /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/in/$directorynumber $!";
            print OUT<<EOF;
fmin sprintf("%.3f", $fstart+$ii*$fstep)
fspan $fspan
Tobs $Tobs
Tcoh $Tcoh
SFToverlap $sftOverlap
ihsfar $ihsfar
ihsfomfar $ihsfomfar
tmplfar $tmplfar
Pmin $Pmin
Pmax $Pmax
dfmin $dfmin
dfmax $dfmax
t0 $t0
blksize $blksize
avesqrtSh $avesqrtSh
minTemplateLength $minTemplateLength
maxTemplateLength $maxTemplateLength
sftDir $sftDir
ephemDir $ephemDir
ephemYear $ephemYear
outdirectory $outdirectory
sftType $sfttype
IFO $ifokey
markBadSFTs
skyRegionFile $skyRegionFile
FFTplanFlag $FFTplanFlag
fastchisqinv
useSSE
lineDetection 1.5
randSeed $startbin+$directorynumber
keepOnlyTopNumIHS 5
EOF
            close(OUT);
            
            print DAG "JOB A$directorynumber /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/condor\n";
            print DAG "VARS A$directorynumber PID=\"$directorynumber\"\n";
            
            $directorynumber++;
            
            $numberindag = 0;
            $firstpointset = 0;
            
         } else {
            print SKYFILE "$1 $2\n";
            $numberindag++;
         }
      }
   }
   if($numberindag>0) {
      close(SKYFILE);
      my $outdirectory = $outdirectory0.$directorynumber;
      
      my $lowestFneeded = ($fstart+$ii*$fstep) - 0.1 - 0.1 - 4e-3;
      my $lowestFinteger = int(floor($lowestFneeded - ($lowestFneeded % 2) + 0.5)) + 1;
      if ($lowestFinteger>$lowestFneeded) $lowestFinteger -= 2;
      my $highestFinteger = $lowestFinteger + 3;
      my $sftDir = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz";
      
      open(OUT,">/home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/in/$directorynumber") or die "Cannot write to /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/in/$directorynumber $!";
      print OUT<<EOF;
fmin sprintf("%.3f", $fstart+$ii*$fstep)
fspan $fspan
Tobs $Tobs
Tcoh $Tcoh
SFToverlap $sftOverlap
ihsfar $ihsfar
ihsfomfar $ihsfomfar
tmplfar $tmplfar
Pmin $Pmin
Pmax $Pmax
dfmin $dfmin
dfmax $dfmax
t0 $t0
blksize $blksize
avesqrtSh $avesqrtSh
minTemplateLength $minTemplateLength
maxTemplateLength $maxTemplateLength
sftDir $sftDir
ephemDir $ephemDir
ephemYear $ephemYear
outdirectory $outdirectory
sftType $sfttype
IFO $ifokey
markBadSFTs
skyRegionFile $skyRegionFile
FFTplanFlag $FFTplanFlag
fastchisqinv
useSSE
lineDetection 1.5
randSeed $startbin+$directorynumber
keepOnlyTopNumIHS 5
EOF
      close(OUT);
      
      print DAG "JOB A$directorynumber /home/egoetz/TwoSpect/S6/$band\_$ifo/$analysisdate/condor\n";
      print DAG "VARS A$directorynumber PID=\"$directorynumber\"\n";      
      $directorynumber++;
      
   }
   
   close(INPUT);
}

close(DAG);


