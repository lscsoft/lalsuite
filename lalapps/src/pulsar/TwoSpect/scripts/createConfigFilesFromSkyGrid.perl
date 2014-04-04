#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
use Getopt::Long;

my $band = '';
my $fstart = 0.0;
my $fstep = 0.0;
my $numberBands = 0;
my $ifo = '';
my $outputPathBase = '';
my $analysisdate = '';
my $sftVersion = 0;
my $Tcoh = 1800.0;
my $Tobs = 40551300;
my $Pmin = 7200.0;
my $Pmax = 8110260.0;
my $dfmin = 0.0002;
my $dfmax = 0.1;
my $ihsfar = 1.0e-14;
my $ihsfomfar = 1.0;
my $tmplfar = 1.0e-18;
my $blksize = 101;
my $minTemplateLength = 1;
my $maxTemplateLength = 500;
my $avesqrtSh = 1.0e-22;
my $sftOverlap = 900;
my $FFTplanFlag = 3;
my $ephemDir = "/home/egoetz/TwoSpect/S6";
my $ephemYear = "08-11-DE405";
my $maxnumperdag = 200;
my $linedetection = -1.0;
my $keepTopIHS = -1;
my $nofastchisq = 0;
my $noSSE = 0;
my $noMarkBadSFTs =0;
my $t0 = 0.0;
#LHO t0 S6 = 931081500
#LLO t0 S6 = 931113900
#Virgo t0 VSR2/3 = 931131900

GetOptions('band=s' => \$band, 'fstart=f' => \$fstart, 'fstep=f' => \$fstep, 'numberBands=i' => \$numberBands, 'ifo=s' => \$ifo, 'outputPathBase=s' => \$outputPathBase, 'analysisDate=s' => \$analysisdate, 'sftVersion=i' => \$sftVersion, 'Tcoh=:f' => \$Tcoh, 'Tobs:f' => \$Tobs, 'Pmin=:f' => \$Pmin, 'Pmax=:f' => \$Pmax, 'dfmin=:f' => \$dfmin, 'dfmax=:f' => \$dfmax, 'ihsfar=:f' => \$ihsfar, 'ihsfomfar:f' => \$ihsfomfar, 'tmplfar:f' => \$tmplfar, 'blksize:i' => \$blksize, 'minTemplateL:i' => \$minTemplateLength, 'maxTemplateL:i' => \$maxTemplateLength, 'avesqrtSh:f' => \$avesqrtSh, 'sftOverlap:f' => \$sftOverlap, 'FFTplan:i' => \$FFTplanFlag, 'maxnumperdag:i' => \$maxnumperdag, 'linedetection:f' => \$linedetection, 'keepTopIHS:i' => \$keepTopIHS, 'nofastchisq' => $nofastchisq, 'noSSE' => \$noSSE, 'noMarkBadSFTs' => \$noMarkBadSFTs, 't0=f' => \$t0);

my $fspan = $fstep;
my $ifokey = "";
if ($ifo eq "LHO" || $ifo eq "H1") {
   $ifo = "LHO";
   $ifokey = "H1";
} elsif ($ifo eq "LLO" || $ifo eq "L1") {
   $ifo = "LLO";
   $ifokey = "L1";
} elsif ($ifo eq "Virgo" || $ifo eq "V1") {
   $ifo = "Virgo";
   $ifokey = "V1";
}
my $sfttype = "";
if ($sftVersion==1) {
   $sfttype = "vladimir";
   $sftVersion = "";
} elsif ($sftVersion==2) {
   $sfttype = "standard";
   $sftVersion = "_v2";
} else {
   die "sftVersion must equal 1 or 2, but input was $sftVersion";
}

my $directorynumber = 0;
my $outdirectory0 = "$outputPathBase/$band\_$ifokey/$analysisdate/output/";

system("mkdir $outputPathBase/$band\_$ifokey/$analysisdate");
die "mkdir $outputPathBase/$band\_$ifokey/$analysisdate failed: $?" if $?;
system("mkdir $outputPathBase/$band\_$ifokey/$analysisdate/in");
die "mkdir $outputPathBase/$band\_$ifokey/$analysisdate/in failed: $?" if $?;
system("mkdir $outputPathBase/$band\_$ifokey/$analysisdate/out");
die "mkdir $outputPathBase/$band\_$ifokey/$analysisdate/out failed: $?" if $?;
system("mkdir $outputPathBase/$band\_$ifokey/$analysisdate/err");
die "mkdir $outputPathBase/$band\_$ifokey/$analysisdate/err failed: $?" if $?;
system("mkdir $outputPathBase/$band\_$ifokey/$analysisdate/output");
die "mkdir $outputPathBase/$band\_$ifokey/$analysisdate/output failed: $?" if $?;

open(CONDORFILE,">>$outputPathBase/$band\_$ifokey/$analysisdate/condor") or die "Cannot write to $outputPathBase/$band\_$ifokey/$analysisdate/condor $!";
print CONDORFILE<<EOF;
universe=standard
executable=/home/egoetz/opt/lscsoft/bin/lalapps_TwoSpect
input=/dev/null
output=$outputPathBase/$band\_$ifokey/$analysisdate/out/out.\$(PID)
error=$outputPathBase/$band\_$ifokey/$analysisdate/err/err.\$(PID)
arguments=--config=$outputPathBase/$band\_$ifokey/$analysisdate/in/\$(PID)
log=/local/user/egoetz/$band\_$ifokey.$analysisdate.log
request_memory = 1550
notification=Never
notify_user=evan.goetz\@aei.mpg.de
queue
EOF
close(CONDORFILE);

open(DAG,">>$outputPathBase/$band\_$ifokey/$analysisdate/dag") or die "Cannot write to $outputPathBase/$band\_$ifokey/$analysisdate/dag $!";
for(my $ii=0; $ii<$numberBands; $ii++) {
   open(INPUT, "$outputPathBase/$band\_$ifokey/skygrid-".sprintf("%.3f", $fstart+$ii*$fstep)."-${fstep}HzBand.dat") or die "Cannot open $outputPathBase/$band\_$ifokey/skygrid-".sprintf("%.3f", $fstart+$ii*$fstep)."-${fstep}HzBand.dat $!";
   my $numberindag = 0;
   my $firstpointset = 0;
   my $startbin = int(($fstart+$ii*$fstep)*$Tcoh + 0.5);
   my $skyregionfile = "";
   
   while(my $line=<INPUT>) {
      if ($line =~ /^(\d+.\d+)\s(-?\d+.\d+)/) {
         if ($firstpointset==0) {
            $skyregionfile = "$outputPathBase/$band\_$ifokey/$analysisdate/in/skygrid.$directorynumber";
            
            open(SKYFILE,">$skyregionfile") or die "Cannot write to $skyregionfile $!";
            print SKYFILE "$1 $2\n";
            $numberindag++;
            $firstpointset = 1;
            
         } elsif ($numberindag==$maxnumperdag-1) {
            print SKYFILE "$1 $2\n";
            $numberindag++;
            close(SKYFILE);
            
            my $lowestFneeded = ($fstart+$ii*$fstep) - 0.1 - 0.1 - 4e-3;
            my $lowestFinteger = int(floor($lowestFneeded - fmod($lowestFneeded, 2.0) + 0.5)) + 1;
            if ($lowestFinteger>$lowestFneeded) {
               $lowestFinteger -= 2;
            }
            my $highestFinteger = $lowestFinteger + 3;
            my $sftDir = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz";
            
            my $fmin = sprintf("%.3f", $fstart+$ii*$fstep);
            my $randseedval = $startbin+$directorynumber;
                        
            my $outdirectory = $outdirectory0.$directorynumber;
            open(OUT,">$outputPathBase/$band\_$ifokey/$analysisdate/in/$directorynumber") or die "Cannot write to $outputPathBase/$band\_$ifokey/$analysisdate/in/$directorynumber $!";
            print OUT<<EOF;
fmin $fmin
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
sftFile $sftDir/*$sftVersion.sft
ephemDir $ephemDir
ephemYear $ephemYear
outdirectory $outdirectory
sftType $sfttype
IFO $ifokey
skyRegionFile $skyregionfile
FFTplanFlag $FFTplanFlag
randSeed $randseedval
EOF

            if ($linedetection>=0.0) {
               print OUT "lineDetection $linedetection\n";
            }
            if ($keepTopIHS>=0) {
               print OUT "keepOnlyTopNumIHS $keepTopIHS\n";
            }
            if ($noMarkBadSFTs==0) {
               print OUT "markBadSFTs\n";
            }
            if ($noSSE==0) {
               print OUT "useSSE\n";
            }
            if ($nofastchisq==0) {
               print OUT "fastchisqinv\n";
            }
            close(OUT);
            
            print DAG "JOB A$directorynumber $outputPathBase/$band\_$ifokey/$analysisdate/condor\n";
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
      my $lowestFinteger = int(floor($lowestFneeded - fmod($lowestFneeded, 2.0) + 0.5)) + 1;
      if ($lowestFinteger>$lowestFneeded) {
         $lowestFinteger -= 2;
      }
      my $highestFinteger = $lowestFinteger + 3;
      my $sftDir = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz";
      
      my $fmin = sprintf("%.3f", $fstart+$ii*$fstep);
      my $randseedval = $startbin+$directorynumber;
      
      open(OUT,">$outputPathBase/$band\_$ifokey/$analysisdate/in/$directorynumber") or die "Cannot write to $outputPathBase/$band\_$ifokey/$analysisdate/in/$directorynumber $!";
      print OUT<<EOF;
fmin $fmin
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
sftFile $sftDir/*$sftVersion.sft
ephemDir $ephemDir
ephemYear $ephemYear
outdirectory $outdirectory
sftType $sfttype
IFO $ifokey
skyRegionFile $skyregionfile
FFTplanFlag $FFTplanFlag
randSeed $randseedval
EOF

      if ($linedetection>=0.0) {
         print OUT "lineDetection $linedetection\n";
      }
      if ($keepTopIHS>=0) {
         print OUT "keepOnlyTopNumIHS $keepTopIHS\n";
      }
      if ($noMarkBadSFTs==0) {
         print OUT "markBadSFTs\n";
      }
      if ($noSSE==0) {
         print OUT "useSSE\n";
      }
      if ($nofastchisq==0) {
         print OUT "fastchisqinv\n";
      }
      close(OUT);
      
      print DAG "JOB A$directorynumber $outputPathBase/$band\_$ifokey/$analysisdate/condor\n";
      print DAG "VARS A$directorynumber PID=\"$directorynumber\"\n";      
      $directorynumber++;
      
   }
   
   close(INPUT);
}

close(DAG);




