#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

my $jobnum = $ARGV[0];
srand(42+$jobnum);

system("mkdir /local/user/egoetz/$$");
die "mkdir failed: $?" if $?;

my $h0ul = 4.5709e-24;
my $Tsft = 1800.0;
my $dur = 40551300.0;
my $skygrid = "/home/egoetz/TwoSpect/mismatchDist/skygrid.dat";

for(my $ii=0; $ii<10; $ii++) {
   my $h0 = $h0ul;
   my $psi = sprintf("%.6f",0.5*pi*rand()-0.25*pi);
   my $phi0 = sprintf("%.6f",2.0*pi*rand());
   my $alpha = sprintf("%.6f",2.0*pi*rand());
   my $delta = sprintf("%.6f",acos(2.0*rand()-1.0)-0.5*pi);
   my $f0 = 401.25 + 0.24*rand();
   my $df = rand()*0.1;
   while ($df-0.5/$Tsft<1.0e-6) {
      $df = rand()*0.1;
   }
   my $P = rand()*0.2*($dur-7200.0)+7200.0;
   while ($P<1.2*2.0*$df*$Tsft*$Tsft) {
      $P = rand()*0.2*($dur-7200.0)+7200.0;
   }
   $f0 = sprintf("%.6f", $f0);
   $df = sprintf("%.6f", $df);
   $P = sprintf("%.6f", $P);
   my $asini = sprintf("%.6f",$df*$P/2.0/pi/$f0);
   my $randseedval = int(rand(1000000));
   
   open(MFDCONFIG,">/local/user/egoetz/$$/mfdconfig") or die "Cannot write to /local/user/egoetz/$$/mfdconfig $!";
   print MFDCONFIG<<EOF;
outSFTbname /local/user/egoetz/$$/testsfts.sft
outSingleSFT TRUE
IFO H1
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
timestampsFile /home/egoetz/TwoSpect/mismatchDist/timestamps.dat
generationMode 0
fmin 401.0
Band 2.9992
Tsft 1800
window Hann
Alpha $alpha
Delta $delta
h0 $h0
cosi 1.0
psi $psi
phi0 $phi0
Freq $f0
orbitasini $asini
orbitEcc 0.0
orbitTpSSBsec 900000000
orbitTpSSBnan 0
orbitPeriod $P
orbitArgp 0.0
f1dot 0.0
refTime 900000000
noiseSqrtSh 3.0e-23
randSeed $randseedval
EOF
   close(MFDCONFIG);

   open(INJECTION, ">>/home/egoetz/TwoSpect/mismatchDist/$jobnum/injections.dat") or die "Cannot write to /home/egoetz/TwoSpect/mismatchDist/$jobnum/injections.dat $!";
   print INJECTION "$alpha $delta $h0 $psi $phi0 $f0 $P $df\n";
   close(INJECTION);
   
   open(TWOSPECTCONFIG, ">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
   print TWOSPECTCONFIG<<EOF;
fmin 401.25
fspan 0.25
Tobs 40551300
Tcoh 1800
SFToverlap 900
ihsfar 1.0
ihsfomfar 1.0
tmplfar 1.0
Pmin $P
Pmax 8110260
dfmin $df
dfmax 0.1
skyRegion ($alpha,$delta)
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 50
maxTemplateLength 500
sftDir /local/user/egoetz/$$
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/mismatchDist/$jobnum
sftType MFD
IFO H1
markBadSFTs
FFTplanFlag 1
ULfmin $f0
ULfspan 0.1
fastchisqinv
useSSE
outfilename logfile_$ii.txt
ULfilename uls_$ii.dat
configCopy input_copy_$ii.conf
keepOnlyTopNumIHS 5
EOF
   close(TWOSPECTCONFIG);

   system("/home/egoetz/opt/lscsoft/bin/lalapps_Makefakedata_v4 \@/local/user/egoetz/$$/mfdconfig");
   die "system lalapps_Makefakedata_v4 failed: $?" if $?;
   
   system("/home/egoetz/TwoSpect/mismatchDist/TwoSpect_templateTest --config=/local/user/egoetz/$$/twospectconfig");
   die "TwoSpect_templateTest failed: $?" if $?;


   my $dist1 = -1.0;
   my $dist2 = -1.0;
   my $dist3 = -1.0;
   my $dist4 = -1.0;
   my @pt1 = (0, 0);
   my @pt2 = (0, 0);
   my @pt3 = (0, 0);
   my @pt4 = (0, 0);
   open(SKYFILE, $skygrid) or die "Cannot open $skygrid $!";
   while(my $line=<SKYFILE>) {
      if($line =~ /^(\d+.\d+) (-?\d+.\d+)/) {

         my $disttest1 = sqrt(($1-$alpha)*($1-$alpha) + ($2-$delta)*($2-$delta));
         my $disttest2 = sqrt((2.0*pi-$1+$alpha)*(2.0*pi-$1+$alpha) + ($2-$delta)*($2-$delta));
         my $dist = $disttest1;
         if($disttest1<=$disttest2) {
            $dist = $disttest1;
         } else {
            $dist = $disttest2;
         }

         if($dist<$dist1 || $dist1<0.0) {
            @pt4 = @pt3;
            @pt3 = @pt2;
            @pt2 = @pt1;
            $dist4 = $dist3;
            $dist3 = $dist2;
            $dist2 = $dist1;
            $pt1[0] = $1;
            $pt1[1] = $2;
            $dist1 = $dist;
         } elsif ($dist<$dist2 || $dist2<0.0) {
            @pt4 = @pt3;
            @pt3 = @pt2;
            $dist4 = $dist3;
            $dist3 = $dist2;
            $pt2[0] = $1;
            $pt2[1] = $2;
            $dist2 = $dist;
         } elsif ($dist<$dist3 || $dist3<0.0) {
            @pt4 = @pt3;
            $dist4 = $dist3;
            $pt3[0] = $1;
            $pt3[1] = $2;
            $dist3 = $dist;
         } elsif ($dist<$dist4 || $dist4<0.0) {
            $pt4[0] = $1;
            $pt4[1] = $2;
            $dist4 = $dist;
         }
      }
   }
   close(SKYFILE);

   open(SKYFILE2,">/local/user/egoetz/$$/skygrid2.dat") or die "Cannot write to /local/user/egoetz/$$/skygrid2.dat $!";
   print SKYFILE2<<EOF;
$pt1[0] $pt1[1]
$pt2[0] $pt2[1]
$pt3[0] $pt3[1]
$pt4[0] $pt4[1]
EOF
   close(SKYFILE2);

   open(TWOSPECTCONFIG,">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
   print TWOSPECTCONFIG<<EOF;
fmin 401.25
fspan 0.25
Tobs 40551300
Tcoh 1800
SFToverlap 900
ihsfar 1.0
ihsfomfar 1.0
tmplfar 1.0
Pmin 7200
Pmax 8110260
dfmin 0.0002
dfmax 0.1
skyRegionFile /local/user/egoetz/$$/skygrid2.dat
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 50
maxTemplateLength 500
sftDir /local/user/egoetz/$$
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/mismatchDist/$jobnum
sftType MFD
IFO H1
markBadSFTs
FFTplanFlag 1
ULfmin $f0
ULfspan 0.1
fastchisqinv
useSSE
outfilename logfile_${ii}_2.txt
ULfilename uls_${ii}_2.dat
configCopy input_copy_${ii}_2.conf
keepOnlyTopNumIHS 5
EOF
   close(TWOSPECTCONFIG);

   system("/home/egoetz/opt/lscsoft/bin/lalapps_TwoSpect --config=/local/user/egoetz/$$/twospectconfig");
   die "system lalapps_TwoSpect failed: $?" if $?;

   system("rm /local/user/egoetz/$$/*.sft");
   die "rm failed: $?" if $?;

}



