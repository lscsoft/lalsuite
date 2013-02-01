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
my $skygrid2 = "/home/egoetz/TwoSpect/mismatchDist/skygrid_L1.dat";
for(my $ii=0; $ii<10; $ii++) {
   my $h0 = 0.5*$h0ul;
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
   my $mfdrandseed = int(rand(1000000));
   
   open(MFDCONFIG,">/local/user/egoetz/$$/mfdconfig") or die "Cannot write to /local/user/egoetz/$$/mfdconfig $!";
   print MFDCONFIG<<EOF;
outSFTbname /local/user/egoetz/$$/testsfts.sft
outSingleSFT TRUE
IFO H1
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
timestampsFile /home/egoetz/TwoSpect/IFOmismatch/timestamps.dat
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
randSeed $mfdrandseed
EOF
   close(MFDCONFIG);
   
   open(INJECTION, ">>/home/egoetz/TwoSpect/IFOmismatch/$jobnum/injections.dat") or die "Cannot write to /home/egoetz/TwoSpect/IFOmismatch/$jobnum/injections.dat $!";
   print INJECTION "$alpha $delta $h0 $psi $phi0 $f0 $P $df\n";
   close(INJECTION);

   my $dist1 = -1.0;
   my $dist2 = -1.0;
   my $dist3 = -1.0;
   my $dist4 = -1.0;
   my $dist5 = -1.0;
   my $dist6 = -1.0;
   my $dist7 = -1.0;
   my $dist8 = -1.0;
   my $dist9 = -1.0;
   my $dist10 = -1.0;
   my $dist11 = -1.0;
   my $dist12 = -1.0;
   my $dist13 = -1.0;
   my $dist14 = -1.0;
   my $dist15 = -1.0;
   my $dist16 = -1.0;
   my @pt1 = (0, 0);
   my @pt2 = (0, 0);
   my @pt3 = (0, 0);
   my @pt4 = (0, 0);
   my @pt5 = (0, 0);
   my @pt6 = (0, 0);
   my @pt7 = (0, 0);
   my @pt8 = (0, 0);
   my @pt9 = (0, 0);
   my @pt10 = (0, 0);
   my @pt11 = (0, 0);
   my @pt12 = (0, 0);
   my @pt13 = (0, 0);
   my @pt14 = (0, 0);
   my @pt15 = (0, 0);
   my @pt16 = (0, 0);
   open(SKYFILE, $skygrid) or die "Cannot open $skygrid $!";
   while(my $line=<SKYFILE>) {
      if($line =~ /^(\d+.\d+) (-?\d+.\d+)/) {

         my $dist = acos(sin(abs($2-0.5*pi))*sin(abs($delta-0.5*pi))*cos($1-$alpha)+cos(abs($2-0.5*pi))*cos(abs($delta-0.5*pi)));

         if($dist<$dist1 || $dist1<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            @pt4 = @pt3;
            @pt3 = @pt2;
            @pt2 = @pt1;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $dist4 = $dist3;
            $dist3 = $dist2;
            $dist2 = $dist1;
            $pt1[0] = $1;
            $pt1[1] = $2;
            $dist1 = $dist;
         } elsif ($dist<$dist2 || $dist2<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            @pt4 = @pt3;
            @pt3 = @pt2;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $dist4 = $dist3;
            $dist3 = $dist2;
            $pt2[0] = $1;
            $pt2[1] = $2;
            $dist2 = $dist;
         } elsif ($dist<$dist3 || $dist3<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            @pt4 = @pt3;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $dist4 = $dist3;
            $pt3[0] = $1;
            $pt3[1] = $2;
            $dist3 = $dist;
         } elsif ($dist<$dist4 || $dist4<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $pt4[0] = $1;
            $pt4[1] = $2;
            $dist4 = $dist;
         } elsif ($dist<$dist5 || $dist5<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $pt5[0] = $1;
            $pt5[1] = $2;
            $dist5 = $dist;
         } elsif ($dist<$dist6 || $dist6<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $pt6[0] = $1;
            $pt6[1] = $2;
            $dist6 = $dist;
         } elsif ($dist<$dist7 || $dist7<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $pt7[0] = $1;
            $pt7[1] = $2;
            $dist7 = $dist;
         } elsif ($dist<$dist8 || $dist8<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $pt8[0] = $1;
            $pt8[1] = $2;
            $dist8 = $dist;
         } elsif ($dist<$dist9 || $dist9<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $pt9[0] = $1;
            $pt9[1] = $2;
            $dist9 = $dist;
         } elsif ($dist<$dist10 || $dist10<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $pt10[0] = $1;
            $pt10[1] = $2;
            $dist10 = $dist;
         } elsif ($dist<$dist11 || $dist11<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $pt11[0] = $1;
            $pt11[1] = $2;
            $dist11 = $dist;
         } elsif ($dist<$dist12 || $dist12<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $pt12[0] = $1;
            $pt12[1] = $2;
            $dist12 = $dist;
         } elsif ($dist<$dist13 || $dist13<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $pt13[0] = $1;
            $pt13[1] = $2;
            $dist13 = $dist;
         } elsif ($dist<$dist14 || $dist14<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $pt14[0] = $1;
            $pt14[1] = $2;
            $dist14 = $dist;
         } elsif ($dist<$dist15 || $dist15<0.0) {
            @pt16 = @pt15;
            $dist16 = $dist15;
            $pt15[0] = $1;
            $pt15[1] = $2;
            $dist15 = $dist;
         } elsif ($dist<$dist16 || $dist16<0.0) {
            $pt16[0] = $1;
            $pt16[1] = $2;
            $dist16 = $dist;
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
$pt5[0] $pt5[1]
$pt6[0] $pt6[1]
$pt7[0] $pt7[1]
$pt8[0] $pt8[1]
$pt9[0] $pt9[1]
$pt10[0] $pt10[1]
$pt11[0] $pt11[1]
$pt12[0] $pt12[1]
$pt13[0] $pt13[1]
$pt14[0] $pt14[1]
$pt15[0] $pt15[1]
$pt16[0] $pt16[1]
EOF
   close(SKYFILE2);
   
   open(TWOSPECTCONFIG, ">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
   print TWOSPECTCONFIG<<EOF;
fmin 401.25
fspan 0.25
Tobs 40551300
Tcoh 1800
SFToverlap 900
ihsfar 1.0e-14
ihsfomfar 1.0
tmplfar 1.0e-18
Pmin 7200
Pmax 8110260
dfmin 0.0002
dfmax 0.1
skyRegionFile /local/user/egoetz/$$/skygrid2.dat
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength 500
sftDir /local/user/egoetz/$$
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/IFOmismatch/$jobnum
sftType standard
IFO H1
FFTplanFlag 1
fastchisqinv
useSSE
outfilename logfile_H1$ii.txt
ULfilename uls_H1$ii.dat
configCopy input_copy_H1$ii.conf
keepOnlyTopNumIHS 5
EOF
   close(TWOSPECTCONFIG);
   
   system("/home/egoetz/opt/lscsoft/bin/lalapps_Makefakedata_v4 @/local/user/egoetz/$$/mfdconfig");
   die "system lalapps_Makefakedata_v4 failed: $?" if $?;
   
   system("/home/egoetz/opt/lscsoft/bin/lalapps_TwoSpect --config=/local/user/egoetz/$$/twospectconfig");
   die "system lalapps_TwoSpect failed: $?" if $?;
   
   system("rm /local/user/egoetz/$$/*.sft");
   die "rm failed: $?" if $?;

   $mfdrandseed = int(rand(1000000));
   
   open(MFDCONFIG,">/local/user/egoetz/$$/mfdconfig") or die "Cannot write to /local/user/egoetz/$$/mfdconfig $!";
   print MFDCONFIG<<EOF;
outSFTbname /local/user/egoetz/$$/testsfts.sft
outSingleSFT TRUE
IFO L1
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
timestampsFile /home/egoetz/TwoSpect/IFOmismatch/timestamps.dat
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
randSeed $mfdrandseed
EOF
   close(MFDCONFIG);

   $dist1 = -1.0;
   $dist2 = -1.0;
   $dist3 = -1.0;
   $dist4 = -1.0;
   $dist5 = -1.0;
   $dist6 = -1.0;
   $dist7 = -1.0;
   $dist8 = -1.0;
   $dist9 = -1.0;
   $dist10 = -1.0;
   $dist11 = -1.0;
   $dist12 = -1.0;
   $dist13 = -1.0;
   $dist14 = -1.0;
   $dist15 = -1.0;
   $dist16 = -1.0;
   @pt1 = (0, 0);
   @pt2 = (0, 0);
   @pt3 = (0, 0);
   @pt4 = (0, 0);
   @pt5 = (0, 0);
   @pt6 = (0, 0);
   @pt7 = (0, 0);
   @pt8 = (0, 0);
   @pt9 = (0, 0);
   @pt10 = (0, 0);
   @pt11 = (0, 0);
   @pt12 = (0, 0);
   @pt13 = (0, 0);
   @pt14 = (0, 0);
   @pt15 = (0, 0);
   @pt16 = (0, 0);
   open(SKYFILE, $skygrid2) or die "Cannot open $skygrid2 $!";
   while(my $line=<SKYFILE>) {
      if($line =~ /^(\d+.\d+) (-?\d+.\d+)/) {

         my $dist = acos(sin(abs($2-0.5*pi))*sin(abs($delta-0.5*pi))*cos($1-$alpha)+cos(abs($2-0.5*pi))*cos(abs($delta-0.5*pi)));

         if($dist<$dist1 || $dist1<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            @pt4 = @pt3;
            @pt3 = @pt2;
            @pt2 = @pt1;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $dist4 = $dist3;
            $dist3 = $dist2;
            $dist2 = $dist1;
            $pt1[0] = $1;
            $pt1[1] = $2;
            $dist1 = $dist;
         } elsif ($dist<$dist2 || $dist2<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            @pt4 = @pt3;
            @pt3 = @pt2;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $dist4 = $dist3;
            $dist3 = $dist2;
            $pt2[0] = $1;
            $pt2[1] = $2;
            $dist2 = $dist;
         } elsif ($dist<$dist3 || $dist3<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            @pt4 = @pt3;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $dist4 = $dist3;
            $pt3[0] = $1;
            $pt3[1] = $2;
            $dist3 = $dist;
         } elsif ($dist<$dist4 || $dist4<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            @pt5 = @pt4;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $dist5 = $dist4;
            $pt4[0] = $1;
            $pt4[1] = $2;
            $dist4 = $dist;
         } elsif ($dist<$dist5 || $dist5<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            @pt6 = @pt5;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $dist6 = $dist5;
            $pt5[0] = $1;
            $pt5[1] = $2;
            $dist5 = $dist;
         } elsif ($dist<$dist6 || $dist6<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            @pt7 = @pt6;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $dist7 = $dist6;
            $pt6[0] = $1;
            $pt6[1] = $2;
            $dist6 = $dist;
         } elsif ($dist<$dist7 || $dist7<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            @pt8 = @pt7;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $dist8 = $dist7;
            $pt7[0] = $1;
            $pt7[1] = $2;
            $dist7 = $dist;
         } elsif ($dist<$dist8 || $dist8<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            @pt9 = @pt8;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $dist9 = $dist8;
            $pt8[0] = $1;
            $pt8[1] = $2;
            $dist8 = $dist;
         } elsif ($dist<$dist9 || $dist9<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            @pt10 = @pt9;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $dist10 = $dist9;
            $pt9[0] = $1;
            $pt9[1] = $2;
            $dist9 = $dist;
         } elsif ($dist<$dist10 || $dist10<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            @pt11 = @pt10;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $dist11 = $dist10;
            $pt10[0] = $1;
            $pt10[1] = $2;
            $dist10 = $dist;
         } elsif ($dist<$dist11 || $dist11<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            @pt12 = @pt11;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $dist12 = $dist11;
            $pt11[0] = $1;
            $pt11[1] = $2;
            $dist11 = $dist;
         } elsif ($dist<$dist12 || $dist12<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            @pt13 = @pt12;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $dist13 = $dist12;
            $pt12[0] = $1;
            $pt12[1] = $2;
            $dist12 = $dist;
         } elsif ($dist<$dist13 || $dist13<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            @pt14 = @pt13;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $dist14 = $dist13;
            $pt13[0] = $1;
            $pt13[1] = $2;
            $dist13 = $dist;
         } elsif ($dist<$dist14 || $dist14<0.0) {
            @pt16 = @pt15;
            @pt15 = @pt14;
            $dist16 = $dist15;
            $dist15 = $dist14;
            $pt14[0] = $1;
            $pt14[1] = $2;
            $dist14 = $dist;
         } elsif ($dist<$dist15 || $dist15<0.0) {
            @pt16 = @pt15;
            $dist16 = $dist15;
            $pt15[0] = $1;
            $pt15[1] = $2;
            $dist15 = $dist;
         } elsif ($dist<$dist16 || $dist16<0.0) {
            $pt16[0] = $1;
            $pt16[1] = $2;
            $dist16 = $dist;
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
$pt5[0] $pt5[1]
$pt6[0] $pt6[1]
$pt7[0] $pt7[1]
$pt8[0] $pt8[1]
$pt9[0] $pt9[1]
$pt10[0] $pt10[1]
$pt11[0] $pt11[1]
$pt12[0] $pt12[1]
$pt13[0] $pt13[1]
$pt14[0] $pt14[1]
$pt15[0] $pt15[1]
$pt16[0] $pt16[1]
EOF
   close(SKYFILE2);
   
   open(TWOSPECTCONFIG, ">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
   print TWOSPECTCONFIG<<EOF;
fmin 401.25
fspan 0.25
Tobs 40551300
Tcoh 1800
SFToverlap 900
ihsfar 1.0e-14
ihsfomfar 1.0
tmplfar 1.0e-18
Pmin 7200
Pmax 8110260
dfmin 0.0002
dfmax 0.1
skyRegionFile /local/user/egoetz/$$/skygrid2.dat
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength 500
sftDir /local/user/egoetz/$$
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/IFOmismatch/$jobnum
sftType standard
IFO L1
FFTplanFlag 1
fastchisqinv
useSSE
outfilename logfile_L1$ii.txt
ULfilename uls_L1$ii.dat
configCopy input_copy_L1$ii.conf
keepOnlyTopNumIHS 5
EOF
   close(TWOSPECTCONFIG);
   
   system("/home/egoetz/opt/lscsoft/bin/lalapps_Makefakedata_v4 @/local/user/egoetz/$$/mfdconfig");
   die "system lalapps_Makefakedata_v4 failed: $?" if $?;
   
   system("/home/egoetz/opt/lscsoft/bin/lalapps_TwoSpect --config=/local/user/egoetz/$$/twospectconfig");
   die "system lalapps_TwoSpect failed: $?" if $?;
   
   system("rm /local/user/egoetz/$$/*.sft");
   die "rm failed: $?" if $?;
   
}

system("rm -rf /local/user/egoetz/$$");
die "rm failed: $?" if $?;

