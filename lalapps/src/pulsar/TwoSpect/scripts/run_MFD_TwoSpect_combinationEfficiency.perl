#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

my $jobnum = $ARGV[0];
srand(42+$jobnum);

system("mkdir /local/user/egoetz/$$");
die "mkdir failed: $?" if $?;

my $h0ul = 4.5709e-24;
my @h0fact = (0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.9, 1.1, 1.5, 2.0);
my $Tsft = 1800.0;
my $dur = 40551300.0;
my $skygrid = "/home/egoetz/TwoSpect/combinationEfficiency/skygrid.dat";
my $skygrid2 = "/home/egoetz/TwoSpect/combinationEfficiency/skygrid2.dat";
for(my $ii=0; $ii<10; $ii++) {
   my $h0 = $h0ul*$h0fact[$ii];
   my $psi = sprintf("%.6f",0.5*pi*rand()-0.25*pi);
   my $phi0 = sprintf("%.6f",2.0*pi*rand());
   my $alpha = sprintf("%.6f",2.0*pi*rand());
   my $delta = sprintf("%.6f",acos(2.0*rand()-1.0)-0.5*pi);
   my $f0 = 200.0 + 0.25*rand();
   my $df = rand()*0.1;
   while ($df-0.5/$Tsft<1.0e-6) {
      $df = rand()*0.1;
   }
   my $P = rand()*0.2*($dur-7200.0)+7200.0;
   while ($P<2.0*$df*$Tsft*$Tsft) {
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
timestampsFile /home/egoetz/TwoSpect/combinationEfficiency/timestamps.dat
generationMode 0
fmin 199.0
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
   
   open(INJECTION, ">>/home/egoetz/TwoSpect/combinationEfficiency/$jobnum/injections.dat") or die "Cannot write to /home/egoetz/TwoSpect/combinationEfficiency/$jobnum/injections.dat $!";
   print INJECTION "$alpha $delta $h0 $psi $phi0 $f0 $P $df\n";
   close(INJECTION);

   open(SKYFILE, $skygrid) or die "Cannot open $skygrid $!";
   my $grid_alpha = 0.0;
   my $grid_delta = 0.0;
   my $best_dist = -1.0;
   while(my $line=<SKYFILE>) {
      if($line =~ /^(\d+.\d+) (-?\d+.\d+)/) {
         my $dist = acos(sin(abs($2-0.5*pi))*sin(abs($delta-0.5*pi))*cos($1-$alpha)+cos(abs($2-0.5*pi))*cos(abs($delta-0.5*pi)));
         if($best_dist<0.0) {
            $best_dist = $dist;
            $grid_alpha = $1;
            $grid_delta = $2;
         } elsif ($dist<$best_dist) {
            $best_dist = $dist;
            $grid_alpha = $1;
            $grid_delta = $2;
         }
      }
   }
   close(SKYFILE);
   
   open(TWOSPECTCONFIG, ">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
   print TWOSPECTCONFIG<<EOF;
fmin 200
fspan 0.25
Tobs 40551300
Tcoh 1800
SFToverlap 900
ihsfar 1e-14
ihsfomfar 1
tmplfar 1e-18
Pmin 7200
Pmax 8110260
dfmin 0.0002
dfmax 0.1
skyRegion ($grid_alpha,$grid_delta)
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength 500
sftDir /local/user/egoetz/$$
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/combinationEfficiency/$jobnum
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
timestampsFile /home/egoetz/TwoSpect/combinationEfficiency/timestamps.dat
generationMode 0
fmin 199.0
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
   
   open(SKYFILE, $skygrid2) or die "Cannot open $skygrid2 $!";
   $grid_alpha = 0.0;
   $grid_delta = 0.0;
   $best_dist = -1.0;
   while(my $line=<SKYFILE>) {
      if($line =~ /^(\d+.\d+) (-?\d+.\d+)/) {
         my $dist = acos(sin(abs($2-0.5*pi))*sin(abs($delta-0.5*pi))*cos($1-$alpha)+cos(abs($2-0.5*pi))*cos(abs($delta-0.5*pi)));
         if($best_dist<0.0) {
            $best_dist = $dist;
            $grid_alpha = $1;
            $grid_delta = $2;
         } elsif ($dist<$best_dist) {
            $best_dist = $dist;
            $grid_alpha = $1;
            $grid_delta = $2;
         }
      }
   }
   close(SKYFILE);

   open(TWOSPECTCONFIG, ">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
   print TWOSPECTCONFIG<<EOF;
fmin 200
fspan 0.25
Tobs 40551300
Tcoh 1800
SFToverlap 900
ihsfar 1e-14
ihsfomfar 1
tmplfar 1e-18
Pmin 7200
Pmax 8110260
dfmin 0.0002
dfmax 0.1
skyRegion ($grid_alpha,$grid_delta)
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength 500
sftDir /local/user/egoetz/$$
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/combinationEfficiency/$jobnum
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

