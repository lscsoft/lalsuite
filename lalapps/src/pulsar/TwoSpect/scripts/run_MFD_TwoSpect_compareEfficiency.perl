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
for(my $ii=0; $ii<10; $ii++) {
   my $h0 = $h0ul*$h0fact[$ii];
   my $psi = 0.384099;
   my $phi0 = 2.153257;
   my $alpha = 0.697969;
   my $delta = 0.155953;
   my $f0 = 401.269467;
   my $df = 0.0856440;
   my $P = 4051874.730676;
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
startTime 931081500
duration 40551300
SFToverlap 900
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
Pmin 4051874.730676
Pmax 8110260
dfmin 0.085644
dfmax 0.1
skyRegion ($alpha,$delta)
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength 500
sftDir /local/user/egoetz/$$
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/efficiency/$jobnum
sftType standard
IFO H1
FFTplanFlag 1
ULfmin 401.269467
ULfspan 0.1
fastchisqinv
useSSE
outfilename logfile_$ii.txt
ULfilename uls_$ii.dat
configCopy input_copy_$ii.conf
keepOnlyTopNumIHS 5
EOF
   close(TWOSPECTCONFIG);
   
   system("/home/egoetz/opt/lscsoft/bin/lalapps_Makefakedata_v4 @/local/user/egoetz/$$/mfdconfig");
   die "system lalapps_Makefakedata_v4 failed: $?" if $?;
   
   system("/home/egoetz/TwoSpect/compareEfficiency/TwoSpect_templateTest --config=/local/user/egoetz/$$/twospectconfig");
   die "system lalapps_TwoSpect failed: $?" if $?;
   
   system("rm /local/user/egoetz/$$/*.sft");
   die "rm failed: $?" if $?;
   
}

system("rm -rf /local/user/egoetz/$$");
die "rm failed: $?" if $?;

