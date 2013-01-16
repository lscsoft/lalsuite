#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

my $jobnum = $ARGV[0];
srand(42+$jobnum);

system("mkdir /local/user/egoetz/$$");
die "mkdir failed: $?" if $?;

chdir("/local/user/egoetz/$$") or die "Cannot change to directory /local/user/egoetz/$$ $!";

my $h0ul = 4.5709e-24;
my $h0 = $h0ul;
my $psi = sprintf("%.6f",0.5*pi*rand()-0.25*pi);
my $phi0 = sprintf("%.6f",2.0*pi*rand());
my $alpha = sprintf("%.6f",2.0*pi*rand());
my $delta = sprintf("%.6f",acos(2.0*rand()-1.0)-0.5*pi);
my $f0 = 401.25 + 0.24*rand();
my $df = rand()*0.1;
$df = 0.01;
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
my @ecc = (0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45);
my $argp = sprintf("%.6f", 2.0*pi*rand());
my $randseedval = int(rand(1000000));

for(my $ii=0; $ii<10; $ii++) {
   open(MFDCONFIG,">/local/user/egoetz/$$/mfdconfig") or die "Cannot write to /local/user/egoetz/$$/mfdconfig $!";
   print MFDCONFIG<<EOF;
outSFTbname /local/user/egoetz/$$/testsfts.sft
outSingleSFT TRUE
IFO H1
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
startTime 931081500
duration 40551300
generationMode 0
fmin 401.0
Band 2.9992
Tsft 1800
SFToverlap 900
window Hann
Alpha $alpha
Delta $delta
h0 $h0
cosi 1.0
psi $psi
phi0 $phi0
Freq $f0
orbitasini $asini
orbitEcc $ecc[$ii]
orbitTpSSBsec 900000000
orbitTpSSBnan 0
orbitPeriod $P
orbitArgp $argp
f1dot 0.0
refTime 900000000
noiseSqrtSh 3.0e-23
randSeed $randseedval
EOF
   close(MFDCONFIG);

   open(INJECTION, ">>/home/egoetz/TwoSpect/circTemplatesEccOrbits/$jobnum/injections.dat") or die "Cannot write to /home/egoetz/TwoSpect/circTemplatesEccOrbits/$jobnum/injections.dat $!";
   print INJECTION "$alpha $delta $h0 $psi $phi0 $f0 $P $df $ecc[$ii] $argp\n";
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
outdirectory /home/egoetz/TwoSpect/circTemplatesEccOrbits/$jobnum
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
EOF
   close(TWOSPECTCONFIG);

   system("/home/egoetz/opt/lscsoft/bin/lalapps_Makefakedata_v4 \@/local/user/egoetz/$$/mfdconfig");
   die "system lalapps_Makefakedata_v4 failed: $?" if $?;
   
   system("/home/egoetz/TwoSpect/circTemplatesEccOrbits/TwoSpect_templateTest --config=/local/user/egoetz/$$/twospectconfig");
   die "TwoSpect_templateTest failed: $?" if $?;

   system("rm /local/user/egoetz/$$/*.sft");
   die "rm failed: $?" if $?;

}



