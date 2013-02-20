#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

my $jobnum = $ARGV[0];
srand(424242+$jobnum);

my $pidvalue = $$;

system("mkdir /local/user/egoetz/$$");
die "mkdir failed: $?" if $?;

my $Tsft = 1800.0;
my $dur = 40551300.0;
my $skygrid = "/home/egoetz/TwoSpect/UL_randPolarizations/skygrid.dat";
for(my $ii=0; $ii<10; $ii++) {
   my $cosi = sprintf("%.6f", 2.0*rand()-1.0);
   my $h0 = sprintf("%.6e",10**(rand()-24.0));
   my $psi = sprintf("%.6f",0.5*pi*rand()-0.25*pi);
   my $phi0 = sprintf("%.6f",2.0*pi*rand());
   my $alpha = sprintf("%.6f",2.0*pi*rand());
   my $delta = sprintf("%.6f",acos(2.0*rand()-1.0)-0.5*pi);
   my $f0 = 401.25 + 0.25*rand();
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
   
   open(MFDCONFIG,">/local/user/egoetz/$pidvalue/mfdconfig") or die "Cannot write to /local/user/egoetz/$pidvalue/mfdconfig $!";
   print MFDCONFIG<<EOF;
outSFTbname /local/user/egoetz/$pidvalue/testsfts.sft
outSingleSFT TRUE
IFO H1
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
timestampsFile /home/egoetz/TwoSpect/UL_randPolarizations/timestamps.dat
generationMode 0
fmin 401.0
Band 2.9992
Tsft 1800
window Hann
Alpha $alpha
Delta $delta
h0 $h0
cosi $cosi
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
   
   open(INJECTION, ">>/home/egoetz/TwoSpect/UL_randPolarizations/$jobnum/injections.dat") or die "Cannot write to /home/egoetz/TwoSpect/UL_randPolarizations/$jobnum/injections.dat $!";
   print INJECTION "$alpha $delta $h0 $psi $phi0 $cosi $f0 $P $df\n";
   close(INJECTION);

   open(SKYFILE, $skygrid) or die "Cannot open $skygrid $!";
   my @distances = ();
   my @ras = ();
   my @decs = ();
   while(my $line=<SKYFILE>) {
      if($line =~ /^(\d+.\d+) (-?\d+.\d+)/) {
         my $dist = acos(sin(abs($2-0.5*pi))*sin(abs($delta-0.5*pi))*cos($1-$alpha)+cos(abs($2-0.5*pi))*cos(abs($delta-0.5*pi)));
         push(@ras, $1);
         push(@decs, $2);
         push(@distances, $dist);
      }
   }
   close(SKYFILE);
   my @sortedindexvalues = sort {$distances[$a] <=> $distances[$b]} 0 .. $#distances;

   open(SKYFILE2,">/local/user/egoetz/$pidvalue/skygrid2.dat") or die "Cannot write to /local/user/egoetz/$pidvalue/skygrid2.dat $!";
   for(my $jj=0; $jj<10; $jj++) {
      print SKYFILE2 "$ras[$sortedindexvalues[$jj]] $decs[$sortedindexvalues[$jj]]\n";
   }
   close(SKYFILE2);
   
   open(TWOSPECTCONFIG, ">/local/user/egoetz/$pidvalue/twospectconfig") or die "Cannot write to /local/user/egoetz/$pidvalue/twospectconfig $!";
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
skyRegionFile /local/user/egoetz/$pidvalue/skygrid2.dat
t0 931081500
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength 500
sftDir /local/user/egoetz/$pidvalue
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory /home/egoetz/TwoSpect/UL_randPolarizations/$jobnum
sftType standard
IFO H1
FFTplanFlag 1
fastchisqinv
useSSE
outfilename logfile_$ii.txt
ULfilename uls_$ii.dat
configCopy input_copy_$ii.conf
IHSonly
EOF
   close(TWOSPECTCONFIG);
   
   system("/home/egoetz/opt/lscsoft/bin/lalapps_Makefakedata_v4 @/local/user/egoetz/$pidvalue/mfdconfig");
   die "system lalapps_Makefakedata_v4 failed: $?" if $?;
   
   system("/home/egoetz/opt/lscsoft/bin/lalapps_TwoSpect --config=/local/user/egoetz/$pidvalue/twospectconfig");
   die "system lalapps_TwoSpect failed: $?" if $?;
   
   system("rm /local/user/egoetz/$pidvalue/*.sft");
   die "rm failed: $?" if $?;
   
}

system("rm -rf /local/user/egoetz/$pidvalue");
die "rm failed: $?" if $?;

