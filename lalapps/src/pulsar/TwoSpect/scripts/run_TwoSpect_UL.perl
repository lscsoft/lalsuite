#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Math::Trig;
use Math::Random;
use POSIX;

my $jobnum = 0;
my $noiseswitch = 0;
my $gapsswitch = 0;
my $randpolswitch = 0;
my $linpolswitch = 0;
my $eccentricityswitch = 0;
my $spindownswitch = 0;
my $ifo = '';
my $fmin = '';
my $h0min = '';
my $h0val = '';
my $skylocations = 0;
my $timestampsfile = '';
my $injskyra = '';
my $injskydec = '';
my $ihsfactor = 5;
my $seedstart = 42;
my $scox1switch = 0;
my $directory = '';
GetOptions('dir=s' => \$directory, 'jobnum=i' => \$jobnum, 'realnoise' => \$noiseswitch, 'gaps' => \$gapsswitch, 'randpol' => \$randpolswitch, 'linpol' => \$linpolswitch, 'eccOrbit' => \$eccentricityswitch, 'spindown' => \$spindownswitch, 'ifo=s' => \$ifo, 'fmin=f' => \$fmin, 'h0min:f' => \$h0min, 'h0val:f' => \$h0val, 'skylocations:i' => \$skylocations, 'timestampsfile:s' => \$timestampsfile, 'injskyra:f' => \$injskyra, 'injskydec:f' => \$injskydec, 'ihsfactor:i' => \$ihsfactor, 'seed:i' => \$seedstart, 'scox1' => \$scox1switch);

die "Must set only one of randpolswitch, linpolswitch, or circpolswitch" if ($randpolswitch==1 && $linpolswitch==1);
die "Must specify one of --h0min or --h0val" if (($h0min ne "" && $h0val ne "") || ($h0min eq "" && $h0val eq ""));

srand($seedstart+$jobnum);

system("mkdir /local/user/egoetz/$$");
die "mkdir failed: $?" if $?;

my $ifokey = "";
my $t0 = 0;
if ($ifo eq "LHO" || $ifo eq "H1") {
   $t0 = 931081500;  #H1 start
   $ifo = "LHO";
   $ifokey = "H1";
} elsif ($ifo eq "LLO" || $ifo eq "L1") {
   $t0 = 931113900;  #L1 start
   $ifo = "LLO";
   $ifokey = "L1";
} elsif ($ifo eq "Virgo" || $ifo eq "V1") {
   $t0 = 931131900;
   $ifo = "Virgo";
   $ifokey = "V1";
}

my $lowestFneeded = $fmin - 0.1 - 0.1 - 4e-3;
my $lowestFinteger = int(floor($lowestFneeded - fmod($lowestFneeded, 2.0) + 0.5)) + 1;
if ($lowestFinteger>$lowestFneeded) {
   $lowestFinteger -= 2;
}
my $highestFinteger = $lowestFinteger + 3;

my $sftFile = "";
my $sftType = "";
if ($lowestFinteger<173) {
   $sftType = "vladimir";
   $sftFile = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz/*.sft"
} else {
   $sftType = "standard";
   $sftFile = "/atlas/user/atlas3/egoetz/twospect/$ifo/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz/*.sft"
}

my $skygridfile = "/local/user/egoetz/$$/skygrid.dat";
if ($skylocations>=1) {
   system("/atlas/user/atlas3/egoetz/lalsuite/lalapps/src/pulsar/TwoSpect/skygridsetup --fmin=$fmin --fspan=0.25 --IFO=$ifokey --outfilename=$skygridfile --ephemDir=/home/egoetz/TwoSpect/S6 --ephemYear=08-11-DE405");
   die "skygridsetup failed: $?" if $?;
}

my $Tsft = 1800.0;
my $dur = 40551300.0;
for(my $ii=0; $ii<10; $ii++) {
   my $h0 = 0.0;
   if ($h0min ne "") {
      $h0 = sprintf("%.6e",10**(1.5*rand())*$h0min);
   } else {
      $h0 = $h0val;
   }
   my $psi = sprintf("%.6f",0.5*pi*rand()-0.25*pi);
   my $phi0 = sprintf("%.6f",2.0*pi*rand());
   my $alpha = 0.0;
   my $delta = 0.0;
   if ($injskyra eq "" && $injskydec eq "") {
      $alpha = sprintf("%.6f",2.0*pi*rand());
      $delta = sprintf("%.6f",acos(2.0*rand()-1.0)-0.5*pi);
   } elsif ($injskyra ne "" && $injskydec ne "") {
      $alpha = $injskyra;
      $delta = $injskydec;
   } else {
      die "Need both --injskyra and --injskydec";
   }
   my $f0 = $fmin + 0.25*rand();

   my $df = 0.0;
   my $P = 0.0;
   my $asini = 0.0;
   if ($scox1switch==0) {
      $df = rand()*0.1;
      while ($df-0.5/$Tsft<1.0e-6) {
         $df = rand()*0.1;
      }
      $P = rand()*0.2*($dur-7200.0)+7200.0;
      while ($P<2.0*$df*$Tsft*$Tsft) {
         $P = rand()*0.2*($dur-7200.0)+7200.0;
      }
      $f0 = sprintf("%.6f", $f0);
      $df = sprintf("%.6f", $df);
      $P = sprintf("%.6f", $P);
      $asini = sprintf("%.6f",$df*$P/2.0/pi/$f0);
   } else {
      $alpha = 4.275699238500;
      $delta = -0.272973858335;
      $P = 68023.7136 + random_normal(1, 0, 0.0432);
      $asini = 1.44 + random_normal(1, 0, 0.18);
      $df = 2.0*pi*$f0*$asini/$P;
   }

   my $mfdrandseed = int(rand(1000000));

   my $cosi = 1.0;
   if ($randpolswitch!=0) {
      $cosi = rand()*2.0 - 1.0;
   } elsif ($linpolswitch!=0) {
      $cosi = 0.0;
   }
   my $ecc = 0.0;
   my $argp = 0.0;
   if ($eccentricityswitch!=0 && $scox1switch==0) {
      $ecc = 10**(rand()*5.0-5.05);
      $argp = 2.0*pi*rand();
   }
   my $f1dot = 0.0;
   if ($spindownswitch!=0 && $scox1switch==0) {
      $f1dot = -10**(rand()*2.2-11.25);
   }

   open(MFDCONFIG,">/local/user/egoetz/$$/mfdconfig") or die "Cannot write to /local/user/egoetz/$$/mfdconfig $!";
   print MFDCONFIG<<EOF;
Alpha $alpha
Delta $delta
h0 $h0
cosi $cosi
psi $psi
phi0 $phi0
Freq $f0
orbitasini $asini
orbitEcc $ecc
orbitTpSSB 900000000
orbitPeriod $P
orbitArgp $argp
f1dot $f1dot
refTime 900000000
EOF
   close(MFDCONFIG);
   
   open(INJECTION, ">>$directory/$jobnum/injections.dat") or die "Cannot write to $directory/$jobnum/injections.dat $!";
   print INJECTION "$alpha $delta $h0 $cosi $psi $phi0 $f0 $asini $ecc $P $argp $f1dot $df\n";
   close(INJECTION);
   
   open(TWOSPECTCONFIG, ">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
   print TWOSPECTCONFIG<<EOF;
fmin $fmin
fspan 0.25
Tobs $dur
Tcoh 1800
SFToverlap 900
ihsfar 1.0
ihsfomfar 1.0
tmplfar 1.0
Pmin 7200
Pmax 8110260
dfmin 0.0002
dfmax 0.1
t0 $t0
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength 500
ephemDir /home/egoetz/TwoSpect/S6
ephemYear 08-11-DE405
outdirectory $directory/$jobnum
IFO $ifokey
FFTplanFlag 1
fastchisqinv
useSSE
outfilename logfile_$ii.txt
ULfilename uls_$ii.dat
configCopy input_copy_$ii.conf
IHSonly
injectionSources \@/local/user/egoetz/$$/mfdconfig
ihsfactor $ihsfactor
EOF

   if ($noiseswitch!=0) {
      print TWOSPECTCONFIG "sftFile $sftFile\n";
      print TWOSPECTCONFIG "sftType $sftType\n";
   } elsif ($noiseswitch==0 && $gapsswitch!=0) {
      print TWOSPECTCONFIG "sftType standard\n";
      print TWOSPECTCONFIG "injRandSeed $mfdrandseed\n";
      print TWOSPECTCONFIG "timestampsFile $timestampsfile\n";
   } else {
      print TWOSPECTCONFIG "sftType standard\n";
      print TWOSPECTCONFIG "injRandSeed $mfdrandseed\n";
   }

   if ($skylocations==0) {
      print TWOSPECTCONFIG "skyRegion ($alpha,$delta)\n";
   } elsif ($skylocations>0) {
      print TWOSPECTCONFIG "skyRegionFile /local/user/egoetz/$$/skygrid2.dat\n";

      open(SKYFILE, $skygridfile) or die "Cannot open $skygridfile $!";
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

      open(SKYFILE2,">/local/user/egoetz/$$/skygrid2.dat") or die "Cannot write to /local/user/egoetz/$$/skygrid2.dat $!";
      for(my $jj=0; $jj<$skylocations; $jj++) {
         print SKYFILE2 "$ras[$sortedindexvalues[$jj]] $decs[$sortedindexvalues[$jj]]\n";
      }
      close(SKYFILE2);
   } else {
      die "Need 1 or more sky locations (though 0 is exact sky location)";
   }

   close(TWOSPECTCONFIG);
   
   system("/home/egoetz/opt/lscsoft-master/bin/lalapps_TwoSpect --config=/local/user/egoetz/$$/twospectconfig");
   die "system lalapps_TwoSpect failed: $?" if $?;
   
}

system("rm -rf /local/user/egoetz/$$");
die "rm failed: $?" if $?;

