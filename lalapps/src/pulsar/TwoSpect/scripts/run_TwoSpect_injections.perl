#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Math::Trig;
use Math::Random;
use POSIX;

=head1 NAME

run_TwoSpect_injections - script to perform injection studies

=head1 SYNOPSIS

run_TwoSpect_injections.perl [options]

=head1 OPTIONS

=item B<--help:>
Show the help information

=item B<--directory:>
Base directory to store the job subdirectories [R]

=item B<--jobnum:>
Job number, also the base directory/subdirectory path [R]

=item B<--SFTnoise:>
Use noise from SFTs

=item B<--gaussianNoiseWithSFTGaps:>
Gaussian noise from the gaps of SFT data

=item B<--timestampsfile:>
File containing start times of SFTs to produce

=item B<--segmentfile:>
File containing <start end> of segments to produce SFTs

=item B<--randpol:>
Polarization of the injected waveform is random

=item B<--linpol:>
Polarization of the injected waveform is linear

=item B<--eccOrbit:>
Source is in a random eccentric orbit

=item B<--spindown:>
Source has spindown

=item B<--ifo:>
Interferometer to use (may be multiple) [R]

=item B<--fmin:>
Minimum frequency of band [R]

=item B<--h0min:>
Minimum strain value to inject

=item B<--h0val:>
Constant strain value to inject

=item B<--skylocations:>
Number of sky locations to search

=item B<--injskyra:>
Right ascension of injection

=item B<--injskydec:>
Declination of injection

=item B<--injPmin:>
Minimum period of injection

=item B<--injPmax:>
Maximum period of injection

=item B<--injDfmin:>
Minimum modulation depth of injection

=item B<--injDfmax:>
Maximum modulation depth of injection

=item B<--ihsfactor:>
IHS folding factor (typically 5 for all-sky)

=item B<--seed:>
Seed value for producing injections and Gaussian noise

=item B<--scox1:>
Inject Sco X-1 signals

=item B<--weightedIHS:>
Use noise-weighted IHS statistic

=item B<--ulonly:>
Only produce ULs from IHS statistic, no follow up

=item B<--templateTest:>
Brute force template test

=item B<--ihsfar:>
IHS false alarm rate

=item B<--ihsfomfar:>
IHS figure of merit false alarm rate

=item B<--tmplfar:>
Template statistic false alarm rate

=item B<--tmplLength:>
Maximum length of a template

=item B<--markBadSFTs:>
Mark and remove bad SFTs

=cut

my $help = 0;

my $jobnum = 0;

my $SFTnoise = 0;
my $gaussianNoiseWithSFTGaps = 0;
my @timestampsfile = ();
my @segmentfile = ();

my $randpolswitch = 0;
my $linpolswitch = 0;
my $eccentricityswitch = 0;
my $spindownswitch = 0;
my $fmin = '';
my $h0min = '';
my $h0val = '';
my $injskyra = '';
my $injskydec = '';
my $injPmin = 7200;
my $injPmax = 0;
my $injDfmin = 0;
my $injDfmax = 0.1;
my $scox1switch = 0;

my @ifo = (); 
my $skylocations = 0;

my $ihsfactor = 5;
my $seedstart = 42;
my $weightedIHS = 0;
my $ulonly = 0;
my $templateTest = 0;
my $ihsfar = 1.0;
my $ihsfomfar = 1.0;
my $tmplfar = 1.0;
my $tmplLength = 500;
my $markBadSFTs = 0;
my $directory = '';

GetOptions('help' => \$help, 
           'dir=s' => \$directory, 
           'jobnum=i' => \$jobnum, 
           'SFTnoise' => \$SFTnoise, 
           'gaussianNoiseWithSFTGaps' => \$gaussianNoiseWithSFTGaps, 
           'randpol' => \$randpolswitch, 
           'linpol' => \$linpolswitch, 
           'eccOrbit' => \$eccentricityswitch, 
           'spindown' => \$spindownswitch, 
           'ifo=s' => \@ifo, 
           'fmin=f' => \$fmin, 
           'h0min:f' => \$h0min, 
           'h0val:f' => \$h0val, 
           'skylocations:i' => \$skylocations, 
           'timestampsfile:s' => \@timestampsfile, 
           'segmentfile:s' => \@segmentfile, 
           'injskyra:f' => \$injskyra, 
           'injskydec:f' => \$injskydec, 
           'injPmin:f' => \$injPmin,
           'injPmax:f' => \$injPmax,
           'injDfmin:f' => \$injDfmin,
           'injDfmax:f' => \$injDfmax,
           'ihsfactor:i' => \$ihsfactor, 
           'seed:i' => \$seedstart, 
           'scox1' => \$scox1switch, 
           'weightedIHS' => \$weightedIHS, 
           'ulonly' => \$ulonly, 
           'templateTest' => \$templateTest,
           'ihsfar:f' => \$ihsfar, 
           'ihsfomfar:f' => \$ihsfomfar, 
           'tmplfar:f' => \$tmplfar, 
           'tmplLength:i' => \$tmplLength, 
           'markBadSFTs' => \$markBadSFTs) or pod2usage(2);
pod2usage(1) if $help;

my $numIFOs = @ifo;
my $numberTimestampFiles = @timestampsfile;
my $numberSegmentFiles = @segmentfile;

die "No more than 3 interferometers can be used" if $numIFOs > 3;
die "Must specify one of --randpol or --linpol" if ($randpolswitch==1 && $linpolswitch==1);
die "Must specify one of --h0min or --h0val" if (($h0min ne "" && $h0val ne "") || ($h0min eq "" && $h0val eq ""));
die "Must specify one of --SFTnoise or --gaussianNoiseWithSFTGaps or --timestampsfile or --segmentfile" if (($SFTnoise!=0 && ($gaussianNoiseWithSFTGaps!=0 || $numberTimestampFiles>0 || $numberSegmentFiles>0)) || ($gaussianNoiseWithSFTGaps!=0 && ($SFTnoise!=0 || $numberTimestampFiles>0 || $numberSegmentFiles>0)) || ($numberTimestampFiles>0 && ($SFTnoise!=0 || $gaussianNoiseWithSFTGaps!=0 || $numberSegmentFiles>0)) || ($numberSegmentFiles>0 && ($SFTnoise!=0 || $gaussianNoiseWithSFTGaps!=0 || $numberTimestampFiles>0)));
die "Need both --injskyra and --injskydec" if (($injskyra ne "" && $injskydec eq "") ||  ($injskyra eq "" && $injskydec ne ""));

srand($seedstart+$jobnum);
my $Tsft = 1800.0;
my $SFToverlap = 900.0;
my $dur = 40551300.0;
my $injPmax = 0.2*$dur;
my $scoX1P = 68023.70;
my $scoX1asini = 1.44;
my $lowestFneeded = $fmin - 0.1 - 0.1 - 4e-3;
my $lowestFinteger = int(floor($lowestFneeded - fmod($lowestFneeded, 2.0) + 0.5)) + 1;
if ($lowestFinteger>$lowestFneeded) {
   $lowestFinteger -= 2;
}
my $highestFinteger = $lowestFinteger + 3;

#system("mkdir /local/user/egoetz/$$");
system("mkdir /Users/evgoet/$$");
die "mkdir failed: $?" if $?;

my @ifokey = ();
my @ifos = ();
my @t0 = ();
my @sftFile = ();
my @sftType = ();
my @skygridfile = ();
foreach my $ifoval (@ifo) {
   if ($ifoval eq "LHO" || $ifoval eq "H1") {
      push(@t0, 931081500);  #H1 start
      push(@ifos, "LHO");
      push(@ifokey, "H1");
   } elsif ($ifoval eq "LLO" || $ifoval eq "L1") {
      push(@t0, 931113900);  #L1 start
      push(@ifos, "LLO");
      push(@ifokey, "L1");
   } elsif ($ifoval eq "Virgo" || $ifoval eq "V1") {
      push(@t0, 931131900);  #V1 start
      push(@ifos, "Virgo");
      push(@ifokey, "V1");
   }

   if ($lowestFinteger<173 && ($ifokey[-1] eq "H1" || $ifokey[-1] eq "L1")) {
      push(@sftFile, "/atlas/user/atlas3/egoetz/twospect/$ifos[-1]/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz/*.sft");
      if ($SFTnoise!=0) {
         push (@sftType, "vladimir");
      } else {
         push (@sftType, "standard");
      }
   } else {
      push(@sftType, "standard");
      push(@sftFile, "/atlas/user/atlas3/egoetz/twospect/$ifos[-1]/1800s_sfts/$lowestFinteger\-${highestFinteger}Hz/*.sft");
   }

   if ($skylocations>=1) {
      push(@skygridfile, "/local/user/egoetz/$$/skygrid-$ifokey[-1].dat");
      system("/atlas/user/atlas3/egoetz/lalsuite-master/lalapps/src/pulsar/TwoSpect/skygridsetup --fmin=$fmin --fspan=0.25 --IFO=$ifokey[-1] --Tcoh=$Tsft --SFToverlap=$SFToverlap --t0=$t0[-1] --Tobs=$dur --v2 --outfilename=$skygridfile[-1]");
   die "skygridsetup failed: $?" if $?;
   }
}

my $Pmin = $injPmin;
my $Pmax = $injPmax;
my $dfmin = $injDfmin;
my $dfmax = $injDfmax;
if ($scox1switch!=0) {
   $Pmin = $dur/(int($dur/$scoX1P+0.5)+1.0);
   $Pmax = $dur/(int($dur/$scoX1P+0.5)-1.0);
   $dfmin = 2*pi*$fmin*($scoX1asini-3.0*0.18)/($scoX1P+3.0*.0432);
   $dfmax = 2*pi*($fmin+0.25)*($scoX1asini+3.0*0.18)/($scoX1P-3.0*.0432);
}

for(my $ii=0; $ii<10; $ii++) {
   my $h0 = 0.0;
   if ($h0min ne "") {
      $h0 = 10**(1.5*rand())*$h0min;
   } else {
      $h0 = $h0val;
   }
   my $psi = 0.5*pi*rand()-0.25*pi;
   my $phi0 = 2.0*pi*rand();
   my $alpha = 0.0;
   my $delta = 0.0;
   if ($injskyra eq "" && $injskydec eq "") {
      $alpha = 2.0*pi*rand();
      $delta = acos(2.0*rand()-1.0)-0.5*pi;
   } else {
      $alpha = $injskyra;
      $delta = $injskydec;
   }
   my $f0 = $fmin + 0.25*rand();

   my $df = 0.0;
   my $P = 0.0;
   my $asini = 0.0;
   my $ecc = 0.0;
   my $argp = 0.0;
   my $f1dot = 0.0;
   if ($scox1switch==0) {
      $df = ($injDfmax-$injDfmin)*rand() + $injDfmin;
      while ($df-0.5/$Tsft<1.0e-6) {
         $df = ($injDfmax-$injDfmin)*rand() + $injDfmin;
      }
      $P = rand()*($injPmax-$injPmin) + $injPmin;
      while ($P<2.0*$df*$Tsft*$Tsft) {
         $P = rand()*($injPmax-$injPmin) + $injPmin;
      }
      $asini = $df*$P/2.0/pi/$f0;
      if ($eccentricityswitch!=0) {
         $ecc = 10**(rand()*5.0-5.05);
         $argp = 2.0*pi*rand();
      }
      if ($spindownswitch!=0) {
         $f1dot = -10**(rand()*2.2-11.25);
      }
   } else {
      $alpha = 4.275699238500;
      $delta = -0.272973858335;
      $P = $scoX1P + random_normal(1, 0, 0.0432);
      $asini = 1.44 + random_normal(1, 0, 0.18);
      $df = 2.0*pi*$f0*$asini/$P;
   }

   my $cosi = 1.0;
   if ($randpolswitch!=0) {
      $cosi = rand()*2.0 - 1.0;
   } elsif ($linpolswitch!=0) {
      $cosi = 0.0;
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

   for (my $jj=0; $jj<$numIFOs; $jj++) {
      open(TWOSPECTCONFIG, ">/local/user/egoetz/$$/twospectconfig") or die "Cannot write to /local/user/egoetz/$$/twospectconfig $!";
      print TWOSPECTCONFIG<<EOF;
fmin $fmin
fspan 0.25
Tobs $dur
Tcoh $Tsft
SFToverlap $SFToverlap
ihsfar $ihsfar
ihsfomfar $ihsfomfar
tmplfar $tmplfar
Pmin $Pmin
Pmax $Pmax
dfmin $dfmin
dfmax $dfmax
t0 $t0[$jj]
blksize 101
avesqrtSh 1.0e-23
minTemplateLength 1
maxTemplateLength $tmplLength
outdirectory $directory/$jobnum
IFO $ifokey[$jj]
FFTplanFlag 1
fastchisqinv
useSSE
outfilename ${ifokey[$jj]}logfile_$ii.txt
ULfilename ${ifokey[$jj]}uls_$ii.dat
configCopy ${ifokey[$jj]}input_copy_$ii.conf
injectionSources \@/local/user/egoetz/$$/mfdconfig
ihsfactor $ihsfactor
sftType $sftType[$jj]
EOF

      if ($SFTnoise!=0 || $gaussianNoiseWithSFTGaps!=0) {
         print TWOSPECTCONFIG "sftFile $sftFile[$jj]\n";
      } elsif ($timestampsfile[$jj] ne "") {
         print TWOSPECTCONFIG "timestampsFile $timestampsfile[$jj]\n";
      } elsif ($segmentfile[$jj] ne "") {
         print TWOSPECTCONFIG "segmentFile $segmentfile[$jj]\n";
      }
      if ($SFTnoise==0) {
         my $mfdrandseed = int(rand(1000000));
         print TWOSPECTCONFIG "injRandSeed $mfdrandseed\n";
      }
      if ($markBadSFTs!=0) {
         print TWOSPECTCONFIG "markBadSFTs\n";
      }
      if ($weightedIHS!=0) {
         print TWOSPECTCONFIG "weightedIHS\n";
      }
      if ($ulonly!=0) {
         print TWOSPECTCONFIG "IHSonly\n";
      }
      if ($skylocations==0) {
         print TWOSPECTCONFIG "skyRegion ($alpha,$delta)\n";
      } elsif ($skylocations>0) {
         print TWOSPECTCONFIG "skyRegionFile /local/user/egoetz/$$/skygrid2.dat\n";

         open(SKYFILE, $skygridfile[$jj]) or die "Cannot open $skygridfile[$jj] $!";
         my @distances = ();
         my @ras = ();
         my @decs = ();
         while (my $line=<SKYFILE>) {
            if ($line =~ /^(\d+.\d+) (-?\d+.\d+)/) {
               my $dist = acos(sin(abs($2-0.5*pi))*sin(abs($delta-0.5*pi))*cos($1-$alpha)+cos(abs($2-0.5*pi))*cos(abs($delta-0.5*pi)));
               push(@ras, $1);
               push(@decs, $2);
               push(@distances, $dist);
            }
         }
         close(SKYFILE);
         my @sortedindexvalues = sort {$distances[$a] <=> $distances[$b]} 0 .. $#distances;

         open(SKYFILE2,">/local/user/egoetz/$$/skygrid2.dat") or die "Cannot write to /local/user/egoetz/$$/skygrid2.dat $!";
         for (my $kk=0; $kk<$skylocations; $kk++) {
            print SKYFILE2 "$ras[$sortedindexvalues[$kk]] $decs[$sortedindexvalues[$kk]]\n";
         }
         close(SKYFILE2);
      }

      if ($templateTest!=0) {
         print TWOSPECTCONFIG "bruteForceTemplateTest\n";
         print TWOSPECTCONFIG "templateTestF $f0\n";
         print TWOSPECTCONFIG "templateTestP $P\n";
         print TWOSPECTCONFIG "templateTestDf $df\n";
      }

      close(TWOSPECTCONFIG);
   
      system("/home/egoetz/opt/lscsoft-master/bin/lalapps_TwoSpect --config=/local/user/egoetz/$$/twospectconfig");
      die "system lalapps_TwoSpect failed: $?" if $?;
   }

}

system("rm -rf /local/user/egoetz/$$");
die "rm failed: $?" if $?;

