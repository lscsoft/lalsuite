#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Math::Trig qw(pi);
use Math::Random;
use POSIX;

=head1 NAME

run_TwoSpect_injections - script to perform injection studies

=head1 SYNOPSIS

run_TwoSpect_injections.perl [options]

=head1 OPTIONS

=item B<help:>
Show the help information

=item B<directory:>
Base directory to store the job subdirectories [R]

=item B<jobnum:>
Job number, also the base directory/subdirectory path [R]

=item B<fmin:>
Minimum frequency in Hz of injection [R]

=item B<fspan:>
Frequency span in Hz of injection band (0.25)

=item B<h0min:>
Minimum strain value to inject

=item B<h0max:>
Maximum strain value to inject

=item B<h0val:>
Constant strain value to inject

=item B<injPol:>
Polarizations, 0 = linear, 1 = random, 2 = circular (2)

=item B<injskyra:>
Right ascension of injection

=item B<injskydec:>
Declination of injection

=item B<injPmin:>
Minimum period of injection

=item B<injPmax:>
Maximum period of injection

=item B<periodDist:>
Orbit period dist. 0 = linear, 1 = log, -1 = inv. log (0)

=item B<injDfmin:>
Minimum modulation depth of injection in Hz (0.0)

=item B<injDfmax:>
Maximum modulation depth of injection (0.1)

=item B<minEcc:>
Minimum of eccentricity of orbit (0.0)

=item B<maxEcc:>
Maximum of eccentricity of orbit (0.0)

=item B<eccDist:>
Ecc. distribution, 0 = linear, 1 = log, -1 = inv. log (0)

=item B<minSpindown:>
Minimum spindown of source (0.0)

=item B<maxSpindown:>
Maximum spindown of source (0.0)

=item B<spindownDist:>
Spindown distribution, 0 = linear (0)

=item B<ifo:>
Interferometer to use (may be multiple) [R]

=item B<t0:>
Start time of search (may be multiple) [R]

=item B<Tobs:>
Observation time of the search in seconds (40551300)

=item B<Tsft:>
Coherence length of each SFT (1800)

=item B<SFToverlap:>
Overlap of each SFT in seconds (900)

=item B<SFTnoise:>
Use noise from SFTs

=item B<gaussianNoiseWithSFTGaps:>
Gaussian noise from the gaps of SFT data

=item B<timestampsfile:>
File containing start times of SFTs to produce

=item B<segmentfile:>
File containing <start end> of segments to produce SFTs

=item B<skylocations:>
Number of sky locations to search, 0 exact location (0)

=item B<ihsfactor:>
IHS folding factor (5)

=item B<seed:>
Seed value for producing injections and Gaussian noise

=item B<scox1:>
Inject Sco X-1 signals (flag, off)

=item B<weightedIHS:>
Use noise-weighted IHS statistic (flag, off)

=item B<ulonly:>
Only produce ULs from IHS statistic, no follow up (flag, off)

=item B<templateTest:>
Brute force template test (flag, off)

=item B<ihsfar:>
IHS false alarm rate (1.0)

=item B<ihsfomfar:>
IHS figure of merit false alarm rate (1.0)

=item B<tmplfar:>
Template statistic false alarm rate (1.0)

=item B<tmplLength:>
Maximum length of a template (500)

=item B<markBadSFTs:>
Mark and remove bad SFTs (flag, off)

=cut

my $help = 0;

my $jobnum = 0;

my @sftFile = ();
my $gaussianNoiseWithSFTGaps = 0;
my @timestampsfile = ();
my @segmentfile = ();

my $dur = 40551300.0;
my $Tsft = 1800.0;
my $SFToverlap = 900.0;

my $injPol = 2;
my $minEcc = 0;
my $maxEcc = 0;
my $eccDist = 0;
my $minSpindown = 0;
my $maxSpindown = 0;
my $spindownDist = 0;
my $fmin = '';
my $fspan = 0.25;
my $h0min = '';
my $h0max = '';
my $h0dist = 1;
my $h0val = '';
my $injskyra = '';
my $injskydec = '';
my $injPmin = 7200;
my $injPmax = 0.2*$dur;
my $periodDist = 0;
my $injDfmin = 0;
my $injDfmax = 0.1;
my $injDfExpansionAllowance = 1.0;
my $scox1switch = 0;

my @ifo = ();
my @t0 = (); 
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

#H1 t0 = 931081500
#L1 t0 = 931113900
#V1 t0 = 931131900

GetOptions('help' => \$help,
           'dir=s' => \$directory,
           'jobnum=i' => \$jobnum,
           'sftFile:s' => \@sftFile,
           'gaussianNoiseWithSFTGaps' => \$gaussianNoiseWithSFTGaps,
           'injPol:i' => \$injPol,
           'minEcc:f' => \$minEcc,
           'maxEcc:f' => \$maxEcc,
           'eccDist:i' => \$eccDist,
           'minSpindown:f' => \$minSpindown, 
           'maxSpindown:f' => \$maxSpindown, 
           'spindownDist:i' => \$spindownDist, 
           'ifo=s' => \@ifo,
           't0=f' => \@t0,
           'Tobs:f' => \$dur,
           'Tsft:f' => \$Tsft,
           'SFToverlap:f' => \$SFToverlap,
           'fmin=f' => \$fmin,
           'fspan:f' => \$fspan,
           'h0min:f' => \$h0min, 
           'h0max:f' => \$h0max, 
           'h0dist:i' => \$h0dist, 
           'h0val:f' => \$h0val, 
           'skylocations:i' => \$skylocations, 
           'timestampsfile:s' => \@timestampsfile, 
           'segmentfile:s' => \@segmentfile, 
           'injskyra:f' => \$injskyra, 
           'injskydec:f' => \$injskydec, 
           'injPmin:f' => \$injPmin,
           'injPmax:f' => \$injPmax, 
           'periodDist:i' => $periodDist,
           'injDfmin:f' => \$injDfmin,
           'injDfmax:f' => \$injDfmax,
           'injDfExpAllow:f' => \$injDfExpansionAllowance,
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
my $numberSFTfiles = @sftFile;
my $numberT0 = @t0;

die "Number of IFOs and concatenated SFT files must be the same" if $numIFOs!=$numberSFTfiles;
die "Must provide SFT file when gaussianNoiseWithSFTGaps is specifed" if $numberSFTfiles<1 && $gaussianNoiseWithSFTGaps!=0;
die "Only choose sftFile OR sftFile and gaussianNoiseWithSFTgaps OR timestampsfile OR segmentfile" if (($numberSFTfiles>0 && ($numberTimestampFiles>0 || $numberSegmentFiles>0)) || ($numberTimestampFiles>0 && $numberSegmentFiles>0));
die "Number of IFOs and timestamp files must be the same" if $numIFOs!=$numberTimestampFiles;
die "Number of IFOs and segment files must be the same" if $numIFOs!=$numberSegmentFiles;
die "injPol must be 0, 1, or 2" if $injPol<0 || $injPol>2;
die "Minimum eccentricity must be 0 or larger" if $minEcc<0.0;
die "Maximum eccentricity must be smaller than 1" if $maxEcc>=1.0;
die "eccDist must be 0, 1, or -1" if $eccDist<-1 || $eccDist>1;
die "spindownDist must be 0" if $spindownDist!=0;
die "Number of IFOs must be 3 or less" if $numIFOs>3;
die "Number of IFOs and t0 must be the same" if $numIFOs!=$numberT0;
die "Must specify both h0min and h0max" if (($h0min ne "" && $h0max eq "") || ($h0max ne "" && $h0min eq ""));
die "h0dist must be 0, 1, or -1" if $h0dist<-1 || $h0dist>1;
die "h0val cannot be specified with h0min or h0max" if ($h0val ne "" && ($h0min ne "" || $h0max ne ""));
die "skylocations must be 0 or greater" if $skylocations<0;
die "Need both injskyra and injskydec" if (($injskyra ne "" && $injskydec eq "") || ($injskyra eq "" && $injskydec ne ""));
die "periodDist must be 0, 1, or -1" if $periodDist<-1 || $periodDist>1;
die "ihsfactor must be 1 or greater" if $ihsfactor<1;
die "injDfExpAllow must be 1 or greater" if $injDfExpansionAllowance<1.0;
die "seed must be greater than 0" if $seedstart<1;
die "ihsfar must be less than or equal to 1" if $ihsfar>1.0;
die "ihsfomfar must be less than or equal to 1" if $ihsfomfar>1.0;
die "tmplfar must be less than or equal to 1" if $tmplfar>1.0;
die "tmplLength must be 1 or greater" if $tmplLength<1;

my @seedval = ($seedstart, $jobnum+1);
random_set_seed(@seedval);

my $scoX1P = 68023.70;
my $scoX1asini = 1.44;

system("mkdir /local/user/egoetz/$$");
die "mkdir failed: $?" if $?;

my @ifokey = ();
my @ifos = ();
my @skygridfile = ();
foreach my $ifoval (@ifo) {
   if ($ifoval eq "LHO" || $ifoval eq "H1") {
      push(@ifos, "LHO");
      push(@ifokey, "H1");
   } elsif ($ifoval eq "LLO" || $ifoval eq "L1") {
      push(@ifos, "LLO");
      push(@ifokey, "L1");
   } elsif ($ifoval eq "Virgo" || $ifoval eq "V1") {
      push(@ifos, "Virgo");
      push(@ifokey, "V1");
   } else { die "Interferometer not recognized!"; }

   if ($skylocations>=1) {
      push(@skygridfile, "/local/user/egoetz/$$/skygrid-$ifokey[-1].dat");
      system("/atlas/user/atlas3/egoetz/lalsuite-master/lalapps/src/pulsar/TwoSpect/skygridsetup --fmin=$fmin --fspan=$fspan --IFO=$ifokey[-1] --Tcoh=$Tsft --SFToverlap=$SFToverlap --t0=$t0[-1] --Tobs=$dur --v2 --outfilename=$skygridfile[-1]");
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
   my $alpha = 0.0;
   my $delta = 0.0;
   my $df = 0.0;
   my $P = 0.0;
   my $asini = 0.0;
   my $ecc = 0.0;
   my $argp = 0.0;
   my $f1dot = 0.0;
   my $cosi = 1.0;

   if ($h0val ne "") { $h0 = $h0val; }
   else {
      if ($h0dist==0) { $h0 = ($h0max - $h0min)*random_uniform() + $h0min; }
      elsif ($h0dist==1) { $h0 = 10**((log10($h0max)-log10($h0min))*random_uniform()) * $h0min; }
      else { $h0 = ($h0max + $h0min) - 10**((log10($h0max)-log10($h0min))*random_uniform()) * $h0min; }
   }

   my $psi = 0.5*pi*random_uniform()-0.25*pi;
   my $phi0 = 2.0*pi*random_uniform();

   if ($injskyra eq "" && $injskydec eq "") {
      $alpha = 2.0*pi*random_uniform();
      $delta = acos(random_uniform(1, -1, 1))-0.5*pi;
   } else {
      $alpha = $injskyra;
      $delta = $injskydec;
   }

   my $f0 = $fmin + $fspan*random_uniform();

   if ($scox1switch==0) {
      if ($periodDist==0) { $P = ($injPmax - $injPmin)*random_uniform() + $injPmin; }
      elsif ($periodDist==1) { $P = 10**((log10($injPmax)-log10($injPmin))*random_uniform()) * $injPmin; }
      else { $P = ($injPmax + $injPmin) - 10**((log10($injPmax)-log10($injPmin))*random_uniform()) * $injPmin; }

      $df = ($injDfmax-$injDfmin)*random_uniform() + $injDfmin;
      while ($df-0.5/$Tsft<1.0e-6 || $df>$injDfExpansionAllowance*$P/(2.0*$Tsft*$Tsft)) { $df = ($injDfmax-$injDfmin)*random_uniform() + $injDfmin; }

      $asini = $df*$P/2.0/pi/$f0;

      if ($maxEcc > 0.0) {
         if ($eccDist==0) { $ecc = ($maxEcc - $minEcc)*random_uniform() + $minEcc; }
         elsif ($eccDist==1) { $ecc = 10**((log10($maxEcc)-log10($minEcc))*random_uniform()) * $minEcc; }
         else { $ecc = ($maxEcc + $minEcc) - 10**((log10($maxEcc)-log10($minEcc))*random_uniform()) * $minEcc; }
         $argp = 2.0*pi*random_uniform();
      }

      if ($minSpindown != 0.0 || $maxSpindown!=0.0) {
         if ($spindownDist==0) { $f1dot = ($maxSpindown - $minSpindown)*random_uniform() + $minSpindown; }
         else { die "Spindowns other than linear distributions are not yet implemented." }
      }
   } else {
      $alpha = 4.275699238500;
      $delta = -0.272973858335;
      $P = $scoX1P + random_normal(1, 0, 0.0432);
      $asini = $scoX1asini + random_normal(1, 0, 0.18);
      $df = 2.0*pi*$f0*$asini/$P;
   }

   if ($injPol==0) { $cosi = 0.0; }
   elsif ($injPol==1) { $cosi = random_uniform(1, -1, 1); }

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
fspan $fspan
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
EOF

      if ($numberSFTfiles>0) { print TWOSPECTCONFIG "sftFile $sftFile[$jj]\n"; }
      elsif ($timestampsfile[$jj] ne "") { print TWOSPECTCONFIG "timestampsFile $timestampsfile[$jj]\n"; } 
      elsif ($segmentfile[$jj] ne "") { print TWOSPECTCONFIG "segmentFile $segmentfile[$jj]\n"; }

      if ($numberSFTfiles==0) {
         my $mfdrandseed = random_uniform_integer(1, 0, 1000000);
         print TWOSPECTCONFIG "injRandSeed $mfdrandseed\n";
      }

      if ($markBadSFTs!=0) { print TWOSPECTCONFIG "markBadSFTs\n"; }
      if ($weightedIHS!=0) { print TWOSPECTCONFIG "weightedIHS\n"; }
      if ($ulonly!=0) { print TWOSPECTCONFIG "IHSonly\n"; }

      if ($skylocations==0) { print TWOSPECTCONFIG "skyRegion ($alpha,$delta)\n"; }
      elsif ($skylocations>0) {
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
         for (my $kk=0; $kk<$skylocations; $kk++) { print SKYFILE2 "$ras[$sortedindexvalues[$kk]] $decs[$sortedindexvalues[$kk]]\n"; }
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

