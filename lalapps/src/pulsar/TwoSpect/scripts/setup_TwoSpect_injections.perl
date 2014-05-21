#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

setup_TwoSpect_injections - script to setup injection studies

=head1 SYNOPSIS

setup_TwoSpect_injections.perl [options]

=head1 OPTIONS

=item B<help:>
Show the help information

=item B<dir:>
Base directory to store the job subdirectories [R]

=item B<jobs:>
Number of jobs to execute [R]

=item B<logfile:>
Name of condor logfile path/filename [R]

=item B<fmin:>
Minimum frequency in Hz of injection [R]

=item B<fspan:>
Frequency span in Hz of injection band (0.25)

=item B<h0min:>
Minimum strain value to inject

=item B<h0max:>
Maximum strain value to inject

=item B<h0dist:>
Distribution of h0, 0 = linear, 1 = log, -1 = inv. log (1)

=item B<h0val:>
Constant strain value to inject

=item B<injPol:>
Polarizations, 0 = linear, 1 = random, 2 = circular (2)

=item B<injskyra:>
Right ascension of injection

=item B<injskydec:>
Declination of injection

=item B<injPmin:>
Minimum period of injection (7200)

=item B<injPmax:>
Maximum period of injection (8110260)

=item B<periodDist:>
Orbit period dist. 0 = linear, 1 = log, -1 = inv. log (0)

=item B<injDfmin:>
Minimum modulation depth of injection in Hz (0.0)

=item B<injDfmax:>
Maximum modulation depth of injection (0.1)

=item B<injDfExpAllow:>
Allow larger Df values by factor (1)

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

=item B<sftFile:>
Noise SFT file (may be multiple)

=item B<sftType:>
SFT type (may be multiple)

=item B<gaussianNoiseWithSFTGaps:>
Gaussian noise from the gaps of SFT data

=item B<timestampsfile:>
File of start times of SFTs to produce (may be multiple)

=item B<segmentfile:>
File of <start end> times to produce SFTs (may be multiple)

=item B<skylocations:>
Number of sky locations to search, 0 exact location (0)

=item B<ihsfactor:>
IHS folding factor (5)

=item B<seed:>
Seed value for producing injections and Gaussian noise (42)

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

my $directory = '';
my $jobs = 1;
my $logfile = '';

my $dur = 40551300.0;
my $Tsft = 1800.0;
my $SFToverlap = 900.0;

my @sftFile = ();
my @sftType = ();
my $gaussianNoiseWithSFTGaps = 0;
my @timestampsfile = ();
my @segmentfile = ();

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

#H1 t0 = 931081500
#L1 t0 = 931113900
#V1 t0 = 931131900

GetOptions('help' => \$help, 
           'dir=s' => \$directory, 
           'jobs=i' => \$jobs,
           'logfile=s' => \$logfile,
           'sftFile:s' => \@sftFile,
           'sftType:s' => \@sftType,
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
           'Tcoh:f' => \$Tsft,
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
my $numberSFTtypes = @sftType;
my $numberT0 = @t0;

die "Jobs must be 1 or more" if $jobs<1;
die "Number of IFOs and concatenated SFT files must be the same" if $numIFOs!=$numberSFTfiles;
die "Must provide SFT file when gaussianNoiseWithSFTGaps is specifed" if $numberSFTfiles<1 && $gaussianNoiseWithSFTGaps!=0;
die "Only choose sftFile OR sftFile and gaussianNoiseWithSFTgaps OR timestampsfile OR segmentfile" if (($numberSFTfiles>0 && ($numberTimestampFiles>0 || $numberSegmentFiles>0)) || ($numberTimestampFiles>0 && $numberSegmentFiles>0));
die "Number of IFOs and timestamp files must be the same" if ($numberTimestampFiles>0 && $numIFOs!=$numberTimestampFiles);
die "Number of IFOs and segment files must be the same" if ($numberSegmentFiles>0 && $numIFOs!=$numberSegmentFiles);
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

if ($dur!=40551300.0 && $injPmax>0.2*$dur) { $injPmax = $0.2*$dur; }

system("mkdir $directory/err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">$directory/dag") or die "Cannot write to $directory/dag $!";
for (my $ii=0; $ii<$jobs; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii $directory/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $directory/$ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">$directory/condor") or die "Cannot write to $directory/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/atlas/user/atlas3/egoetz/lalsuite-master/lalapps/src/pulsar/TwoSpect/scripts/run_TwoSpect_injections.perl
input=/dev/null
output=/dev/null
error=$directory/err/err.\$(JOBNUM)
log=$logfile
request_memory=2500
notification=Never
EOF

print CONDORFILE "arguments=\"--dir=$directory --jobnum=\$(JOBNUM) --Tobs=$dur --Tcoh=$Tsft --SFToverlap=$SFToverlap --fmin=$fmin --fspan=$fspan --injPol=$injPol";

if ($h0val ne "") { print CONDORFILE " --h0val=$h0val"; }
else { print CONDORFILE " --h0min=$h0min --h0max=$h0max --h0dist=$h0dist"; }

if ($scox1switch==1) { print CONDORFILE " --scox1"; }
else {
   print CONDORFILE " --injPmin=$injPmin --injPmax=$injPmax --periodDist=$periodDist --injDfmin=$injDfmin --injDfmax=$injDfmax";
   if ($injskyra ne "" && $injskydec ne "") { print CONDORFILE " --injskyra=$injskyra --injskydec=$injskydec"; }
   if ($maxEcc > 0.0) { print CONDORFILE " --minEcc=$minEcc --maxEcc=$maxEcc --eccDist=$eccDist"; }
   if ($minSpindown != 0.0 || $maxSpindown!=0.0) { print CONDORFILE " --minSpindown=$minSpindown --maxSpindown=$maxSpindown --spindownDist=$spindownDist"; }
}

for (my $ii=0; $ii<$numIFOs; $ii++) {
   print CONDORFILE " --ifo=$ifo[$ii] --t0=$t0[$ii]";
   if ($numberSFTfiles>0) {
      print CONDORFILE " --sftFile=$sftFile[$ii]";
      print CONDORFILE " --sftType=$sftType[$ii]";
   } elsif ($numberTimestampFiles>0) { print CONDORFILE " --timestampsfile=$timestampsfile[$ii]"; }
   elsif ($numberSegmentFiles>0) { print CONDORFILE " --segmentfile=$segmentfile[$ii]"; }
}

if ($gaussianNoiseWithSFTGaps!=0) { print CONDORFILE " --gaussianNoiseWithSFTGaps"; }

if ($skylocations!=0) { print CONDORFILE " --skylocations=$skylocations"; }

if ($ihsfactor!=5) { print CONDORFILE " --ihsfactor=$ihsfactor"; }

if ($seedstart!=42) { print CONDORFILE " --seed=$seedstart"; }

if ($weightedIHS!=0) { print CONDORFILE " --weightedIHS"; }

if ($markBadSFTs!=0) { print CONDORFILE " --markBadSFTs"; }

if ($ulonly!=0) { print CONDORFILE " --ulonly"; }

if ($templateTest!=0) { print CONDORFILE " --templateTest"; }

if ($ihsfar!=1.0) { print CONDORFILE " --ihsfar=$ihsfar"; }
if ($ihsfomfar!=1.0) { print CONDORFILE " --ihsfomfar=$ihsfomfar"; }
if ($tmplfar!=1.0) { print CONDORFILE " --tmplfar=$tmplfar"; }
if ($tmplLength!=500) { print CONDORFILE " --tmplLength=$tmplLength"; }
if ($injDfExpansionAllowance!=1.0) { print CONDORFILE " --injDfExpAllow=$injDfExpansionAllowance"; }

print CONDORFILE "\"\n";

print CONDORFILE "queue\n";
close(CONDORFILE);
