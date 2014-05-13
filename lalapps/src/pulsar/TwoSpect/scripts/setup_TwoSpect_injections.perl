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
Name of logfile to be stored at /local/user/egoetz [R]

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

my $directory = '';
my $jobs = 0;
my $logfile = '';

my $SFTnoise = 0;
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
my $injPmax = 0;
my $periodDist = 0;
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

GetOptions('help' => \$help, 
           'dir=s' => \$directory, 
           'jobs=i' => \$jobs,
           'logfile=s' => \$logfile,
           'SFTnoise' => \$SFTnoise, 
           'gaussianNoiseWithSFTGaps' => \$gaussianNoiseWithSFTGaps, 
           'injPol:i' => \$injPol, 
           'minEcc:f' => \$minEcc, 
           'maxEcc:f' => \$maxEcc, 
           'eccDist:i' => \$eccDist, 
           'minSpindown:f' => \$minSpindown, 
           'maxSpindown:f' => \$maxSpindown, 
           'spindownDist:i' => \$spindownDist, 
           'ifo=s' => \@ifo, 
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

system("mkdir $directory/out");
die "mkdir failed: $?" if $?;
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
output=$directory/out/out.\$(JOBNUM)
error=$directory/err/err.\$(JOBNUM)
log=/local/user/egoetz/$logfile
request_memory=2500
notification=Never
EOF

print CONDORFILE "arguments=\"--dir=$directory --jobnum=\$(JOBNUM) --fmin=$fmin --fspan=$fspan --injPol=$injPol";

if ($h0val ne "") { print CONDORFILE " --h0val=$h0val"; }
else { print CONDORFILE " --h0min=$h0min --h0max=$h0max --h0dist=$h0dist"; }

if ($scox1switch==1) { print CONDORFILE " --scox1"; }
else {
   print CONDORFILE " --injPmin=$injPmin --injPmax=$injPmax --periodDist=$periodDist --injDfmin=$injDfmin --injDfmax=$injDfmax";
   if ($injskyra ne "" && $injskydec ne "") { print CONDORFILE " --injskyra=$injskyra --injskydec=$injskydec"; }
   if ($maxEcc > 0.0) { print CONDORFILE " --minEcc=$minEcc --maxEcc=$maxEcc --eccDist=$eccDist"; }
   if ($minSpindown != 0.0 || $maxSpindown!=0.0) { print CONDORFILE " --minSpindown=$minSpindown --maxSpindown=$maxSpindown --spindownDist=$spindownDist"; }
}

for (my $ii=0; $ii<$numIFOs; $ii++) { print CONDORFILE " --ifo=$ifo[$ii]"; }

if ($SFTnoise!=0) { print CONDORFILE " --SFTnoise"; }
elsif ($gaussianNoiseWithSFTGaps!=0) { print CONDORFILE " --gaussianNoiseWithSFTGaps"; }
elsif ($numberTimestampFiles>0) {
   for (my $ii=0; $ii<$numberTimestampFiles; $ii++) { print CONDORFILE " --timestampsfile=$timestampsfile[$ii]"; }
} elsif ($numberSegmentFiles>0) {
   for (my $ii=0; $ii<$numberSegmentFiles; $ii++) { print CONDORFILE " --segmentfile=$segmentfile[$ii]"; }
}

if ($skylocations!=0) { print CONDORFILE " --skylocations=$skylocations"; }

if ($ihsfactor!=5) { print CONDORFILE " --ihsfactor=$ihsfactor"; }

if ($seedstart!=42) { print CONDORFILE " --seed=$seedstart"; }

if ($weightedIHS!=0) { print CONDORFILE " --weightedIHS"; }

if ($markBadSFTs!=0) { print CONDORFILE " --markBadSFTs"; }

if ($ulonly!=0) { print CONDORFILE " --IHSonly"; }

if ($templateTest!=0) { print CONDORFILE " --templateTest"; }

if ($ihsfar!=1.0) { print CONDORFILE " --ihsfar=$ihsfar"; }
if ($ihsfomfar!=1.0) { print CONDORFILE " --ihsfomfar=$ihsfomfar"; }
if ($tmplfar!=1.0) { print CONDORFILE " --tmplfar=$tmplfar"; }
if ($tmplLength!=500) { print CONDORFILE " --tmplLength=$tmplLength"; }

print CONDORFILE "\"\n";

print CONDORFILE "queue\n";
close(CONDORFILE);
