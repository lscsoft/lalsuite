#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $directory = '';
my $jobs = 0;
my $noiseswitch = 0;
my $gapsswitch = 0;
my $randpolswitch = 0;
my $linpolswitch = 0;
my $eccentricityswitch = 0;
my $spindownswitch = 0;
my $ifo = '';
my $fmin = '';
my $h0min = '';
my $skylocations = 0;
my $timestampsfile = '';
my $injskyra = '';
my $injskydec = '';
my $ihsfactor = 5;
my $seedstart = 42;
my $scox1switch = 0;
GetOptions('dir=s' => \$directory, 'jobs=i' => \$jobs, 'realnoise' => \$noiseswitch, 'gaps' => \$gapsswitch, 'randpol' => \$randpolswitch, 'linpol' => \$linpolswitch, 'eccOrbit' => \$eccentricityswitch, 'spindown' => \$spindownswitch, 'ifo=s' => \$ifo, 'fmin=f' => \$fmin, 'h0min=f' => \$h0min, 'skylocations:i' => \$skylocations, 'timestampsfile:s' => \$timestampsfile, 'injskyra:f' => \$injskyra, 'injskydec:f' => \$injskydec, 'ihsfactor:i' => \$ihsfactor, 'seed:i' => \$seedstart, 'scox1' => \$scox1switch);

system("mkdir $directory/out");
die "mkdir failed: $?" if $?;
system("mkdir $directory/err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">$directory/dag") or die "Cannot write to $directory/dag $!";
for(my $ii=0; $ii<$jobs; $ii++) {
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
executable=$directory/run_TwoSpect_UL.perl
input=/dev/null
output=$directory/out/out.\$(JOBNUM)
error=$directory/err/err.\$(JOBNUM)
log=/local/user/egoetz/UL_test.log
request_memory=2500
notification=Never
EOF

print CONDORFILE "arguments=\"--dir=$directory --jobnum=\$(JOBNUM) --ifo=H1 --fmin=$fmin --h0min=$h0min";
if ($noiseswitch!=0) {
   print CONDORFILE " --realnoise";
}
if ($gapsswitch!=0) {
   print CONDORFILE " --gaps";
}
if ($randpolswitch!=0) {
   print CONDORFILE " --randpol";
}
if ($linpolswitch!=0) {
   print CONDORFILE " --linpol";
}
if ($eccentricityswitch!=0) {
   print CONDORFILE " --eccOrbit";
}
if ($spindownswitch!=0) {
   print CONDORFILE " --spindown";
}
if ($skylocations!=0) {
   print CONDORFILE " --skylocations=$skylocations";
}
if ($timestampsfile~='') {
   print CONDORFILE " --timestampsfile=$timestampsfile";
}
if ($injskyra!='' && $injskydec!='') {
   print CONDORFILE " --injskyra=$injskyra --injskydec=$injskydec";
}
if ($ihsfactor!=5) {
   print CONDORFILE " --ihsfactor=$ihsfactor";
}
if ($seedstart!=42) {
   print CONDORFILE " --seed=$seedstart";
}
if ($scox1switch!=0) {
   print CONDORFILE " --scox1";
}
print CONDORFILE "\"\n";

print CONDORFILE "queue\n";
close(CONDORFILE);
