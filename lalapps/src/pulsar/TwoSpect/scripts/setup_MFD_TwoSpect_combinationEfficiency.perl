#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/combinationEfficiency/dag") or die "Cannot write to /home/egoetz/TwoSpect/combinationEfficiency/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/combinationEfficiency/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/combinationEfficiency/condor") or die "Cannot write to /home/egoetz/TwoSpect/combinationEfficiency/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/combinationEfficiency/run_MFD_TwoSpect_combinationEfficiency.perl
input=/dev/null
output=/home/egoetz/TwoSpect/combinationEfficiency/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/combinationEfficiency/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/combinationEfficiency.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
