#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/compareEfficiency/dag") or die "Cannot write to /home/egoetz/TwoSpect/compareEfficiency/dag $!";
for(my $ii=0; $ii<200; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/compareEfficiency/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/compareEfficiency/condor") or die "Cannot write to /home/egoetz/TwoSpect/compareEfficiency/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/compareEfficiency/run_MFD_TwoSpect_compareEfficiency.perl
input=/dev/null
output=/home/egoetz/TwoSpect/compareEfficiency/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/compareEfficiency/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/compareEfficiency.log
request_memory=3000
notification=Never
queue
EOF
close(CONDORFILE);
