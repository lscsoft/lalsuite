#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/efficiency/dag") or die "Cannot write to /home/egoetz/TwoSpect/efficiency/dag $!";
for(my $ii=0; $ii<300; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/efficiency/condor
VARS A$ii JOBNUM="$ii"
RETRY A$ii 5
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/efficiency/condor") or die "Cannot write to /home/egoetz/TwoSpect/efficiency/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/efficiency/run_MFD_TwoSpect_efficiency.perl
input=/dev/null
output=/home/egoetz/TwoSpect/efficiency/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/efficiency/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/efficiency.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
