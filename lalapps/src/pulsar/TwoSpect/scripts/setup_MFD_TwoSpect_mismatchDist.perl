#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/mismatchDist/dag") or die "Cannot write to /home/egoetz/TwoSpect/mismatchDist/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/mismatchDist/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/mismatchDist/condor") or die "Cannot write to /home/egoetz/TwoSpect/mismatchDist/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/mismatchDist/run_MFD_TwoSpect_mismatchDist.perl
input=/dev/null
output=/home/egoetz/TwoSpect/mismatchDist/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/mismatchDist/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/mismatchDist.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
