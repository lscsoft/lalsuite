#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/UL_randPolarizations/dag") or die "Cannot write to /home/egoetz/TwoSpect/UL_randPolarizations/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/UL_randPolarizations/condor
VARS A$ii JOBNUM="$ii"
RETRY A$ii 5
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/UL_randPolarizations/condor") or die "Cannot write to /home/egoetz/TwoSpect/UL_randPolarizations/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/UL_randPolarizations/run_MFD_TwoSpect_UL_randPolarizations.perl
input=/dev/null
output=/home/egoetz/TwoSpect/UL_randPolarizations/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/UL_randPolarizations/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/UL_randPolarizations.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
