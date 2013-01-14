#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/IFOmismatch/dag") or die "Cannot write to /home/egoetz/TwoSpect/IFOmismatch/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/IFOmismatch/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/IFOmismatch/condor") or die "Cannot write to /home/egoetz/TwoSpect/IFOmismatch/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/IFOmismatch/run_MFD_TwoSpect_IFOmismatch.perl
input=/dev/null
output=/home/egoetz/TwoSpect/IFOmismatch/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/IFOmismatch/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/IFOmismatch.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
