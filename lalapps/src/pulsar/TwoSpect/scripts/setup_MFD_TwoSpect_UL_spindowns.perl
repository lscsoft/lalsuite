#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/UL_spindowns/dag") or die "Cannot write to /home/egoetz/TwoSpect/UL_spindowns/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/UL_spindowns/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/UL_spindowns/condor") or die "Cannot write to /home/egoetz/TwoSpect/UL_spindowns/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/UL_spindowns/run_MFD_TwoSpect_UL_spindowns.perl
input=/dev/null
output=/home/egoetz/TwoSpect/UL_spindowns/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/UL_spindowns/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/UL_spindowns.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
