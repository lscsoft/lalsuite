#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/confirmUL/dag") or die "Cannot write to /home/egoetz/TwoSpect/confirmUL/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/confirmUL/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/confirmUL/condor") or die "Cannot write to /home/egoetz/TwoSpect/confirmUL/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/confirmUL/run_MFD_TwoSpect_confirmUL.perl
input=/dev/null
output=/home/egoetz/TwoSpect/confirmUL/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/confirmUL/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/confirmUL.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
