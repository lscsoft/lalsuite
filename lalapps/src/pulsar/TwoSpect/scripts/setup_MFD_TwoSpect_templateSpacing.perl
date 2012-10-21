#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/templateSpacing/dag") or die "Cannot write to /home/egoetz/TwoSpect/templateSpacing/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/templateSpacing/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/templateSpacing/condor") or die "Cannot write to /home/egoetz/TwoSpect/templateSpacing/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/templateSpacing/run_MFD_TwoSpect_templateSpacing.perl
input=/dev/null
output=/home/egoetz/TwoSpect/templateSpacing/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/templateSpacing/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/templateSpacing.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
