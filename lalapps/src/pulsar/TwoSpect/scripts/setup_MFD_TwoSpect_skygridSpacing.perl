#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/skygridSpacing/dag") or die "Cannot write to /home/egoetz/TwoSpect/skygridSpacing/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/skygridSpacing/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/skygridSpacing/condor") or die "Cannot write to /home/egoetz/TwoSpect/skygridSpacing/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/skygridSpacing/run_MFD_TwoSpect_skygridSpacing.perl
input=/dev/null
output=/home/egoetz/TwoSpect/skygridSpacing/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/skygridSpacing/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/skygridSpacing.log
request_memory=2500
notification=Never
queue
EOF
close(CONDORFILE);
