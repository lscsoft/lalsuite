#!/usr/bin/perl

use strict;
use warnings;

system("mkdir out");
die "mkdir failed: $?" if $?;
system("mkdir err");
die "mkdir failed: $?" if $?;

open(DAGFILE,">/home/egoetz/TwoSpect/circTemplatesEccOrbits/dag") or die "Cannot write to /home/egoetz/TwoSpect/circTemplatesEccOrbits/dag $!";
for(my $ii=0; $ii<100; $ii++) {
   print DAGFILE<<EOF;
JOB A$ii /home/egoetz/TwoSpect/circTemplatesEccOrbits/condor
VARS A$ii JOBNUM="$ii"
EOF
   
   system("mkdir $ii");
   die "mkdir failed: $?" if $?;
}
close(DAGFILE);

open(CONDORFILE,">/home/egoetz/TwoSpect/circTemplatesEccOrbits/condor") or die "Cannot write to /home/egoetz/TwoSpect/circTemplatesEccOrbits/condor $!";
print CONDORFILE<<EOF;
universe=vanilla
executable=/home/egoetz/TwoSpect/circTemplatesEccOrbits/run_MFD_TwoSpect_confirmUL.perl
input=/dev/null
output=/home/egoetz/TwoSpect/circTemplatesEccOrbits/out/out.\$(JOBNUM)
error=/home/egoetz/TwoSpect/circTemplatesEccOrbits/err/err.\$(JOBNUM)
arguments=\$(JOBNUM)
log=/local/user/egoetz/circTemplatesEccOrbits.log
request_memory=3200
notification=Never
queue
EOF
close(CONDORFILE);
