#!/usr/bin/perl

use strict;
use warnings;

for(my $ii=0; $ii<100; $ii++) {
  for(my $jj=0; $jj<10; $jj++) {
    system("cat /home/egoetz/TwoSpect/templateSpacing/$ii/templatespacingout_$jj.dat >> /home/egoetz/TwoSpect/templateSpacing/results.dat");
    die "cat failed: $?" if $?;
  }
}

close(RESULTS);
