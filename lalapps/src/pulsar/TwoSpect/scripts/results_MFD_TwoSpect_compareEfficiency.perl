#!/usr/bin/perl

use strict;
use warnings;

open(RESULTS, ">/home/egoetz/TwoSpect/compareEfficiency/results.dat") or die "Cannot write to /home/egoetz/TwoSpect/compareEfficiency/results.dat $!";

for(my $ii=0; $ii<200; $ii++) {
   for(my $jj=0; $jj<5; $jj++) {
      open(RECOVERED, "/home/egoetz/TwoSpect/compareEfficiency/$ii/logfile_$jj.txt") or die "Cannot open /home/egoetz/TwoSpect/compareEfficiency/$ii/logfile_$jj.txt $!";
      my @twospectout = reverse <RECOVERED>;
      
      my $kk = 1;
      while($twospectout[$kk] =~ /^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (-?\d+.\d+), h0 = (0|\d+.\d+e-\d+), Prob = (-?\d+.\d+), TF norm = (\d+.\d+e\+\d+)/) {
         print RESULTS "$jj $ii $1 $2 $3 $6 $7 $8\n";
         $kk++;
      }
      
      close(RECOVERED);
   }
}

close(RESULTS);
