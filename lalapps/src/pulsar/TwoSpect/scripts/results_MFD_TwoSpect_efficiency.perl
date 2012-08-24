#!/usr/bin/perl

use strict;
use warnings;

open(RESULTS, ">/home/egoetz/TwoSpect/efficiency/results.dat") or die "Cannot write to /home/egoetz/TwoSpect/efficiency/results.dat $!";

for(my $ii=0; $ii<300; $ii++) {
   open(INJECTEDVALS, "/home/egoetz/TwoSpect/efficiency/$ii/injections.dat") or die "Cannot open /home/egoetz/TwoSpect/efficiency/$ii/injections.dat $!";
   my @injections = reverse <INJECTEDVALS>;
   my @injections2;
   my $jj;
   for($jj=0; $jj<10; $jj++) {
      push(@injections2, $injections[$jj]);
   }
   close(INJECTEDVALS);
   @injections2 = reverse @injections2;
   
   $jj = 0;
   foreach my $injection (@injections2) {
      chomp($injection);
      
      open(RECOVERED, "/home/egoetz/TwoSpect/efficiency/$ii/logfile_$jj.txt") or die "Cannot open /home/egoetz/TwoSpect/efficiency/$ii/logfile_$jj.txt $!";
      my @twospectout = reverse <RECOVERED>;
      
      my $kk = 1;
      while($twospectout[$kk] =~ /^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)/) {
         print RESULTS "$jj $ii $injection $1 $2 $3 $6 $7 $8\n";
         $kk++;
      }
      
      close(RECOVERED);
      $jj++;
   }
}

close(RESULTS);
