#!/usr/bin/perl

use strict;
use warnings;

open(RESULTS, ">/home/egoetz/TwoSpect/skygridSpacing/results.dat") or die "Cannot write to /home/egoetz/TwoSpect/skygridSpacing/results.dat $!";

for(my $ii=0; $ii<100; $ii++) {
   open(INJECTEDVALS, "/home/egoetz/TwoSpect/skygridSpacing/$ii/injections.dat") or die "Cannot open /home/egoetz/TwoSpect/skygridSpacing/$ii/injections.dat $!";
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
      
      open(H1RECOVERED, "/home/egoetz/TwoSpect/skygridSpacing/$ii/logfile_H1$jj.txt") or die "Cannot open /home/egoetz/TwoSpect/skygridSpacing/$ii/logfile_H1$jj.txt $!";
      my @twospectout = reverse <H1RECOVERED>;
      
      my $kk = 1;
      while($twospectout[$kk] =~ /^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)/) {
         print RESULTS "$jj $ii $injection $4 $5 $1 $2 $3 $6 $7 $8 H1\n";
         $kk++;
      }
      
      close(H1RECOVERED);
      
      $jj++;
   }
}

close(RESULTS);
