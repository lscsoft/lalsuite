#!/usr/bin/perl

use strict;
use warnings;

open(RESULTS, ">/home/egoetz/TwoSpect/UL_spindowns/results.dat") or die "Cannot write to /home/egoetz/TwoSpect/UL_spindowns/results.dat $!";

for(my $ii=0; $ii<100; $ii++) {
   open(INJECTEDVALS, "/home/egoetz/TwoSpect/UL_spindowns/$ii/injections.dat") or die "Cannot open /home/egoetz/TwoSpect/UL_spindowns/$ii/injections.dat $!";
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
      
      open(ULFILE, "/home/egoetz/TwoSpect/UL_spindowns/$ii/uls_$jj.dat") or die "Cannot open /home/egoetz/TwoSpect/UL_spindowns/$ii/uls_$jj.dat $!";
      while(my $line = <ULFILE>) {
         chomp($line);
         print RESULTS "$injection $line\n";
      }
     
      close(ULFILE);
      $jj++;
   }
}

close(RESULTS);
