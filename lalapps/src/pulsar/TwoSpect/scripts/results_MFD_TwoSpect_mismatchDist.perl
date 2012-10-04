#!/usr/bin/perl

use strict;
use warnings;

open(RESULTS, ">/home/egoetz/TwoSpect/mismatchDist/results.dat") or die "Cannot write to /home/egoetz/TwoSpect/mismatchDist/results.dat $!";

for(my $ii=0; $ii<100; $ii++) {
   open(INJECTEDVALS, "/home/egoetz/TwoSpect/mismatchDist/$ii/injections.dat") or die "Cannot open /home/egoetz/TwoSpect/mismatchDist/$ii/injections.dat $!";
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
      
      open(LOGFILE, "/home/egoetz/TwoSpect/mismatchDist/$ii/logfile_$jj.txt") or die "Cannot open /home/egoetz/TwoSpect/mismatchDist/$ii/logfile_$jj.txt $!";
      my @twospectout = reverse <LOGFILE>;
      my $kk = 1;
      while($twospectout[$kk] =~ /^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)/) {
         print RESULTS "$jj $ii $1 $2 $3 $6 $7 $8\n";
         $kk++;
      }
      close(LOGFILE);

      open(LOGFILE2, "/home/egoetz/TwoSpect/mismatchDist/$ii/logfile_$jj_2.txt") or die "Cannot open /home/egoetz/TwoSpect/mismatchDist/$ii/logfile_$jj_2.txt $!";
      @twospectout = reverse <LOGFILE2>;
      my $kk = 1;
      while($twospectout[$kk] =~ /^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)/) {
         print RESULTS "$jj $ii $1 $2 $3 $6 $7 $8\n";
         $kk++;
      }
      close(LOGFILE2);
      
      $jj++;
   }
}

close(RESULTS);
