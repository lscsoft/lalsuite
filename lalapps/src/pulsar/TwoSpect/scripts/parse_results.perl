#!/usr/bin/perl

use strict;
use warnings;

my $directory = $ARGV[0];
my $analysisdate = $ARGV[1];

open(MISSING, ">$directory/$analysisdate/missingbands.dat") or die "Cannot write to $directory/$analysisdate/missingbands.dat $!";
open(CANDIDATES, ">$directory/$analysisdate/candidates.dat") or die "Cannot write to $directory/$analysisdate/candidates.dat $!";

for(my $ii=0; $ii<; $ii++) {
   my $ulfilename = "$directory/$analysisdate/output/$ii/uls.dat";
   if (-e $ulfilename) {
      system("cat $directory/$analysisdate/output/$ii/uls.dat >> $directory/$analysisdate/ULresults.dat");
      die "cat failed: $?" if $?;
      
      open(LOGFILE, "$directory/$analysisdate/output/$ii/logfile.txt") or die "Cannot open $directory/$analysisdate/output/$ii/logfile.txt $!";
      @lines = reverse LOGFILE;
      my $jj = 1;
      while ($lines[$jj] =~ /^fsig/) {
         print CANDIDATES $lines[$jj];
         $jj++
      }
      close(LOGFILE);
      
   } else {
      open(CONFIG, "$directory/$analysisdate/output/$ii/input_values.conf") or die "Cannot open $directory/$analysisdate/output/$ii/input_values.conf $!";
      while (my $line=<CONFIG>) {
         my $ulfmin = 0.0;
         my $ulfspan = 0.0;
         if ($line =~ /^ULfmin="(\d+.?\d*)"/) {
            $ulfmin = $1; 
         } elsif ($line =~ /^ULfspan="(\d+.?\d*)"/) {
            $ulfspan = $1;
         }
         if ($ulfmin != 0.0 && $ulfspan != 0.0) {
            print MISSING "$ii $ulfmin $ulfspan\n";
            last;
         }
      }
      close(CONFIG);
   }
}

close(MISSING);
close(CANDIDATES);
