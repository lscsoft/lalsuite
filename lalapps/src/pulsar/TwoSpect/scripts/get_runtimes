#!/usr/bin/perl

use strict;
use warnings;
use Date::Parse;

my $directory = $ARGV[0];
my $analysisdate = $ARGV[1];
my @jobfolders = <$directory/$analysisdate/output/*>;
my $numberofjobs = @jobfolders;

open(RUNTIMES, ">$directory/$analysisdate/runtimes.dat") or die "Cannot write to $directory/$analysisdate/runtimes.dat $!";

for(my $ii=0; $ii<$numberofjobs; $ii++) {
   system("head -n 1 $directory/$analysisdate/output/$ii/logfile.txt > $directory/$analysisdate/output/$ii/starttime");
   die "head failed: $?" if $?;
   system("tail -n 1 $directory/$analysisdate/output/$ii/logfile.txt > $directory/$analysisdate/output/$ii/endtime");
   die "tail failed: $?" if $?;
   
   open(STARTTIME, "$directory/$analysisdate/output/$ii/starttime") or die "Cannot open $directory/$analysisdate/output/$ii/starttime $!";
   open(ENDTIME, "$directory/$analysisdate/output/$ii/endtime") or die "Cannot open $directory/$analysisdate/output/$ii/endtime $!";
   my $start = 0;
   my $end = 0;
   while (my $line = <STARTTIME>) {
      if ($line =~ /^Program lalapps_TwoSpect 1.1.24 executed on (.*)/) {
         $start = str2time($1);
      }
   }
   while (my $line = <ENDTIME>) {
      if ($line =~ /^Program finished on (.*)/) {
         $end = str2time($1);
      }
   }
   close(STARTTIME);
   close(ENDTIME);
   
   if ($start != 0 && $end != 0) {
      my $difference = $end-$start;
      print RUNTIMES "$ii $difference\n";
   }
   
   system("rm $directory/$analysisdate/output/$ii/starttime");
   die "rm failed: $?" if $?;
   system("rm $directory/$analysisdate/output/$ii/endtime");
   die "rm failed: $?" if $?;
   
}

close(RUNTIMES);
