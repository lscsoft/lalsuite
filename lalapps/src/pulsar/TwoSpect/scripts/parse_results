#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $directory = '';

GetOptions('dir=s' => \$directory);

my @jobfolders = <$directory/output/*>;
my $numberofjobs = @jobfolders;

open(MISSING, ">$directory/missingbands.dat") or die "Cannot write to $directory/missingbands.dat $!";
open(CANDIDATES, ">$directory/candidates.dat") or die "Cannot write to $directory/candidates.dat $!";

for(my $ii=0; $ii<$numberofjobs; $ii++) {
   my $ulfilename = "$directory/output/$ii/uls.dat";
   if (-e $ulfilename) {
      system("cat $ulfilename >> $directory/ULresults.dat");
      die "cat failed: $?" if $?;
      
      open(LOGFILE, "$directory/output/$ii/logfile.txt") or die "Cannot open $directory/output/$ii/logfile.txt $!";
      if (-z "$directory/output/$ii/logfile.txt") {
         print STDERR "Missing $ii!\n";
         system("echo JOB A$ii $directory/condor >> $directory/dag4missing");
         die "echo failed: $?" if $?;
         system("echo VARS A$ii PID=\\\"$ii\\\" >> $directory/dag4missing");
         die "echo failed: $?" if $?;
      } else {
         my @lines = reverse <LOGFILE>;
         if ( $lines[0] !~ /^Program finished/ ) {
            print STDERR "Missing $ii!\n";
            system("echo JOB A$ii $directory/condor >> $directory/dag4missing");
            die "echo failed: $?" if $?;
            system("echo VARS A$ii PID=\\\"$ii\\\" >> $directory/dag4missing");
            die "echo failed: $?" if $?;
         } else {
            my $jj = 1;
            while ($lines[$jj] =~ /^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)/) {
               print CANDIDATES "$1 $2 $3 $4 $5 $6 $7 $8 $9 $ii\n";
               $jj++
            }
         }
      }
      close(LOGFILE);

   } else {
      open(CONFIG, "$directory/output/$ii/input_values.conf") or die "Cannot open $directory/output/$ii/input_values.conf $!";
      my $ulfmin = 0.0;
      my $ulfspan = 0.0;
      while (my $line=<CONFIG>) {
         if ($line =~ /^fmin="(.*?)"/) {
            $ulfmin = $1; 
         } elsif ($line =~ /^fspan="(.*?)"/) {
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
