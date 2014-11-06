#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $directory = '';
my $jobs = 0;
my @ifos;
my $ULoutput = '';
my $candidateOutput = '';
GetOptions('dir=s' => \$directory,
           'jobs=i' => \$jobs,
           'ifo=s' => \@ifos,
           'ULoutput:s' => \$ULoutput,
           'candidateOutput:s' => \$candidateOutput);

if ($ULoutput ne "") { open(ULRESULTS, ">$directory/$ULoutput") or die "Cannot write to $directory/$ULoutput $!"; }
if ($candidateOutput ne "") { open(CANDIDATERESULTS, ">$directory/$candidateOutput") or die "Cannot write to $directory/$candidateOutput $!"; }

for(my $ii=0; $ii<$jobs; $ii++) {
   open(INJECTEDVALS, "$directory/$ii/injections.dat") or die "Cannot open $directory/$ii/injections.dat $!";
   my @injections = reverse <INJECTEDVALS>;
   my $injectionLength = @injections;
   next if $injectionLength<10;

   my @injections2;
   my $jj = 0;
   while ($jj<10 && $jj<$injectionLength) {
      push(@injections2, $injections[$jj]);
      $jj++;
   }
   close(INJECTEDVALS);
   @injections2 = reverse @injections2;
   $injectionLength = @injections2;
   next if $injectionLength<10;

   my %repeated;
   my $numberRepeated = 0;
   $jj = 0;
   while ($jj<10 && $jj<$injectionLength) {
      $jj++;
      my $injection = $injections2[$jj-1];
      next unless $repeated{$injection}++;
      $numberRepeated++;
   }
   next if $injectionLength-$numberRepeated<10;

   $jj = 0;
   foreach my $injection (@injections2) {
      chomp($injection);
      foreach my $ifo (@ifos) {
         my $ifoval = '';
         if ($ifo eq "H1" || $ifo eq "LHO") { $ifoval = "0"; }
         elsif ($ifo eq "L1" || $ifo eq "LLO") { $ifoval = "1"; }
         elsif ($ifo eq "V1" || $ifo eq "Virgo") { $ifoval = "2"; }

         if ($ULoutput ne "") {
            open(ULFILE, "$directory/$ii/uls_$jj.dat") or die "Cannot open $directory/$ii/uls_$jj.dat $!";
            while (my $line = <ULFILE>) {
               chomp($line);
               print ULRESULTS "$injection $line $ifoval\n";
            }
            close(ULFILE);
         }

         if ($candidateOutput ne "") {
            open(TWOSPECTOUT, "$directory/$ii/logfile_$jj.txt") or die "Cannot open $directory/$ii/logfile_$jj.txt $!";
            my @twospectoutput = reverse <TWOSPECTOUT>;
            my $kk = 1;
            my $foundoutlier = 0;
            while ($twospectoutput[$kk] =~ /^fsig = (\d+.\d+), period = (\d+.\d+), df = (\d+.\d+), RA = (\d+.\d+), DEC = (-?\d+.\d+), R = (\d+.\d+), h0 = (\d+.\d+e-\d+), Prob = (-\d+.\d+), TF norm = (\d+.\d+e\+\d+)/) {
               print CANDIDATERESULTS "$jj $ii $injection $1 $2 $3 $4 $5 $6 $7 $8 $ifoval\n";
               $kk++;
               $foundoutlier = 1;
            }
            if ($foundoutlier==0 && $twospectoutput[0] =~ /^system lalapps_TwoSpect failed/) { print CANDIDATERESULTS "$jj $ii $injection -1 -1 -1 -1 -1 -1 -1 $ifoval\n"; }
            elsif ($foundoutlier==0) { print CANDIDATERESULTS "$jj $ii $injection NaN NaN NaN NaN NaN NaN NaN NaN $ifoval\n"; }
            close(TWOSPECTOUT);
         }
      }
      $jj++;
   }
}

if ($ULoutput ne "") { close(ULRESULTS); }
if ($candidateOutput ne "") { close(CANDIDATERESULTS); }
