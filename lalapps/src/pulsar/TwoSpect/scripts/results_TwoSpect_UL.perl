#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $directory = '';
my $jobs = 0;
my $output = '';
GetOptions('dir=s' => \$directory, 'jobs=i' => \$jobs, 'output=s' => \$output);

open(RESULTS, ">$directory/$output") or die "Cannot write to $directory/$output $!";

for(my $ii=0; $ii<$jobs; $ii++) {
   open(INJECTEDVALS, "$directory/$ii/injections.dat") or die "Cannot open $directory/$ii/injections.dat $!";
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
      
      open(ULFILE, "$directory/$ii/uls_$jj.dat") or die "Cannot open $directory/$ii/uls_$jj.dat $!";
      while(my $line = <ULFILE>) {
         chomp($line);
         print RESULTS "$injection $line\n";
      }
      
      close(ULFILE);
      $jj++;
   }
}

close(RESULTS);
