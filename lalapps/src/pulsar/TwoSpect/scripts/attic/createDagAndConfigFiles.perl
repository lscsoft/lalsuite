#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig qw(pi);
use POSIX;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

createDagAndConfigFiles.perl - create the dag and configuration files for a TwoSpect run

=head1 SYNOPSIS

createDagAndConfigFiles.perl [options]

=head1 OPTIONS

=item B<help:>
Show the help information

=item B<outputDir:>
Output directory path. Detector name will be appended [R]

=item B<logfile:>
Output path/filename of log file, detector name appended [R]

=item B<program:>
Path and file name of executable [R]

=item B<memRequest:>
Amount of memory to request of the nodes in MB (1550)

=item B<fstart:>
Frequency in Hz to start making config files [R]

=item B<fspan:>
Frequency span in Hz of each band being searched [R]

=item B<numberBands:>
Number of bands to produce >= 1 (1)

=item B<ifo:>
Interferometer to use, multiple (order specific) [R]

=item B<sourceSFTs:>
Path/file to concatenated SFT file, multiple (ordered) [R]

=item B<sourceSFTtypes:>
SFTs are 'vladimir' or 'standard', multiple (ordered) [R]

=item B<t0:>
Start time of the search, multiple (ordered) [R]

=item B<Tcoh:>
Coherence lengthe of the SFTs in seconds (1800)

=item B<Tobs:>
Observation time of the search in seconds (40551300)

=item B<Pmin:>
Minimum period of the search in seconds (7200)

=item B<Pmax:>
Maximum period of the search in seconds (8110260)

=item B<dfmin:>
Minimum modulation depth of the search in Hz (0.000277)

=item B<dfmax:>
Maximum modulation depth of the search in Hz (0.1)

=item B<ihsfar:>
IHS false alarm rate (1e-14)

=item B<ihsfomfar:>
IHS figure of merit false alarm rate (1.0)

=item B<tmplfar:>
Template false alarm rate (1e-18)

=item B<blksize:>
Block size of running median (101)

=item B<minTemplateL:>
Minimum template length in pixels (1)

=item B<maxTemplateL:>
Maximum template length in pixels (500)

=item B<avesqrtSh:>
Expected noise background of SFTs (1e-22)

=item B<SFToverlap:>
Overlap of the SFTs in seconds (900)

=item B<FFTplan:>
Plan flag for the SFTs: 0, 1, 2, or 3 (3)

=item B<maxSkyLocPerJob:>
Maximum sky locations per job (200)

=item B<linedetection:>
Line detection threshold setting

=item B<keepTopIHS:>
Keep only the top N IHS candidates

=item B<fastchisqinv:>
Use fast chi squared inversion methods

=item B<useSSE:>
Use SSE intrinsics

=item B<useAVX:>
Use AVX intrinsics, implies use SSE intrinsics

=item B<markBadSFTs:>
Mark the non-Gaussian SFTs and remove them

=item B<scoX1:>
Perform the restricted, Sco X-1 search

=item B<ihsfactor:>
IHS folding factor (5)

=cut

my $help = 0;
my $outputDir = '';
my $logfile = '';
my $executable = '';
my $memRequest = 1550;
my $fstart = 0.0;
my $fspan = 0.0;
my $numberBands = 1;
my @ifos;
my @sourceSFTs;
my @t0;
my $Tcoh = 1800;
my $SFToverlap = 900;
my $Tobs = 40551300;
my $Pmin = 7200;
my $Pmax = 8110260;
my $dfmin = 0.0;
my $dfmax = 0.1;
my $ihsfar = 1.0e-14;
my $ihsfomfar = 1.0;
my $tmplfar = 1.0e-18;
my $blksize = 101;
my $minTemplateLength = 1;
my $maxTemplateLength = 500;
my $avesqrtSh = 1.0e-22;
my $FFTplanFlag = 3;
my $maxSkyLocationsPerJob = 200;
my $linedetection = -1.0;
my $keepTopIHS = -1;
my $fastchisqinv = 0;
my $useSSE = 0;
my $useAVX = 0;
my $markBadSFTs = 0;
my $scoX1 = 0;
my $ihsfactor = 5;
#LHO t0 S6 = 931081500
#LLO t0 S6 = 931113900
#Virgo t0 VSR2/3 = 931131900

GetOptions('help' => \$help,
           'outputDir=s' => \$outputDir,
           'logfile=s' => \$logfile,
           'program=s' => \$executable,
           'memRequest:i' => $memRequest,
           'fstart=f' => \$fstart,
           'fspan=f' => \$fspan,
           'numberBands:i' => \$numberBands,
           'ifo=s' => \@ifos,
           'sourceSFTs=s' => \@sourceSFTs,
           't0=f' => \@t0,
           'Tcoh:f' => \$Tcoh,
           'Tobs:f' => \$Tobs,
           'Pmin:f' => \$Pmin,
           'Pmax:f' => \$Pmax,
           'dfmin:f' => \$dfmin,
           'dfmax:f' => \$dfmax,
           'ihsfar:f' => \$ihsfar,
           'ihsfomfar:f' => \$ihsfomfar,
           'tmplfar:f' => \$tmplfar,
           'blksize:i' => \$blksize,
           'minTemplateL:i' => \$minTemplateLength,
           'maxTemplateL:i' => \$maxTemplateLength,
           'avesqrtSh:f' => \$avesqrtSh,
           'SFToverlap:f' => \$SFToverlap,
           'FFTplan:i' => \$FFTplanFlag,
           'maxSkyLocPerJob:i' => \$maxSkyLocationsPerJob,
           'linedetection:f' => \$linedetection,
           'keepTopIHS:i' => \$keepTopIHS,
           'fastchisqinv' => $fastchisqinv,
           'useSSE' => \$useSSE,
           'useAVX' => \$useAVX,
           'markBadSFTs' => \$markBadSFTs,
           'scoX1' => \$scoX1,
           'ihsfactor:i' => \$ihsfactor) or pod2usage(2);
pod2usage(1) if $help;

my $numIFOs = @ifos;
my @ifokey;
my @outdirectory0;
for (my $ii=0; $ii<$numIFOs; $ii++) {
   if ($ifos[$ii] eq "LHO" || $ifos[$ii] eq "H1") {
      $ifos[$ii] = "LHO";
      push(@ifokey, "H1");
   } elsif ($ifos[$ii] eq "LLO" || $ifos[$ii] eq "L1") {
      $ifos[$ii] = "LLO";
      push(@ifokey, "L1");
   } elsif ($ifos[$ii] eq "Virgo" || $ifos[$ii] eq "V1") {
      $ifos[$ii] = "Virgo";
      push(@ifokey, "V1");
   }
   push(@outdirectory0, "$outputDir\_$ifokey[$ii]/output/");
   unless(-e "$outputDir\_$ifokey[$ii]" || mkdir "$outputDir\_$ifokey[$ii]") { die "mkdir $outputDir\_$ifokey[$ii] failed: $?"; }
   unless(-e "$outputDir\_$ifokey[$ii]/in" || mkdir "$outputDir\_$ifokey[$ii]/in") { die "mkdir $outputDir\_$ifokey[$ii]/in failed: $?"; }
   unless(-e "$outputDir\_$ifokey[$ii]/err" || mkdir "$outputDir\_$ifokey[$ii]/err") { die "mkdir $outputDir\_$ifokey[$ii]/err failed: $?"; }
   unless(-e "$outputDir\_$ifokey[$ii]/output" || mkdir "$outputDir\_$ifokey[$ii]/output") { die "mkdir $outputDir\_$ifokey[$ii]/output failed: $?"; }

   open(CONDORFILE,">$outputDir\_$ifokey[$ii]/condor") or die "Cannot write to $outputDir\_$ifokey[$ii]/condor $!";
   print CONDORFILE<<EOF;
universe=standard
executable=$executable
input=/dev/null
output=/dev/null
error=$outputDir\_$ifokey[$ii]/err/err.\$(PID)
arguments=--config=$outputDir\_$ifokey[$ii]/in/\$(PID)
log=$logfile\_$ifokey[$ii]
request_memory = $memRequest
notification=Never
queue
EOF
   close(CONDORFILE);

   my @skyfilesexisting = <$outputDir\_$ifokey[$ii]/in/skygrid.*>;
   my $numberofjobssofar = @skyfilesexisting;
   my $directorynumber = 0;
   if ($numberofjobssofar > 0) { $directorynumber = $numberofjobssofar; }

   open(DAG,">>$outputDir\_$ifokey[$ii]/dag") or die "Cannot write to $outputDir\_$ifokey[$ii]/dag $!";
   for (my $jj=0; $jj<$numberBands; $jj++) {
      my $fmin = sprintf("%.3f", $fstart+$jj*$fspan);

      die "You already are going to analyze $fmin Hz" if -e "$outputDir\_$ifokey[$ii]/skygrid\-${fmin}Hz\-${fspan}Hz.dat";

      if ($scoX1==0) {
         system("../skygridsetup --fmin=$fmin --fspan=$fspan --Tcoh=$Tcoh --IFO=$ifokey[$ii] --SFToverlap=$SFToverlap --t0=$t0[$ii] --Tobs=$Tobs --v2 --outfilename=$outputDir\_$ifokey[$ii]/skygrid\-${fmin}Hz\-${fspan}Hz.dat");
         die "skygridsetup failed: $?" if $?;
      } else {
         open(SKYPOINT, ">$outputDir\_$ifokey[$ii]/skygrid\-${fmin}Hz\-${fspan}Hz.dat") or die "Cannot write to $outputDir\_$ifokey[$ii]/skygrid\-${fmin}Hz\-${fspan}Hz.dat $!";
         print SKYPOINT "4.275699238500 -0.272973858335";
         close(SKYPOINT);

         my $scoX1P = 68023.70;
         my $scoX1asini = 1.44;
         $Pmin = $Tobs/(int($Tobs/$scoX1P+0.5)+1.0);
         $Pmax = $Tobs/(int($Tobs/$scoX1P+0.5)-1.0);
         $dfmin = 2*pi*$fmin*($scoX1asini-3.0*0.18)/($scoX1P+3.0*.0432);
         $dfmax = 2*pi*($fmin+$fspan)*($scoX1asini+3.0*0.18)/($scoX1P-3.0*.0432);
      }

      open(INPUT, "$outputDir\_$ifokey[$ii]/skygrid\-${fmin}Hz\-${fspan}Hz.dat") or die "Cannot open $outputDir\_$ifokey[$ii]/skygrid\-${fmin}Hz\-${fspan}Hz.dat $!";
      my $numberSkyLocations = 0;
      my $firstpointset = 0;
      my $startbin = int(($fstart+$jj*$fspan)*$Tcoh + 0.5);
      my $skyregionfile = '';
   
      while (my $line=<INPUT>) {
         chomp($line);
         if ($firstpointset==0 && !eof(INPUT)) {
            $skyregionfile = "$outputDir\_$ifokey[$ii]/in/skygrid.$directorynumber";

            open(SKYFILE,">$skyregionfile") or die "Cannot write to $skyregionfile $!";
            print SKYFILE "$line\n";
            $numberSkyLocations++;
            $firstpointset = 1;

         } elsif ($numberSkyLocations==$maxSkyLocationsPerJob-1 || eof(INPUT)) {
            if ($firstpointset==0) {
               $skyregionfile = "$outputDir\_$ifokey[$ii]/in/skygrid.$directorynumber";
               open(SKYFILE,">$skyregionfile") or die "Cannot write to $skyregionfile $!";
               $firstpointset = 1;
            }
            print SKYFILE "$line\n";
            $numberSkyLocations++;
            close(SKYFILE);

            my $randseedval = $startbin+$directorynumber;

            my $outdirectory = $outdirectory0[$ii].$directorynumber;
            open(OUT,">$outputDir\_$ifokey[$ii]/in/$directorynumber") or die "Cannot write to $outputDir\_$ifokey[$ii]/in/$directorynumber $!";
            print OUT<<EOF;
fmin $fmin
fspan $fspan
Tobs $Tobs
Tcoh $Tcoh
SFToverlap $SFToverlap
ihsfar $ihsfar
ihsfomfar $ihsfomfar
tmplfar $tmplfar
Pmin $Pmin
Pmax $Pmax
dfmin $dfmin
dfmax $dfmax
t0 $t0[$ii]
blksize $blksize
avesqrtSh $avesqrtSh
minTemplateLength $minTemplateLength
maxTemplateLength $maxTemplateLength
sftFile $sourceSFTs[$ii]
outdirectory $outdirectory
IFO $ifokey[$ii]
skyRegionFile $skyregionfile
FFTplanFlag $FFTplanFlag
randSeed $randseedval
ihsfactor $ihsfactor
EOF

            if ($linedetection>=0.0) { print OUT "lineDetection $linedetection\n"; }
            if ($keepTopIHS>=0) { print OUT "keepOnlyTopNumIHS $keepTopIHS\n"; }
            if ($markBadSFTs!=0) { print OUT "markBadSFTs\n"; }
            if ($useSSE!=0) { print OUT "useSSE\n"; }
            elsif ($useAVX!=0) { print OUT "useAVX\n"; }
            if ($fastchisqinv!=0) { print OUT "fastchisqinv\n"; }
            close(OUT);

            print DAG "JOB A$directorynumber $outputDir\_$ifokey[$ii]/condor\n";
            print DAG "VARS A$directorynumber PID=\"$directorynumber\"\n";

            $directorynumber++;

            $numberSkyLocations = 0;
            $firstpointset = 0;

         } else {
            print SKYFILE "$line\n";
            $numberSkyLocations++;
         }
      } #ends while loop over sky locations
      close(INPUT);
   } #ends for loop over bands
   close(DAG);
} #ends loop over IFOs

