#!/usr/bin/perl -w

$fstart = 50.0;
$totalband = "50-60Hz-0.05HzBands";
$ifo = "H1";
$fstep = 0.25;
@fmins = ();
for($ii=0; $ii<40; $ii++) {
   push(@fmins, sprintf("%.3f", $fstart+$ii*$fstep));
}

$fspan = "0.25";
$tcoh = 1800;
$ephemdir = "/Users/evgoet";
$ephemyear = "08-11-DE405";
$sky = "allSky";
$sftoverlap = 900;

for($ii=0; $ii<40; $ii++) {
   $outfile = "/Users/evgoet/S6/band".$totalband."/skygrid-".$fmins[$ii]."-".$fspan."HzBand.dat";
   system("./skygridsetup --fmin=".$fmins[$ii]." --fspan=".$fspan." --Tcoh=".$tcoh." --IFO=".$ifo." --ephemDir=".$ephemdir." --ephemYear=".$ephemyear." --skyRegion=".$sky." --SFToverlap=".$sftoverlap." --outfilename=".$outfile);
}


