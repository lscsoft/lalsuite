December 19, 2017: Add updated code for using Python for plotting and comb finding.

Before running, set this location to get working Python and Matplotlib:

. /home/detchar/opt/gwpysoft-2.7/bin/activate

Here are some notes on the Fscans.

1. The main documentation for Fscans is on this page:

https://wiki.ligo.org/viewauth/CW/FscanNotesAndInvestigations

This is from S6 but gives an overview and how the daily Fscans were run
then.

2. More links to Fscan pages are here,

https://wiki.ligo.org/CW/InvestigationsOverview

under "Detector characterisation".

3. These talks give information about the Fscans:

Fscan Update January 2014: https://dcc.ligo.org/LIGO-G1400026

Fscan Update March 2012: https://dcc.ligo.org/LIGO-G1200176

Fscan Update March 2010: https://dcc.ligo.org/LIGO-G1000242

Fscan Update December 2009: https://dcc.ligo.org/LIGO-G0901064

Fscan Discussion: https://dcc.ligo.org/LIGO-G0900102

4. The code is in lalapps here:

https://ligo-vcs.phys.uwm.edu/cgit/lalsuite/tree/lalapps/src/pulsar/fscan

5. Example Fscans can be found here:

i. Daily Fscans of the primary channels:

https://ldas-jobs.ligo-wa.caltech.edu/~pulsar/fscan/H1_DUAL_ARM/H1_DUAL_ARM_HANN/fscanNavigation.html

https://ldas-jobs.ligo-la.caltech.edu/~pulsar/fscan/L1_DUAL_ARM/L1_DUAL_ARM_DCREADOUT_HANN/fscanNavigation.html

Scroll to the bottom of the calendar to get to ER8/O1 times. Click on a
date then a channel, then on "old or new" to get the plots.

(Old for the old arrangement from low frequency to high frequency, new
for the new arrangement from high frequency to low. Probably need to
rename the buttons. Both arrangements are useful.)

The Fscans also appear on the daily summary pages, e.g., here:

https://ldas-jobs.ligo-wa.caltech.edu/~detchar/summary/day/20151012/detchar/fscan/

https://ldas-jobs.ligo-la.caltech.edu/~detchar/summary/day/20151012/detchar/fscan/

6. Steps for running the Fscan code:

i. Login ldas-grid at LHO:

$ gsissh ldas-grid.ligo-wa.caltech.edu

ii. Change to the fscan working directory:

For example:

$ cd /home/[YOUR USER NAME]/public_html/test

iii. Make a directory under the test directory, which will be the top
level directory for the output:

$ mkdir myOutput

iii. If you want to run on "locked" times, you need to query the segment
database to find the times.

For example, to find the times for the past 24 hrs up to 6 am PDT today,
you would need to run while logged into the cluster:

$ ligolw_segment_query_dqsegdb --segment-url=https://segments.ligo.org
--query-segments --include-segments H1:DMT-ANALYSIS_READY:1
--gps-start-time 1128776417 --gps-end-time 1128862817 | ligolw_print -t
segment:table -c start_time -c end_time -d " "

This will return a list like this:

1128776417 1128782611
1128811907 1128837796
1128858099 1128862817

You can put these segments into a segment file, e.g.,

$HOME/public_html/test/mysegs.txt

iv. Edit the .rsc (resource) file as needed:

An example .rsc files is here:

https://ldas-jobs.ligo-wa.caltech.edu/~pulsar/fscan/H1_DUAL_ARM/H1_DUAL_ARM_HANN/autoFscanGeneratorHann_H1_DUAL_ARM_HANN.rsc

If you make a file,

myFscans.rsc

which is a copy of the above file then

$ vi myFscans.rsc

The edit the channel list at the top of the .rsc file.

Edit the start and end times, e.g., to be this:

set startTime 1128776417
set endTime 1128862817

You will also need this line:

set parentOutputDirectory "myOutput";

The Fscans will put directories based on date and channel name and the
fscans under this directory.

For example, the daily Fscans I pointed to above, for H1, that ran
today, ran from:

/home/pulsar/public_html/fscan/H1_DUAL_ARM/H1_DUAL_ARM_HANN

And from that location the Fscan output for today's STRAIN data will go
under:

H1_DUAL_ARM_HANN/fscans_2015_10_14_06_00_02_PDT_Wed/H1_GDS-CALIB_STRAIN/

v. You will need these env variables set:

$ export LSC_DATAFIND_PATH=/usr/bin
$ export MAKESFTS_PATH=/archive/home/pulsar/searchcode/bin
$ export SPECAVG_PATH=/archive/home/pulsar/searchcode/bin
$ export
PLOTSPECAVGOUTPUT_PATH=/archive/home/pulsar/searchcode/src/lalsuite/lalapps/src/pulsar/fscan

Also update your PATH:

$ export PATH=$PATH:/home/pulsar/searchcode/bin

This is no longer needed: add the location of ligotools to your env:

$ eval `/ligotools/bin/use_ligotools`

However, once you have ligotools in your env, you can run things like FrChannels,
FrDump, and also tconvert to convert between GPS and standard time.

For example:

$ tconvert Oct 13 2015 06:00:00 PDT
1128776417

$ tconvert Oct 14 2015 06:00:00 PDT
1128862817

vi. Start the Fscans:

$ cd $HOME/public_html/test

$ mkdir myOutput

$ ./multiFscanGenerator.tcl myFscans.rsc $HOME/public_html/test/myseg.txt -R

You will see lots of messages scroll by about the script running
MakeSFTDAG, fscanDriver.py, and condor_submit_dag, and hopefully the
word "Succeeded" a lot!

The jobs will run under condor.

Once the messages stop scrolling by, run "condor_q pulsar" to monitor
progress.

vii. Make the web pages by running,

$ ./generateMultiFscanHTML.tcl myOutput.

Your results will be here:

https://ldas-jobs.ligo-wa.caltech.edu/~[YOUR USER
NAME]/test/fscanNavigation.html

There are a lot steps done when the scripts are run that go on under the
hood. 

See the following older instructions for running the fscanDriver.py script, which gets
run above by the multiFscanGenerator.tcl script. 

Older instructions:

1. Run ./fscanDriver.py -h to get help

2. Example Trial Run That Generates A Dag To Make H1 DARM_CTRL SFTs in /archive/home/pulsar/fscan/sfts:

../fscanDriver.py -s 842705233 -L 36000 -G exampleTest -d RDS_R_L1 -k 40 -T 1800 -p /archive/home/pulsar/searchcode/src/lalapps/src/pulsar/fscan/test/sfts -N H1:LSC-DARM_CTRL -F 100 -B 20 -b 5 -X fscanTest -o /usr1/pulsar -C

3. Example That Generates And Runs A Dag To Make H1 DARM_CTRL SFTs in /archive/home/pulsar/fscan/sfts:

../fscanDriver.py -s 842705233 -L 36000 -G exampleTest -d RDS_R_L1 -k 40 -T 1800 -p /archive/home/pulsar/searchcode/src/lalapps/src/pulsar/fscan/test/sfts -N H1:LSC-DARM_CTRL -F 100 -B 20 -b 5 -X fscanTest -o /usr1/pulsar -C --run 

4. Example That Generates And Runs A Dag To Make H1 DARM_CTRL SFTs in /archive/home/pulsar/fscan/sfts and then run matlab driver script to output plots:

../fscanDriver.py -s 842705233 -L 36000 -G exampleTest -d RDS_R_L1 -k 40 -T 1800 -p /archive/home/pulsar/searchcode/src/lalapps/src/pulsar/fscan/test/sfts -N H1:LSC-DARM_CTRL -F 100 -B 20 -b 5 -X fscanTest -o /usr1/pulsar -O . -C --run 

5. Obsolete instructions for using Matlab (replaced by plotSpecAvgOutput.py). 

How to compile plotSpecAvgOutput.m, e.g., as user pulsar on ldas-grid:

Make sure /archive/home/pulsar/.usematlab_r2008a exists, otherwise touch this file, and logout and login so that "which matlab" returns:
/ldcg/matlab_r2008a/bin/matlab 

Then run these commands:

$ cd /archive/home/pulsar/searchcode/src/lalsuite/lalapps/src/pulsar/fscan
$ source MatlabSetup_R2008a_glnxa64.sh
$ matlab -nodisplay -nodesktop -nojvm
>> mbuild -setup
>> mcc -mv plotSpecAvgOutput.m
>> exit

You many need to run the code once by hand before it will work on the cluster (it does not matter if this fails, this just makes Matlab set up the libraries under ~/.mcr_cache_v78): 

./plotSpecAvgOutput S4/spec_50.00_100.00_H1_793181266_795677345 /archive/home/pulsar/public_html/fscan/test/spec_50.00_100.00_H1_793181266_795677345.pdf H1:hoft 793181266 795677345 50 100 10 5 10

After this, logout, and login to unset the environment set by MatlabSetup_R2008a_glnxa64.sh.

You do not want to source MatlabSetup_R2008a_glnxa64.sh when running plotSpecAvgOutput!!!

Instead, the fscan code runs run_plotSpecAvgOutput.sh /ldcg/matlab_r2008a, which is a wrapper script the matlab mcc commands generated that runs plotSpecAvgOutput, and it will take care of setting up the environment.   

