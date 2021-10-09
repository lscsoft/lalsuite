1. Fscan Information:

i. The code is in lalapps here:

https://git.ligo.org/lscsoft/lalsuite/tree/master/lalapps/src/pulsar/Fscan

The latest production code is at LHO under,

/home/pulsar/public_html/fscan/H1/scripts/

The production code runs via cron jobs on ldas-grid as user pulsar from,

/home/pulsar/public_html/fscan/H1/daily/H1Fscan_coherence
/home/pulsar/public_html/fscan/H1/monthly/H1Fscan_coherence
/home/pulsar/public_html/fscan/H1/weekly/H1Fscan_coherence

at LHO, and from,

/home/pulsar/public_html/fscan/L1/daily/L1Fscan_coherence
/home/pulsar/public_html/fscan/L1/monthly/L1Fscan_coherence
/home/pulsar/public_html/fscan/L1/weekly/L1Fscan_coherence

at LLO.

A copy of the code is kept at Caltech under,

/home/gregory.mendell/scripts/fscan/newFscanCode2017/

ii. This Fscan README file is in git here:

https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalapps/src/pulsar/Fscan/README.txt

iii. See these talks,

https://dcc.ligo.org/LIGO-G1400026

https://dcc.ligo.org/LIGO-G2000677

to get an overview of what the Fscans do.

iv. Links to Fscan information is here:

https://wiki.ligo.org/DetChar/O3LinesCombsInvestigations

v. Links to older Fscan investigations are also here,

https://wiki.ligo.org/CW/InvestigationsOverview

vi. Older documentation for Fscans is on this page:

https://wiki.ligo.org/viewauth/CW/FscanNotesAndInvestigations

2. Example Fscans:

The Fscans appear on the daily summary pages, e.g., here:

https://ldas-jobs.ligo-wa.caltech.edu/~detchar/summary/day/20200212/detchar/fscan/

https://ldas-jobs.ligo-la.caltech.edu/~detchar/summary/day/20200212/detchar/fscan/

Click on green button, e.g., Coherence Daily, to get to the Fscan calendars.

Click on a date then a channel and then on: Fscan Plots or Fscans Plots
(Reverse Order) or Coherence Plots.

Note that there links to the data in the plots and additional
information below or to the side of each plot.

3. Steps for running the Fscan code:

Note that the examples below are for running at LHO from

$HOME/public_html/fscan/test

If running at another site or a different directory, you will need to
change URLs and paths accordingly.

i. Login ldas-grid at LHO:

$ ssh <albert.einstein@>ldas-grid.ligo-wa.caltech.edu

or

$ ssh <albert.einstein@>ssh.ligo.org

and follow the menu to ldas-grid at LHO.

ii. Make an fscan working directory:

Make sure you're in your home directory, e.g.,

$ cd /home/albert.einstein

and run,

$ mkdir -p public_html/fscan/test

iii. Make a directories under the test directory, which will be the top
level directory for the output:

$ cd $HOME/public_html/fscan/test

$ mkdir output

$ mkdir output/comparisonFscans

iv. For now, don't worry about locked times.

But, if you want to run on "locked" times, you need to query the segment
database to find the times.

For example, while logged into the cluster, run:

$ ligo-proxy-init albert.einstein

You will be prompted for your ligo.org password. For example:

$ ligo-proxy-init gregory.mendell
Your identity: gregory.mendell@LIGO.ORG
Enter pass phrase for this identity:
Creating proxy .................................... Done
Your proxy is valid until: Jan 27 09:19:48 2021 GMT

Then, to get the analysis ready times, e.g., run:

$ ligolw_segment_query_dqsegdb --segment-url=https://segments.ligo.org
--query-segments --include-segments H1:DMT-ANALYSIS_READY:1
--gps-start-time 1265500818 --gps-end-time 1265587218 | ligolw_print -t
segment:table -c start_time -c end_time -d " "

This should output this list:

1265500818 1265520839
1265521077 1265522807
1265527149 1265572718
1265575316 1265583488

You can put the list of segments into a segment file, e.g.,

$HOME/public_html/fscan/test/mysegs.txt

However, the Fscan code can also get these segments for you, as you will
see below.

v. Get the Fscan code:

$ cd $HOME/public_html/fscan/test

$ cp /home/pulsar/public_html/fscan/H1/scripts/*.tcl .

$ cp /home/pulsar/public_html/fscan/H1/scripts/*.py .

or see item 1 above for the location of the code.

vi. Make a resource file called myFscans.rsc, like the following,

https://ldas-jobs.ligo-wa.caltech.edu/~gmendell/fscan/test/myFscans.rsc

For comparision, see also this production .rsc file:

https://ldas-jobs.ligo-wa.caltech.edu/~pulsar/fscan/H1/daily/H1Fscan_coherence/autoFscanGeneratorHann_H1FscanCoherence.rsc.sav22April2021

You can edit the start and end times, e.g., change the lines like these,

set startTime 1265500818

set endTime 1265587218

to the start/end times you want.

You can use tconvert to get times.

For example, to get the above times I ran:

$ tconvert Feb 11 2020 16:00:00 PST
1265500818

$ tconvert Feb 12 2020 16:00:00 PST
1265587218

*** WARNING!!! ***

In the .rsc file, this setting,

set startTimeFile "lastTimeUsedByH1FscanCoherenceGeneratorAuto.txt";

gives a file with the last time the Fscan ran. If this file exists, the
time in this file, INSTEAD of the startTime will be used. Thus, if you
start the Fscans, jobs will be launched all the way back to the last
time the Fscan ran! If the Fscans have not run for a while, you will
want to update this file to the time from which you actually want to start.

Also note this line:

set parentOutputDirectory "output";

The Fscans will put directories based on date and channel name and the
fscans under this directory.

vii. You need to set env variables. Make env.txt with

export S6_SEGMENT_SERVER=https://segments.ligo.org
export LSC_DATAFIND_PATH=/usr/bin
export MAKESFTS_PATH=/usr/bin
export SPECAVG_PATH=/home/pulsar/searchcode/bin
export PLOTSPECAVGOUTPUT_PATH=/home/[YOUR USER NAME] /public_html/fscan/test
export CHECK_COHERENCE_PATH=/home/[YOUR USER NAME] /public_html/fscan/test
export COHERENCE_FROM_SFTS_PATH=/home/[YOUR USER
NAME]/public_html/fscan/test
export PATH=$PATH:/home/pulsar/searchcode/bin

but with [YOUR USER NAME] replaced with your albert.einstein user name.

*** ToDo: Check whether the Fscans work with lalapps_spec_avg and
MakeSFTDAG under /usr/bin rather than using the older versions under
/home/pulsar/searchcode/bin ***

viii. Before running the Fscans, run

$ source env.txt

and to get the correct version of Python, run,

$ . /home/detchar/opt/gwpysoft-2.7/bin/activate

*** ToDo: Update the Fscans to use a newer version of Python ***

(Note that the "." is the same as "source" in the syntax above. Either
should work, but "." has been used when activating gwpysoft-2.7.)

ix. Update the accountingGroupUser:

In

fscanDriver.py

you must change,

accountingGroupUser = "gregory.mendell"

to

accountingGroupUser = "[YOUR USER NAME]"

There is also a command line option for setting this with fscanDriver.py,

-U, --accounting-group-user  (optional) accounting group albert.einstein
username to be added to the condor submit files.

However, that option is not set by multiFscanGenerator.tcl.

*** ToDo: accountingGroupUser should be added as a configuration option
for multiFscanGenerator.tcl and then  accountingGroupUser could be set
in the .rsc file.***

x. Start the Fscans:

$ ./multiFscanGenerator.tcl myFscans.rsc -R

You will see lots of messages scroll by about the script running
MakeSFTDAG, fscanDriver.py, and condor_submit_dag, and hopefully the
word "Succeeded" a lot!

The jobs will run under condor.

Once the messages stop scrolling by, run "condor_q -nobatch" to monitor
progress.

xi. Make the web pages by running,

$ ./generateMultiFscanHTML.tcl output

Your results will be here:

https://ldas-jobs.ligo-wa.caltech.edu/~[YOUR USER
NAME]/fscans/test/fscanNavigation.html

It will take 1-2 hrs for the jobs to run (or longer if the cluster is busy).

On the filesystem, there should be a list of channels and
fscanChannels.html under under, e.g., something like,

output/fscans_2020_02_11_16_00_00_PST_Tue/H1_GDS-CALIB_STRAIN

To understand the details of how the jobs are set up, go to the working
directory, e.g.,

$ cd output/fscans_2020_02_12_16_00_02_PST_Wed/H1_GDS-CALIB_STRAIN/

and study the .sub and .dag files.

But don't worry if you don't see any output for a while.

However, if things fail:

Look in

output/fscans_2020_02_11_16_00_00_PST_Tue/H1_GDS-CALIB_STRAIN/sft/tmp

for .sft files.

Check for errors in files like,

output/fscans_2020_02_11_16_00_00_PST_Tue/H1_GDS-CALIB_STRAIN/SUPERDAGH1_GDS-CALIB_STRAIN_H1_1265414418_1265500818tmp.dag.dagman.out

and the logs directory, e.g, under,

output/fscans_2020_02_11_16_00_00_PST_Tue/H1_GDS-CALIB_STRAIN/logs

xii. To retry the test, you will need to remove a few things first:

$ rm lastTimeUsedByH1FscanCoherenceGeneratorAuto.txt
$ cd output
$ rm -rf fscans_2020_02_12_16_00_02_PST_Wed

xiii. If you already have mysegs.txt, you can also run on the segments
in mysegs.txt by changing in myFscans.rsc, this line from,

set useLSCsegFind 1

to

set useLSCsegFind 0

and running,

$ ./multiFscanGenerator.tcl myFscans.rsc mysegs.txt -R

However, in most cases during a run, you'll want to keep these lines in
in the .rsc file:

set useLSCsegFind 1

set typeLSCsegFind "H1:DMT-ANALYSIS_READY:1",

or to run on locked times instead of analysis read times, set in the
.rsc file,

set typeLSCsegFind "H1:DMT-DC_READOUT_LOCKED:1"

Also, whenever you generate a new set of Fscans, update the html to the
results by running:

$ ./generateMultiFscanHTML.tcl output

xiv. To see how to automate fscan running, at LHO see:

/home/pulsar/public_html/fscan/H1/daily/H1Fscan_coherence/runDailyH1FscanCoherence.sh

/home/pulsar/public_html/fscan/H1/daily/H1Fscan_coherence/genDailyH1FscanCoherenceHTML.sh

and

/home/pulsar/.bash_profile

Also, in myFscans.rsc, change,

set startTime "1265414418";
set endTime "1265500818";
set timeLag "0";
#set endTime "now";
#set timeLag "7800";
set startTimeFile "lastTimeUsedByH1FscanCoherenceGeneratorAuto.txt";

to,

set startTime "1265414418";
#set endTime "1265500818";
#set timeLag "0";
set endTime "now";
set timeLag "7800";
set startTimeFile "lastTimeUsedByH1FscanCoherenceGeneratorAuto.txt";

to automatically pick up the endTime as "now" with a lag of 7800
seconds. Note that the startTime will always get set by
lastTimeUsedByH1FscanCoherenceGeneratorAuto.txt, if that file exists.
This file is updated each time multiFscanGenerator.tcl is run. (See the
WARNING about updating this file if the Fscans haven't run for awhile
above.)

Finally, see the cron jobs run by user pulsar,

10 19 * * *
/home/pulsar/public_html/fscan/H1/daily/H1Fscan_coherence/runDailyH1FscanCoherence.sh
20 21 * * *
/home/pulsar/public_html/fscan/H1/daily/H1Fscan_coherence/genDailyH1FscanCoherenceHTML.sh

