# Copyright (C) 2013 Lindy Blackburn, Reed Essick and Ruslan Vaulin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import os
import sys
from laldetchar.idq import idq
#from laldetchar.idq import ovl ### already imported as part of idq module, reference from there
                                ### this lets us use 'ovl' as a boolean flag
import math
import time
from laldetchar.idq import event
import logging
import tempfile
import optparse
import subprocess
import ConfigParser

from laldetchar import git_version

#===================================================================================================

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
    """ This program performs real-time identification and classification of glitch events. Currently it operates on equal length strides (e.g. 32 seconds). Within the iDQ pipeline this script plays the central role by scheduling all tasks including periodic training and generation of summary html pages."""

#===================================================================================================

def compare_segments(
    gstlal_good,
    gstlal_covered,
    dmt_good,
    dmt_covered,
    logger,
    ):
    """
....Temporary function written to compare gstlal and dmt segments.
...."""

    # first determine coverage

    if gstlal_covered:
        gstlal_coverage = 1
    else:
        gstlal_coverage = 0
    if dmt_covered:
        dmt_coverage = 1
    else:
        dmt_coverage = 0

    # compare coverage

    if gstlal_covered != dmt_covered:
        logger.info('WARNING: gstlal covered and dmt covered segments are different'
                    )
        logger.info('gstlal covered segments: ' + str(gstlal_covered))
        logger.info('dmt covered segments: ' + str(gmt_covered))

    if gstlal_coverage == 1 and dmt_coverage == 1:
        logger.info('both gstlal and dmt data coverage is OK')
    elif gstlal_coverage == 0 and dmt_coverage == 1:
        logger.info('No gstlal coverage, dmt coverage is OK')
    elif gstlal_coverage == 1 and dmt_coverage == 0:
        logger.info('gstlal coverage is OK, no dmt coverage')
    else:
        logger.info('no gstlal or dmt coverage')

    # compare science segments

    if gstlal_coverage == 1 and dmt_coverage == 1:
        if gstlal_good != dmt_good:
            logger.info('WARNING: gstlal science and dmt science segments are different'
                        )
            logger.info('gstlal science segments: ' + str(gstlal_good))
            logger.info('dmt science segments: ' + str(gmt_good))

    if gstlal_coverage == 1:

        # use gstlal segments

        logger.info('Using gstlal science segments')
        return (gstlal_good, gstlal_covered)
    else:

        # use dmt segments

        logger.info('Using dmt science segments')
        return (dmt_good, dmt_covered)

#===================================================================================================

parser = optparse.OptionParser(version='Name: %%prog\n%s'
                               % git_version.verbose_msg,
                               usage='%prog [options]',
                               description=description)
parser.add_option(
    '-c',
    '--config-file',
    dest='config_file',
    help='configuration file',
    metavar='FILE',
    default='idq.ini',
    )
parser.add_option(
    '-s',
    '--start-time',
    dest='startgps',
    type='int',
    help='first stride begins at or after GPS time',
    metavar='GPS',
    default=idq.nowgps(),
    )
parser.add_option(
    '-e',
    '--end-time',
    dest='endgps',
    type='int',
    help='last stride ends at or before GPS time',
    metavar='GPS',
    default=1e10,
    )
parser.add_option(
    '-k',
    '--lock-file',
    dest='lockfile',
    help='use custom lockfile',
    metavar='FILE',
    default='.idq_realtime.lock',
    )
parser.add_option('-l', '--log-file', default='idq_realtime.log',
                  type='string', help='log file')

parser.add_option('', '--ignore-science-segments', 
    default=False, action="store_true", 
    help='analyze strides regardless of the science segement content. \
    This is NOT passed on to training job, which always uses science segments'
    )
parser.add_option('', '--ignore-science-segments-training', 
    default=False, action="store_true", 
    help="forces training jobs to ignore science segments as well. \
    This should NOT be used to generate actual DQ information, but may be useful when debugging the pipeline."
    )

parser.add_option("", "--no-robot-cert",
    default=False, action="store_true",
    help="do not use robot cert with segment query in training jobs"
    )

(options, args) = parser.parse_args()

#===================================================================================================
### setup logger to record processes
logger = logging.getLogger('idq_logger')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(message)s')
hdlr1 = logging.StreamHandler(sys.stdout)
hdlr1.setFormatter(formatter)
hdlr1.setLevel(logging.INFO)
logger.addHandler(hdlr1)
hdlr2 = logging.FileHandler(options.log_file)
hdlr2.setFormatter(formatter)
hdlr2.setLevel(logging.INFO)
logger.addHandler(hdlr2)

### redirect stdout and stderr into logger
sys.stdout = idq.LogFile(logger)
sys.err = idq.LogFile(logger)

#=================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(options.config_file)

### prevent multiple copies from running
idq.dieiflocked(config.get('general', 'idqdir') + '/'
                + options.lockfile)

#=================================================
### pull out general stuff
usertag = config.get('general', 'usertag')

### get the main directory where idq pipeline is going to be running.
mainidqdir = config.get('general', 'idqdir')

### define the gravitational wave channel/what is a glitch
gwchannel = config.get('general', 'gwchannel')
gwthreshold = config.getfloat('general', 'gw_kwsignif_thr')

### figure out which classifiers are present
classifiers = config.get('general', 'classifiers').split(' ')
ovl = "ovl" in classifiers
mvsc = "mvsc" in classifiers
svm = "svm" in classifiers
mla = mvsc or svm  # True if machine learning algorithms are present

### pick up which ifo this is
ifo = config.get('general', 'ifo')

### load settings for accessing gstlal segments
# gstlal_segments_location = config.get("lldq_science_segments", "xmlurl")
# lldq_name = ifo + ":" + config.get("lldq_science_segments", "lldq_name")
# hoft_name = ifo + ":" + config.get("lldq_science_segments", "hoft_name")

### load settings for accessing dmt segment files
dmt_segments_location = config.get('get_science_segments', 'xmlurl')
dmtdq_name = config.get('get_science_segments', 'include').split(':')[1]

### kleineWelle config
kwconfig = idq.loadkwconfig(config.get('general', 'kwconfig'))
kwbasename = kwconfig['basename']

#=================================================
### pull out things specific to realtime job
myconf = dict(config.items('idq_realtime'))

### sampling frequency for time-series files
ts_fs = float(myconf["sampling_rate"]) 

### output directory for realtime predictions
realtimedir = config.get('general', 'realtimedir')
if not os.path.exists(realtimedir):
    os.makedirs(realtimedir)

### slave the realtime job to the kw stride
stride = int(float(kwconfig['stride'])) ### this is given as a decimal, so we must cast as float first
delay = config.getint('idq_realtime', 'delay') ### buffer to let jobs finish

### headers for dat files
samples_header = config.get('idq_realtime', 'dat_columns').split()

#=================================================
### required information about summary jobs
summarydir = config.get('general', 'summarydir')

summary_stride = config.getint('idq_summary', 'stride')
idq_summary_script = config.get('condor', 'idq_summary')

#=================================================
### required information about training jobs

training_stride = config.getint('idq_train', 'stride')
idq_train_script = config.get('condor', 'idq_train')

# set up lookback time for initial and regular training jobs
training_lookback = config.get('idq_train', 'lookback')
if training_lookback == 'infinity':
    lookback = -1
else:
    lookback = int(training_lookback)

initial_training_lookback = config.get("idq_train","initial-lookback") ### if used, this is passed on directly and we don't need to manipulate it

#=================================================
###figure out if we need to do initial training or initial summary jobs
### training jobs allow us to compute ranks from sets of auxiliary glitches
### summary jobs allow us to map those ranks into FAP measurements

initial_training = False
initial_summary = False

if ovl:
    ### ovl training cache
    ovl_train_cache = config.get('general', 'ovl_train_cache')
    if not os.path.exists(ovl_train_cache):
        initial_training = True
    else: ### make sure these files are not empty
        file_obj = open(ovl_train_cache, 'r')
        initial_training += len(file_obj.readlines()) < 1
        file_obj.close()

    ### ovl summary cache
    ovl_uroc_cache = summarydir + '/ovl_uroc.cache'
    if not os.path.exists(ovl_uroc_cache):
        initial_summary = True
    else:
        file_obj = open(ovl_uroc_cache, 'r')
        initial_summary += len(file_obj.readlines()) < 1
        file_obj.close()

    ### ovl cache times
    ### used to determine if there is new information present
    ovl_train_cache_modification_time = 0
    ovl_uroc_cache_modification_time = 0

if mvsc:
    ### mvsc training cache
    mvsc_train_cache = config.get('general', 'mvsc_train_cache')
    if not os.path.exists(mvsc_train_cache):
        initial_training = True
    else:
        file_obj = open(mvsc_train_cache, 'r')
        initial_training += len(file_obj.readlines()) < 1
        file_obj.close()

    ### mvsc summary cache
    mvsc_uroc_cache = summarydir + '/mvsc_uroc.cache'
    if not os.path.exists(mvsc_uroc_cache):
        initial_summary = True
    else:
        file_obj = open(mvsc_uroc_cache, 'r')
        initial_summary += len(file_obj.readlines()) < 1
        file_obj.close()

    ### mvsc cache times
    mvsc_train_cache_modification_time = 0
    mvsc_uroc_cache_modification_time = 0

if svm:
    ### svm training cache
    svm_train_cache = config.get('general', 'svm_train_cache')
    if not os.path.exists(svm_train_cache):
        initial_training = True
    else:
        file_obj = open(svm_train_cache, 'r')
        initial_training += len(file_obj.readlines()) < 1
        file_obj.close()

    ### svm summary cache
    svm_uroc_cache = summarydir + '/svm_uroc.cache'
    if not os.path.exists(svm_uroc_cache):
        initial_summary = True
    else:
        file_obj = open(svm_uroc_cache, 'r')
        initial_summary += len(file_obj.readlines()) < 1
        file_obj.close()

    ### svm cache times
    svm_train_cache_modification_time = 0
    svm_uroc_cache_modification_time = 0

#===================================================================================================
### figure out where we are in time, and launch initial jobs in needed

t = int(math.ceil(options.startgps / stride)) * stride  # require times to be integer multiples of stride
global_tstart = t

#=================================================
# Initial training 
#=================================================
if initial_training:
    logger.info('Initial training of classifiers')
    logger.info('Begin: initial training')

    ### set up options for training script
    training_start_time = t - training_stride
    training_command = [
        idq_train_script,
        '--config',
        options.config_file,
        '--gps-start',
        str(training_start_time),
        '--gps-stop',
        str(t),
        '--lookback',
        str(initial_training_lookback),
	'--force'
        ]

    if options.ignore_science_segments_training:
        training_command.append( "--ignore-science-segments" )

    elif options.no_robot_cert:
        training_command.append( "--no-robot-cert" )

    logger.info('Submiting idq_train script with the following options:')
    logger.info(training_command)

    ### launch training script, wait for it to finish
    it_exit_code = idq.submit_command(training_command,
            process_name='initial_idq_training', dir=os.getcwd(),
            verbose=True)
    if it_exit_code != 0:
        logger.info('ERROR: idq_train script for initial training failed.')
        raise StandardError, "initial training job failed"

    logger.info('Done: initial training')

else:
    logger.info('All required training cache files already exist, skipping initial training.')

#=================================================
# Initial summary
#=================================================
if initial_summary:
    logger.info('Initial summary job for classifiers')
    logger.info('Begin: initial summary')
    logger.info('WARNING: the *.uroc files generated by this job do not have much meaning')

    summary_start_time = str(t / summary_stride * summary_stride - summary_stride)
    summary_stop_time = str(t / summary_stride * summary_stride)

    summary_command = [
        idq_summary_script,
        '--config',
        options.config_file,
        '--gps-start',
        summary_start_time,
        '--gps-stop',
        summary_stop_time,
        ]

    logger.info('Submiting idq_summary script with the following options:')
    logger.info(summary_command)

    summary_process = subprocess.Popen(summary_command,
            stdout=open('idq_summary.out', 'a'),
            stderr=open('idq_summary.err', 'a'))
    if summary_process.wait() != 0:
        logger.info('ERROR: idq_summary script for intitial summary failed' )
        raise StandardError, "initial summary job failed"

    logger.info('Done: initial summary')

else:
    logger.info('All required *uroc.cache files already exist, skipping initial summary.' )

#===================================================================================================
#
# MAIN realtime jobs
#
#===================================================================================================

logger.info('Begin real-time evaluation')

while t + stride < options.endgps:

    #===============================================================================================
    # manage training and summary cadences
    #===============================================================================================

    #=============================================
    # Check if next training cycle is reached 
    # if it is, launch training
    # 
    # training cycle launched if the remainder from divison of the time of the current evaluation step
    # by the training stride is positive and less than a single evaluation stride.
    # To avoid premature training launch, we also require that the difference between
    # the time of the current evaluation step and the start time exceeds a single training stride.

    t_remainder = t - (t / training_stride) * training_stride
    if (t - options.startgps > training_stride) and (t_remainder > 0) and (t_remainder <= stride):
        logger.info('****************************************************')

        training_start_time = str((t / training_stride) * training_stride - training_stride)
        training_stop_time = str((t / training_stride) * training_stride)

        if training_lookback == 'infinity': ### increment so we always look back to the first stride
            lookback += 1

        logger.info('Begin: launching training script for period: %s - %s'%(training_start_time,training_stop_time))

        training_command = [
            idq_train_script,
            '--config',
            options.config_file,
            '--gps-start',
            training_start_time,
            '--gps-stop',
            training_stop_time,
            '--lookback',
            str(lookback)
            ]

        if options.ignore_science_segments_training:
            training_command.append( "--ignore-science-segments" )

        elif options.no_robot_cert:
            training_command.append( "--no-robot-cert" )

        logger.info('Submiting idq_train script with the following options:')
        logger.info(training_command)

        ### launch idq_train
        # check if we are in a catch-up mode

        if abs(int(idq.nowgps()) - t) > 3600: ### 3600 is a magic number. 
                                              ### If training is launched more than 1 hour after it should have, we're in "catch-up" mode

            ### catching-up, wait for training job to complete before proceeding with evaluation
            logger.info('Running in catch-up mode, will wait until training is complete before proceeding with evaluation')
            subprocess.call(training_command,
                            stdout=open('idq_train.out', 'a'),
                            stderr=open('idq_train.err', 'a')) ### call blocks until the process to finishes
        else:
            ### running in realtime, don't wait for training job to complete so we can keep the latency low

            train_pid = subprocess.Popen(training_command,
                    stdout=open('idq_train.out', 'a'),
                    stderr=open('idq_train.err', 'a')).pid ### only remember the pid, and let the process float

        logger.info('Done: launching training script for period: %s - %s'%(training_start_time,training_stop_time))
        logger.info('****************************************************')

    #=============================================
    # Check if next summary cycle is reached 
    # if it is, launch idq_summary
    #
    # summary cycle launched if the remainder from divison of the time of the current evaluation step
    # by the summary stride is positive and less than a single evaluation stride.
    # To avoid premature summary launch, we also require that the difference between
    # the time of the current evaluation step and the start time exceeds a single summary stride.

    t_remainder = t - (t / summary_stride) * summary_stride
    if (t - options.startgps > summary_stride) and (t_remainder > 0) and (t_remainder <= stride):
        logger.info('****************************************************')

        summary_start_time = str((t / summary_stride) * summary_stride - summary_stride)
        summary_stop_time = str((t / summary_stride) * summary_stride)

        logger.info('Begin: launching summary script for period: %s - %s'%(summary_start_time,summary_stop_time))

        summary_command = [
            idq_summary_script,
            '--config',
            options.config_file,
            '--gps-start',
            summary_start_time,
            '--gps-stop',
            summary_stop_time,
            ]

        logger.info('Submiting idq_summary script with the following options:')
        logger.info(summary_command)

        summary_pid = subprocess.Popen(summary_command,
                stdout=open('idq_summary.out', 'a'),
                stderr=open('idq_summary.err', 'a')).pid ### only remember the pid, let the process float

        logger.info('Done: launching summary script for period: %s - %s'%(summary_start_time,summary_stop_time))
        logger.info('****************************************************')

    #===============================================================================================
    # realtime evaluation
    #===============================================================================================
    logger.info('----------------------------------------------------')

    logger.info('Begin: stride %d-%d' % (t, t + stride))

    ### make sure end time of analysis stride is not ahead of now-delay
    wait = t + stride + delay - int(idq.nowgps())
    if wait > 0:
        logger.info('  waiting %.1f seconds before processing' % wait)
        time.sleep(wait)

    #=============================================
    # get science segments (still a mess)
    #=============================================
    if options.ignore_science_segments:
        logger.info('analyzing data regardless of science segements')
        
    else:
        logger.info('Begin: querrying science segments')

        ### get gstlal segments
        ### not yet ready
#       try: ### get gstlal segments
#           gstlalxmlfiles = idq.get_all_files_in_range(gstlal_segments_location, t, t+stride, pad=0, suffix='.xml.gz')
#           (gstlal_good, gstlal_covered) = idq.extract_lldq_segments(gstlalxmlfiles, lldq_name, hoft_name)
#       except Exception:
#           e = sys.exc_info()[0]
#           logger.warning("error from gstlal segment query: %s" % e)
#           (gstlal_good, gstlal_covered) = ([], [])

        ### get DMT xml segments (changing the query)
        try:
            dmtfiles = idq.get_all_files_in_range(dmt_segments_location, t, t + stride, pad=0, suffix='.xml')
            (good, covered) = idq.extract_dmt_segments(dmtfiles, dmtdq_name)
        except Exception:
            e = sys.exc_info()
            logger.warning('error from DMT segment query: %s, %s, %s' % (e[0], e[1], e[2]))
            (good, covered) = ([], [])

        logger.info('Done: querrying science segments')

        ### handle incomplete science coverage (incomplete DMT DQ info)
        if event.livetime(covered) < stride:
            if wait > -delay:
                logger.info('  unknown science coverage, waiting additional %d seconds' % delay)
                time.sleep(delay)

                ### decrease wait time to account for this delay
                wait -= delay

                ### retry segment query
                logger.info('  retrying science segment query')
                try:
                    dmtfiles = idq.get_all_files_in_range(dmt_segments_location, t, t + stride, pad=0, suffix='.xml')
                    (good, covered) = idq.extract_dmt_segments(dmtfiles, dmtdq_name)
                except Exception:
                    e = sys.exc_info()
                    logger.warning('error from DMT segment query: %s, %s, %s' % (e[0], e[1], e[2]))
                    (good, covered) = ([], [])

                ### is science coverage still unknown?
                if event.livetime(covered) < stride:
                    logger.warning('unknown science coverage, skipping')
                    t += stride
                    continue
            else:
                logger.warning('unknown science coverage, skipping')
                t += stride
                continue

        ### skip if not in science mode for entire stride
        if event.livetime(good) < stride:
            logger.info('incomplete science coverage, skipping')
            t += stride
            continue

        logger.info('complete science coverage')

    #=============================================
    # find and read KW.trg files
    #=============================================
    logger.info('Begin: aggregating triggers')

    ### check if KW file is there
    ### THIS IS HIGHLY DEPENDENT ON KW NAMING CONVENTION
    ### WE MAY NEED TO RE-THINK THIS LOGIC IF WE USE ANOTHER ETG
    kwfilename = '%s/%s/%s-%d/%s-%d-%d.trg' % (
        config.get('general', 'gdsdir'),
        kwbasename,
        kwbasename,
        t / 1e5,
        kwbasename,
        t,
        stride,
        )

    if not os.path.exists(kwfilename):
        if wait > -stride: ### kw files may not have appeared yet
            logger.info('  missing KW triggers, waiting additional %d seconds' % stride)
            time.sleep(stride)

            if not os.path.exists(kwfilename): ### kw file still doesn't exit
                logger.warning('  still missing KW triggers, skipping')
                t += stride
                continue

        else: ### kw files should have appeared
            logger.warning('  missing KW triggers, skipping')
            t += stride
            continue

    logger.info('  loading KW triggers')
    kwdict = event.loadkwm(kwfilename)

    #=============================================
    # set up samples for classifiers
    #=============================================
    ### require GW channel to be in trigger dict
    ### if the GW channel did not contain any glitches, it will not be in the kwdict
    if gwchannel not in kwdict:
        kwdict[gwchannel] = []

    ### pull out only those GW glitches with sufficiently high significance
    gw_trig = [a for a in kwdict[gwchannel] if a[-1] >= gwthreshold]

    ### retrieve GW trigger central times
    gw_gps = [a[2] for a in gw_trig]

    ### generate clean times as a poisson time series
    clean_gps = event.randomrate(float(myconf['clean_rate']), [[t, t + stride]])

    ### define segemnts that are too close to glitch times to be really "clean"
    ### we use all GW triggers that are above "clean_threshold" to define dirty times
    dirtyseg = event.vetosegs(kwdict[gwchannel], float(myconf['clean_window']), float(myconf['clean_threshold']))

    ### set up the "samples" list that is passed to classifiers
    samples = []
    for trig in gw_trig: 
        ### check whether the trigger lives in the dirty segments
        unclean = 1 ### by definition, this must be true
        samples.append([trig[2], 1, unclean, trig[-1], math.sqrt(trig[-3] - trig[-2])]) ### magic numbers correspond to KW trigger parameters
                                                                                        ### FIXME: MAGIC NUMBERS ARE EVIL

    for gps in clean_gps:
        ### check whether the trigger lives in the dirty segments
        unclean = len(event.include([[gps]], dirtyseg, tcent=0))  # will be 1 if inside dirtyseg
        ### only include clean samples if they are truly clean
        ### THIS IS IMPORTANT FOR TRAINING JOBS THAT RELY ON PAT FILES GENERATED IN THIS SCRIPT
        if unclean == 0: ### truly clean
            samples.append([gps, 0, unclean, 0, 0])


    ### check whether we need to include auxiliary triggers from neighboring strides
    if len(gw_gps + clean_gps) > 0:
        mintime = min(gw_gps + clean_gps) ### earliest trigger present
        maxtime = max(gw_gps + clean_gps) ### latest trigger present
    else:
         ### do not mess up if no triggers or clean samples
        (mintime, maxtime) = (t + stride / 2, t + stride / 2) ### nothing to do, so both are the midpoint of the stride

    ### look for neighboring triggers if needed
    ### do not abort if missing
    padding = float(myconf['padding'])

    if mintime < t + padding: ### previous triggers must be included
        logger.info('  previous stride KW triggers are needed')
        kwprevious = '%s/%s/%s-%d/%s-%d-%d.trg' % (
            config.get('general', 'gdsdir'),
            kwbasename,
            kwbasename,
            (t - stride) / 1e5,
            kwbasename,
            t - stride,
            stride,
            )

        if os.path.exists(kwprevious):
            kwdict = event.loadkwm(kwprevious, trigs_dict=kwdict) ### add these to the existing kwdict
            logger.info('  adding previous stride KW triggers')
        else:
            logger.info('  missing previous stride KW triggers')

    if maxtime > t + stride - padding: ### triggers from the following stride must be included
        logger.info('  next stride KW triggers are needed')
        kwnext = '%s/%s/%s-%d/%s-%d-%d.trg' % (
            config.get('general', 'gdsdir'),
            kwbasename,
            kwbasename,
            (t + stride) / 1e5,
            kwbasename,
            t + stride,
            stride,
            )

        if os.path.exists(kwnext): ### load next triggers since they already exist
            kwdict = event.loadkwm(kwnext, trigs_dict=kwdict)
            logger.info('  adding next stride KW triggers')

        else:
            if wait > -2 * stride: ### we wait up to one additional stride for these KW triggers to appear
                logger.info('  missing next KW triggers, waiting additional %d seconds' % stride)
                time.sleep(stride)

                if os.path.exists(kwnext): ### file now exists
                    kwdict = event.loadkwm(kwnext, trigs_dict=kwdict) ### add it into the existing kwdict
                else:
                    logger.info('  still missing next KW triggers')
            else: ### file should have been there, but it isn't
                logger.info('  missing next KW triggers')

    kwdict.resort() ### make sure the kwdict's trigger lists are stored in the correct order

    ### get rid of unwanted AUX triggers from previous or next strides
    kwdict.include([[min(mintime - padding, t), max(maxtime + padding, t + stride)]])

    ### get rid of unwanted GW triggers from previous or next strides
    ### these will be evaluated in their own strides
    kwdict.include([[t, t + stride]], channels=[gwchannel])

    logger.info('Done: aggregating triggers')

    logger.info('N_GW: %d, N_GW > thr: %d, N_clean: %d, N_total: %d'
                % (len(kwdict[gwchannel]), len(gw_trig),
                len(clean_gps), len(gw_gps) + len(clean_gps)))

    #=============================================
    # generate patfiles for mla classifiers
    #=============================================
    #### output dir for evaluation jobs
    this_output_dir = '%s/%s-%d/' % (realtimedir, kwbasename, t / 1e5)
    if not os.path.exists(this_output_dir):
        os.makedirs(this_output_dir)

    if mla:  # only build patfiles if machine-learning algorithms are present
        logger.info('Begin: building auxmvc feature vectors ...')

        ### reading parameters from config file
        auxmvc_coinc_window = config.getfloat('build_auxmvc_vectors','time-window')
        auxmc_gw_signif_thr = config.getfloat('build_auxmvc_vectors','signif-threshold')
        auxmvc_selected_channels = config.get('general','selected-channels')
        auxmvc_unsafe_channels = config.get('general','unsafe-channels')

        auxmvc_pat_file_name = "%s%s_mvsc-%d-%d.pat"%(this_output_dir, kwbasename, t, stride)

        # generating auxmvc vector samples. result is saved into pat file
        # FIXME: depending how padding is done we  should adjust behavior of  build_auxmvc_vectors
        # Currently it keeps gw trigger from [t, t + stride] and uses time_window to pad this segment for auxiliary triggers
        auxmvc_vectors = idq.build_auxmvc_vectors(
            kwdict,
            config.get('general', 'gwchannel'),
            auxmvc_coinc_window,
            auxmc_gw_signif_thr,
            auxmvc_pat_file_name,
            gps_start_time=t,
            gps_end_time=t + stride,
            channels=auxmvc_selected_channels,
            unsafe_channels=auxmvc_unsafe_channels,
            clean_times=clean_gps,
            clean_window=float(myconf['clean_window']),
            )

        logger.info('Done: building auxmvc feature vectors.')

    #=============================================
    # mvsc prediction
    #=============================================
    if mvsc:
        logger.info('Begin: MVSC cycle ...')
        logger.info('  Begin: executing MVSC evaluate job ...')

        ### get lastest trained forest configuration
        if mvsc_train_cache_modification_time != int(os.path.getmtime(mvsc_train_cache)): ### file has been modified since we last checked
            try:
                trainedforest = open(mvsc_train_cache, 'r').readlines()[-1].strip('\n')
            except IndexError:
                raise StandardError, "%s appears to be empty. Make sure you train your classifiers before evaluation"%mvsc_train_cache

            mvsc_train_cache_modification_time = int(os.path.getmtime(mvsc_train_cache))
            logger.info('    Use the newest MVSC model :' + trainedforest)

            ### getting start and end of training period from the name of the trained forest file
            mvsc_start_training_period = int(trainedforest.split('-')[-2])
            mvsc_end_training_period = mvsc_start_training_period + int(trainedforest.split('-')[-1].split('.')[0])
            mvsc_trained_range = "%d_%d"%(mvsc_start_training_period, mvsc_end_training_period)

        ### output *dat file
        ranked_file = "%smvsc_%s_%s-%d-%d.dat"%(auxmvc_pat_file_name.split('mvsc')[0], mvsc_trained_range, usertag, t, stride)

        if len(auxmvc_vectors) > 0: ### run mvsc prediction jobs
            (exit_status, _) = idq.execute_forest_evaluate(
                auxmvc_pat_file_name,
                trainedforest,
                ranked_file,
                config,
                gps_start_time=t,
                gps_end_time=t + stride,
                dir=this_output_dir,
                )
        else: ### nothing to classify, write empty dat file
            nonaux_vars = [
                'index',
                'i',
                'w',
                'GPS_s',
                'GPS_ms',
                'signif',
                'SNR',
                'unclean',
                'glitch-rank',
                ]
            varline = 'GPS i w unclean signif SNR rank ' + ' '.join([var for var in auxmvc_vectors.dtype.names if not var in nonaux_vars]) + '\n'
            file = open(ranked_file, 'w')
            file.write(varline)
            file.close()
            exit_status = 0

        logger.info('  Done: executing MVSC evaluate job.')

        ### check if process has been execited correctly
        if exit_status != 0:
            logger.info('  WARNING: Forest predict failed')
            logger.info('  skipping Forest timeseries')

        else: ### job completed sucessfully, so we continue with output
            ### create timeseries from datfile and save as gzip numpy array
            ### there are boundary effects here
            logger.info('  Begin: generating MVSC rank timeseries.')
            ts = idq.datfile2timeseries(ranked_file, window=0.1, fs=ts_fs)

            ### save timeseries to file
            rank_ts_filename = "%s/%s_idq_mvsc_%s_rank_%s-%d-%d.npy.gz"%(this_output_dir, ifo, mvsc_trained_range, usertag, t, stride)
            idq.numpy.save(event.gzopen(rank_ts_filename, 'w'), ts)

            logger.info('  Done: generating MVSC rank timeseries.')

            ### create FAP timeseries
            logger.info('  Begin: generating MVSC FAP time-series')
            ### load most recent mvsc*.uroc file from cache
            if mvsc_uroc_cache_modification_time != int(os.path.getmtime(mvsc_uroc_cache)):
                logger.info('    Updating newest mvsc*uroc filename')
                try:
                    mvsc_uroc_filename = open(mvsc_uroc_cache, 'r').readlines()[-1].strip('\n')
                except IndexError:
                    raise StandardError, "%s appears to be empty. Make sure you run a summary job before evaluation"%mvsc_summary_cache
                mvsc_uroc_cache_modification_time = int(os.path.getmtime(mvsc_uroc_cache))

                logger.info('    Using mvsc_Effmap and mvsc_FAPmap to estimate likelihood map')
                (mvsc_Effmap, mvsc_FAPmap) = idq.rank_to_EffFAPmap(mvsc_uroc_filename)

            ### convert rank time-series into FAP time-series
            fap_ts = mvsc_FAPmap(ts)
            fap_ts_filename = "%s/%s_idq_mvsc_%s_fap_%s-%d-%d.npy.gz"%(this_output_dir, ifo, mvsc_trained_range, usertag, t, stride)
            fap_ts_filename = this_output_dir + '/' + ifo \
            
            ### save file to disk
            idq.numpy.save(event.gzopen(fap_ts_filename, 'w'), fap_ts)
            logger.info('  Done: generating MVSC FAP time-series')

            ### convert datfiles to xml tables
            logger.info('  Begin: converting MVSC dat file into xml files')
            logger.info('    converting %s to xml tables' % ranked_file)

            ### read dafile -> xml docs
            (gchxml_doc, clnxml_doc) = idq.datfilename_to_xmldocs(
                ranked_file,
                ifo,
                mvsc_FAPmap,
                Lmap=False,
                Effmap=mvsc_Effmap,
                classifier='mvsc',
                )

            ### write documents
            gchxml_filename = "%s/%s_idq_mvsc_%s_glitche_%s_%d_%d.xml.gz"%(this_output_dir, ifo, mvsc_trained_range, usertag, t, stride)
            logger.info('    --> writing ' + gchxml_filename)
            idq.ligolw_utils.write_filename(gchxml_doc, gchxml_filename, gz=gchxml_filename.endswith('.gz'))

            clnxml_filename = "%s/%s_idq_mvsc_%s_clean_%s-%d-%d.xml.gz"%(this_output_dir, ifo, mvsc_trained_range, usertag, t, stride)
            logger.info('    --> writing ' + clnxml_filename)
            idq.ligolw_utils.write_filename(clnxml_doc, clnxml_filename, gz=clnxml_filename.endswith('.gz'))

            logger.info('  Done: converting MVSC dat file into xml files.')

            logger.info('Done: MVSC cycle.')

    #=============================================
    # ovl prediction
    #=============================================
    if ovl:
        logger.info('Begin: OVL cycle ...')
        logger.info('	Begin: executing OVL evaluate job ...')

        if ovl_train_cache_modification_time != int(os.path.getmtime(ovl_train_cache)): ### file has been modified since we last checked
            try:
                vetolist = open(ovl_train_cache, 'r').readlines()[-1].strip('\n')
            except IndexError:
                raise StandardError, "%s appears to be empty. Make sure you train your classifiers before evaluation"%ovl_train_cache

            logger.info('    Use the newest OVL vetolist ' + vetolist)
            ovl_train_cache_modification_time = int(os.path.getmtime(ovl_train_cache))

        ### determine time period corresponding to this vetolist file
        ovl_trained_range = [v for v in vetolist.split('/') if v != ''][-3]  # we assume a particular directory structure

        ### ovl prediction filename
        ovl_filename = "%s_ovl_%s_%s-%d-%d.dat"%(kwbasename, ovl_trained_range, usertag, t, stride)

        if len(samples) > 0: ### launch evaluate job
            try:
                (gw_predict, ovl_filename) = idq.ovl_evaluate(
                    vetolist,
                    GPStimes=samples,
                    GPS_headers=samples_header,
                    allvtrg=kwdict,
                    filename=ovl_filename,
                    output_dir=this_output_dir,
                    )

                logger.info('  Done: executing OVL evaluate job ...')
                exit_status = 0
            except:
                e = sys.exc_info()
                logger.warning('WARNING: OVL evaluate job failed with the following error message: %s, %s, %s' % (e[0], e[1], e[2]))
                exit_status = 1
        else: ### nothing to classify, write empty dat file

            ovl_dat_vars = ['rank', 'eff/dt', 'vchan', 'vthr', 'vwin']
            varline = ' '.join(samples_header) + ' ' + ' '.join(ovl_dat_vars) + '\n'
            file = open(this_output_dir + '/' + ovl_filename, 'w')
            file.write(varline)
            file.close()
            exit_status = 0

        logger.info('  Done: executing OVL evaluate job.')

        if exit_status != 0:
            logger.info('  WARNING: OVL predict failed')
            logger.info('  skipping OVL timeseries')

        else: ### create timeseries and save as gzip numpy array with standard filename

            ### generate OVL rank timeseries
            logger.info('  Begin: generating OVL rank timeseries ...')

            ### generate timeseries
            ts = idq.ovl.vetolist2timeseries(vetolist, kwdict, [t, t + stride], fs=ts_fs)

            ### write file to disk
            rank_ts_filename = "%s/%s_idq_ovl_%s_rank_%s-%d-%d.npy.gz"%(this_output_dir, ifo, ovl_trained_range, usertag, t, stride)
            idq.numpy.save(event.gzopen(rank_ts_filename, 'w'), ts)

            logger.info('	  --> generated ' + this_output_dir + '/' + ovl_filename)
#            logger.info('	  --> ovl made predictions about %d events' % len(gw_predict))
            logger.info('  Done: generating OVL rank timeseries ...')

            ### generate OVL FAP timeseries
            logger.info('  Begin: generating OVL FAP timeseries ...')

            ### load most recent ovl*.uroc file from cache
            if ovl_uroc_cache_modification_time != int(os.path.getmtime(ovl_uroc_cache)): ### file was touched since we last checked
                logger.info('	  Updating newest ovl*uroc filename')
                try:
                    ovl_uroc_filename = open(ovl_uroc_cache, 'r').readlines()[-1].strip('\n')
                except IndexError:
                    raise StandardError, "%s appears to be empty. Make sure you run a summary job before evaluation"%ovl_uroc_cache

                ovl_uroc_cache_modification_time = int(os.path.getmtime(ovl_uroc_cache))

                ### compute rank->FAP maps
                logger.info('	  Using ovl_Effmap and ovl_FAPmap to estimate likelihood map')
                (ovl_Effmap, ovl_FAPmap) = idq.rank_to_EffFAPmap(ovl_uroc_filename)

            ### convert rank time-series into FAP time-series
            fap_ts = ovl_FAPmap(ts)

            ### write fap timeseries to disk
            fap_ts_filename = "%s/%s_idq_ovl_%s_fap_%s-%d-%d.npy.gz"%(this_output_dir, ifo, ovl_trained_range, usertag, t, stride)
            idq.numpy.save(event.gzopen(fap_ts_filename, 'w'), fap_ts)

            logger.info('  Done: generating OVL FAP timeseries.')

            ### convert dat files to xml tables
            logger.info('  Begin: converting OVL dat file into xml files ...')
            logger.info('	  converting %s to xml tables' % (this_output_dir + '/' + ovl_filename))

            ### read dafile -> xml docs
            (gchxml_doc, clnxml_doc) = idq.datfilename_to_xmldocs(
                this_output_dir + '/' + ovl_filename,
                ifo,
                ovl_FAPmap,
                Lmap=False,
                Effmap=ovl_Effmap,
                classifier='ovl',
                )

            ### write documents
            gchxml_filename = "%s/%s_idq_ovl_%s_glitch_%s-%d-%d.xml.gz"%(this_output_dir, ifo, ovl_trained_range, usertag, t, stride)
            logger.info('    --> writing ' + gchxml_filename)
            idq.ligolw_utils.write_filename(gchxml_doc, gchxml_filename, gz=gchxml_filename.endswith('.gz'))

            clnxml_filename = "%s/%s_idq_ovl_%s_clean_%s-%d-%d.xml.gz"%(this_output_dir, ifo, ovl_trained_range, usertag, t, stride)
            logger.info('	  --> writing ' + clnxml_filename)
            idq.ligolw_utils.write_filename(clnxml_doc,
                    clnxml_filename, gz=clnxml_filename.endswith('.gz'))

            logger.info('  Done: converting OVL dat file into xml files ...' )

            logger.info('Done: OVL cycle.')

    #=============================================
    # svm prediction
    #=============================================

    if svm:
        logger.info('Begin: SVM cycle ...')
        logger.info('  Begin: executing SVM evaluate job ...')

        ### get lastest trained SVM configuration
        if svm_train_cache_modification_time != int(os.path.getmtime(svm_train_cache)): ### file was touched since we last checked
            try:
                [svm_model, svm_range_file] = open(svm_train_cache, 'r').readlines()[-2:]
            except IndexError:
                raise StandardError, "%s appears to be empty. Make sure you train your classifiers before evalutation"%svm_train_cache

            svm_model = svm_model.strip('\n')
            svm_range_file = svm_range_file.strip('\n')
            svm_train_cache_modification_time = int(os.path.getmtime(svm_train_cache))

            logger.info('    Use the newest SVM model ' + svm_model)

            ## start and end of training period from the name of the svm model file
            svm_start_training_period = int(svm_model.split('-')[-2])
            svm_end_training_period = svm_start_training_period + int(svm_model.split('-')[-1].split('.')[0])
            svm_trained_range = "%d_%d"%(svm_start_training_period,svm_end_training_period)

        svm_predict_file = "%ssvm_%s_%s-%d-%d.dat"%(auxmvc_pat_file_name.split('mvsc')[0] , svm_trained_range, usertag, t, stride)

        if len(auxmvc_vectors) > 0: ### run svm evaluate job
            exit_label = idq.execute_svm_evaluate(
                config,
                auxmvc_pat_file_name,
                svm_range_file,
                svm_model,
                svm_predict_file,
                dir=this_output_dir,
                )
        else: ### nothing to classify, write empty dat file
            nonaux_vars = [
                'index',
                'i',
                'w',
                'GPS_s',
                'GPS_ms',
                'signif',
                'SNR',
                'unclean',
                'glitch-rank',
                ]
            varline = 'GPS i w unclean signif SNR rank ' + ' '.join([var for var in auxmvc_vectors.dtype.names if not var in nonaux_vars]) + '\n'
            file = open(svm_predict_file, 'w')
            file.write(varline)
            file.close()
            exit_label = 0

        logger.info('  Done: executing SVM evaluate job ...')

        if exit_label != 0:
            logger.info('  WARNING: SVM predict failed')
            logger.info('  skipping SVM timeseries')

        else: ### create timeseries from datfile and save as gzip numpy array, there are boundary effects here
            logger.info('  Begin: generating SVM rank timeseries ...')

            ### generate rank timeseries
            ts = idq.datfile2timeseries(svm_predict_file, window=0.1, fs=ts_fs)

            ### write to file
            rank_ts_filename = "%s/%s_idq_svm_%s_rank_%s-%d-%d.npy.gz"%(this_output_dir, ifo, svm_trained_range, usertag, t, stride)
            idq.numpy.save(event.gzopen(rank_ts_filename, 'w'), ts)

            logger.info('  Done: generating SVM rank timeseries.')

            ### generate FAP timeseries
            logger.info('  Begin: generating SVM FAP timeseries ...')

            ### load most recent svm*.uroc file from cache
            if svm_uroc_cache_modification_time != int(os.path.getmtime(svm_uroc_cache)):
                logger.info('    Updating newest svm*uroc filename')
                try:
                    svm_uroc_filename = open(svm_uroc_cache, 'r').readlines()[-1].strip('\n')
                except IndexError:
                    raise StandardError, "%s appears to be empty. Make sure you run a summary job before evaluation"%svm_uroc_cache

                svm_uroc_cache_modification_time = int(os.path.getmtime(svm_uroc_cache))

                logger.info('    Using svm_Effmap and svm_FAPmap to estimate likelihood map')
                (svm_Effmap, svm_FAPmap) = idq.rank_to_EffFAPmap(svm_uroc_filename)

            ### convert rank time-series into FAP time-series
            logger.info('  Begin: generating SVM FAP time-series ...')
            fap_ts = svm_FAPmap(ts)

            ### save to disk
            fap_ts_filename = "%s/%s_idq_svm_%s_fap_%s-%d-%d.npy.gz"%(this_output_dir, ifo, svm_trained_range, usertag, t, stride)
            idq.numpy.save(event.gzopen(fap_ts_filename, 'w'), fap_ts)

            logger.info('  Done: generating SVM FAP timeseries.')

            ### convert dat files to xml tables
            logger.info('  Begin: converting SVM dat file into xml files ...')
            logger.info('    converting %s to xml tables' % svm_predict_file)

            ### read dafile -> xml docs
            (gchxml_doc, clnxml_doc) = idq.datfilename_to_xmldocs(
                svm_predict_file,
                ifo,
                svm_FAPmap,
                Lmap=False,
                Effmap=svm_Effmap,
                classifier='svm',
                )

            ### write documents
            gchxml_filename = "%s/%s_idq_svm_%s_glitch_%s-%d-%d.xml.gz"%(this_output_dir, ifo, svm_trained_range, usertag, t, stride)
            logger.info('    --> writing ' + gchxml_filename)
            idq.ligolw_utils.write_filename(gchxml_doc, gchxml_filename, gz=gchxml_filename.endswith('.gz'))

            clnxml_filename = "%s/%s_idq_svm_%s_clean_%s-%d-%d.xml.gz"%(this_output_dir, ifo, svm_trained_range, usertag, t, stride)
            logger.info('    --> writing ' + clnxml_filename)
            idq.ligolw_utils.write_filename(clnxml_doc, clnxml_filename, gz=clnxml_filename.endswith('.gz'))

            logger.info('  Done: converting SVM dat file into xml files.')

            logger.info('Done: SVM cycle.')

    ### proceed to next epoch
    logger.info('Done: stride %d-%d' % (t, t + stride))

    t += stride

### loop exited, which means we must have run past options.gpsend
logger.info('End real-time evaluation')
logger.info('t + stride = %d > %d = endgps'%(t+stride, options.endgps))


