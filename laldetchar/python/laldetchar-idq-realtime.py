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

import subprocess

import time
import logging
import tempfile
import ConfigParser

import numpy
import math

#from laldetchar.idq import idq
from laldetchar.idq import reed
from laldetchar.idq import event

from laldetchar import git_version

from optparse import OptionParser

#===================================================================================================

__prog__ = 'laldetchar-idq-realtime'

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
""" 
This program performs real-time identification and classification of glitch events. 
Currently it operates on equal length strides (e.g. 32 seconds). 
Within the iDQ pipeline this script plays the central role by scheduling all tasks including periodic training and generation of summary html pages.
"""

#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'%git_version.verbose_msg, 
                                usage='%prog [options]', 
                                description=description)
parser.add_option('-c', '--config-file', dest='config_file',
    help='configuration file', metavar='FILE', default='idq.ini' )

parser.add_option('-s', '--start-time', dest='startgps',
    type='int', help='first stride begins at or after GPS time',
    metavar='GPS', default=reed.nowgps(), )
parser.add_option('-e', '--end-time', dest='endgps',
    type='int', help='last stride ends at or before GPS time',
    metavar='GPS', default=1e10 )

parser.add_option('-k', '--lock-file', dest='lockfile',
    help='use custom lockfile', metavar='FILE', default='.idq_realtime.lock',
 )
parser.add_option('-l', '--log-file', default='idq_realtime.log',
                  type='string', help='log file')

parser.add_option('', '--ignore-science-segments', 
    default=False, action="store_true", 
    help='analyze strides regardless of the science segment content. \
    This is NOT passed on to training job, which always uses science segments'
    )
parser.add_option('', '--ignore-science-segments-training', 
    default=False, action="store_true", 
    help="forces training jobs to ignore science segments as well. \
    This should NOT be used to generate actual DQ information, but may be useful when debugging the pipeline."
    )
parser.add_option('', '--ignore-science-segments-summary', 
    default=False, action="store_true", 
    help="forces summary jobs to ignore science segments as well. \
    This should NOT be used to generate actual DQ information, but may be useful when debugging the pipeline."
    )
parser.add_option('', '--ignore-science-segments-calibration',
    default=False, action="store_true",
    help="forces calibration jobs to ignore science segments as well. \
    this should NOT be used to generate actual DQ information, but may be useful when debugging the pipeline."
    )

parser.add_option("", "--no-robot-cert",
    default=False, action="store_true",
    help="do not use robot cert with segment query in training jobs"
    )
parser.add_option("", "--offline",
    default=False, action="store_true",
    help="wait for training and summary jobs to finish before continuing. Only meant for offline analysis"
    )

(opts, args) = parser.parse_args()

cwd = os.getcwd()

#===================================================================================================
### setup logger to record processes
logger = reed.setup_logger('idq_logger', opts.log_file, sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = reed.LogFile(logger)
sys.stderr = reed.LogFile(logger)

#===================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config_file)

mainidqdir = config.get('general', 'idqdir') ### get the main directory where idq pipeline is going to be running.

reed.dieiflocked('%s/%s'%(config.get('general', 'idqdir'),opts.lockfile) ) ### prevent multiple copies from running

gchtag = "_glitch" ### used for xml filenames
clntag = "_clean"
rnktag = "_rank"
faptag = "_fap"
usertag = config.get('general', 'usertag')
if usertag:
    usertag = "_%s"%usertag
    gchtag = "%s_%s"%(gchtag, usertag)
    clntag = "%s_%s"%(clntag, usertag)
    rnktag = "%s_%s"%(rnktag, usertag)
    faptag = "%s_%s"%(faptag, usertag)

ifo = config.get('general', 'ifo')

#========================
# which classifiers
#========================
classifiers = sorted(set(config.get('general', 'classifiers').split()))
if not classifiers:
    raise ValueError("no classifiers in general section of %s"%opts.config_file)

### ensure we have a section for each classifier and fill out dictionary of options
classifiersD, mla, ovl = reed.config_to_classifiersD( config )

if mla:
    ### reading parameters from config file needed for mla
    auxmvc_coinc_window = config.getfloat('build_auxmvc_vectors','time-window')
    auxmc_gw_signif_thr = config.getfloat('build_auxmvc_vectors','signif-threshold')
    auxmvc_selected_channels = config.get('general','selected-channels')
    auxmvc_unsafe_channels = config.get('general','unsafe-channels')

#========================
# data discovery
#========================
if not opts.ignore_science_segments:
    ### load settings for accessing dmt segment files
    dmt_segments_location = config.get('get_science_segments', 'xmlurl')
    dmtdq_name = config.get('get_science_segments', 'include')
#    dmtdq_name = config.get('get_science_segments', 'include').split(':')[1]

### define the gravitational wave channel/what is a glitch
gwchannel = config.get('general', 'gwchannel')
gwthreshold = config.getfloat('general', 'gw_kwsignif_thr')

### kleineWelle config
GWkwconfig = reed.loadkwconfig(config.get('data_discovery', 'GWkwconfig'))
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')

AUXkwconfig = reed.loadkwconfig(config.get('data_discovery', 'AUXkwconfig'))
AUXkwbasename = AUXkwconfig['basename']
AUXgdsdir = config.get('data_discovery', 'AUXgdsdir')

identical_trgfile = (GWgdsdir == AUXgdsdir) and (GWkwbasename == AUXkwbasename) ### used to avoid re-loading the same file for both GW and AUX triggers

#========================
# realtime job
#========================
ts_fs = config.getfloat('realtime','sampling_rate') ### sampling frequency for time-series files

clean_rate = config.getfloat('realtime', 'clean_rate')
clean_window = config.getfloat('realtime', 'clean_window')
clean_threshold = config.getfloat('realtime', 'clean_threshold')

padding = config.getfloat('realtime', 'padding')

realtimedir = config.get('general', 'realtimedir')### output directory for realtime predictions
if not os.path.exists(realtimedir):
    os.makedirs(realtimedir)

samples_header = config.get('realtime', 'dat_columns').split()### headers for dat files

### slave the realtime job to the kw stride
stride = int(float(GWkwconfig['stride'])) ### this is given as a decimal, so we must cast as float first
delay = config.getint('realtime', 'delay') ### buffer to let jobs finish

#========================
# summary jobs
#========================
summarydir = config.get('general', 'summarydir')

summary_stride = config.getint('summary', 'stride')
summary_delay = config.getint('summary', 'delay')
summary_lookback = config.get('summary', 'lookback')
if summary_lookback != "infinity":
    summary_lookback = int(summary_lookback)
summary_trending = config.get('summary', 'trending')
if summary_trending != 'infinity':
    summary_trending = int(summary_trending)

summary_script = config.get('condor', 'summary')

#summary_cache = dict( (classifier, reed.Cachefile(reed.cache(summarydir, classifier, tag='_summary%s'%usertag))) for classifier in classifiers ) ### not needed?
#for cache in summary_cache.values():
#    cache.time = 0

summary_log = "%s/summary%s.log"%(mainidqdir, usertag)
summary_out = "%s/summary%s.out"%(mainidqdir, usertag)
summary_err = "%s/summary%s.out"%(mainidqdir, usertag)

#========================
# training jobs
#========================
traindir = config.get('general', 'traindir')
#if ovl: ### need snglchandir
#   snglchndir = config.get('general', 'snglchndir') 

train_stride = config.getint('train', 'stride')
train_delay = config.getint('train', 'delay')
train_lookback = config.get('train', 'lookback')
if train_lookback != "infinity":
    train_lookback = int(train_lookback)

train_script = config.get('condor', 'train')

train_cache = dict( (classifier, reed.Cachefile(reed.cache(traindir, classifier, tag='_train%s'%usertag))) for classifier in classifiers )
for cache in train_cache.values():
    cache.time = 0

train_log = "%s/train%s.log"%(mainidqdir, usertag)
train_out = "%s/train%s.out"%(mainidqdir, usertag)
train_err = "%s/train%s.err"%(mainidqdir, usertag)

#========================
# calibration jobs
#========================
calibrationdir = config.get('general', 'calibrationdir')

calibration_stride = config.getint('calibration', 'stride')
calibration_delay = config.getint('calibration', 'delay')
calibration_lookback = config.get('calibration', 'lookback')
if calibration_lookback != "infinity":
    calibration_lookback = int(calibration_lookback)

calibration_FAPthrs = [float(l) for l in config.get('calibration','FAP').split()]

calibration_cache = dict( (classifier, reed.Cachefile(reed.cache(calibrationdir, classifier, tag='_calibration%s'%usertag))) for classifier in classifiers )
for cache in calibration_cache.values():
    cache.time = 0

calibration_script = config.get('condor', 'calibration')

calibration_log = "%s/calibration%s.log"%(mainidqdir, usertag)
calibration_out = "%s/calibration%s.out"%(mainidqdir, usertag)
calibration_err = "%s/calibration%s.err"%(mainidqdir, usertag)

#===================================================================================================
### figure out where we are in time, and launch initial jobs in needed
t = int(math.ceil(opts.startgps / stride)) * stride  # require times to be integer multiples of stride
global_tstart = t

#===================================================================================================
#
# INITIAL TRAINING / CALIBRATION
#
#===================================================================================================
### figure out if we need to run initial training or calibration jobs
### training jobs allow us to compute ranks from sets of auxiliary glitches
### calibration jobs allow us to map those ranks into FAP measurements

initial_training = False
initial_calibration = False

for classifier in classifiers:
    if not initial_training:
        initial_training =  train_cache[classifier].is_empty()

    if not initial_calibration:		
        initial_calibration = calibration_cache[classifier].is_empty()

#=================================================
# Initial training 
#=================================================
if initial_training:
    logger.info('At least one training cache file was empty or did not exist. Launching initial training of classifiers.')
    logger.info('Begin: initial training')

    ### set up options for training script
    initial_training_lookback = config.get("train", "initial-lookback")

    train_start_time = t - train_stride
    train_command = "%s --config %s -l %s --gps-start %d --gps-stop %d --lookback %s --force"%(train_script, opts.config_file, train_log, train_start_time, t, initial_training_lookback)

    if opts.ignore_science_segments_training:
        train_command += " --ignore-science-segments" 
    elif opts.no_robot_cert:
        train_command += " --no-robot-cert"

    logger.info('Submiting training job with the following options:')
    logger.info(train_command)

    ### launch training script, wait for it to finish
    proc = subprocess.Popen( train_command.split() , stdout=open(train_out, 'a'), stderr=open(train_err, 'a'), cwd=cwd )
    returncode = proc.wait() ### block
    if returncode:
        logger.info('ERROR: initial training failed.')
        raise StandardError, "initial training failed"

    logger.info('Done: initial training')

else:
    logger.info('All required training cache files already exist and are not empty.')

#=================================================
# Initial calibration
#=================================================
if initial_calibration:
    logger.info('At least one calibration cache file was empty or did not exist. Launching initial calibration job for classifiers.')
    logger.info('Begin: initial calibration')
    logger.warning('WARNING: the *.uroc files generated by this job do not have much meaning')

    ### set up options for calibration script
    initial_calibration_lookback = config.get("calibration", "initial-lookback")

    calibration_start_time = t - calibration_stride
    calibration_command = "%s --config %s -l %s --gps-start %d --gps-stop %d --lookback %s --force"%(calibration_script, opts.config_file, calibration_log, calibration_start_time, t, initial_calibration_lookback)
    ### we skip validation (no FAPthr arguments) for the initial job

    if opts.ignore_science_segments_calibration:
        calibration_command += " --ignore-science-segments" 
    elif opts.no_robot_cert:
        calibration_command += " --no-robot-cert" 

    logger.info('Submiting calibration job with the following options:')
    logger.info(calibration_command) ### block

    ### launch calibration script, wait for it to finish
    proc = subprocess.Popen( calibration_command.split() , stdout=open(calibration_out, 'a'), stderr=open(calibration_err, 'a'), cwd=cwd )
    returncode = proc.wait()
    if returncode:
        logger.info('ERROR: initial calibration failed.')
        raise StandardError, "initial calibration failed"

    logger.info('Done: initial calibration')

else:
    logger.info('All required calibration cache files already exist and are not empty.')

#===================================================================================================
#
# PERSISTENT LOOP
#
#===================================================================================================
logger.info('Begin: realtime evaluation')

#========================
# set up template commands
#========================
train_template = "%s --config %s -l %s --gps-start %d --gps-stop %d --lookback %d"
if opts.ignore_science_segments_training:
    train_template += " --ignore-science-segments"

summary_template = "%s --config %s -l %s --gps-start %d --gps-stop %d --lookback %d --trending %d"
if opts.ignore_science_segments_summary:
    summary_template += " --ignore-science-segments"

calibration_template = "%s --config %s -l %s --gps-start %d --gps-stop %d --lookback %d"
if opts.ignore_science_segments_calibration:
    calibration_template += " --ignore-science-segments"

if opts.no_robot_cert:
    train_template += " --no-robot-cert"
    summary_template += " --no-robot-cert"
    calibration_template += " --no-robot-cert"

#========================
# begin loop
#========================
global_start = t

while t  < opts.endgps:

    #====================
    # TRAIN
    # training cycle launched if:
    #    FIGURE OUT LOGIC AND PUT IT HERE

    ### RUSLAN'S OLD DESCRIPTION
    # the remainder from divison of the time of the current evaluation step
    # by the training stride is positive and less than a single evaluation stride.
    # To avoid premature training launch, we also require that the difference between
    # the time of the current evaluation step and the start time exceeds a single training stride.
    #====================    
    t_remainder = t - (t / train_stride) * train_stride
    if (2*(t - opts.startgps) > train_stride) and (t_remainder > 0) and (t_remainder <= stride):

        if train_lookback == "infinity":
            lookback = t-global_start
        else:
            lookback = train_lookback

        logger.info('****************************************************')

        train_stop_time = (t / train_stride) * train_stride
        train_start_time = train_stop_time - train_stride

        logger.info('Begin: launching train script for period: %s - %s'%(train_start_time,train_stop_time))

        train_command = train_template%(train_script, opts.config_file, train_log, train_start_time, train_stop_time, lookback)

        logger.info('Submiting train script with the following opts:')
        logger.info(train_command)

        if opts.offline: ### running in "offline mode"
            logger.info('Running in OFFLINE mode, will wait until train is complete before proceeding with evaluation.')
            proc = subprocess.Popen( train_command.split() , stdout=open(train_out, 'a'), stderr=open(train_err, 'a') , cwd=cwd )
            proc.wait() ### blocks!
        else:
            ### running in realtime, don't wait for training job to complete so we can keep the latency low
            train_pid = subprocess.Popen(train_command, stdout=open(train_out, 'a'), stderr=open(train_err, 'a'), cwd=cwd ).pid ### only remember the pid, and let the process float

        logger.info('Done: launching train script for period: %d - %d'%(train_start_time, train_stop_time))
        logger.info('****************************************************')

    #====================
    # CALIBRATION
    # calibration cycle launched if:
    #    FIGURE OUT LOGIC AND PUT IT HERE

    ### OLD DESCRIPTION
    # the remainder from divison of the time of the current evaluation step
    # by the summary stride is positive and less than a single evaluation stride.
    # To avoid premature summary launch, we also require that the difference between
    # the time of the current evaluation step and the start time exceeds a single summary stride.
    #====================
    t_remainder = t - (t / calibration_stride) * calibration_stride
    if (2*(t - opts.startgps) > calibration_stride) and (t_remainder > 0) and (t_remainder <= stride):
        logger.info('****************************************************')

        if calibration_lookback == "infinity":
            lookback = t - global_start
        else:
            lookback = calibration_lookback

        calibration_stop_time = (t / calibration_stride) * calibration_stride
        calibration_start_time = calibration_stop_time - calibration_stride

        logger.info('Begin: launching calibration script for period: %s - %s'%(calibration_start_time, calibration_stop_time))

        calibration_command = calibration_template%(calibration_script, opts.config_file, calibration_log, calibration_start_time, calibration_stop_time, lookback)

        logger.info('Submiting calibration script with the following options:')
        logger.info(calibration_command)

        if opts.offline:
            logger.info("Running in OFFLINE mode, will wait until calibration is complete before proceeding with evaluation")
            proc = subprocess.Popen(calibration_command.split(), stdout=open(calibration_out, "a"), stderr=open(calibration_err, "a"), cwd=cwd)
            proc.wait() # block!
        else:
            summary_pid = subprocess.Popen(calibration_command.split(), stdout=open(calibration_out, 'a'), stderr=open(calibration_err, 'a'), cwd=cwd).pid ### only remember the pid, let the process float

        logger.info('Done: launching calibration script for period: %s - %s'%(calibration_start_time,calibration_stop_time))
        logger.info('****************************************************')

    #====================
    # SUMMARY 
    # summary cycle launched if:
    #    FIGURE OUT LOGIC AND PUT IT HERE

    ### OLD DESCRIPTION
    # the remainder from divison of the time of the current evaluation step
    # by the summary stride is positive and less than a single evaluation stride.
    # To avoid premature summary launch, we also require that the difference between
    # the time of the current evaluation step and the start time exceeds a single summary stride.
    #====================
    t_remainder = t - (t / summary_stride) * summary_stride
    if (2*(t - opts.startgps) > summary_stride) and (t_remainder > 0) and (t_remainder <= stride):
        logger.info('****************************************************')

        if summary_lookback == "infinity" and summary_trending == "infinity":
            lookback = trending = t-global_start
        elif summary_lookback == "infinity":
            lookback = t-global_start
            trending = summary_trending
        elif summary_trending == "infinity":
            lookback = summary_lookback
            trending = t-global_start
        else:
            lookback = summary_lookback
            trending = summary_trending

        summary_stop_time = (t / summary_stride) * summary_stride
        summary_start_time = summary_stop_time - summary_stride

        logger.info('Begin: launching summary script for period: %s - %s'%(summary_start_time,summary_stop_time))

        summary_command = summary_template%(summary_script, opts.config_file, summary_log, summary_start_time, summary_stop_time, lookback, trending)

        logger.info('Submiting summary script with the following options:')
        logger.info(summary_command)

        if opts.offline:
            logger.info("Running in OFFLINE mode, will wait until summary is complete before proceeding with evaluation")
            proc = subprocess.Popen(summary_command.split(), stdout=open(summary_out, "a"), stderr=open(summary_err, "a"), cwd=cwd)
            proc.wait() # block!
        else:                                              
            summary_pid = subprocess.Popen(summary_command.split(), stdout=open(summary_out, 'a'), stderr=open(summary_err, 'a'), cwd=cwd).pid ### only remember the pid, let the process float

        logger.info('Done: launching summary script for period: %s - %s'%(summary_start_time,summary_stop_time))
        logger.info('****************************************************')

    #===============================================================================================
    # REALTIME EVALUATION
    #===============================================================================================
    logger.info('----------------------------------------------------')

    logger.info('Begin: stride %d-%d' % (t, t + stride))

    ### make sure end time of analysis stride is not ahead of now-delay
    wait = t + stride + delay - int(reed.nowgps())
    if wait > 0:
        logger.info('  waiting %.1f seconds before processing' % wait)
        time.sleep(wait)

    #====================
    # science segments
    #====================
    if opts.ignore_science_segments:
        logger.info('analyzing data regardless of science segments')
        
    else:
        logger.info('Begin: querrying science segments')

        ### get DMT xml segments from disk
        ### logic about waiting and re-trying is built into retrieve_scisegs
        good, covered = reed.retrieve_scisegs(dmt_segments_location, dmtdq_name, t, stride, pad=0, sleep=delay, nretry=1, logger=logger)

        logger.info('Done: querrying science segments')

        if event.livetime(covered) < stride:
            logger.warning('unknown science coverage, skipping')
            t += stride
            continue
        elif event.livetime(good) < stride:
            logger.info('incomplete science coverage, skipping')
            t += stride
            continue

        logger.info('complete science coverage')

    #===================
    # trigger discovery
    #===================
    logger.info('Begin: aggregating triggers')

    trgdict = reed.retrieve_kwtrig(GWgdsdir, GWkwbasename, t, stride, sleep=(t+stride-reed.nowgps()) + delay, ntrials=2, logger=logger)
    if trgdict == None: ### no file found, so we skip this stride
        logger.warning('  skipping this stride because no gwchannel kwtrg found')
        t += stride
        continue

    #====================
    # set up samples for classifiers
    #====================
    samples = []
    clean_gps = []

    #====================
    # glitches
    #====================
    ### require GW channel to be in trigger dict
    if gwchannel not in trgdict:
        trgdict[gwchannel] = []

    gw_trig = [a for a in trgdict[gwchannel] if a[-1] >= gwthreshold] ### pull out only those GW glitches with sufficiently high significance   
    gw_gps = [a[2] for a in gw_trig] ### retrieve GW trigger central times
                                     ### MAGIC NUMBERS ARE EVIL

    unclean = 1 ### by definition, this must be true for glitches
    samples += [ [trg[2], 1, unclean, trg[-1], math.sqrt(trg[-3] - trg[-2])] for trg in gw_trig] ### format gw_trgs for samples
                                                                                           ### FIXME: MAGIC NUMBERS ARE EVIL
    #====================
    # cleans
    #====================
    ### define segments that are too close to glitch times to be really "clean"
    ### we use all GW triggers that are above "clean_threshold" to define dirty times
    dirtyseg = event.vetosegs(trgdict[gwchannel], clean_window, clean_threshold)

    clean_gps = sorted(event.randomrate(clean_rate, [[t, t + stride]])) ### generate random clean times as a poisson time series
    clean_gps = [ l[0] for l in event.exclude( [[gps] for gps in clean_gps], dirtyseg, tcent=0)] ### keep only those gps times that are outside of dirtyseg

    unclean = 0 ### axiomatic, because we threw out everything within dirtyseg
    samples += [ [gps, 0, unclean, 0, 0] for gps in clean_gps] 

    #====================
    # aux triggers
    #====================
    ### load auxiliary triggers from current stride
    if not identical_trgfile: ### current AUX file is different from current GW file
        aux_trgdict = reed.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t, stride, sleep=(t+stride-reed.nowgps()) + delay, ntrials=2, logger=logger)
        if aux_trgdict == None:
            logger.warning('  no auxiliary triggers were found')
            ### we do not skip, although we might want to?
        else:
            trgdict.add( aux_trgdict )

    ### check whether we need to include auxiliary triggers from neighboring strides
    if len(gw_gps + clean_gps) > 0:
        mintime = min(gw_gps + clean_gps) ### earliest trigger present
        maxtime = max(gw_gps + clean_gps) ### latest trigger present
    else:
        ### do not mess up if no triggers or clean samples
        mintime = maxtime = t + stride/2 ### nothing to do, so both are the midpoint of the stride

    if mintime < t + padding: ### previous triggers must be included
        logger.info('  previous stride KW triggers are needed')
        prev_trgdict = reed.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t-stride, stride, sleep=0, ntrials=1, logger=logger)
        if prev_trgdict == None:
            logger.warning('  missing previous stride KW triggers')
            ### we do not skip this stride when previous kw stride is required and missing!
        else:
            trgdict.add( prev_trgdict ) ### we add the previous triggers into the current trgdict

    if maxtime > t + stride - padding: ### triggers from the following stride must be included
        logger.info('  next stride KW triggers are needed')
        next_trgdict = reed.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t+stride, stride, sleep=(t+stride+stride-reed.nowgps()) + delay, ntrials=2, logger=logger)
        if next_trgdict == None:
            logger.warning('  missing next KW triggers')
            ### we do not skip this stride when next kw stride is required and missing!
        else:
            trgdict.add( next_trgdict ) ### we add the next triggers to the current trgdict

    trgdict.resort() ### make sure the trgdict's trigger lists are stored in the correct order

    ### get rid of unwanted AUX triggers
    trgdict.include([[min(mintime - padding, t), max(maxtime + padding, t + stride)]])

    ### get rid of unwanted GW triggers from previous or next strides; these will be evaluated in their own strides
    trgdict.include([[t, t + stride]], channels=[gwchannel])

    logger.info('Done: aggregating triggers')

    #====================
    # report total numbers
    #====================
    Ntotal = len(samples)
    logger.info('N_GW: %d, N_GW > thr: %d, N_clean: %d, N_total: %d' % (len(trgdict[gwchannel]), len(gw_trig), len(clean_gps), Ntotal))

    #### output dir for evaluation jobs
    this_output_dir = '%s/%s-%d/' % (realtimedir, GWkwbasename, t / 1e5)
    if not os.path.exists(this_output_dir):
        os.makedirs(this_output_dir)

    #====================
    # generate patfiles for mla classifiers
    #====================
    if mla:  # only build patfiles if machine-learning algorithms are present
        logger.info('Begin: building auxmvc feature vectors ...')

        pat = reed.pat(this_output_dir, GWkwbasename, t, stride)

        # generating auxmvc vector samples. result is saved into pat file
        # FIXME: depending how padding is done we should adjust behavior of build_auxmvc_vectors
        # Currently it keeps gw trigger from [t, t + stride] and uses time_window to pad this segment for auxiliary triggers
		# we do not filter out unclean beacuse it is already done when clean_gps times are formed
        auxmvc_vectors = reed.build_auxmvc_vectors(trgdict, gwchannel, auxmvc_coinc_window, auxmc_gw_signif_thr, pat, gps_start_time=t,
                                gps_end_time=t + stride,  channels=auxmvc_selected_channels, unsafe_channels=auxmvc_unsafe_channels, clean_times=clean_gps,
                                clean_window=clean_window, filter_out_unclean=False )

        logger.info('Done: building auxmvc feature vectors.')

    #=============================================
    # predictions
    #=============================================
    for classifier in classifiers:
        flavor = classifiersD[classifier]['flavor']

        logger.info('Begin: %s cycle -> flavor=%s'%(classifier, flavor))

        ### check for new training data
        cache = train_cache[classifier]
        if cache.was_modified():
            lines = cache.tail(nlines=reed.traincache_nlines[flavor])
            trained_range = reed.extract_trained_range(lines, flavor)
            
            classifiersD[classifier]['lines'] = lines
            classifiersD[classifier]['trained_range'] = trained_range

            cache.set_time()
            logger.info('using newest %s train data : %s'%(classifier, ",".join(lines)))

        else:
            lines = classifiersD[classifier]['lines']
            trained_range = classifiersD[classifier]['trained_range']

        ### check for new calibration data
        cache = calibration_cache[classifier]
        if cache.was_modified():
            uroc = cache.tail(nlines=1)[0]
            calib_range = reed.extract_calib_range(uroc)
            r, c, g, _, _ = reed.file_to_rcg(uroc)
            effmap, fapmap = reed.rcg_to_EffFAPmap(r, c, g)
            
            classifiersD[classifier]['uroc'] = uroc
            classifiersD[classifier]['calib_range'] = calib_range
            classifiersD[classifier]['effmap'] = effmap
            classifiersD[classifier]['fapmap'] = fapmap

            cache.set_time()
            logger.info('using newest %s calibration data : %s'%(classifier, uroc))
        else:
            calib_range = classifiersD[classifier]['calib_range']
            effmap = classifiersD[classifier]['effmap']
            fapmap = classifiersD[classifier]['fapmap']

        ### compute filenames
        dat = reed.dat(this_output_dir, classifier, ifo, trained_range, usertag, t, stride)

        gchxml = reed.xml(this_output_dir, classifier, ifo, trained_range, calib_range, gchtag, t, stride)
        clnxml = reed.xml(this_output_dir, classifier, ifo, trained_range, calib_range, clntag, t, stride)

        rnknp = reed.timeseries(this_output_dir, classifier, ifo, trained_range, calib_range, rnktag, t, stride)
        fapnp = reed.timeseries(this_output_dir, classifier, ifo, trained_range, calib_range, faptag, t, stride)

        ### perform evalutation
        logger.info('  Begin %s evaluation -> %s'%(classifier, dat) )
        if classifiersD[classifier]['mla']:
            returncode = reed.evaluate(flavor, lines, dat, config, gps_start_time=t, gps_end_time=t+stride, dir=this_output_dir, trgdict=pat, auxmvc_vectors=auxmvc_vectors)
        else:
            returncode = reed.evaluate(flavor, lines, dat, config, gps_start_time=t, gps_end_time=t+stride, dir=this_output_dir, trgdict=trgdict, samples=samples, samples_header=samples_header)

        ### check if process has been execited correctly
        if returncode:
            logger.warning('  WARNING: %s predict failed'%classifier)
            logger.info('  skipping %s timeseries'%classifier)
            continue
        logger.info('  Done: %s evaluation '%classifier)

        ### create timeseries from datfile and save as gzip numpy array, there are boundary effects here
        logger.info('  Begin: generating %s rank timeseries -> %s'%(classifier, rnknp))
        ts = reed.datfile2timeseries(flavor, dat, lines, trgdict=trgdict, start=t, stride=stride, window=0.1, fs=256)
        numpy.save(event.gzopen(rnknp, 'w'), ts)
        logger.info('  Done: generating %s rank timeseries'%classifier)

        ### create FAP timeseries
        logger.info('  Begin: generating %s FAP time-series -> %s'%(classifier, fapnp))

        fap_ts = fapmap(ts) ### MLE estimate of FAP (interpolated)
        fapUL_ts = fapmap.ul(ts, conf=0.99) ### upper limit with FAPconf (around interpolated points)
                                            ### we hard-code this intentionally to make it difficult to change/mess up
        numpy.save(event.gzopen(fapnp, 'w'), (fap_ts, fapUL_ts) )

        logger.info('  Done: generating %s FAP time-series'%classifier)

        ### convert datfiles to xml tables
        logger.info('  Begin: converting %s dat file into xml files'%classifier)
        logger.info('    converting %s to xml tables' % dat)

        ### read dafile -> xml docs
        (gchxml_doc, clnxml_doc) = reed.datfile2xmldocs(dat, ifo, fapmap, Lmap=False, Effmap=effmap, flavor=flavor,
                                         gwchan=gwchannel, gwtrigs=gw_trig, prog=__prog__, options=opts.__dict__, version=__version__ )

        ### write documents
        logger.info('    --> writing ' + gchxml)
        reed.ligolw_utils.write_filename(gchxml_doc, gchxml, gz=gchxml.endswith('.gz'))

        logger.info('    --> writing ' + clnxml)
        reed.ligolw_utils.write_filename(clnxml_doc, clnxml, gz=clnxml.endswith('.gz'))

        logger.info('  Done: converting %s dat file into xml files.'%classifier)

        ### done with this classifier
        logger.info('Done: %s cycle.'%classifier)

    #=============================================
    # proceed to next epoch
    #=============================================
    logger.info('Done: stride %d-%d' % (t, t + stride))

    t += stride

#===================================================================================================
### loop exited, which means we must have run past options.gpsend
logger.info('End real-time evaluation')
logger.info('t + stride = %d > %d = endgps'%(t+stride, opts.endgps))

