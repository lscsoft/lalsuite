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

import time
import logging
import tempfile
import ConfigParser

import numpy
import math

from laldetchar.idq import idq
#from laldetchar.idq import reed as idq

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

parser.add_option('-s', '--gpsstart', dest='startgps',
    type='int', help='first stride begins at or after GPS time',
    metavar='GPS', default=idq.nowgps(), )
parser.add_option('-e', '--gpsstop', dest='endgps',
    type='int', help='last stride ends at or before GPS time',
    metavar='GPS', default=1e10 )

parser.add_option('-k', '--lock-file', dest='lockfile',
    help='use custom lockfile', metavar='FILE', default='',
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
parser.add_option('', '--lock-training', default=False, action="store_true", help="forces training jobs to use a lock file. Useful when preventing jobs from stacking up.")

parser.add_option('', '--ignore-science-segments-summary', 
    default=False, action="store_true", 
    help="forces summary jobs to ignore science segments as well. \
    This should NOT be used to generate actual DQ information, but may be useful when debugging the pipeline."
    )

parser.add_option("", "--dont-cluster-summary", default=False, action="store_true")
parser.add_option('', '--lock-summary', default=False, action="store_true", help="forces summary jobs to use a lock file. Useful when preventing jobs from stacking up.")

parser.add_option('', '--ignore-science-segments-calibration',
    default=False, action="store_true",
    help="forces calibration jobs to ignore science segments as well. \
    this should NOT be used to generate actual DQ information, but may be useful when debugging the pipeline."
    )

parser.add_option("", "--dont-cluster-calibration", default=False, action="store_true")
parser.add_option('', '--lock-calibration', default=False, action="store_true", help="forces calibration jobs to use a lock file. Useful when preventing jobs from stacking up.")

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
logger = idq.setup_logger('idq_logger', opts.log_file, sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = idq.LogFile(logger)
sys.stderr = idq.LogFile(logger)

#===================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config_file)

mainidqdir = config.get('general', 'idqdir') ### get the main directory where idq pipeline is going to be running.
if not opts.lockfile:
    opts.lockfile = "%s/.idq_realtime.lock"%mainidqdir

idq.dieiflocked( opts.lockfile ) ### prevent multiple copies from running

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
### ensure we have a section for each classifier and fill out dictionary of options
classifiersD, mla, ovl = idq.config_to_classifiersD( config )
classifiers = sorted(classifiersD.keys())

### get combiners information and DO NOT add this to classifiersD, we treat them separately!
combinersD, referenced_classifiers = idq.config_to_combinersD( config )
combiners = sorted(combinersD.keys())

### compute channel names to be stored in frames
channameD = dict( (name, {'rank':idq.channame(ifo, name, "%s_rank"%usertag), 
                          'fap':idq.channame(ifo, name, "%s_fap"%usertag),
                          'fapUL':idq.channame(ifo, name, "%s_fapUL"%usertag)}) for name in classifiers+combiners )

if mla:
    ### reading parameters from config file needed for mla
#    auxmvc_coinc_window = config.getfloat('build_auxmvc_vectors','time-window')
#    auxmc_gw_signif_thr = config.getfloat('build_auxmvc_vectors','signif-threshold')
    auxmvc_coinc_window = config.getfloat('realtime', 'padding')
    auxmc_gw_signif_thr = config.getfloat('general', 'gw_kwsignif_thr')
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
GWkwconfig = idq.loadkwconfig(config.get('data_discovery', 'GWkwconfig'))
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')

AUXkwconfig = idq.loadkwconfig(config.get('data_discovery', 'AUXkwconfig'))
AUXkwbasename = AUXkwconfig['basename']
AUXgdsdir = config.get('data_discovery', 'AUXgdsdir')

identical_trgfile = (GWgdsdir == AUXgdsdir) and (GWkwbasename == AUXkwbasename) ### used to avoid re-loading the same file for both GW and AUX triggers

#========================
# realtime job
#========================
ts_fs = config.getfloat('realtime','sampling_rate') ### sampling frequency for time-series files
ts_dt = 1.0/ts_fs

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
summary_trending = config.get('summary', 'trending').split()
for ind, trend in enumerate(summary_trending):
    if trend != 'infinity':
        summary_trending[ind] = int(trend)

summary_FAPthrs = [float(l) for l in config.get('summary','FAP').split()]

summary_script = config.get('summary', 'executable')

#summary_cache = dict( (classifier, idq.Cachefile(idq.cache(summarydir, classifier, tag='_summary%s'%usertag))) for classifier in classifiers ) ### not needed?
#for cache in summary_cache.values():
#    cache.time = 0

summary_log = "%s/summary%s.log"%(mainidqdir, usertag)
summary_out = "%s/summary%s.out"%(mainidqdir, usertag)
summary_err = "%s/summary%s.err"%(mainidqdir, usertag)

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

train_script = config.get('train', 'executable')

train_cache = dict( (classifier, idq.Cachefile(idq.cache(traindir, classifier, tag='_train%s'%usertag))) for classifier in classifiers )
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

calibration_mode = config.get('calibration','mode')

### cache for uroc files
calibration_cache = dict( (classifier, idq.Cachefile(idq.cache(calibrationdir, classifier, tag='_calibration%s'%usertag))) for classifier in classifiers+combiners )
for cache in calibration_cache.values():
    cache.time = 0

### cache for kde estimates
kde_cache = dict( (classifier, idq.Cachefile(idq.cache(calibrationdir, classifier, tag='_calibration-kde%s'%usertag))) for classifier in classifiers )
for cache in kde_cache.values():
    cache.time = 0

kdeD = dict( (classifier,{}) for classifier in classifiers ) ### used to store kde information read in from files within kde_cache

calibration_script = config.get('calibration', 'executable')

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
for classifier in classifiers:
    if not initial_training:
        initial_training =  train_cache[classifier].is_empty()

initial_calibration = False
for classifier in classifiers+combiners:
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
    train_command = "%s --config %s -l %s -s %d -e %d --lookback %s --force"%(train_script, opts.config_file, train_log, train_start_time, t, initial_training_lookback)

    if opts.ignore_science_segments_training:
        train_command += " --ignore-science-segments" 
    if opts.no_robot_cert:
        train_command += " --no-robot-cert"
    if opts.lock_training:
        train_command += " --lock-file %s/.idq_train.lock"%mainidqdir

    logger.info('Submiting training job with the following options:')
    logger.info(train_command)

    ### launch training script, wait for it to finish
    train_out_file = open(train_out, 'a')
    train_err_file = open(train_err, 'a')
    proc = idq.fork( train_command.split() , stdout=train_out_file, stderr=train_err_file, cwd=cwd )
    train_out_file.close()
    train_err_file.close()
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
    calibration_command = "%s --config %s -l %s -s %d -e %d --lookback %s --mode %s --force"%(calibration_script, opts.config_file, calibration_log, calibration_start_time, t, initial_calibration_lookback, calibration_mode)
    ### we skip validation (no FAPthr arguments) for the initial job

    if opts.ignore_science_segments_calibration:
        calibration_command += " --ignore-science-segments" 
    if opts.no_robot_cert:
        calibration_command += " --no-robot-cert" 
    if opts.lock_calibration:
        calibration_command += " --lock-file %s/.idq_calibration.lock"%mainidqdir

    logger.info('Submiting calibration job with the following options:')
    logger.info(calibration_command) ### block

    ### launch calibration script, wait for it to finish
    calibration_out_file = open(calibration_out, 'a')
    calibration_err_file = open(calibration_err, 'a')
    proc = idq.fork( calibration_command.split() , stdout=calibration_out_file, stderr=calibration_err_file, cwd=cwd )
    calibration_out_file.close()
    calibration_err_file.close()
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
train_template = "%s --config %s -l %s -s %d -e %d --lookback %d"
if opts.ignore_science_segments_training:
    train_template += " --ignore-science-segments"
if opts.lock_training:
    train_template += " --lock-file %s/.idq_train.lock"%mainidqdir

summary_template = "%s --config %s -l %s -s %d -e %d --lookback %d %s"
for fap in summary_FAPthrs:
    summary_template += " --FAPthr %.6e"%fap
if opts.ignore_science_segments_summary:
    summary_template += " --ignore-science-segments"
if opts.dont_cluster_summary:
    summary_template += " --dont-cluster"
if opts.lock_summary:
    summary_template += " --lock-file %s/.idq_summary.lock"%mainidqdir

calibration_template = "%s --config %s -l %s -s %d -e %d --lookback %d --mode %s" 
for fap in calibration_FAPthrs:
    calibration_template += " --FAPthr %.6e"%fap
if opts.ignore_science_segments_calibration:
    calibration_template += " --ignore-science-segments"
if opts.dont_cluster_calibration:
    calibration_template += " --dont-cluster"
if opts.lock_calibration:
    calibration_template += " --lock-file %s/.idq_calibration.lock"%mainidqdir

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

        train_out_file = open(train_out, 'a')
        train_err_file = open(train_err, 'a')
        if opts.offline: ### running in "offline mode"
            logger.info('Running in OFFLINE mode, will wait until train is complete before proceeding with evaluation.')
            proc = idq.fork( train_command.split() , stdout=train_out_file, stderr=train_err_file , cwd=cwd )
            train_out_file.close()
            train_err_file.close()
            proc.wait() ### blocks!
        else:
            ### running in realtime, don't wait for training job to complete so we can keep the latency low
            idq.double_fork(train_command.split(), stdout=train_out_file, stderr=train_err_file, cwd=cwd )
            train_out_file.close()
            train_err_file.close()

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

        calibration_command = calibration_template%(calibration_script, opts.config_file, calibration_log, calibration_start_time, calibration_stop_time, lookback, calibration_mode)

        logger.info('Submiting calibration script with the following options:')
        logger.info(calibration_command)

        calibration_out_file = open(calibration_out, 'a')
        calibration_err_file = open(calibration_err, 'a')
        if opts.offline:
            logger.info("Running in OFFLINE mode, will wait until calibration is complete before proceeding with evaluation")
            proc = idq.fork(calibration_command.split(), stdout=calibration_out_file, stderr=calibration_err_file, cwd=cwd)
            calibration_out_file.close()
            calibration_err_file.close()
            proc.wait() # block!
        else:
            idq.double_fork(calibration_command.split(), stdout=calibration_out_file, stderr=calibration_err_file, cwd=cwd)
            calibration_out_file.close()
            calibration_err_file.close()

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

        if summary_lookback == "infinity":
            lookback = t-global_start
        else:
            lookback = summary_lookback

        trend_str = ""
        for trend in summary_trending:
            if trend == "infinity":
                trend_str += " --trending %d"%(t-global_start)
            else:
                trend_str += " --trending %d"%(trend)
        
        summary_stop_time = (t / summary_stride) * summary_stride
        summary_start_time = summary_stop_time - summary_stride

        logger.info('Begin: launching summary script for period: %s - %s'%(summary_start_time,summary_stop_time))

        summary_command = summary_template%(summary_script, opts.config_file, summary_log, summary_start_time, summary_stop_time, lookback, trend_str)

        logger.info('Submiting summary script with the following options:')
        logger.info(summary_command)

        summary_out_file = open(summary_out, 'a')
        summary_err_file = open(summary_err, 'a')
        if opts.offline:
            logger.info("Running in OFFLINE mode, will wait until summary is complete before proceeding with evaluation")
            proc = idq.fork(summary_command.split(), stdout=summary_out_file, stderr=summary_err_file, cwd=cwd)
            summary_out_file.close()
            summary_err_file.close()
            proc.wait() # block!
        else:                                              
            idq.double_fork(summary_command.split(), stdout=summary_out_file, stderr=summary_err_file, cwd=cwd)
            summary_out_file.close()
            summary_err_file.close()

        logger.info('Done: launching summary script for period: %s - %s'%(summary_start_time,summary_stop_time))
        logger.info('****************************************************')

    #===============================================================================================
    # REALTIME EVALUATION
    #===============================================================================================
    logger.info('----------------------------------------------------')

    logger.info('Begin: stride %d-%d' % (t, t + stride))

    ### make sure end time of analysis stride is not ahead of now-delay
    wait = t + stride + delay - int(idq.nowgps())
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
        good, covered = idq.retrieve_scisegs(dmt_segments_location, dmtdq_name, t, stride, pad=0, sleep=delay, nretry=1, logger=logger)

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

    trgdict = idq.retrieve_kwtrig(GWgdsdir, GWkwbasename, t, stride, sleep=(t+stride-idq.nowgps()) + delay, ntrials=2, logger=logger)
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
        aux_trgdict = idq.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t, stride, sleep=(t+stride-idq.nowgps()) + delay, ntrials=2, logger=logger)
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
        prev_trgdict = idq.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t-stride, stride, sleep=0, ntrials=1, logger=logger)
        if prev_trgdict == None:
            logger.warning('  missing previous stride KW triggers')
            ### we do not skip this stride when previous kw stride is required and missing!
        else:
            trgdict.add( prev_trgdict ) ### we add the previous triggers into the current trgdict

    if maxtime > t + stride - padding: ### triggers from the following stride must be included
        logger.info('  next stride KW triggers are needed')
        next_trgdict = idq.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t+stride, stride, sleep=(t+stride+stride-idq.nowgps()) + delay, ntrials=2, logger=logger)
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

        pat = idq.pat(this_output_dir, ifo, usertag, t, stride)

        # generating auxmvc vector samples. result is saved into pat file
        # FIXME: depending how padding is done we should adjust behavior of build_auxmvc_vectors
        # Currently it keeps gw trigger from [t, t + stride] and uses time_window to pad this segment for auxiliary triggers
		# we do not filter out unclean beacuse it is already done when clean_gps times are formed
        auxmvc_vectors = idq.build_auxmvc_vectors(trgdict, gwchannel, auxmvc_coinc_window, auxmc_gw_signif_thr, pat, gps_start_time=t,
                                gps_end_time=t + stride,  channels=auxmvc_selected_channels, unsafe_channels=auxmvc_unsafe_channels, clean_times=clean_gps,
                                clean_window=clean_window, filter_out_unclean=False )

        logger.info('Done: building auxmvc feature vectors.')

    #=============================================
    # predictions
    #=============================================
    dats = {}
    rnks = {}
    for classifier in classifiers:
        flavor = classifiersD[classifier]['flavor']

        logger.info('Begin: %s cycle -> flavor=%s'%(classifier, flavor))

        ### check for new training data
        cache = train_cache[classifier]
        if cache.was_modified():
            lines = cache.tail(nlines=idq.traincache_nlines[flavor])
            trained_range = idq.extract_trained_range(lines, flavor)
            
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
            calib_range = idq.extract_calib_range(uroc)
            r, c, g, _, _ = idq.file_to_rcg(uroc)
            effmap, fapmap = idq.rcg_to_EffFAPmap(r, c, g)
            
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
        dat = idq.dat(this_output_dir, classifier, ifo, trained_range, usertag, t, stride)

        gchxml = idq.xml(this_output_dir, classifier, ifo, trained_range, calib_range, gchtag, t, stride)
        clnxml = idq.xml(this_output_dir, classifier, ifo, trained_range, calib_range, clntag, t, stride)

        rnknp = idq.timeseries(this_output_dir, classifier, ifo, trained_range, calib_range, rnktag, t, stride)
        rnkfr = idq.timeseriesgwf(this_output_dir, classifier, ifo, trained_range, calib_range, rnktag, t, stride)

        fapnp = idq.timeseries(this_output_dir, classifier, ifo, trained_range, calib_range, faptag, t, stride)
        fapfr = idq.timeseriesgwf(this_output_dir, classifier, ifo, trained_range, calib_range, faptag, t, stride)

        ### store relevant filenames for combiners
        dats[classifier] = dat
        rnks[classifier] = rnknp

        ### perform evalutation
        logger.info('  Begin %s evaluation -> %s'%(classifier, dat) )
        miniconfig = classifiersD[classifier]['config']
        if classifiersD[classifier]['mla']:
            returncode = idq.evaluate(flavor, lines, dat, miniconfig, gps_start_time=t, gps_end_time=t+stride, dir=this_output_dir, trgdict=pat, auxmvc_vectors=auxmvc_vectors)
        else:
            returncode = idq.evaluate(flavor, lines, dat, miniconfig, gps_start_time=t, gps_end_time=t+stride, dir=this_output_dir, trgdict=trgdict, samples=samples, samples_header=samples_header)

        ### check if process has been execited correctly
        if returncode:
            logger.warning('  WARNING: %s predict failed'%classifier)
            logger.info('  skipping %s timeseries'%classifier)
            continue
        logger.info('  Done: %s evaluation '%classifier)

        ### create timeseries from datfile and save as gzip numpy array, there are boundary effects here
        logger.info('  Begin: generating %s rank timeseries -> %s'%(classifier, rnknp))
        ts = idq.datfile2timeseries(flavor, dat, lines, trgdict=trgdict, start=t, stride=stride, window=0.1, fs=ts_fs)

        rnknp_file = event.gzopen(rnknp, 'w')
        numpy.save(rnknp_file, ts)
        rnknp_file.close()

        logger.info('    writing %s'%rnkfr)
        idq.timeseries2frame( rnkfr, {channameD[classifier]['rank']:ts}, t, ts_dt )

        logger.info('  Done: generating %s rank timeseries'%classifier)

        ### create FAP timeseries
        logger.info('  Begin: generating %s FAP time-series -> %s'%(classifier, fapnp))

        fap_ts = fapmap(ts) ### MLE estimate of FAP (interpolated)
        fapUL_ts = fapmap.ul(ts, conf=0.99) ### upper limit with FAPconf (around interpolated points)
                                            ### we hard-code this intentionally to make it difficult to change/mess up
        fapnp_file = event.gzopen(fapnp, 'w')
        numpy.save(fapnp_file, (fap_ts, fapUL_ts) )
        fapnp_file.close()

        logger.info('    writing %s'%fapfr)
        idq.timeseries2frame( fapfr, {channameD[classifier]['fap']:fap_ts, channameD[classifier]['fapUL']:fapUL_ts}, t, ts_dt )

        logger.info('  Done: generating %s FAP time-series'%classifier)

        ### convert datfiles to xml tables
        logger.info('  Begin: converting %s dat file into xml files'%classifier)
        logger.info('    converting %s to xml tables' % dat)

        ### read dafile -> xml docs
        (gchxml_doc, clnxml_doc) = idq.datfile2xmldocs(dat, ifo, fapmap, Lmap=False, Effmap=effmap, flavor=flavor,
                                         gwchan=gwchannel, gwtrigs=gw_trig, prog=__prog__, options=opts.__dict__, version=__version__ )

        ### write documents
        logger.info('    --> writing ' + gchxml)
        idq.ligolw_utils.write_filename(gchxml_doc, gchxml, gz=gchxml.endswith('.gz'))

        logger.info('    --> writing ' + clnxml)
        idq.ligolw_utils.write_filename(clnxml_doc, clnxml, gz=clnxml.endswith('.gz'))

        logger.info('  Done: converting %s dat file into xml files.'%classifier)

        ### done with this classifier
        logger.info('Done: %s cycle.'%classifier)

    #=============================================
    # cycle through combiners
    #=============================================
    if combiners:
        ### load in dat files only once!
        outputs = {}
        ts = {}
        for classifier in referenced_classifiers:
            dat = dats[classifier]
            logger.info('reading in data from %s'%dat)
            ### we sort according to GPS so everything is ordered in the same way!
            outputs[classifier] = idq.sort_output( idq.slim_load_datfile( dat, skip_lines=0, columns=samples_header+['rank'] ), 'GPS' )

            rnk = rnks[classifier]
            logger.info('reading in data from %s'%rnk)
            rnk_file = event.gzopen(rnk, 'r')
            ts[classifier] = numpy.load(rnk_file)
            rnk_file.close()

            ### check kde_caches
            cache = kde_cache[classifier]
            if cache.was_modified():
                kde_cln_name, kde_gch_name = cache.tail(nlines=2)

                kde_range = idq.extract_kde_range( kde_cln_name )

                kde_cln_file = event.gzopen(kde_cln_name, 'r')
                kde, kde_cln = numpy.load(kde_cln_file)
                kde_cln_file.close()

                kde_gch_file = event.gzopen(kde_gch_name, 'r')
                _  , kde_gch = numpy.load(kde_gch_file)
                kde_gch_file.close()

                ### store kdes
                kdeD[classifier]['kde'] = kde

                kdeD[classifier]['kde_cln'] = kde_cln
                kdeD[classifier]['ckde_cln'] = idq.kde_to_ckde( kde_cln )
                kdeD[classifier]['kde_gch'] = kde_gch
                kdeD[classifier]['ckde_gch'] = idq.kde_to_ckde( kde_gch )

                kdeD[classifier]['kde_cln_name'] = kde_cln_name
                kdeD[classifier]['kde_gch_name'] = kde_gch_name

                kdeD[classifier]['kde_range'] = kde_range

                logger.info('using newest %s kde data : %s, %s'%(classifier, kde_cln_name, kde_gch_name))
                cache.set_time()

        n_samples = len(outputs[classifier]['rank']) ### assumes all outputs have the same length, which they should

        ### iterate over each combiner
        for combiner in combiners:
            logger.info('Begin: %s cycle'%combiner)

            flavor = combinersD[combiner]['flavor'] ### should always be None, but we pull the reference nonetheless

            these_classifiers = combinersD[combiner]['classifiers'] ### the classifiers this combiner combines

            ### check for new calibration data
            cache = calibration_cache[combiner]
            if cache.was_modified():
                uroc = cache.tail(nlines=1)[0]
                calib_range = idq.extract_calib_range(uroc)
                r, c, g, _, _ = idq.file_to_rcg(uroc)
                effmap, fapmap = idq.rcg_to_EffFAPmap(r, c, g)

                combinersD[combiner]['uroc'] = uroc
                combinersD[combiner]['calib_range'] = calib_range
                combinersD[combiner]['effmap'] = effmap
                combinersD[combiner]['fapmap'] = fapmap

                cache.set_time()
                logger.info('using newest %s calibration data : %s'%(combiner, uroc))
            else:
                calib_range = combinersD[combiner]['calib_range']
                effmap = combinersD[combiner]['effmap']
                fapmap = combinersD[combiner]['fapmap']

            ### compute filenames
            trained_range = "__".join( ["%s-%s"%(classifier, kdeD[classifier]['kde_range']) for classifier in these_classifiers] ) ### string together all kde ranges
             
            dat = idq.dat(this_output_dir, combiner, ifo, trained_range, usertag, t, stride)

            gchxml = idq.xml(this_output_dir, combiner, ifo, trained_range, calib_range, gchtag, t, stride)
            clnxml = idq.xml(this_output_dir, combiner, ifo, trained_range, calib_range, clntag, t, stride)

            rnknp = idq.timeseries(this_output_dir, combiner, ifo, trained_range, calib_range, rnktag, t, stride)
            rnkfr = idq.timeseriesgwf(this_output_dir, combiner, ifo, trained_range, calib_range, rnktag, t, stride)

            fapnp = idq.timeseries(this_output_dir, combiner, ifo, trained_range, calib_range, faptag, t, stride)
            fapfr = idq.timeseriesgwf(this_output_dir, combiner, ifo, trained_range, calib_range, faptag, t, stride)

            #===========================================================================================
            #
            # perform combiners evaluation here!
            #
            #===========================================================================================

            c_pofg = combinersD[combiner]['joint_p(g)'].endswith("cdf") ### check to see whether we use pdf or cdf
            c_pofc = combinersD[combiner]['joint_p(c)'].endswith("cdf")

            combinersD[combiner]['columns']
            output = dict( [(c,[]) for c in samples_header+combinersD[combiner]['columns']] ) 

            pofg_joint = [] ### store the p(g) data from dat files
            pofc_joint = []

            ts_pofg_joint = [] ### store the p(g) data from rank files
            ts_pofc_joint = []

            for classifier in these_classifiers:
                logger.info('  extracting data from %s'%classifier)

                ### pull out kde data for this classifier
                kde = kdeD[classifier]['kde']
                if c_pofc:
                    kde_cln = kdeD[classifier]['ckde_cln']
                else:   
                    kde_cln = kdeD[classifier]['kde_cln']
                if c_pofg:
                    kde_gch = kdeD[classifier]['ckde_gch']
                else:
                    kde_gch = kdeD[classifier]['kde_gch']

                #========
                # dat files
                #========
                ### pull out ranks
                ranks = [outputs[classifier]['rank'][i] for i in xrange(n_samples)]
                
                ### map ranks into p(g), p(c)
                pofg = numpy.interp( ranks, kde, kde_gch )
                pofc = numpy.interp( ranks, kde, kde_cln )

                ### fill out output dictionary
                output["%s_rank"%classifier] = ranks
                output["%s_p(g)"%classifier] = pofg
                output["%s_p(c)"%classifier] = pofc

                ### add to joint probabilities
                pofg_joint.append( pofg )
                pofc_joint.append( pofc )

                ### fill in samples_header if not already there
                for c in samples_header:
                    if not output[c]: ### only do this once!
                        output[c] = outputs[classifier][c]

                #========
                # rnk files
                #========
                ranks = ts[classifier]
                pofg = numpy.interp( ranks, kde, kde_gch )
                pofc = numpy.interp( ranks, kde, kde_cln )

                ts_pofg_joint.append( pofg )
                ts_pofc_joint.append( pofc )

            ### compute joint probabilities
            logger.info('  combining classifier probability distributions')
            if combinersD[combiner]['joint_p(g)'][:3] == "max":
                pofg_joint = numpy.max( pofg_joint, axis=0 )
                ts_pofg_joint = numpy.max( ts_pofg_joint, axis=0 )
            elif combinersD[combiner]['joint_p(g)'][:4] == "prod":
                pofg_joint = numpy.prod( pofg_joint, axis=0 )
                ts_pofg_joint = numpy.prod( ts_pofg_joint, axis=0 )
            else:
                raise ValueError("combiner=%s joint_p(g)=%s not understood"%(combiner, combinersD[combiner]['joint_p(g)']))

            if combinersD[combiner]['joint_p(c)'][:3] == "max":
                pofc_joint = numpy.max( pofc_joint, axis=0 )
                ts_pofc_joint = numpy.max( ts_pofc_joint, axis=0 )
            elif combinersD[combiner]['joint_p(c)'][:4] == "prod":
                pofc_joint = numpy.prod( pofc_joint, axis=0 )
                ts_pofc_joint = numpy.prod( ts_pofc_joint, axis=0 )
            else:
                raise ValueError("combiner=%s joint_p(c)=%s not understood"%(combiner, combinersD[combiner]['joint_p(c)']))

            ### compute combined statistics
            logger.info('  computing approximation to joint likelihood ratio')
            L_joint = pofg_joint / pofc_joint ### compute likelihood
            r_joint = L_joint / ( 1 + L_joint ) ### compute rank

            ### put them into output
            output['rank'] = r_joint
            output['Likelihood'] = L_joint

            ### write datfile
            logger.info('  writing %s'%dat)
            idq.output_to_datfile( output, dat )
           
            ### create timeseries from datfile and save as gzip numpy array, there are boundary effects here
            logger.info('  writing %s'%(rnknp))
            ts_L_joint = ts_pofg_joint / ts_pofc_joint
            ts_r_joint = ts_L_joint / (1 + ts_L_joint)

            rnknp_file = event.gzopen(rnknp, 'w')
            numpy.save(rnknp_file, ts_r_joint)
            rnknp_file.close()

            logger.info('    writing %s'%rnkfr)
            idq.timeseries2frame( rnkfr, {channameD[combiner]['rank']:ts_r_joint}, t, ts_dt )

            ### create FAP timeseries
            logger.info('  Begin: generating %s FAP time-series -> %s'%(combiner, fapnp))

            fap_ts = fapmap(ts_r_joint) ### MLE estimate of FAP (interpolated)
            fapUL_ts = fapmap.ul(ts_r_joint, conf=0.99) ### upper limit with FAPconf (around interpolated points)
                                                        ### we hard-code this intentionally to make it difficult to change/mess up
            fapnp_file = event.gzopen(fapnp, 'w')
            numpy.save(fapnp_file, (fap_ts, fapUL_ts) )
            fapnp_file.close()

            logger.info('    writing %s'%fapfr)
            idq.timeseries2frame( fapfr, {channameD[combiner]['fap']:fap_ts, channameD[combiner]['fapUL']:fapUL_ts}, t, ts_dt )

            logger.info('  Done: generating %s FAP time-series'%combiner)

            ### convert datfiles to xml tables
            logger.info('  Begin: converting %s dat file into xml files'%combiner)
            logger.info('    converting %s to xml tables' % dat)

            ### read dafile -> xml docs
            (gchxml_doc, clnxml_doc) = idq.datfile2xmldocs(dat, ifo, fapmap, Lmap=False, Effmap=effmap, flavor=flavor,
                                         gwchan=gwchannel, gwtrigs=gw_trig, prog=__prog__, options=opts.__dict__, version=__version__ )

            ### write documents
            logger.info('    --> writing ' + gchxml)
            idq.ligolw_utils.write_filename(gchxml_doc, gchxml, gz=gchxml.endswith('.gz'))

            logger.info('    --> writing ' + clnxml)
            idq.ligolw_utils.write_filename(clnxml_doc, clnxml, gz=clnxml.endswith('.gz'))

            logger.info('  Done: converting %s dat file into xml files.'%combiner)

            ### done with this classifier
            logger.info('Done: %s cycle.'%combiner)

    #=============================================
    # proceed to next epoch
    #=============================================
    logger.info('Done: stride %d-%d' % (t, t + stride))

    t += stride

#===================================================================================================
### loop exited, which means we must have run past options.gpsend
logger.info('End real-time evaluation')
logger.info('t + stride = %d > %d = endgps'%(t+stride, opts.endgps))

### unlock file
idq.release(lockfile)
