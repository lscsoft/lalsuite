# Copyright (C) 2013 Lindy Blackburn, Reed Essick and Ruslan Vaulin
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

import time
import sys
import os
import traceback
import logging

import numpy

import ConfigParser
from optparse import OptionParser

import multiprocessing as mp

from laldetchar.idq import idq
#from laldetchar.idq import reed as idq
from laldetchar.idq import event
from laldetchar.idq import auxmvc_utils

from pylal import frutils

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table

from laldetchar import git_version

#===================================================================================================

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
    """This program performs periodic training of the classifiers used in iDQ pipeline. Some intensive training jobs are submitted to cluster using condor. """

#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('-c', '--config', default='idq.ini', type='string', help='configuration file')

parser.add_option('-k', '--lock-file', dest='lockfile', help='use custom lockfile', metavar='FILE', default=None )

parser.add_option('-s', '--gpsstart', dest="gpsstart", default=False, type='int', help='a GPS start time for the analysis. If default, gpsstart is calculated from the current time.')
parser.add_option('-e', '--gpsstop', dest="gpsstop", default=False, type='int', help='a GPS stop time for the analysis. If default, gpsstop is calculated from the current time.')
parser.add_option('-b', '--lookback', default='0', type='string', help="Number of seconds to look back and get data for training. Default is zero.\
	Can be either positive integer or 'infinity'. In the latter case, the lookback will be incremented at every stride and all data after --gps-start will be used in every training.")

parser.add_option('-l', '--log-file', default='idq_train.log', type='string', help='log file')

parser.add_option('-f','--force',default=False, action='store_true', help="forces classifiers to be trained, \
	and if errors are recovered it passes them along. Intended for use when we must have a trained classifier or else we should raise an error. \
	This may produce an infinite loop if the condor dags get lost (for mla classifiers). Use with caution.")

parser.add_option("", "--no-robot-cert", default=False, action="store_true", help="does not set up robot certificate. This can cause jobs to fail because of expired certificates")

parser.add_option("", "--ignore-science-segments", default=False, action="store_true", help="ignores science segments and trains over all available data within range")

parser.add_option("", "--sngl_chan-xml", default=False, action="store_true", help="convert single channel trg files to xml tables")

(opts, args) = parser.parse_args()

if opts.lookback != "infinity":
    opts.lookback = int(opts.lookback)

cwd = os.getcwd()

#===================================================================================================
### setup logger to record processes
logger = idq.setup_logger('idq_logger', opts.log_file, sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = idq.LogFile(logger)
sys.stderr = idq.LogFile(logger)

#===================================================================================================
### check lockfile
if opts.lockfile:
    idq.dieiflocked( opts.lockfile )

#===================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

ifo = config.get('general', 'ifo')

usertag = config.get('general', 'usertag')
if usertag:
    usertag = "_%s"%usertag

#========================
# which classifiers
#========================
### ensure we have a section for each classifier and fill out dictionary of options
classifiersD, mla, ovl = idq.config_to_classifiersD( config )

classifiers = sorted(classifiersD.keys())

if mla:
    ### reading parameters from config file needed for mla
#    auxmvc_coinc_window = config.getfloat('build_auxmvc_vectors','time-window')
#    auxmc_gw_signif_thr = config.getfloat('build_auxmvc_vectors','signif-threshold')
    auxmvc_coinc_window = config.getfloat('realtime', 'padding')
    auxmc_gw_signif_thr = config.getfloat('general', 'gw_kwsignif_thr')

auxmvc_selected_channels = config.get('general','selected-channels')
auxmvc_unsafe_channels = config.get('general','unsafe-channels')

#min_samples = config.getint('train', 'min_samples') ### minimum number of samples a training set should have
#min_svm_samples = config.getint('idq_train', 'min_svm_samples')
#min_ovl_samples = config.getint('ovl_train', 'min_num_glitches')

dag_classifiers = []
blk_classifiers = []
for classifier in classifiers:
    if classifiersD[classifier]['flavor'] in idq.train_with_dag:
        dag_classifiers.append( classifier )
    else:
        blk_classifiers.append( classifier )

#========================
# realtime
#========================
realtimedir = config.get('general', 'realtimedir')

padding = config.getfloat('realtime', 'padding')

clean_rate = config.getfloat('realtime', 'clean_rate')
clean_window = config.getfloat('realtime', 'clean_window')
clean_threshold = config.getfloat('realtime', 'clean_threshold')

#========================
# train
#========================
traindir = config.get('general', 'traindir')
if ovl: ### need snglchandir
   snglchndir = config.get('general', 'snglchndir') 

stride = config.getint('train', 'stride')
delay = config.getint('train', 'delay')

#train_script = config.get('condor', 'train')

train_cache = dict( (classifier, idq.Cachefile(idq.cache(traindir, classifier, tag='_train%s'%usertag))) for classifier in classifiers )

build_auxmvc_vectors = mla and (not os.path.exists(realtimedir)) ### if realtimedir does not exist, we cannot rely on patfiles from the realtime job
                                                                 ### we need to build our own auxmvc_vectors

max_gch_samples = config.getint("train", "max-glitch-samples")
max_cln_samples = config.getint("train", "max-clean-samples")

#========================
# data discovery
#========================
if not opts.ignore_science_segments:
    ### load settings for accessing dmt segment files
#    dmt_segments_location = config.get('get_science_segments', 'xmlurl')
    dq_name = config.get('get_science_segments', 'include')
#    dq_name = config.get('get_science_segments', 'include').split(':')[1]
    segdb_url = config.get('get_science_segments', 'segdb')
else:
    dq_name = None ### needs to be defined

### set up content handler for xml files
lsctables.use_in(ligolw.LIGOLWContentHandler)

### define the gravitational wave channel/what is a glitch
gwchannel = config.get('general', 'gwchannel')
gwthreshold = config.getfloat('general', 'gw_kwsignif_thr')

### kleineWelle config
GWkwconfigpath = config.get('data_discovery', 'GWkwconfig')
GWkwconfig = idq.loadkwconfig(GWkwconfigpath)
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')
GWkwstride = int(float(GWkwconfig['stride']))

GWkwtrgdir = "%s/%s"%(GWgdsdir, GWkwbasename)

AUXkwconfigpath = config.get('data_discovery', 'AUXkwconfig')
AUXkwconfig = idq.loadkwconfig(AUXkwconfigpath)
AUXkwbasename = AUXkwconfig['basename']
AUXgdsdir = config.get('data_discovery', 'AUXgdsdir')
AUXkwstride = int(float(AUXkwconfig['stride']))

AUXkwtrgdir = "%s/%s"%(AUXgdsdir, AUXkwbasename)

identical_trgfile = (GWgdsdir == AUXgdsdir) and (GWkwbasename == AUXkwbasename) ### used to avoid re-loading the same file for both GW and AUX triggers

#==================================================
### set up ROBOT certificates
### IF ligolw_segment_query FAILS, THIS IS A LIKELY CAUSE
if opts.no_robot_cert:
    logger.warning("Warning: running without a robot certificate. Your personal certificate may expire and this job may fail")
else:
    ### unset ligo-proxy just in case
    if os.environ.has_key("X509_USER_PROXY"):
        del os.environ['X509_USER_PROXY']

    ### get cert and key from ini file
    robot_cert = config.get('ldg_certificate', 'robot_certificate')
    robot_key = config.get('ldg_certificate', 'robot_key')

    ### set cert and key
    os.environ['X509_USER_CERT'] = robot_cert
    os.environ['X509_USER_KEY'] = robot_key

#=================================================
### data storage during training jobs

### create condorlogs directory if needed
condorlogs = config.get('train', 'condorlogs')
if not os.path.exists(condorlogs):
    os.makedirs(condorlogs)

#==================================================
### current time and boundaries

t = int(idq.nowgps())

gpsstop = opts.gpsstop
if not gpsstop: ### stop time of this analysis
    logger.info('computing gpsstop from current time')
    gpsstop = t ### We do not require boundaries to be integer multiples of stride

gpsstart = opts.gpsstart
if not gpsstart:
    logger.info('computing gpsstart from gpsstop')
    gpsstart = gpsstop - stride

#===================================================================================================
#
# LOOP
#
#===================================================================================================
logger.info('Begin: training')

n_train = 0 ### this is used to track the "original starting place" via incrementation within the persistent loop

### wait until all jobs are finished
wait = gpsstart + stride + delay - t
if wait > 0:
    logger.info('----------------------------------------------------')
    logger.info('waiting %.1f seconds to reach gpsstart+stride+delay=%d' % (wait, gpsstart+stride+delay))
    time.sleep(wait)

global_start = gpsstart 

### iterate over all ranges
while gpsstart < gpsstop:

    logger.info('----------------------------------------------------')

    wait = gpsstart + stride + delay - idq.nowgps()
    if wait > 0:
        logger.info('waiting %.1f seconds to reach gpsstart+stride+delay=%d' %(wait, gpsstart+stride+delay))
        time.sleep(wait)

    launch_gps_time = idq.nowgps()

    ### increment lookback time if it is set to infinity
    ### this forces the job to pick up all data since opts.gps_start
    if opts.lookback == 'infinity':
        lookback = gpsstart - global_start
    else:
        lookback = opts.lookback

    logger.info('processing data in [%d, %d]' % (gpsstart - lookback , gpsstart + stride))

    ### directory into which we write data
    output_dir = "%s/%d_%d/"%(traindir, gpsstart, gpsstart + stride)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #=============================================
    # veto segments (external to idq)
    #=============================================
    ### NOTHING implemented for vetosegs
    ### these are only used by ovl
    vetosegs = False

    #=============================================
    # science segments
    #=============================================
    if opts.ignore_science_segments:
        logger.info('analyzing data regardless of science segments')
        seg_file = None
        if ovl:
            ovlsegs = False ### required for OVL input

    else:
        logger.info('Begin: querrying science segments')

        try: 
           ### this returns a string
            seg_xml_file = idq.segment_query(config, gpsstart - lookback , gpsstart + stride, url=segdb_url)

            ### load xml document
            xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

            ### science segments xml filename
            seg_file = idq.segxml(output_dir, "_%s"%dq_name, gpsstart - lookback , lookback+stride)

            logger.info('writing science segments to file : '+seg_file)
            ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))

        except Exception as e:
            traceback.print_exc()
            logger.info('ERROR: segment generation failed. Skipping this training period.')

            if opts.force: ### we are require successful training or else we want errors
                logger.info(traceback.print_exc())
                raise e
            else: ### we don't care if any particular training job fails
                gpsstart += stride
                continue

        if ovl:
            if (not mla) or build_auxmvc_vectors or (not os.path.exists(realtimedir)): ### mla will use science segments, so we need to write those for ovl
                                                                         ### if realtimedir doesn't exits, we need to use queried scisegs
                try:
                    (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name) ### read in segments from xml file

                    ### write segments to ascii list
                    sciseg_path = idq.segascii(output_dir, "_%s"%dq_name, gpsstart-lookback, lookback+stride)
                    logger.info('writing science segments to file : '+sciseg_path)
                    f = open(sciseg_path, 'w')
                    for line in scisegs:
                        print >> f, line[0], line[1]
                    f.close()
                    ovlsegs = sciseg_path

                except Exception as e:
                    traceback.print_exc()
                    logger.warning('WARNING: conversion from xml to ASCII segment file failed.')
    
                    if opts.force:
                        logger.info(traceback.print_exc())
                        raise e
                    else:
                        gpsstart += stride 
                        continue 

            ### if we aren't building auxmvc vectors, we re-use pat files from realtime job
            ### this requires us to redefine the 'science-segments' as the intersection of scisegs with realtime segs
            ### we call this intersection "idq_segs"
            else: ### we're re-using pat files!
                try:
                    ### determine segments from realtime filenames
                    realtime_segs = idq.get_idq_segments(realtimedir, gpsstart - lookback, gpsstart + stride, suffix='.pat')

                    ### read in science segments
                    (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name)

                    ### take the intersection of these segments
                    idq_segs = event.andsegments([scisegs, realtime_segs])

                    ### write segment file
                    idqseg_path = idq.idqsegascii(output_dir, '_%s'%dq_name, gpsstart - lookback, lookback+stride)
                    f = open(idqseg_path, 'w')
                    for seg in idq_segs:
                        print >> f, seg[0], seg[1]
                    f.close()
 ### we may want to remove the unsafe channels, but this could be tricky and we don't want to throw away GW channels accidentally
                    ovlsegs = idqseg_path

                except Exception as e:
                    traceback.print_exc()
                    logger.info('WARNING: failed to generate iDQ segments from realtime output.')

                    if opts.force:
                        raise e
                    else:
                        gpsstart += stride
                        continue

            logger.info('Done.')

    #===============================================================================================
    # preparing auxmvc training samples
    #===============================================================================================
    if mla:
        logger.info('preparing training auxmvc samples')

        ### output file for training samples
        pat = idq.pat(output_dir, ifo, usertag, gpsstart-lookback, lookback+stride)

        if not build_auxmvc_vectors: ### we cat together pat files instead of building vectors from scratch
            ### run job that prepares training samples
            (ptas_exit_status, _) = idq.execute_prepare_training_auxmvc_samples(output_dir, realtimedir, config, gpsstart - lookback, gpsstart + stride, pat, dq_segments=seg_file, dq_segments_name=dq_name )
            os.chdir(cwd) ### go back to starting directory

        if build_auxmvc_vectors or ptas_exit_status!=0: ### we need to build vectors
            if build_auxmvc_vectors: ### no realtime directory...
                logger.warning('WARNING: building auxmvc vectors, this should be necessary only at the very first training cycle')
            else:
                logger.warning('WARNING: patfile generation failed for some reason. Attempt to build auxmvc vectors from scratch')
                if ovl and (not opts.ignore_science_segments): ### need to reset sciseg pointer!
                    ### write segments to ascii list
                    sciseg_path = idq.segascii(output_dir, "_%s"%dq_name, gpsstart-lookback, lookback+stride)
                    logger.info('writing science segments to file : '+sciseg_path)
                    f = open(sciseg_path, 'w')
                    for line in scisegs:
                        print >> f, line[0], line[1]
                    f.close()
                    ovlsegs = sciseg_path

            ### build auxmvc_vectors by hand!
            ### get sciseg info
            if not opts.ignore_science_segments:
                (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name)
            else:
                scisegs = None

            ### get triggers
            logger.info('looking for triggers')
            trigger_dict = idq.retrieve_kwtrigs(GWgdsdir, GWkwbasename, gpsstart-lookback, lookback+stride, GWkwstride, sleep=0, ntrials=1, logger=logger, segments=scisegs) ### go find GW kwtrgs
            if gwchannel not in trigger_dict:
                trigger_dict[gwchannel] = []

            if not identical_trgfile: ### add AUX triggers
                logger.info('looking for additional AUX triggers')
                aux_trgdict = idq.retrieve_kwtrigs(AUXgdsdir, AUXkwbasename, gpsstart-lookback-padding, lookback+stride+padding, AUXkwstride, sleep=0, ntrials=1, logger=logger, segments=scisegs) ### find AUX kwtrgs
                if aux_trgdict == None:
                    logger.warning('  no auxiliary triggers were found')
                    ### we do not skip, although we might want to?
                else:
                    triggger_dict.add( aux_trgdict )
                    trigger_dict.resort()

            trigger_dict.include([[gpsstart-lookback-padding, gpsstart + stride+padding]])
            trigger_dict.include([[gpsstart-lookback, gpsstart + stride]], channels=[gwchannel])
                
            ### define cleans
            logger.info('  constructing cleans')
            dirtyseg = event.vetosegs(trigger_dict[gwchannel], clean_window, clean_threshold)

            if not opts.ignore_science_segments:
                clean_gps = sorted(event.randomrate(clean_rate, scisegs)) ### generate random clean times as a poisson time series within scisegs
            else:
                clean_gps = sorted(event.randomrate(clean_rate, [[gpsstart-lookback, gpsstart + stride]])) ### generate random clean times as a poisson time series within analysis range
            clean_gps = [ l[0] for l in event.exclude( [[gps] for gps in clean_gps], dirtyseg, tcent=0)] ### keep only those gps times that are outside of dirtyseg

            ### keep only times that are within science time
            if not opts.ignore_science_segments:
                logger.info('  filtering trigger_dict through scisegs')
                trigger_dict.include(scisegs) ### already loaded into memory above here

            ### build vectors, also writes them into pat
            logger.info('  writting %s'%pat)
            idq.build_auxmvc_vectors(trigger_dict, gwchannel, auxmvc_coinc_window, auxmc_gw_signif_thr, pat, gps_start_time=gpsstart-lookback,
                                gps_end_time=gpsstart + stride,  channels=auxmvc_selected_channels, unsafe_channels=auxmvc_unsafe_channels, clean_times=clean_gps,
                                clean_window=clean_window, filter_out_unclean=False, max_glitch_samples=max_gch_samples, max_clean_samples=max_cln_samples ,
                                science_segments=None ) ### we handle scisegs in this script rather than delegating to idq.build_auxmvc_vectors, so science_segments=None is appropriate

            ptas_exit_status = 0 ### used to check for success

#            (ptas_exit_status, _) = idq.execute_build_auxmvc_vectors( config, output_dir, AUXkwtrgdir, gwchannel, pat, gpsstart - lookback, gpsstart + stride, channels=auxmvc_selected_channels, unsafe_channels=auxmvc_unsafe_channels, dq_segments=seg_file, dq_segments_name=dq_name )
#            os.chdir(cwd) ### go back to starting directory
 
        # check if process has been executed correctly
        if ptas_exit_status != 0: ### check that process executed correctly
            logger.warning('WARNING: Preparing training auxmvc samples failed')
            if opts.force:
                raise StandardError, "auxmvc samples required for successful training"
            else:
                logger.warning('WARNING: skipping re-training the MLA classifiers')
        else:
            ### figure out training set size
            ### load auxmvc vector samples
            auxmvc_samples = auxmvc_utils.ReadMVSCTriggers([pat], Classified=False)

            ### get numbers of samples
            N_clean = len(auxmvc_samples[numpy.nonzero(auxmvc_samples['i']==0)[0], :])
            N_glitch = len(auxmvc_samples[numpy.nonzero(auxmvc_samples['i']==1)[0], :])

            del auxmvc_samples

    logger.info('Done: preparing training auxmvc samples')

    #===============================================================================================
    # launch training jobs
    #===============================================================================================
    dags = {} ### dictionary that holds the dags submitted for each classifier

    #=============================================
    # training with a dag
    #=============================================
    for classifier in dag_classifiers: ### these are trained with a dag
        classD = classifiersD[classifier]
        flavor = classD['flavor']

        if flavor in idq.mla_flavors and ptas_exit_status:
            logger.warning("WARNING: mla training samples could not be built. skipping %s training"%classifier)
            if opts.force:
                raise StandardError("mla training samples could not be built.")
            continue

        min_num_cln = float(classD['min_num_cln'])
        min_num_gch = float(classD['min_num_gch'])

        if opts.force or ((N_clean >= min_num_cln) and (N_glitch >= min_num_gch)):

            logger.info('submitting %s training dag'%classifier)

            train_dir = "%s/%s/"%(output_dir, classifier)
            if not os.path.exists(train_dir):
                os.makedirs(train_dir)
            else:
                os.system("rm %s/*.dag*"%train_dir) ### remove existing condor files!
      
            ### submit training job
            miniconfig = classD['config']
            (submit_dag_exit_status, dag_file) = idq.dag_train(flavor, pat,  train_cache[classifier], miniconfig, train_dir, cwd)

            if submit_dag_exit_status:
                logger.warning("WARNING: was not able to submit %s training dag"%classifier)
                dags[dag_file] = "submit failed"
                if opts.force:
                    raise StandardError, "submission of %s training dag failed"%classifier
            else:
                dags[dag_file] = "incomplete"
                logger.info('dag file %s/%s'%(train_dir, dag_file))

        elif (N_clean >= min_num_cln):
            logger.warning("WARNING: not enough glitches in training set. skipping %s training"%classifier)
        elif (N_glitch >= min_num_gch):
            logger.warning("WARNING: not enough cleans in training set. skipping %s training"%classifier)
        else:
            logger.warning("WARNING: neither enough cleans nor enough glitches in training set. skipping %s training"%classifier)
            
    #=============================================
    # sngl_chan collection
    #=============================================
    ### OVL reads triggers from single channel summary files. these are not produced in low latency, so we build them now
    ### Currently, OVL is the only method that doesn't train with a dag, so we can put this here and not slow down the other methods' training
    ### HOWEVER, this is fragile and may break in the future.
    ### a possible work around is to define yet another group of flavors to distinguish "blk_train" from "sngchn_train"
    if ovl:
        logger.info('generating single-channel summary files')

        ### pull out only the channels we want to move   
        file_obj = open(auxmvc_selected_channels, "r")
        channels = [line.strip() for line in file_obj.readlines() if line.strip()] ### we may want to remove the unsafe channels, but this could be tricky and we don't want to throw away GW channels accidentally
        file_obj.close()

        ### pull out scisegs from ovlsegs. Will already contain any sciseg info and idq_seg info
        if ovlsegs:
            file_obj = open(ovlsegs, "r")
            ovl_segments = [ [float(l) for l in line.strip().split()] for line in file_obj.readlines() ]
            file_obj.close()
        else:
            ovl_segments = None
        new_dirs = idq.collect_sngl_chan_kw( gpsstart, gpsstart + stride, GWkwconfigpath, width=stride, source_dir=GWkwtrgdir, output_dir=snglchndir, chans=channels, scisegs=ovl_segments )
        if not identical_trgfile:
            new_dirs += idq.collect_sngl_chan_kw( gpsstart, gpsstart+stride, AUXkwconfigpath, width=stride, source_dir=AUXkwtrgdir, output_dir=snglchndir, chans=channels, scisegs=ovl_segments )

    #=============================================
    # training on submit node
    #=============================================
    blk_procs = [] ### for parallelization
    for classifier in blk_classifiers: ### these are trained via routines on the head node
                                       ### right now, this is just OVL, so we're safe putting sngl_chan routine ahead of this 
        classD = classifiersD[classifier]
        flavor = classD['flavor']

        train_dir = "%s/%s/"%(output_dir, classifier)
        if not os.path.exists(train_dir):
            os.makedirs(train_dir)

        logger.info('launching %s training job'%classifier)

        min_num_cln = float(classD['min_num_cln'])
        min_num_gch = float(classD['min_num_gch'])

        ### parallelize through multiprocessing module!
        conn1, conn2 = mp.Pipe()

        proc = mp.Process(target=idq.blk_train, args=(flavor, config, classifiersD[classifier], gpsstart-lookback, gpsstart+stride, ovlsegs, vetosegs, train_dir, cwd, opts.force, train_cache[classifier], min_num_gch, min_num_cln, padding, conn2) )
	proc.start()
        conn2.close()
        blk_procs.append( (proc, classifier, conn1) )

        ### run in series
#        exit_status, _ = idq.blk_train(flavor, config, classifiersD[classifier], gpsstart-lookback, gpsstart+stride, ovlsegs=ovlsegs, vetosegs=vetosegs, train_dir=train_dir, cwd=cwd, force=opts.force, cache=train_cache[classifier], min_num_gch=min_num_gch, min_num_cln=min_num_cln, padding=padding)

    ### loop over processes and wait...
    while blk_procs:
        proc, classifier, conn = blk_procs.pop(0)
        proc.join()

        logger.info('%s training finished'%classifier)

        exit_status, _ = conn.recv()
        conn.close()
        if exit_status:
            logger.info('WARNING: number glitches in the training set is less than %s'%classifiersD[classifier]['min_num_gch'])
            logger.info('WARNING: Not enough glitches in the training set. Skip %s training.'%classifier)

    logger.info('Done.')

    #===============================================================================================
    # jobs completion checkpoint (for dags)
    #===============================================================================================
    if dags: ### dags is not empty

        logger.info('Begin: checking on dags')

        list_of_dags = sorted(dags.keys())

        ### loop until either all dags are complete or gpsstop+stride time is reached
        incompleted_dags = []

        ### if opts.force, we wait for all jobs to finish
        ### otherwise, we wait for some amount of time before proceeding
        wait = 1 ### amount of time to wait between checking-dags epochs
        message = True 
        while opts.force or (idq.nowgps() < gpsstop + stride) or (idq.nowgps() < launch_gps_time + stride):

            for dag in list_of_dags: ### check dag status
                dag_status = idq.get_condor_dag_status(dag)
                dags[dag] = dag_status

                if dag_status == 0:
                    logger.info('%s completed' % dag)
                else:
                    incompleted_dags.append(dag)

            list_of_dags = incompleted_dags
            incompleted_dags = []

            if 'incomplete' in dags.values(): ### check if any of the dags are still incomplete
                if message:
                    logger.info('WARNING: Some dags have not been completed. Waiting for dags to complete')
                    message = False
                time.sleep(wait) 
            else: ### no incomplete dags, break the loop
                break

        ### print warnings for dags that failed or are incomplete
        for dag in list_of_dags:
            if dags[dag] != 0:
                if dags[dag] == 'incomplete':
                    logger.warning('WARNING: %s was  not complete' % dag)
                else:
                    logger.warning('WARNING: %s failed with exit status %s'% (dag, str(dags[dag])))
                    if opts.force:
                       raise StandardError, "%s failed with exit status %s"%(dag, str(dags[dag]))

    #===============================================================================================
    # conversion of snglchn files to xml
    #===============================================================================================
    ### we collect the single channel trg files into xml files
    ### this shouldn't be necessary for OVL training and i doubt these files
    ### will be used by anyone else...
    ### FIXME: WE MAY WANT TO REMOVE THIS STEP
    if opts.sngl_chan_xml:
        logger.info('launching conversion from .trg to .xml files')
        for dir in new_dirs:
#            trg_to_xml_exit_code = idq.submit_command([config.get('condor', 'convertkwtosb'), dir], process_name='convertkwtosb', dir=dir)
            proc = idq.fork([config.get('convertkwtosb','executable'), dir], cwd=dir)
            proc.wait()
            trg_to_xml_exit_code = proc.returncode
            os.chdir(cwd)
            if trg_to_xml_exit_code != 0:
                logger.info('WARNING: Conversion from single-channel KW trig files to xml failed in '+ dir)

        logger.info('Done.')

    #===============================================================================================

    ### continue onto the next stride
    gpsstart += stride

#===================================================================================================
if opts.lockfile:
    idq.release(opts.lockfile) ### unlock lockfile
