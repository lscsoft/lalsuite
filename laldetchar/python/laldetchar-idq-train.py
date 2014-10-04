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
from laldetchar.idq import idq
from laldetchar.idq import event
import subprocess
import ConfigParser
from optparse import *
import traceback
import logging
from glue.ligolw import table, lsctables, utils
from laldetchar.idq import auxmvc_utils
from pylal import frutils
import numpy

# prevent multiple copies from running
# idq.dieiflocked('/home/detchar/idq/.idq_train.lock')

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

parser.add_option('-s', '--gps-start', default=False, type='int', help='a GPS start time for the analysis. If default, gpsstart is calculated from the current time.')
parser.add_option('-e', '--gps-stop', default=False, type='int', help='a GPS stop time for the analysis. If default, gpsstop is calculated from the current time.')
parser.add_option('-b', '--lookback', default='0', type='string', help="Number of strides to look back and get data for training. Default is zero.\
	Can be either positive integer or 'infinity'. In the latter case, the lookback will be incremented at every stride and all data after --gps-start will be used in every training.")

parser.add_option('-l', '--log-file', default='idq_train.log', type='string', help='log file')

parser.add_option('-f','--force',default=False, action='store_true', help="forces classifiers to be trained, \
	and if errors are recovered it passes them along. Intended for use when we must have a trained classifier or else we should raise an error. \
	This may produce an infinite loop if the condor dags get lost (for mla classifiers). Use with caution.")

parser.add_option("", "--no-robot-cert", default=False, action="store_true", help="does not set up robot certificate. This can cause jobs to fail because of expired certificates")

parser.add_option("", "--ignore-science-segments", default=False, action="store_true", help="ignores science segments and trains over all available data within range")

parser.add_option("", "--no-sngl_chan-xml", default=False, action="store_true", help="skip conversion of single channel trg files to xml tables. This is just here for speed when testing")

(opts, args) = parser.parse_args()

#===================================================================================================
### setup logger to record processes
logger = logging.getLogger('idq_logger')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s %(message)s')
hdlr1 = logging.StreamHandler(sys.stdout)
hdlr1.setFormatter(formatter)
hdlr1.setLevel(logging.INFO)
logger.addHandler(hdlr1)
hdlr2 = logging.FileHandler(opts.log_file)
hdlr2.setFormatter(formatter)
hdlr2.setLevel(logging.INFO)
logger.addHandler(hdlr2)

### redirect stdout and stderr into logger
sys.stdout = idq.LogFile(logger)
sys.err = idq.LogFile(logger)

#===================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

ifo = config.get('general', 'ifo')

### figure out which classifiers are present
classifiers = config.get('general', 'classifiers').split(' ')
ovl = "ovl" in classifiers
mvsc = "mvsc" in classifiers
svm = "svm" in classifiers

mla = mvsc or svm  ### label whether machine-learning algorithms are present

#=================================================
### setting paths for input and output data

kw_config_path = config.get('general', 'kwconfig') ### config file for kw pipeline
kwtrgdir = config.get('general', 'kwtrgdir') ### kw triggers directory
realtimedir = config.get('general', 'realtimedir') ### realtime predictions directory
traindir = config.get('general', 'traindir') ### training directory
snglchndir = config.get('general', 'snglchndir') ### single channels trigger directory
                                                 ### required for ovl training

#=================================================
### figure out if we need to train the classifiers before begining evaluation

# training cache files
if ovl:
    ovl_cache = config.get('general', 'ovl_train_cache') ### cache file storing locations of trained classifier files
    if not os.path.exists(ovl_cache):
        if not os.path.exists(os.path.split(ovl_cache)[0]): ### make directory if needed
            os.makedirs(os.path.split(ovl_cache)[0])
        os.mknod(ovl_cache, 0644)

build_auxmvc_vectors = mla and (not os.path.exists(realtimedir)) ### if realtimedir does not exist, we cannot rely on patfiles from the realtime job
                                                                 ### we need to build our own auxmvc_vectors
if mvsc:
    mvsc_cache = config.get('general', 'mvsc_train_cache')
    if not os.path.exists(mvsc_cache):
        if not os.path.exists(os.path.split(mvsc_cache)[0]):
            os.makedirs(os.path.split(mvsc_cache)[0])
        os.mknod(mvsc_cache, 0644)
        build_auxmvc_vectors = True ### we need to build auxmvc vectors because we cannot rely on existence of patfiles

if svm:
    svm_cache = config.get('general', 'svm_train_cache')
    if not os.path.exists(svm_cache):
        if not os.path.exists(os.path.split(svm_cache)[0]):
            os.makedirs(os.path.split(svm_cache)[0])
        os.mknod(svm_cache, 0644)
        build_auxmvc_vectors = True

#=================================================
### parameters for training jobs

min_mvsc_samples = config.getint('idq_train', 'min_mvsc_samples') ### minimum number of samples a training set should have
min_svm_samples = config.getint('idq_train', 'min_svm_samples')
min_ovl_samples = config.getint('ovl_train', 'min_num_glitches')

stride = int(config.get('idq_train', 'stride')) ### stride for training job
delay = int(config.get('idq_train', 'delay')) ### delay for training job
                                              ### provides a buffer for prediction jobs to finish, etc

#==================================================
### current time and boundaries

t = int(idq.nowgps())
if not opts.gps_stop: ### stop time of this analysis
    logger.info('computing gpsstop from current time')
    gpsstop = ( (t - delay) / stride ) * stride  # require boundaries to be integer multiples of stride

else:
    gpsstop = ( opts.gps_stop / stride ) * stride

if not opts.gps_start:
    logger.info('computing gpsstart from gpsstop')
    gpsstart = gpsstop - stride
else:
    gpsstart = ( opts.gps_start / stride ) * stride # require boundaries to be integer multiples of stride

### setting look-back time
### this is how many strides back we look in time to gather triggers
### useful if glitch rate is low but training cadence is high
if opts.lookback == 'infinity': ### we always look back to the beginning of the job (~opts.gps_start)
    lookback = -1
else:
    lookback = int(opts.lookback)

### create condorlogs directory if needed
if not os.path.exists(config.get('idq_train', 'condorlogs')):
    os.makedirs(config.get('idq_train', 'condorlogs'))

#=================================================
### set up ROBOT certificates
### IF ligolw_segement_query FAILS, THIS IS A LIKELY CAUSE
if not opts.no_robot_cert:
    ### unset ligo-proxy just in case
    del os.environ['X509_USER_PROXY']

    ### get cert and key from ini file
    robot_cert = config.get('ldg_certificate', 'robot_certificate')
    robot_key = config.get('ldg_certificate', 'robot_key')

    ### set cert and key
    os.environ['X509_USER_CERT'] = robot_cert
    os.environ['X509_USER_KEY'] = robot_key

#===================================================================================================
#
# MAIN
#
#===================================================================================================

### wait until all jobs are finished
wait = gpsstop + delay - t
if wait > 0:
    logger.info('----------------------------------------------------')
    logger.info('waiting %.1f seconds for jobs to finish' % wait)
    time.sleep(wait)

### iterate over all ranges
while gpsstart < gpsstop:
    logger.info('----------------------------------------------------')

    launch_gps_time = idq.nowgps()

    ### increment lookback time if it is set to infinity
    ### this forces the job to pick up all data since opts.gps_start
    if opts.lookback == 'infinity':
        lookback += 1

    logger.info('\nprocessing data in [%d, %d]' % (gpsstart - lookback * stride, gpsstart + stride))

    ### directory into which we write data
    output_dir = traindir + '/' + str(gpsstart) + '_' + str(gpsstart
            + stride) + '/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #===============================================================================================
    ### generate science times/CAT1 vetoes
    #===============================================================================================
    if opts.ignore_science_segments:
        logger.info('ignoring science segemnts')
        if ovl:
            ovlsegs = False  ### required for OVL input
            vetosegs = False

    else:
        logger.info('generating science segments')
        try: ### query database
            ### this returns a string
            seg_xml_file = idq.segment_query(config, gpsstart - lookback * stride, gpsstart + stride,
                    url=config.get('get_science_segments', 'segdb'))

            ### load xml document
            ### converts string to an object
            xmldoc = utils.load_fileobj(seg_xml_file)[0]

            ### science segments xml filename
            seg_file = "%sscience_segments-%d-%d.xml.gz"%(output_dir, int(gpsstart - lookback * stride), int((lookback+1)*stride))

            logger.info('writing science segments to file : '+seg_file)
            utils.write_filename(xmldoc, seg_file, gz=True)

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
            ### FIXME: writing ascii segment file for OVL, in the future will switch to xml file?
            if build_auxmvc_vectors or (not os.path.exists(realtimedir)): ### mla will use science segments, so we need to write those for ovl
                                                                          ### if realtimedir doesn't exits, we need to use queried scisegs
                try:
                    ### read in segments from xml file
                    (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, config.get('get_science_segments', 'include'))

                    ### write segments to ascii list
                    sciseg_path = "%sscience_segments-%d-%d.seg"%(output_dir, int(gpsstart - lookback * stride), int((lookback+1)*stride))
                    logger.info('writing science segments to file : '+sciseg_path)
                    f = open(sciseg_path, 'w')
                    for line in scisegs:
                        print >> f, line[0], line[1]
                    f.close()
                    ovlsegs = sciseg_path

                except Exception as e:
                    traceback.print_exc()
                    logger.info('WARNING: conversion from xml to ASCII segment file failed.')
    
                    if opts.force:
                        raise e

             ### if we aren't building auxmvc vectors, we re-use pat files from realtime job
            ### this requires us to redefine the 'science-segments' as the intersection of scisegs with realtime segs
            ### currently, we assume realtime only runs when completely in science time, so the intersection is identical to realtime segs
            else:
                try:
                    ### determine segments from realtime filenames
                    realtime_segs = idq.get_idq_segments(realtimedir, gpsstart - lookback * stride, gpsstart + stride, suffix='.pat')

                    ### write segment file
                    idqseg_path = "%s/idq_segements-%d-%d.seg"%(output_dir, int(gpsstart - lookback * stride), int((lookback+1) * stride))
                    f = open(idqseg_path, 'w')
                    for seg in realtime_segs:
                        print >> f, seg[0], seg[1]
                    f.close()

                    ovlsegs = idqseg_path

                except Exception as e:
                    traceback.print_exc()
                    logger.info('WARNING: failed to generate iDQ segments from realtime output.')

                    if opts.force:
                        raise e

            ### NOTHING implemented for vetosegs
            ### these are only used by ovl
            vetosegs = False
            logger.info('Done.')

    #===============================================================================================
    # preparing auxmvc training samples
    #===============================================================================================

    if mla:
        logger.info('preparing training auxmvc samples')

    ### check if auxmvc vectors needs to be build first
    ### should be necessary only at the very first training cycle
    if build_auxmvc_vectors:
        logger.info('WARNING: building auxmvc vectors, this should be necessary only at the very first training cycle')

        ptas_execute_dir = output_dir ### execution directory

        ### output file for training samples
        ptas_output_file = "%s%s-auxmvc_training_samples-%d-%d.pat"%(output_dir, ifo, int(gpsstart - lookback*stride), int((lookback+1)*stride))

        current_dir = os.getcwd() ### remember current directory

        ### pull out pointers to lists of channels
        selected_channels = config.get('general', 'selected-channels')
        unsafe_channels = config.get('general', 'unsafe-channels')

        ### launch job that builds auxmvc_vectors
        (ptas_exit_status, training_samples_file) = \
            idq.execute_build_auxmvc_vectors(
            config,
            ptas_execute_dir,
            kwtrgdir,
            config.get('general', 'gwchannel'),
            ptas_output_file,
            gpsstart - lookback * stride,
            gpsstart + stride,
            channels=selected_channels,
            unsafe_channels=unsafe_channels,
            dq_segments=seg_file,
            dq_segments_name=config.get('get_science_segments',
                    'include'),
            )

        os.chdir(current_dir) ### go back to starting directory

        ### In the next strides we will use samples built by realtime process.
        ### we leave build_auxmvc_vectors=True so that the training jobs can run without an instance of the realtime job
        ### this is NOT the expected use case, but this won't interfere with how the job is launched from the realtime process
#        build_auxmvc_vectors = False 

    elif mla: ### we cat together pat files instead of building vectors from scratch
        ptas_execute_dir = output_dir ### excecution directory

        ### output file for training samples
        ptas_output_file = "%s%s-auxmvc_training_samples-%d-%d.pat"%(output_dir, ifo, int(gpsstart - lookback*strid), int((lookback+1)*stride))

        current_dir = os.getcwd() ### remember current directory

        ### run job that prepares training samples
        (ptas_exit_status, training_samples_file) = \
            idq.execute_prepare_training_auxmvc_samples(
            ptas_execute_dir,
            realtimedir,
            config,
            gpsstart - lookback * stride,
            gpsstart + stride,
            ptas_output_file,
            dq_segments=seg_file,
            dq_segments_name=config.get('get_science_segments',
                    'include'),
            )

        os.chdir(current_dir) ### go back to starting directory

    else:
        ptas_exit_status = 0  # never generated, so it was successful....

    # check if process has been executed correctly

    if ptas_exit_status != 0: ### check that process executed correctly
        logger.info('WARNING: Preparing training auxmvc samples failed')
        if opts.force:
            raise StandardError, "auxmvc samples required for successful training"
        else:
            logger.info('WARNING: skipping re-training the MLA classifiers')

    logger.info('Done.')

    if (ptas_exit_status==0) and mla: ### only load this if it's going to be used
        ### load auxmvc vector samples
        auxmvc_samples = auxmvc_utils.ReadMVSCTriggers([training_samples_file], Classified=False)

        ### get numbers of samples
#        random_samples = auxmvc_samples[numpy.nonzero(auxmvc_samples['i']==0)[0], :])
#        N_clean = len(auxmvc_utils.get_clean_samples(random_samples))

        N_clean = len(auxmvc_samples[numpy.nonzero(auxmvc_samples['i']==0)[0], :])
        N_glitch = len(auxmvc_samples[numpy.nonzero(auxmvc_samples['i']==1)[0], :])


    #===============================================================================================
    # launch training jobs
    #===============================================================================================

    dags = {} ### dictionary that holds the dags submitted for each classifier

    #=============================================
    # mvsc training
    #=============================================

    if mvsc and (ptas_exit_status==0):
        if opts.force or ((N_clean >= min_mvsc_samples) and (N_glitch >= min_mvsc_samples)):
            logger.info('submitting  MVSC training dag')

            mvsc_train_dir = output_dir + 'mvsc/' ### execution directory
            if not os.path.exists(mvsc_train_dir):
                os.makedirs(mvsc_train_dir)

            current_dir = os.getcwd() ### remember current directory

            ### submit mvsc training job
            (mvsc_submit_dag_exit_status, mvsc_dag_file, trained_mvsc_file) = \
                idq.execute_forest_train(training_samples_file,
                    mvsc_cache, config, mvsc_train_dir)

            os.chdir(current_dir) ### switch back to the current directory

            ### check whether dag was submitted succesfully and add it to dictionary of dags
            if mvsc_submit_dag_exit_status != 0:
                logger.info('WARNING: Was not able to submit MVSC training dag')
                dags[mvsc_train_dir + mvsc_dag_file] = 'submit failed'
                if opts.force:
                    raise StandardError, "submission of MVSC training dag failed"

            else:
                dags[mvsc_train_dir + mvsc_dag_file] = 'incomplete'
                logger.info('dag file %s' % (mvsc_train_dir + mvsc_dag_file))

            logger.info('Done.')

        else:
            logger.info('WARNING: Either number of clean or glitch samples in the training set is less than %d'%min_mvsc_samples)
            logger.info('WARNING: Not enough samples in the training set. Skip MVSC training.')
          

    #=============================================
    # SVM training
    #=============================================
    if svm and (ptas_exit_status==0):
        if opts.force or ((N_clean >= min_svm_samples) and (N_glitch >= min_svm_samples)):
            logger.info('submitting SVM training dag')

            svm_train_dir = output_dir + 'svm/' ### excecution directory
            if not os.path.exists(svm_train_dir):
                os.makedirs(svm_train_dir)

            current_dir = os.getcwd() ### remember current directory

            # auxmvc_pat_file_name = training_samples_file   # not check the generated data for svm training
            # auxmvc_pat_file_name = "/home/yingsheng.ji/development/svm_test/train_data/svm_train.pat"   # train data file
            # svm_model_dir = config.get("svm_evaluate", "svm_model_dir")

            ### auxiliary files required for SVM training
            svm_range_file = "%s/%s.range"%(svm_train_dir, os.path.split(training_samples_file)[1])
            svm_model_file = "%s/%s.model"%(svm_train_dir, os.path.split(training_samples_file)[1])
           
            ### submit SVM training job
            (svm_submit_dag_exit_status, svm_dag_file,
             svm_output_files) = idq.execute_svm_train(
                config,
                training_samples_file,
                svm_range_file,
                svm_model_file,
                svm_cache,
                svm_train_dir,
                )

            os.chdir(current_dir) ### switch back to the current directory

            ### check whether dag was submitted succesfully and add it to dictionary of dags
            if svm_submit_dag_exit_status != 0:
                logger.info('WARNING: Was not able to submit SVM training dag')
                dags[svm_train_dir + svm_dag_file] = 'submit failed'
                if opts.force:
                    raise StandardError, "submission of SVM training dag failed"

            else:
                dags[svm_train_dir + svm_dag_file] = 'incomplete'
                logger.info('dag file %s' % (svm_train_dir + svm_dag_file))

            logger.info('Done.')

        else:
            logger.info('WARNING: Either number of clean or glitch samples in the training set is less than %d'%min_svm_samples)
            logger.info('WARNING: Not enough samples in the training set. Skip SVM training.')

    #=============================================
    # OVL training
    #=============================================
    ### we launch this last because it runs on the head node (not through condor)
    ### therefore, all MLA algorithms should be launched before we start the OVL job

    if ovl: 

        #=========================================
        # prepare trigger files so OVL can read them
        #=========================================

        ### OVL reads triggers from single channel summary files
        ### these are not produced in low latency, so we build them now
        logger.info('generating single-channel summary files')
        new_dirs = idq.collect_sngl_chan_kw(
            gpsstart,
            gpsstart + stride,
            kw_config_path,
            width=stride,
            source_dir=kwtrgdir,
            output_dir=snglchndir,
            )

        ### we collect the single channel trg files into xml files
        ### this shouldn't be necessary for OVL training and i doubt these files
        ### will be used by anyone else...
        ### FIXME: WE MAY WANT TO REMOVE THIS STEP
        if not opts.no_sngl_chan_xml:
            logger.info('launching conversion from .trg to .xml files')
            for dir in new_dirs:
                current_dir = os.getcwd()
                trg_to_xml_exit_code = \
                    idq.submit_command([config.get('condor', 'convertkwtosb'
                                      ), dir], process_name='convertkwtosb'
                                      , dir=dir)
                os.chdir(current_dir)
                if trg_to_xml_exit_code != 0:
                    logger.info('WARNING: Conversion from single-channel KW trig files to xml failed in '+ dir)

        logger.info('Done.')

        #=========================================
        # actual ovl training job
        #=========================================
        logger.info('launching ovl training job')

        ### launch ovl training job
        if ovlsegs: ### silly formatting thing
          ovlsegs = [ovlsegs]

        logger.info('idq.ovl_train( %d, %d, %s, scisegs=%s, vetosegs=False, output_dir=%s)'%(gpsstart - lookback * stride, gpsstart + stride, opts.config, str(ovlsegs), output_dir))

        vetolists = idq.ovl_train(
            gpsstart - lookback * stride,
            gpsstart + stride,
            config,
            scisegs=ovlsegs,
            vetosegs=vetosegs,
            output_dir=output_dir,
            )
        logger.info('Done.')

        vetolist = vetolists[0] ### there should only be one vetolist

        ### append new vetolist to training cache
	if opts.force: ### force us to append
		go_ahead_and_append = True
	else:
	        go_ahead_and_append = False
        	f = open(vetolists[0], 'r')
	        for line in f:
        	    if line[0] != '#': ### the first un-commented line will contain the total number of gwtriggers
                	# require certain number of glitches for this to be a good training set
	                go_ahead_and_append = float(line.strip().split()[idq.ovl.vD['#gwtrg']]) >= min_ovl_samples
	                break
        	f.close()
	
        if go_ahead_and_append:
            ### make sure the path doesn't contain any "//" elements............
            vetolist = vetolist.replace("//","/")

            ### append file
            logger.info('appending ovl_train.cache')

            f = open(ovl_cache, 'a')
            print >> f, vetolist
            f.close()

        else:
            logger.info('WARNING: number glitches in the training set is less than ' + str(min_ovl_samples))
            logger.info('WARNING: Not enough glitches in the training set. Skip OVL training.')

        logger.info('Done.')

    #===============================================================================================
    # jobs completion checkpoint (for dags)
    #===============================================================================================

    if len(dags): ### dags is not empty

        print "dag"

        list_of_dags = dags.keys()
        list_of_dags.sort()

        ### loop until either all dags are complete or gpsstop+stride time is reached
        incompleted_dags = []

        ### if opts.force, we wait for all jobs to finish
        ### otherwise, we wait for some amount of time before proceeding
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
                logger.info('WARNING: Some dags have not been completed.')
                time.sleep(300) ### wait 600 seconds before continue the loop
                                ### 600 is a magic number, although 10 min. should be much shorter than stride
            else: ### no incomplete dags, break the loop
                break

        ### print warnings for dags that failed or are incomplete
        for dag in list_of_dags:
            if dags[dag] != 0:
                if dags[dag] == 'incomplete':
                    logger.info('WARNING: %s was  not complete' % dag)
                else:
                    logger.info('WARNING: %s failed with exit status %s'% (dag, str(dags[dag])))
                    if opts.force:
                       raise StandardError, "%s failed with exit status %s"%(dag, str(dags[dag]))

    ### continue onto the next stride
    gpsstart += stride

