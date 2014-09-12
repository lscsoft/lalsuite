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

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
    """This program performs periodic training of the classifiers used in iDQ pipeline. Some intensive training jobs are submitted to cluster using condor. """

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('-c', '--config', default='idq.ini', type='string',
                  help='configuration file')
parser.add_option('-s', '--gps-start', default=False, type='int',
                  help='a GPS start time for the analysis. If default, gpsstart is calculated from the current time.'
                  )
parser.add_option('-e', '--gps-stop', default=False, type='int',
                  help='a GPS stop time for the analysis. If default, gpsstop is calculated from the current time.'
                  )
parser.add_option('-b', '--lookback', default='0', type='string',
                  help="Number of strides to look back and get data for training. Default is zero.\
	Can be either positive integer or 'infinity'. In the latter case, the lookback will be incremented at every stride and all data will be used in training."
                  )
parser.add_option('-l', '--log-file', default='idq_train.log',
                  type='string', help='log file')
(opts, args) = parser.parse_args()

# logging setup

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

# redirect stdout and stderr into logger

sys.stdout = idq.LogFile(logger)
sys.err = idq.LogFile(logger)

# read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

ifo = config.get('general', 'ifo')
classifiers = config.get('general', 'classifiers').split(' ')
mla = 'mvsc' in classifiers or 'svm' in classifiers  # label whether machine-learning algorithms are present

# setting paths for input and output data

# config file for kw pipeline

kw_config_path = config.get('general', 'kwconfig')

# kw triggers directory

kwtrgdir = config.get('general', 'kwtrgdir')

# realtime predictions directory

realtimedir = config.get('general', 'realtimedir')

# training directory

traindir = config.get('general', 'traindir')

# single channels trigger directory

snglchndir = config.get('general', 'snglchndir')

build_auxmvc_vectors = False

# training cache files

if 'ovl' in classifiers:
    ovl_cache = config.get('general', 'ovl_train_cache')
    if not os.path.exists(ovl_cache):
        if not os.path.exists(os.path.split(ovl_cache)[0]):
            os.makedirs(os.path.split(ovl_cache)[0])
        os.mknod(ovl_cache, 0644)
if 'mvsc' in classifiers:
    mvsc_cache = config.get('general', 'mvsc_train_cache')
    if not os.path.exists(mvsc_cache):
        if not os.path.exists(os.path.split(mvsc_cache)[0]):
            os.makedirs(os.path.split(mvsc_cache)[0])
        os.mknod(mvsc_cache, 0644)
        build_auxmvc_vectors = True

if 'svm' in classifiers:
    svm_cache = config.get('general', 'svm_train_cache')
    if not os.path.exists(svm_cache):
        if not os.path.exists(os.path.split(svm_cache)[0]):
            os.makedirs(os.path.split(svm_cache)[0])
        os.mknod(svm_cache, 0644)
        build_auxmvc_vectors = True

# get minimum number of samples a training set should have

min_mvsc_samples = config.getint('idq_train', 'min_mvsc_samples')
min_svm_samples = config.getint('idq_train', 'min_svm_samples')
min_ovl_samples = config.getint('ovl_train', 'min_num_glitches')

# setting stride and delay

stride = int(config.get('idq_train', 'stride'))
delay = int(config.get('idq_train', 'delay'))

# current time and boundaries

t = int(idq.nowgps())
if not opts.gps_stop:
    print 'computing gpsstop from current time'
    gpsstop = (t - delay) / stride * stride  # gpsstop time for this analysis
else:
    gpsstop = opts.gps_stop / stride * stride
if not opts.gps_start:
    print 'computing gpsstart from gpsstop'
    gpsstart = gpsstop - stride
else:
    gpsstart = opts.gps_start / stride * stride

# setting look-back time

if opts.lookback == 'infinity':
    lookback = -1
else:
    lookback = int(opts.lookback)

# create condor logs  directory if it does not exist

if not os.path.exists(config.get('idq_train', 'condorlogs')):
    os.makedirs(config.get('idq_train', 'condorlogs'))

#######
# MAIN
#######

logger.info('----------------------------------------------------')

# logger.info('beginning stride %d-%d' % (t, t+stride))

# wait until all jobs are finished

wait = gpsstop + delay - t
if __name__ == '__main__' and wait > 0:
    logger.info('waiting %.1f seconds for jobs to finish' % wait)
    time.sleep(wait)

# iterate over all ranges

while __name__ == '__main__' and gpsstart < gpsstop:

    # increment lookback time if it is set to infinity

    if opts.lookback == 'infinity':
        lookback += 1
    logger.info('\nprocessing data in [%d, %d]' % (gpsstart - lookback
                * stride, gpsstart + stride))
    launch_gps_time = idq.nowgps()

    # where to output traning data

    output_dir = traindir + '/' + str(gpsstart) + '_' + str(gpsstart
            + stride) + '/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # define dictionary to hold dags

    dags = {}

    # #######################
    # generate science times/CAT1 vetoes

    logger.info('generating science segments')

    # We are switching to robot certificates no need for finding and validating proxy certificate
    # We only neet to set environment to use the robot certificate.

    # unset ligo-proxy just in case

    del os.environ['X509_USER_PROXY']

    # get cert and key from ini file

    robot_cert = config.get('ldg_certificate', 'robot_certificate')
    robot_key = config.get('ldg_certificate', 'robot_key')

    # set cert and key

    os.environ['X509_USER_CERT'] = robot_cert
    os.environ['X509_USER_KEY'] = robot_key

    # find a valid LDG grid certificate and replace with it the current if it is expired
    # try:
    # ....valid_cert = idq.findCredential()[0]
    # ....os.putenv('X509_USER_PROXY', valid_cert)
    # ....os.environ['X509_USER_PROXY'] = valid_cert
    # except:
    # ....traceback.print_exc()
    # ....logger.info('ERROR: Could not find a valid LDG grid certificate. Skipping this training period.')
    # ....gpsstart += stride
    # ....continue

    try:

        # scisegs = get_science_segments(gpsstart, gpsstart+stride)

        seg_xml_file = idq.segment_query(config, gpsstart - lookback
                * stride, gpsstart + stride,
                url=config.get('get_science_segments', 'segdb'))

        # load xml document

        xmldoc = utils.load_fileobj(seg_xml_file)[0]

        # set the name for science segments xml file

        seg_file = output_dir + 'science_segments-' + str(int(gpsstart
                - lookback * stride)) + '-' + str(int(stride + lookback
                * stride)) + '.xml.gz'

        logger.info('writing science segments to file')
        utils.write_filename(xmldoc, seg_file, gz=True)
    except:
        traceback.print_exc()
        logger.info('ERROR: segment generation failed. Skipping this training period.'
                    )
        gpsstart += stride
        continue

    # FIX ME: writing ascii segment file for OVL, in the future will switch to xml file?
    # scisegs = get_science_segments(gpsstart, gpsstart+stride)

    try:
        (scisegs, coveredseg) = idq.extract_dq_segments(seg_file,
                config.get('get_science_segments', 'include'))
        sciseg_path = output_dir + 'science_segments-' \
            + str(int(gpsstart - lookback * stride)) + '-' \
            + str(int(stride + lookback * stride)) + '.seg'
        logger.info('writing idq segments to file')
        f = open(sciseg_path, 'w')
        for line in scisegs:
            print >> f, line[0], line[1]
        f.close()
        ovlsegs = sciseg_path
    except:
        traceback.print_exc()
        logger.info('WARNING: conversion from xml to ASCII segment file failed.'
                    )

    # generate iDQ segments from realtime filenames (for OVL)
    # run it only if pat files already exist
    # FIXME: this step is troubling, it relies on existence of *pat files.
    # FIXME: Also we may need to rewise definition of segments on which OVL should be trained.
    # FIXME: Why not to use science segments from the previous step?

    if not build_auxmvc_vectors:
        try:
            realtime_segs = idq.get_idq_segments(realtimedir, gpsstart
                    - lookback * stride, gpsstart + stride,
                    suffix='.pat')
            idqseg_path = output_dir + '/idq_segments-' \
                + str(int(gpsstart - lookback * stride)) + '-' \
                + str(int(stride + lookback * stride)) + '.seg'
            f = open(idqseg_path, 'w')
            for seg in realtime_segs:
                print >> f, seg[0], seg[1]
            f.close()

            # use idq segments for ovl

            ovlsegs = idqseg_path
        except:
            traceback.print_exc()
            logger.info('WARNING: failed to generate iDQ segments from realtime output.'
                        )

    # NOTHING implemented for vetosegs

    vetosegs = False
    logger.info('Done.')

    # ######################

    # #############################
    #   preparing training auxmvc samples
    # ########################

    logger.info('preparing training auxmvc samples')

    # check if auxmvc vectors needs to be build first
    # should be necessary only at the very first training cycle

    if build_auxmvc_vectors:
        logger.info('WARNING: building auxmvc vectors, this should be necessary only at the very first training cycle'
                    )

        # set directory where the job will be executed

        ptas_execute_dir = output_dir

        # set output file for training samples

        ptas_output_file = output_dir + ifo + '-' \
            + 'auxmvc_training_samples' + '-' + str(int(gpsstart
                - lookback * stride)) + '-' + str(int(stride + lookback
                * stride)) + '.pat'

        # get current directory

        current_dir = os.getcwd()

        selected_channels = config.get('general', 'selected-channels')
        unsafe_channels = config.get('general', 'unsafe-channels')

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

        # set build_auxmvc_vectors to False. In the next strides we will use samples built by realtime process.

        build_auxmvc_vectors = False

        # switch back to the current directory

        os.chdir(current_dir)
    elif mla:

              # check that machine-learning algorithms are present
        # set directory where the job will be executed

        ptas_execute_dir = output_dir

        # set output file for training samples

        ptas_output_file = output_dir + ifo + '-' \
            + 'auxmvc_training_samples' + '-' + str(int(gpsstart
                - lookback * stride)) + '-' + str(int(stride + lookback
                * stride)) + '.pat'

        # get current directory

        current_dir = os.getcwd()

        # run job that prepares training samples

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

        # switch back to the current directory

        os.chdir(current_dir)
    else:

        ptas_exit_status = 0  # never generated, so it was successful....

    # check if process has been executed correctly

    if ptas_exit_status != 0:
        logger.info('WARNING: Preparing training auxmvc samples failed')
        logger.info('WARNING: skipping re-training the MLA classifiers')
    logger.info('Done.')

    if ptas_exit_status == 0:

        # load auxmvc vector samples

        auxmvc_samples = \
            auxmvc_utils.ReadMVSCTriggers([training_samples_file],
                Classified=False)

        # get clean and glitch samples

        random_samples = auxmvc_samples[numpy.nonzero(auxmvc_samples['i'
                ] == 0)[0], :]
        N_clean = len(auxmvc_utils.get_clean_samples(random_samples))
        N_glitch = len(auxmvc_samples[numpy.nonzero(auxmvc_samples['i']
                       == 1)[0], :])

    # ##############
    # mvsc training
    # #############

    if 'mvsc' in classifiers and ptas_exit_status == 0:
        if N_clean >= min_mvsc_samples and N_glitch >= min_mvsc_samples:
            logger.info('submitting  MVSC training dag')

            # set directory where the mvsc training job will be executed

            mvsc_train_dir = output_dir + 'mvsc/'

            # create this directory if it does not exist

            if not os.path.exists(mvsc_train_dir):
                os.makedirs(mvsc_train_dir)

            # get current directory

            current_dir = os.getcwd()

            # submit mvsc training job

            (mvsc_submit_dag_exit_status, mvsc_dag_file,
             trained_mvsc_file) = \
                idq.execute_forest_train(training_samples_file,
                    mvsc_cache, config, mvsc_train_dir)

            # switch back to the current directory

            os.chdir(current_dir)

            # check if dag was submitted succesfully and add it to dictionary of dags

            if mvsc_submit_dag_exit_status != 0:
                logger.info('WARNING: Was not able to submit MVSC training dag'
                            )
                dags[mvsc_train_dir + mvsc_dag_file] = 'submit failed'
            else:
                dags[mvsc_train_dir + mvsc_dag_file] = 'incomplete'
                logger.info('dag file %s' % (mvsc_train_dir
                            + mvsc_dag_file))
            logger.info('Done.')
        else:
            logger.info('WARNING: Either number of clean or glitch samples in the training set is less than '
                         + str(min_mvsc_samples))
            logger.info('WARNING: Not enough samples in the training set. Skip MVSC training.'
                        )

    # #################
    # SVM training
    # ##################

    if 'svm' in classifiers and ptas_exit_status == 0:
        if N_clean >= min_svm_samples and N_glitch >= min_svm_samples:
            logger.info('submitting SVM training dag')

            # set directory where the mvsc training job will be executed

            svm_train_dir = output_dir + 'svm/'

            # create this directory if it does not exist

            if not os.path.exists(svm_train_dir):
                os.makedirs(svm_train_dir)

            # get current directory

            current_dir = os.getcwd()

            # auxmvc_pat_file_name = training_samples_file   # not check the generated data for svm training
            # auxmvc_pat_file_name = "/home/yingsheng.ji/development/svm_test/train_data/svm_train.pat"   # train data file
            # svm_model_dir = config.get("svm_evaluate", "svm_model_dir")

            svm_range_file = svm_train_dir + '/' \
                + os.path.split(training_samples_file)[1] + '.range'
            svm_model_file = svm_train_dir + '/' \
                + os.path.split(training_samples_file)[1] + '.model'
            (svm_submit_dag_exit_status, svm_dag_file,
             svm_output_files) = idq.execute_svm_train(
                config,
                training_samples_file,
                svm_range_file,
                svm_model_file,
                svm_cache,
                svm_train_dir,
                )

            # switch back to the current directory

            os.chdir(current_dir)

            # check if dag was submitted succesfully and add it to dictionary of dags

            if svm_submit_dag_exit_status != 0:
                logger.info('WARNING: Was not able to submit SVM training dag'
                            )
                dags[svm_train_dir + svm_dag_file] = 'submit failed'
            else:
                dags[svm_train_dir + svm_dag_file] = 'incomplete'
                logger.info('dag file %s' % (svm_train_dir
                            + svm_dag_file))
            logger.info('Done.')
        else:
            logger.info('WARNING: Either number of clean or glitch samples in the training set is less than '
                         + str(min_svm_samples))
            logger.info('WARNING: Not enough samples in the training set. Skip SVM training.'
                        )

    # ######################
    # OVL training
    # ######################

    if 'ovl' in classifiers:

        # single channel summary files

        logger.info('generating single-channel summary files')
        new_dirs = idq.collect_sngl_chan_kw(
            gpsstart,
            gpsstart + stride,
            kw_config_path,
            width=stride,
            source_dir=kwtrgdir,
            output_dir=snglchndir,
            )

        # collection into .xml files

        logger.info('launching conversion from .trg to .xml files')
        for dir in new_dirs:
            current_dir = os.getcwd()
            trg_to_xml_exit_code = \
                idq.submit_command([config.get('condor', 'convertkwtosb'
                                   ), dir], process_name='convertkwtosb'
                                   , dir=dir)
            os.chdir(current_dir)
            if trg_to_xml_exit_code != 0:
                logger.info('WARNING: Conversion from single-channel KW trig files to xml failed in '
                             + dir)
        logger.info('Done.')

        # ################################################

        # ovl training job

        logger.info('launching ovl training job')

        vetolists = idq.ovl_train(
            gpsstart - lookback * stride,
            gpsstart + stride,
            config,
            scisegs=[ovlsegs],
            vetosegs=vetosegs,
            output_dir=output_dir,
            )
        logger.info('Done.')

        # append training cache

        go_ahead_and_append = False
        f = open(vetolists[0], 'r')
        for line in f:
            if line[0] != '#':

                # require certain number of glitches for this to be a good training set

                go_ahead_and_append = \
                    float(line.strip().split()[idq.ovl.vD['#gwtrg']]) \
                    >= min_ovl_samples
                break
        f.close()
        if go_ahead_and_append:

            # make sure the path doesn't contain any "//" elements............

            v = [l for l in vetolists[-1].split('/') if l != '']
            vetolist = '/'
            for V in v[:-1]:
                vetolist += V + '/'
            vetolist += v[-1]

            # append file

            f = open(ovl_cache, 'a')
            print >> f, vetolist
            f.close()
            logger.info('appending ovl_train.cache')
        else:
            logger.info('WARNING: number glitches in the training set is less than '
                         + str(min_ovl_samples))
            logger.info('WARNING: Not enough glitches in the training set. Skip OVL training.'
                        )
        logger.info('Done.')

        # #################################################

    # ####################################################
    # jobs completion checkpoint
    # #####################################################

    if dags:
        list_of_dags = dags.keys()
        list_of_dags.sort()

        # loop until either all dags are complete or gpsstop+stride time is reached

        incompleted_dags = []
        while idq.nowgps() < gpsstop + stride or idq.nowgps() \
            < launch_gps_time + stride:
            for dag in list_of_dags:

                # check dag status

                dag_status = idq.get_condor_dag_status(dag)
                dags[dag] = dag_status
                if dag_status == 0:
                    logger.info('%s completed' % dag)
                else:
                    incompleted_dags.append(dag)
            list_of_dags = incompleted_dags
            incompleted_dags = []

            # check if any of the dags are still incomplete

            if 'incomplete' in dags.values():

                # wait 600 seconds before continue the loop

                time.sleep(300)
            else:

                # no incomplete dags, break the loop

                break

        # check if some dags are still incomplete

        if 'incomplete' in dags.values():
            logger.info('WARNING: Some dags have not been completed.')

        # print warnings for dags that failed or are incomplete

        for dag in list_of_dags:
            if not dags[dag] == 0:
                if dags[dag] == 'incomplete':
                    logger.info('WARNING: %s was  not complete' % dag)
                else:
                    logger.info('WARNING: %s failed with exit status %s'
                                 % (dag, str(dags[dag])))

    gpsstart += stride
