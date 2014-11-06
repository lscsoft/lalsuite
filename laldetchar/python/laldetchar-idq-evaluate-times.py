# Copyright (C) 2014 Reed Essick
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
from laldetchar.idq import ovl
import math
import time
from laldetchar.idq import event
import logging
import tempfile
import optparse
import subprocess
import ConfigParser

from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date


####################################################################################################

def find_best_trained_classifier(gps, train_cache, classifier):
    """ determines which trained classifier in train_cache is best to use for gps and returns it """

    dists = []
    if classifier == 'ovl':
        cache = open(train_cache, 'r')
        for line in cache:
            vetolist = line.strip('\n')
            v_range = [v for v in vetolist.split('/') if v != ''][-3]
            (start, stop) = [int(l) for l in v_range.split('_')]
            if start <= gps and gps <= stop:  # training includes gps
                dists.append((vetolist, 0, start, stop - start))
            elif start > gps:

                                      # starts after gps

                dists.append((vetolist, start - gps, start, stop
                             - start))
            else:

                          # stop < gps

                dists.append((vetolist, gps - stop, start, stop
                             - start))
        cache.close()
    elif classifier == 'mvsc':
        cache = open(train_cache, 'r')
        for line in cache:
            trainedforest = line.strip('\n')
            start = int(trainedforest.split('-')[-2])
            stop = start + int(trainedforest.split('-')[-1].split('.'
                               )[0])
            if start <= gps and gps <= stop:  # training includes gps
                dists.append((trainedforest, 0, start, stop - start))
            elif start > gps:

                                          # starts after gps

                dists.append((trainedforest, start - gps, start, stop
                             - start))
            else:

                              # stop < gps

                dists.append((trainedforest, gps - stop, start, stop
                             - start))
        cache.close()
    elif classifier == 'svm':
        cache = open(train_cache, 'r')
        lines = [l.strip('\n') for l in cache.readlines()]
        cache.close()
        ind = 0
        while ind + 1 < len(lines):
            svm_model = lines[ind]
            svm_range_file = lines[ind + 1]
            start = int(svm_model.split('-')[-2])
            stop = start + int(svm_model.split('-')[-1].split('.')[0])
            if start <= gps and gps <= stop:  # training includes gps
                dists.append(((svm_model, svm_range_file), 0, start,
                             stop - start))
            elif start > gps:

                                          # starts after gps

                dists.append(((svm_model, svm_range_file), start - gps,
                             start, stop - start))
            else:

                              # stop < gps

                dists.append(((svm_model, svm_range_file), gps - stop,
                             start, stop - start))
    else:

        raise ValueError('classifier = %s not understood in find_best_trained_classifier'
                          % classifier)

    # ##

    if not len(dists):
        raise StandardError('train_cache=%s does not contain any trained classifiers!'
                             % train_cache)

    # ## determine which trained classifier is best using proximity to gps in time

    dists.sort(key=lambda l: l[2])  # sort by start times, earlier times at the top (favors causal training)
    dists.sort(key=lambda l: l[3], reverse=True)  # sort by durations, longer durations at the top
    dists.sort(key=lambda l: l[1])  # sort by distance, smaller distances at the top

    return (dists[0][0], dists[0][2], dists[0][3])


##################################################

def find_best_uroc(gps, uroc_cache):
    """ determines which uroc in uroc_cache is best to use for gps and returns it """

    dists = []
    cache = open(uroc_cache, 'r')
    for line in cache:
        uroc = line.strip('\n')
        (start, dur) = uroc.strip('.uroc').split('-')[-2:]
        start = int(start)
        stop = start + int(dur)
        if start <= gps and gps <= stop:  # training includes gps
            dists.append((uroc, 0, start, stop - start))
        elif start > gps:

                                  # starts after gps

            dists.append((uroc, start - gps, start, stop - start))
        else:

                      # stop < gps

            dists.append((uroc, gps - stop, start, stop - start))
    cache.close()

        # ##

    if not len(dists):
        raise StandardError('uroc_cache=%s does not contain any uroc files!'
                             % uroc_cache)

        # ## determine which trained classifier is best using proximity to gps in time

    dists.sort(key=lambda l: l[2])  # sort by start times, earlier times at the top
    dists.sort(key=lambda l: l[3], reverse=True)  # sort by durations, longer durations at the top
    dists.sort(key=lambda l: l[1])  # sort by distance, smaller distances at the top

    return dists[0][0]


####################################################################################################

description = \
    """ This program runs iDQ algorithms and classifies (glitch or not) given gps time(s)."""

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
    '-t',
    '--gps-time',
    dest='gps',
    type='float',
    help='the GPS time(s) to be classified',
    metavar='GPS',
    default=0.0,
    )

parser.add_option(
    '',
    '--output-dir',
    dest='outdir',
    type='string',
    help='the output directory',
    default='./',
    )

(opts, args) = parser.parse_args()

if opts.outdir[-1] != '/':
    opts.outdir += '/'

####################################################################################################

# read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config_file)
myconf = dict(config.items('idq_realtime'))
gwchannel = config.get('general', 'gwchannel')
gwthreshold = config.getfloat('general', 'gw_kwsignif_thr')

classifiers = config.get('general', 'classifiers').split(' ')
mla = 'mvsc' in classifiers or 'svm' in classifiers  # True if machine learning algorithms are present

ifo = config.get('general', 'ifo')
usertag = config.get('general', 'usertag')

# kleineWelle config

kwconfig = idq.loadkwconfig(config.get('general', 'kwconfig'))
kwbasename = kwconfig['basename']
delay = config.getint('idq_realtime', 'delay')

summarydir = config.get('general', 'summarydir')

# training cache files

if 'ovl' in classifiers:
    ovl_train_cache = config.get('general', 'ovl_train_cache')

    # check if cache file already exists

    if not os.path.exists(ovl_train_cache):
        raise StandardError('ovl_train_cache file %s does not exist'
                            % ovl_train_cache)

    ovl_uroc_cache = summarydir + '/ovl_uroc.cache'
    if not os.path.exists(ovl_uroc_cache):
        raise StandardError('ovl_uroc_cache file %s does not exist'
                            % ovl_uroc_cache)

if 'mvsc' in classifiers:
    mvsc_train_cache = config.get('general', 'mvsc_train_cache')

    # check if cache file already exists

    if not os.path.exists(mvsc_train_cache):
        raise StandardError('mvsc_train_cache file %s does not exist'
                            % mvsc_train_cache)

    mvsc_uroc_cache = summarydir + '/mvsc_uroc.cache'
    if not os.path.exists(mvsc_uroc_cache):
        raise StandardError('mvsc_uroc_cache file %s does not exist'
                            % mvsc_uroc_cache)

if 'svm' in classifiers:
    svm_train_cache = config.get('general', 'svm_train_cache')

    # check if cache file already exists

    if not os.path.exists(svm_train_cache):
        raise StandardError('svm_train_cache file %s does not exist'
                            % svm_train_cache)

    svm_uroc_cache = summarydir + '/svm_uroc.cache'
    if not os.path.exists(svm_uroc_cache):
        raise StandardError('svm_uroc_cache file %s does not exist'
                            % svm_uroc_cache)

# column headings for datfile output

samples_header = config.get('idq_realtime', 'dat_columns').split()

# ===================================================================================================
#
#                                             MAIN
#
# ===================================================================================================

print '----------------------------------------------------'
print 'beginning analysis for %.6f' % opts.gps

# make sure end time of analysis time is not ahead of now-delay

wait = opts.gps + delay - int(idq.nowgps())
if wait > 0:
    print 'waiting %.1f seconds before processing' % wait
    time.sleep(wait)

# =================================================
#
# finding KW.trg files
#
# =================================================

padding = max(float(myconf['padding']),
              float(config.get('build_auxmvc_vectors', 'time-window')))

kwdirname = '%s/%s/' % (config.get('general', 'gdsdir'), kwbasename)
kwfilenames = idq.get_all_files_in_range(kwdirname, opts.gps - padding,
        opts.gps + padding, pad=0, suffix='.trg')

if not len(kwfilenames):
    raise StandardError('no KW files found!')

print 'loading KW triggers:\n', kwfilenames[0]
kwdict = event.loadkwm(kwfilenames[0])
for kwfilename in kwfilenames[1:]:
    print kwfilename
    kwdict = event.loadkwm(kwfilename, trigs_dict=kwdict)

# require GW channel in trigger dict

if gwchannel not in kwdict:
    kwdict[gwchannel] = []

kwdict.resort()

samples = [[opts.gps, 1, 1.0, 0.0, 0.0]]  # fake sample with fields filled in for kw data

# =================================================
#
# patfile generation
#
# =================================================

if mla:  # only build patfiles if machine-learning algorithms are present
    print 'Building auxmvc feature vectors and saving them to *.pat file ...'

    # reading parameters from config file

    auxmvc_coinc_window = float(config.get('build_auxmvc_vectors',
                                'time-window'))
    auxmc_gw_signif_thr = float(config.get('build_auxmvc_vectors',
                                'signif-threshold'))
    auxmvc_selected_channels = config.get('general', 'selected-channels'
            )
    auxmvc_unsafe_channels = config.get('general', 'unsafe-channels')
    patfile_output_path = opts.outdir + '/'
    auxmvc_pat_file_name = patfile_output_path + kwbasename + '_mvsc-' \
        + str(opts.gps) + '.pat'

    # generating auxmvc vector samples. result is saved into pat file
    # FIXME: depending how padding is done we  should adjust behavior of  build_auxmvc_vectors
    # Currently it keeps gw trigger from [t, t + stride] and uses time_window to pad this segment for auxiliary triggers

    idq.build_auxmvc_vector_at_gps(
        opts.gps,
        kwdict,
        config.get('general', 'gwchannel'),
        padding,
        auxmvc_pat_file_name,
        channels=auxmvc_selected_channels,
        unsafe_channels=auxmvc_unsafe_channels,
        )

    print auxmvc_pat_file_name
    print 'Done.'

# =================================================
#
# predictions
#
# =================================================

for classifier in classifiers:

        # ## ovl prediction

    if classifier == 'ovl':
        print 'Executing OVL evaluate job ...'

        (vetolist, train_start, train_stop) = \
            find_best_trained_classifier(opts.gps, ovl_train_cache,
                'ovl')
        v_range = '%d_%d' % (train_start, train_stop)
        print vetolist

                # ovl prediction filename

        ovl_filename = kwbasename + '_ovl_' + v_range + '-' + usertag \
            + '-' + str(opts.gps).replace('.', 'd') + '.dat'

                # launch evaluate job

        (gw_predict, ovl_eval_file) = idq.ovl_evaluate(  # kw_trgfiles=False, gwchan=config.get("general", "gwchannel"), filename=ovl_filename, output_dir=this_output_dir, patfiles=False, skip_lines=1)
            vetolist,
            GPStimes=samples,
            GPS_headers=samples_header,
            allvtrg=kwdict,
            filename=ovl_filename,
            output_dir=opts.outdir,
            )

                # create timeseries and save as gzip numpy array with standard filename

        print ' --> generated ' + opts.outdir + '/' + ovl_eval_file

        ovl_uroc_filename = find_best_uroc(opts.gps, ovl_uroc_cache)
        print ovl_uroc_filename

        print 'WARNING: using ovl_Effmap and ovl_FAPmap to estimate likelihood map'
        (ovl_Effmap, ovl_FAPmap) = \
            idq.rank_to_EffFAPmap(ovl_uroc_filename)

                # convert ovl datfile --> xml table

        print 'converting %s to xml tables' % (opts.outdir + '/'
                + ovl_eval_file)

                # read dafile -> xml docs

        (gchxml_doc, _) = idq.datfilename_to_xmldocs(
            opts.outdir + '/' + ovl_eval_file,
            ifo,
            ovl_FAPmap,
            Lmap=False,
            Effmap=ovl_Effmap,
            classifier='ovl',
            )

                # write documents

        gchxml_filename = opts.outdir + '/' + ifo + '_idq_ovl_' \
            + v_range + '-glitch-' + usertag + '-' \
            + str(opts.gps).replace('.', 'd') + '.xml.gz'
        print ' --> writing ' + gchxml_filename
        idq.ligolw_utils.write_filename(gchxml_doc, gchxml_filename,
                gz=gchxml_filename.endswith('.gz'))  # write file

        print 'Done.'
    elif classifier == 'mvsc':

    # ## mvsc prediction

        print 'Executing MVSC evaluate job ...'

        # get lastest trained forest configuration

        (trainedforest, train_start, train_stop) = \
            find_best_trained_classifier(opts.gps, mvsc_train_cache,
                'mvsc')
        mvsc_trained_range = '%d_%d' % (train_start, train_stop)
        print trainedforest

        # set output file and directory where executables will be running

        ranked_file = auxmvc_pat_file_name.split('mvsc')[0] + 'mvsc_' \
            + mvsc_trained_range + '-' + usertag + '-' \
            + str(opts.gps).replace('.', 'd') + '.dat'
        execute_forest_evaluate_dir = opts.outdir

        # run mvsc prediction jobs

        (exit_status, dat_file) = idq.execute_forest_evaluate(
            auxmvc_pat_file_name,
            trainedforest,
            ranked_file,
            config,
            gps_start_time=opts.gps - padding,
            gps_end_time=opts.gps + padding,
            dir=execute_forest_evaluate_dir,
            )

        # check if process has been execited correctly

        if exit_status != 0:
            print 'WARNING: Forest predict failed'

        mvsc_uroc_filename = find_best_uroc(opts.gps, mvsc_uroc_cache)
        print mvsc_uroc_filename

        print 'WARNING: using mvsc_Effmap and mvsc_FAPmap to estimate likelihood map'
        (mvsc_Effmap, mvsc_FAPmap) = \
            idq.rank_to_EffFAPmap(mvsc_uroc_filename)

        # convert mvsc datfile --> xml table

        print 'converting %s to xml tables' % ranked_file

        # read dafile -> xml docs

        (gchxml_doc, _) = idq.datfilename_to_xmldocs(
            ranked_file,
            ifo,
            mvsc_FAPmap,
            Lmap=False,
            Effmap=mvsc_Effmap,
            classifier='mvsc',
            )

        # write documents

        gchxml_filename = opts.outdir + '/' + ifo + '_idq_mvsc_' \
            + mvsc_trained_range + '-glitch-' + usertag + '-' \
            + str(opts.gps).replace('.', 'd') + '.xml.gz'
        print ' --> writing ' + gchxml_filename
        idq.ligolw_utils.write_filename(gchxml_doc, gchxml_filename,
                gz=gchxml_filename.endswith('.gz'))  # write file

        print 'Done.'
    elif classifier == 'svm':

    # ## svm prediction

        print 'Executing SVM evaluate job ...'

        ((svm_model, svm_range_file), train_start, train_stop) = \
            find_best_trained_classifier(opts.gps, mvsc_train_cache,
                'mvsc')
        svm_trained_range = '%d_%d' % (train_start, train_stop)

        svm_evaluate_dir = opts.outdir
        svm_predict_file = auxmvc_pat_file_name.split('mvsc')[0] \
            + 'svm_' + svm_trained_range + '-' + usertag + '-' \
            + str(opts.gps).replace('.', 'd') + '.dat'

        exit_label = idq.execute_svm_evaluate(
            config,
            auxmvc_pat_file_name,
            svm_range_file,
            svm_model,
            svm_predict_file,
            dir=svm_evaluate_dir,
            )
        if exit_label != 0:
            print 'WARNING: SVM predict failed'

        svm_uroc_filename = find_best_uroc(opts.gps, svm_uroc_cache)
        print svm_uroc_filename

        print 'WARNING: using svm_Effmap and svm_FAPmap to estimate likelihood map'
        (svn_Effmap, svm_FAPmap) = \
            idq.rank_to_EffFAPmap(svm_uroc_filename)

            # convert svm datfile --> xml table

        print 'converting %s to xml tables' % svm_predict_file

        # read dafile -> xml docs

        (gchxml_doc, _) = idq.datfilename_to_xmldocs(
            svm_predict_file,
            ifo,
            svm_FAPmap,
            Lmap=False,
            Effmap=svm_Effmap,
            classifier='svm',
            )

            # write documents

        gchxml_filename = opts.outdir + '/' + ifo + '_idq_svm_' \
            + svm_trained_range + '-glitch-' + usertag + '-' \
            + str(opts.gps).replace('.', 'd') + '.xml.gz'
        print ' --> writing ' + gchxml_filename
        idq.ligolw_utils.write_filename(gchxml_doc, gchxml_filename,
                gz=gchxml_filename.endswith('.gz'))  # write file

        print 'Done.'
    else:

        raise ValueError('classifier = %s not understood' % classifier)

