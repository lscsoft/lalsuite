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

#===================================================================================================

import os
import sys
import logging
import tempfile
import time

import math

import subprocess

import ConfigParser
from optparse import OptionParser

from laldetchar.idq import event
#from laldetchar.idq import idq
from laldetchar.idq import reed
#from laldetchar.idq import ovl

from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

#====================================================================================================

description = \
    """ This program runs iDQ algorithms and classifies (glitch or not) given gps time(s)."""

parser = OptionParser(version='Name: %%prog\n%s'
                               % git_version.verbose_msg,
                               usage='%prog [options]',
                               description=description)

parser.add_option( '-c', '--config-file', dest='config_file', help='configuration file', metavar='FILE', default='idq.ini' )

parser.add_option( '-t', '--gps-time', dest='gps', type='float', help='the GPS time(s) to be classified', metavar='GPS', default=[], action="append" )

parser.add_option('', '--use-science-segments', default=False, action="store_true", help='analyze times only if they are within the science segements.' )
parser.add_option('', '--padding', default=1.0, type="float", help="the padding around the requested gps time for which we request science segments as needed. \
                                                                    If this is less than the auxmvc_coinc window and we use mla classifiers, \
                                                                    it will be extended to the auxmvc_coinc window.")

parser.add_option( '', '--output-dir', dest='outdir', type='string', help='the output directory', default='.' )
parser.add_option( '', '--usertag', default="", type='string')

(opts, args) = parser.parse_args()

if not opts.gps:
    print "not gps times supplied. Nothing to do..."
    sys.exit(0)

if opts.padding < 1.0:
    raise ValueError("--padding must be >= 1.0")

if not os.path.exists(opts.outdir):
    os.makedirs(opts.outdir)

if opts.usertag:
    opts.usertag = "_%s"%opts.usertag

cwd = os.getcwd()

#====================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config_file)

mainidqdir = config.get('general', 'idqdir') ### get the main directory where idq pipeline is going to be running.

gchtag = "_glitch" ### used for xml filenames
clntag = "_clean"
usertag = "%s%s"%(config.get('general', 'usertag'), opts.usertag)
if usertag:
    usertag = "_%s"%usertag
    gchtag = "%s_%s"%(gchtag, usertag)
    clntag = "%s_%s"%(clntag, usertag)

ifo = config.get('general', 'ifo')

#========================
# which classifiers
#========================
classifiers = sorted(set(config.get('general', 'classifiers').split()))

### ensure we have a section for each classifier and fill out dictionary of options
classifiersD, mla, ovl = reed.config_to_classifiersD( config )

if mla:
    ### reading parameters from config file needed for mla
    auxmvc_coinc_window = config.getfloat('build_auxmvc_vectors','time-window')
    auxmc_gw_signif_thr = config.getfloat('build_auxmvc_vectors','signif-threshold')
    auxmvc_selected_channels = config.get('general','selected-channels')
    auxmvc_unsafe_channels = config.get('general','unsafe-channels')

    ### inflate padding for science time if needed
    if auxmvc_coinc_window > opts.padding:
        print "extending opts.padding to fill the auxmvc_coinc_window"
        opts.padding = auxmvc_coinc_window

twopadding = 2*padding ### useful to have

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

### kleineWelle config
GWkwconfig = reed.loadkwconfig(config.get('data_discovery', 'GWkwconfig'))
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')

AUXkwconfig = reed.loadkwconfig(config.get('data_discovery', 'AUXkwconfig'))
AUXkwbasename = AUXkwconfig['basename']
AUXgdsdir = config.get('data_discovery', 'AUXgdsdir')

#========================
# realtime job
#========================
realtimedir = config.get('general', 'realtimedir')### output directory for realtime predictions
if not os.path.exists(realtimedir):
    os.makedirs(realtimedir)

samples_header = config.get('realtime', 'dat_columns').split()### headers for dat files

### slave the realtime job to the kw stride
stride = int(float(kwconfig['stride'])) ### this is given as a decimal, so we must cast as float first
delay = config.getint('realtime', 'delay') ### buffer to let jobs finish

#========================
# train cache
#========================
traindir = config.get('general', 'traindir')
train_cache = dict( (classifier, reed.Cachefile(reed.cache(traindir, classifier, tag='train'))) for classifier in classifiers )
for cache in train_cache.values():
    cache.time = 0

#========================
# calibration cache
#========================
calibrationdir = config.get('general', 'calibrationdir')
calibration_cache = dict( (classifier, reed.Cachefile(reed.cache(calibrationdir, classifier, tag='calibration'))) for classifier in classifiers )
for cache in calibration_cache.values():
    cache.time = 0

#===================================================================================================
#
# EVALUATION
#
#===================================================================================================

for gps in sorted(opts.gps):
    print '----------------------------------------------------'
    print 'beginning analysis for %.6f' % gps

    t = int(gps - (gps%stride))

    gps_start = int(gps-padding)
    gps_padd = gps_start + twopadding

    ### make sure end time of analysis time is not ahead of now-delay
    wait = opts.gps + delay - int(idq.nowgps())
    if wait > 0:
        print 'waiting %.1f seconds before processing' % wait
        time.sleep(wait)

    #====================
    # science segments
    #====================
    if opts.use_science_segments:
        print 'querrying science segments'

        ### get DMT xml segments from disk
        ### logic about waiting and re-trying is built into retrieve_scisegs
        good, covered = reed.retrieve_scisegs(dmt_segments_location, dmtdq_name, gps_start, twopadding, pad=0, sleep=delay, nretry=1, logger=logger)

        if event.livetime(covered) < twopadding:
            raise Warning('unknown science coverage, skipping')
            continue

        elif event.livetime(good) < twopadding:
            raise Warning('incomplete science coverage, skipping')
            continue

        print 'complete science coverage'

    else:
        print 'analyzing data regardless of science segements'

    #====================
    # samples -> just the supplied gps time
    #====================
    samples = [[gps, 1, 1.0, 0.0, 0.0]]  # fake sample with fields filled in for kw data

    #====================
    # aux triggers
    #====================
    print "aggregating auxiliary triggers"
    trgdict = event.trigdict()

    ### load auxiliary triggers from current stride
    aux_trgdict = reed.retrieve_kwtrg(AUXgdsdir, AUXkwbasename, t, stride, sleep=(t+stride-reed.nowgps()) + delay, ntrials=2, logger=logger)
    if aux_trgdict == None:
        raise Warning('  no auxiliary triggers were found')
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
        print '  previous stride KW triggers are needed'
        prev_trgdict = reed.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t-stride, stride, sleep=0, ntrials=1, logger=logger)
        if prev_trgdict == None:
            raise Warning('  missing previous stride KW triggers')
            ### we do not skip this stride when previous kw stride is required and missing!
        else:
            trgdict.add( prev_trgdict ) ### we add the previous triggers into the current trgdict

    if maxtime > t + stride - padding: ### triggers from the following stride must be included
        print '  next stride KW triggers are needed'
        next_trgdict = reed.retrieve_kwtrg(AUXgdsdir, AUXkwbasename, t+stride, stride, sleep=(t+stride+stride-reed.nowgps()) + delay, ntrials=2, logger=logger)
        if next_trgdict == None:
            raise Warning('  missing next KW triggers')
            ### we do not skip this stride when next kw stride is required and missing!
        else:
            trgdict.add( next_trgdict ) ### we add the next triggers to the current trgdict

    trgdict.resort() ### make sure the trgdict's trigger lists are stored in the correct order

    ### get rid of unwanted AUX triggers
    trgdict.include([[min(mintime - padding, t), max(maxtime + padding, t + stride)]])

    ### get rid of unwanted GW triggers from previous or next strides; these will be evaluated in their own strides
    trgdict.include([[t, t + stride]], channels=[gwchannel])

    #====================
    # generate patfiles for mla classifiers
    #====================
    if mla:  # only build patfiles if machine-learning algorithms are present
        print 'building auxmvc feature vectors ...'

        pat = reed.pat(opts.outdir, ifo, usertag, gps_start, twopadding) #"%s/%s_%d-%d.pat"%(opts.outdir, ifo, gps_start, twopadding)

        # generating auxmvc vector samples. result is saved into pat file
        # FIXME: depending how padding is done we should adjust behavior of build_auxmvc_vectors
        # Currently it keeps gw trigger from [t, t + stride] and uses time_window to pad this segment for auxiliary triggers
                # we do not filter out unclean beacuse it is already done when clean_gps times are formed
        auxmvc_vectors = reed.build_auxmvc_vectors(trgdict, gwchannel, auxmvc_coinc_window, auxmc_gw_signif_thr, pat, gps_start_time=gps_start,
                                gps_end_time=gps_padd,  channels=auxmvc_selected_channels, unsafe_channels=auxmvc_unsafe_channels, clean_times=clean_gps,
                                clean_window=clean_window, filter_out_unclean=False )

    #=============================================
    # predictions
    #=============================================
    for classifier in classifiers:
        flavor = classifiersD[classifier]['flavor']

        print '%s cycle -> flavor=%s'%(classifier, flavor)

        ### find best training data
        cache = train_cache[classifier]

        ### the following is a stupid fudge to allow for multi-line cache file formats
        lines = []
        line = []
        for l in cache.readlines():
            line.append( l )
            if len(line) == reed.traincache_nlines[flavor]:
                lines.append( line )
                line = []
        lines.append( line )
        trained_ranges = dict( zip( reed.extract_trained_ranges( lines, flavor ), lines ) )
        trained_range = reed.best_range( gps, trained_ranges.keys() )
        lines = trained_ranges[trained_range]

        ### find best calibration data
        cache = calibration_cache[classifier]
        lines = cache.readlines() ### calibration caches are forced to have a simpler format that train caches
        calib_ranges = dict( zip(reed.extract_calibration_ranges( lines ), lines ) )
        calib_range = reed.best_range( gps, calib_ranges.keys() )
        effmap, fapmap = reed.rank_to_EffFAPmap( calib_ranges[calib_range] )

        ### compute filenames
        dat = reed.dat(this_output_dir, classifier, ifo, trained_range, usertag, gps_start, twopadding)
        gchxml = reed.xml(directory, classifier, ifo, trained_range, calib_range, gchtag, gps_start, twopadding)




        ### FIXME : check whether these already exist. If so, we need to append instead of overwrite!
        raise Warning("We may over-write output files because of the naming convention! This may need some careful thought if we don't just tack on a random string...")




        ### perform evalutation
        miniconfig = classifiersD[classifier]['config']
        if classifiersD[classifier]['mla']:
            returncode = reed.execute(flavor, lines, dat, miniconfig, gps_start_time=gps_start, gps_end_time=gps_padd, dir=this_output_dir, trgdict=pat, auxmvc_vectors=auxmvc_vectors)
        else:
            returncode = reed.execute(flavor, lines, dat, miniconfig, gps_start_time=gps_start, gps_end_time=gps_padd, dir=this_output_dir, trgdict=trgdict, samples=samples, samples_header=samples_header)

        ### check if process has been execited correctly
        if returncode:
            raise Warning('%s predict failed'%classifier)
            raise Warning('  skipping %s timeseries'%classifier)
            continue
        print '  Done: %s evaluation '%classifier

        ### convert datfiles to xml tables
        print '  converting %s to xml tables'%dat

        ### read dafile -> xml docs
        (gchxml_doc, clnxml_doc) = reed.datfile2xmldocs(dat, ifo, fapmap, Lmap=False, Effmap=effmap, classifier=classifier,
                                         gwchan=gwchannel, gwtrigs=[gps], prog=__prog__, options=opts.__dict__, version=__version__ )

        ### write documents
        print '    --> writing ' + gchxml
        reed.ligolw_utils.write_filename(gchxml_doc, gchxml, gz=gchxml.endswith('.gz'))


