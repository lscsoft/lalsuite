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

import ConfigParser
from optparse import OptionParser

from laldetchar.idq import event
from laldetchar.idq import idq
#from laldetchar.idq import reed as idq

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
classifiersD, mla, ovl = idq.config_to_classifiersD( config )
classifiers = sorted(classifiersD.keys())

### get combiners information and DO NOT add this to classifiersD, we treat them separately!
combinersD, referenced_classifiers = idq.config_to_combinersD( config )
combiners = sorted(combinersD.keys())


if mla:
    ### reading parameters from config file needed for mla
#    auxmvc_coinc_window = config.getfloat('build_auxmvc_vectors','time-window')
#    auxmc_gw_signif_thr = config.getfloat('build_auxmvc_vectors','signif-threshold')
    auxmvc_coinc_window = config.getfloat('realtime', 'padding')
    auxmc_gw_signif_thr = config.getfloat('general', 'gw_kwsignif_thr')
    auxmvc_selected_channels = config.get('general','selected-channels')
    auxmvc_unsafe_channels = config.get('general','unsafe-channels')

    ### inflate padding for science time if needed
    if auxmvc_coinc_window > opts.padding:
        print "extending opts.padding to fill the auxmvc_coinc_window"

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
GWkwconfig = idq.loadkwconfig(config.get('data_discovery', 'GWkwconfig'))
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')

AUXkwconfig = idq.loadkwconfig(config.get('data_discovery', 'AUXkwconfig'))
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
train_cache = dict( (classifier, idq.Cachefile(idq.cache(traindir, classifier, tag='train'))) for classifier in classifiers )
for cache in train_cache.values():
    cache.time = 0

#========================
# calibration cache
#========================
calibrationdir = config.get('general', 'calibrationdir')
calibration_cache = dict( (classifier, idq.Cachefile(idq.cache(calibrationdir, classifier, tag='calibration'))) for classifier in classifiers )
for cache in calibration_cache.values():
    cache.time = 0

### cache for kde estimates
kde_cache = dict( (classifier, idq.Cachefile(idq.cache(calibrationdir, classifier, tag='_calibration-kde%s'%usertag))) for classifier in classifiers )
for cache in kde_cache.values():
    cache.time = 0

kdeD = dict( (classifier,{}) for classifier in classifiers ) ### used to store kde information read in from files within kde_cache

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
        good, covered = idq.retrieve_scisegs(dmt_segments_location, dmtdq_name, gps_start, twopadding, pad=0, sleep=delay, nretry=1, logger=logger)

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
    clean_gps = [] ### we don't need this, but it's faster to define this than to clean up the following code...

    #====================
    # aux triggers
    #====================
    print "aggregating auxiliary triggers"
    trgdict = event.trigdict()

    ### load auxiliary triggers from current stride
    aux_trgdict = idq.retrieve_kwtrg(AUXgdsdir, AUXkwbasename, t, stride, sleep=(t+stride-idq.nowgps()) + delay, ntrials=2, logger=logger)
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
        prev_trgdict = idq.retrieve_kwtrig(AUXgdsdir, AUXkwbasename, t-stride, stride, sleep=0, ntrials=1, logger=logger)
        if prev_trgdict == None:
            raise Warning('  missing previous stride KW triggers')
            ### we do not skip this stride when previous kw stride is required and missing!
        else:
            trgdict.add( prev_trgdict ) ### we add the previous triggers into the current trgdict

    if maxtime > t + stride - padding: ### triggers from the following stride must be included
        print '  next stride KW triggers are needed'
        next_trgdict = idq.retrieve_kwtrg(AUXgdsdir, AUXkwbasename, t+stride, stride, sleep=(t+stride+stride-idq.nowgps()) + delay, ntrials=2, logger=logger)
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

        pat = idq.pat(opts.outdir, ifo, usertag, gps_start, twopadding) #"%s/%s_%d-%d.pat"%(opts.outdir, ifo, gps_start, twopadding)

        # generating auxmvc vector samples. result is saved into pat file
        # FIXME: depending how padding is done we should adjust behavior of build_auxmvc_vectors
        # Currently it keeps gw trigger from [t, t + stride] and uses time_window to pad this segment for auxiliary triggers
                # we do not filter out unclean beacuse it is already done when clean_gps times are formed
        auxmvc_vectors = idq.build_auxmvc_vectors(trgdict, gwchannel, auxmvc_coinc_window, auxmc_gw_signif_thr, pat, gps_start_time=gps_start,
                                gps_end_time=gps_padd,  channels=auxmvc_selected_channels, unsafe_channels=auxmvc_unsafe_channels, clean_times=clean_gps,
                                clean_window=clean_window, filter_out_unclean=False )

    #=============================================
    # predictions
    #=============================================
    dats = {}
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
            if len(line) == idq.traincache_nlines[flavor]:
                lines.append( line )
                line = []
        lines.append( line )
        trained_ranges = dict( zip( idq.extract_trained_ranges( lines, flavor ), lines ) )
        trained_range = idq.best_range( gps, trained_ranges.keys() )
        lines = trained_ranges[trained_range]

        ### find best calibration data
        cache = calibration_cache[classifier]
        lines = cache.readlines() ### calibration caches are forced to have a simpler format that train caches
        calib_ranges = dict( zip(idq.extract_calibration_ranges( lines ), lines ) )
        calib_range = idq.best_range( gps, calib_ranges.keys() )
        effmap, fapmap = idq.rank_to_EffFAPmap( calib_ranges[calib_range] )

        ### compute filenames
        dat = idq.dat(this_output_dir, classifier, ifo, trained_range, usertag, gps_start, twopadding)
        gchxml = idq.xml(directory, classifier, ifo, trained_range, calib_range, gchtag, gps_start, twopadding)

        dats[classifier] = dat ### store for combiners

        ### perform evalutation
        miniconfig = classifiersD[classifier]['config']
        if classifiersD[classifier]['mla']:
            returncode = idq.execute(flavor, lines, dat, miniconfig, gps_start_time=gps_start, gps_end_time=gps_padd, dir=this_output_dir, trgdict=pat, auxmvc_vectors=auxmvc_vectors)
        else:
            returncode = idq.execute(flavor, lines, dat, miniconfig, gps_start_time=gps_start, gps_end_time=gps_padd, dir=this_output_dir, trgdict=trgdict, samples=samples, samples_header=samples_header)

        ### check if process has been execited correctly
        if returncode:
            raise Warning('%s predict failed'%classifier)
            raise Warning('  skipping %s timeseries'%classifier)
            continue
        print '  Done: %s evaluation '%classifier

        ### convert datfiles to xml tables
        print '  converting %s to xml tables'%dat

        ### read dafile -> xml docs
        (gchxml_doc, clnxml_doc) = idq.datfile2xmldocs(dat, ifo, fapmap, Lmap=False, Effmap=effmap, classifier=classifier,
                                         gwchan=gwchannel, gwtrigs=[gps], prog=__prog__, options=opts.__dict__, version=__version__ )

        ### write documents
        print '    --> writing ' + gchxml
        idq.ligolw_utils.write_filename(gchxml_doc, gchxml, gz=gchxml.endswith('.gz'))

    ### iterate through combiners
    if combiners:

        ### load in dat files only once!
        outputs = {}
        for classifier in referenced_classifiers:
            dat = dats[classifier]
            logger.info('reading in data from %s'%dat)
            ### we sort according to GPS so everything is ordered in the same way!
            outputs[classifier] = idq.sort_output( idq.slim_load_datfile( dat, skip_lines=0, columns=samples_header+['rank'] ), 'GPS' )

            ### check kde_caches
            cache = kde_cache[classifier]

            ### need to find best kde information for each classifier
            lines = cache.readlines()
            if len(lines)%2:
                raise ValueError('there must be an even number of lines in kde_cache for %s'%classifier)

            kde_ranges = {}
            for ind in xrange(len(lines)/2):
                kde_cln_name = lines[2*ind]
                kde_gch_name = lines[2*ind+1]
                kde_range = idq.extract_kde_range( kde_cln_name )

                kde_ranges[kde_range] = (kde_cln_name, kde_gch_name)

            kde_range = idq.best_range( gps, kde_ranges.keys() )
            kde_cln_name, kde_gch_name = kde_ranges[kde_range]

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

        n_samples = len(outputs[classifier]['rank']) ### assumes all outputs have the same length, which they should

        for combiner in combiners:

            flavor = combinersD[combiner]['flavor'] ### should always be None, but we pull the reference nonetheless

            these_classifiers = combinersD[combiner]['classifiers'] ### the classifiers this combiner combines

            ### find best calibration data
            cache = calibration_cache[combiner]
            lines = cache.readlines() ### calibration caches are forced to have a simpler format that train caches
            calib_ranges = dict( zip(idq.extract_calibration_ranges( lines ), lines ) )
            calib_range = idq.best_range( gps, calib_ranges.keys() )
            effmap, fapmap = idq.rank_to_EffFAPmap( calib_ranges[calib_range] )

            ### compute filenames
            trained_range = "__".join( ["%s-%s"%(classifier, kdeD[classifier]['kde_range']) for classifier in these_classifiers] ) ### string together all kde ranges

            dat = idq.dat(this_output_dir, classifier, ifo, trained_range, usertag, gps_start, twopadding)
            gchxml = idq.xml(directory, classifier, ifo, trained_range, calib_range, gchtag, gps_start, twopadding)

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


            ### compute joint probabilities
            logger.info('  combining classifier probability distributions')
            if combinersD[combiner]['joint_p(g)'][:3] == "max":
                pofg_joint = numpy.max( pofg_joint, axis=0 )
            elif combinersD[combiner]['joint_p(g)'][:4] == "prod":
                pofg_joint = numpy.prod( pofg_joint, axis=0 )
            else:
                raise ValueError("combiner=%s joint_p(g)=%s not understood"%(combiner, combinersD[combiner]['joint_p(g)']))

            if combinersD[combiner]['joint_p(c)'][:3] == "max":
                pofc_joint = numpy.max( pofc_joint, axis=0 )
            elif combinersD[combiner]['joint_p(c)'][:4] == "prod":
                pofc_joint = numpy.prod( pofc_joint, axis=0 )
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

