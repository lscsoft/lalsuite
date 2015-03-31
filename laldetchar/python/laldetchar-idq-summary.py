# Copyright (C) 2013 Lindy Blackburn, Reed Essick, Ruslan Vaulin
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


import sys
import os
import glob
import time

import traceback
import logging
import ConfigParser
from optparse import OptionParser

import numpy as np

from collections import defaultdict

#from laldetchar.idq import idq
from laldetchar.idq import reed
#from laldetchar.idq import idq_summary_plots as idq_s_p
from laldetchar.idq import reed_summary_plots as rsp 
from laldetchar.idq import event

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
    """This program generates summary html pages from iDQ pipeline output. The summmary pages provide variety of diagnostic and interpretational plots and data."""

#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)

parser.add_option('-c', '--config', default='idq.ini', type='string', help='configuration file')
parser.add_option('-s', '--gps-start', dest="gpsstart", default=False, type='int', help='a GPS start time for the analysis. If default, gpsstart is calculated from the current time.')
parser.add_option('-e', '--gps-stop', dest="gpsstop", default=False, type='int', help='a GPS stop time for the analysis. If default, gpsstop is calculated from the current time.')

parser.add_option('-b', '--lookback', default='0', type='string', help="Number of seconds to look back and get data for training. Default is zero.\
        Can be either positive integer or 'infinity'. In the latter case, the lookback will be incremented at every stride and all data after --gps-start will be used in every training.")

parser.add_option('-t', '--trending', default=[], action='append', type='string', help='Number of seconds to look back for trending plots')

parser.add_option('-l', '--log-file', default='idq_summary.log', type='string', help='log file')

parser.add_option("", "--ignore-science-segments", default=False, action="store_true")

parser.add_option('-f','--force',default=False, action='store_true', help="forces *uroc cache file to be updated, even if we have no data. Use with caution.")

parser.add_option("", "--no-robot-cert", default=False, action="store_true")

(opts, args) = parser.parse_args()

if opts.lookback != "infinity":
    lookback = int(opts.lookback)

for ind, trend in enumerate(opts.trending):
    if trend != "infinity":
        opts.trending[ind] = int(trend)

#===================================================================================================
### setup logger to record processes
logger = reed.setup_logger('idq_logger', opts.log_file, sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = reed.LogFile(logger)
sys.stderr = reed.LogFile(logger)

#===================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

mainidqdir = config.get('general', 'idqdir') ### get the main directory where idq pipeline is going to be running.

usertag = config.get('general', 'usertag')
if usertag:
    usertag = "_%s"%usertag

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

colors = dict( zip( classifiers, 'r b g c m k y'.split() ) )
if len(colors) < len(classifiers):
    raise ValueError("ran out of colors for classifiers!")

#========================
# realtime 
#========================
realtimedir = config.get('general', 'realtimedir')### output directory for realtime predictions

dat_columns = config.get('realtime', 'dat_columns').split()
columns = list(set( dat_columns + "GPS i rank".split() ) )

#========================
# calibration jobs
#========================
calibrationdir = config.get('general', 'calibrationdir')

#========================
# summary jobs
#========================
summarydir = config.get('general', 'summarydir')

stride = config.getint('summary', 'stride')
delay = config.getint('summary', 'delay')

emaillist = config.get('warnings', 'summary')
errorthr = config.getfloat('warnings', 'summary_errorthr')

#summary_script = config.get('condor', 'summary')

#summary_cache = dict( (classifier, reed.Cachefile(reed.cache(summarydir, classifier, tag='summary'))) for classifier in classifiers ) ### not needed?

kde_nsamples = config.getint('summary', 'kde_num_samples')
kde = np.linspace(0, 1, kde_nsamples) ### uniformly spaced ranks used to sample kde

faircoin = kde[:] ### used for ROC plots

cln_linestyle='dashed' ### FIXME: pull from config file?
gch_linestyle='solid' 

FAPthrs = [float(l) for l in config.get('summary','FAPthrs').split()]

#========================
# data discovery
#========================
if not opts.ignore_science_segments:
    ### load settings for accessing dmt segment files
#    dmt_segments_location = config.get('get_science_segments', 'xmlurl')
    dq_name = config.get('get_science_segments', 'include')
#    dq_name = config.get('get_science_segments', 'include').split(':')[1]
    segdb_url = config.get('get_science_segments', 'segdb')

#==================================================
### set up ROBOT certificates
### IF ligolw_segement_query FAILS, THIS IS A LIKELY CAUSE
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

#==================================================
### current time and boundaries

t = int(reed.nowgps())

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
logger.info('Begin: summary')

### wait until all jobs are finished
wait = gpsstart + stride + delay - t
if wait > 0:
    logger.info('----------------------------------------------------')
    logger.info('waiting %.1f seconds to reach gpsstop+delay=%d' % (wait, delay))
    time.sleep(wait)

global_start = gpsstart

### iterate over all ranges
while gpsstart < gpsstop:

    logger.info('----------------------------------------------------')

    wait = gpsstart + stride + delay - t
    if wait > 0:
        logger.info('waiting %.1f seconds to reach gpsstart+stride+delay=%d' %(wait, gpsstart+stride+delay))
        time.sleep(wait)

    logger.info('Begin: stride [%d, %d]'%(gpsstart, gpsstart+stride))

    ### directory into which we write data
    output_dir = "%s/%d_%d/"%(summarydir, gpsstart, gpsstart + stride)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if opts.lookback=="infinity":
        lookback = gpsstart - global_start

    for ind, trend in enumerate(opts.trending):
        if trend == "infinity":
            opts.trending[ind] = gpsstart - global_start

    #===============================================================================================
    # science segments
    # we query the segdb right now, although that latency may be an issue...
    #===============================================================================================
    if opts.ignore_science_segments:
        logger.info('analyzing data regardless of science segements')
        scisegs = [[gpsstart-lookback, gpsstart+stride]] ### set segs to be this stride range
        coveredsegs = [[gpsstart-lookback, gpsstart+stride]] ### set segs to be this stride range

    else:
        logger.info('Begin: querrying science segments')

        try:
            ### this returns a string
            seg_xml_file = reed.segment_query(config, gpsstart - lookback , gpsstart + stride, url=segdb_url)

            ### write seg_xml_file to disk
            lsctables.use_in(ligolw.LIGOLWContentHandler)
            xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

            ### science segments xml filename
            seg_file = reed.segxml(output_dir, "_%s"%dq_name, gpsstart - lookback , lookback+stride)

            logger.info('writing science segments to file : '+seg_file)
            ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))

            (scisegs, coveredseg) = reed.extract_dq_segments(seg_file, dq_name) ### read in segments from xml file

        except Exception as e:
            traceback.print_exc()
            logger.info('ERROR: segment generation failed. Skipping this training period.')

            if opts.force: ### we are require successful training or else we want errors
                logger.info(traceback.print_exc())
                raise e
            else: ### we don't care if any particular training job fails
                gpsstart += stride
                continue

    logger.info('finding idq segments')
    idqsegs = reed.get_idq_segments(realtimedir, gpsstart-lookback, gpsstart+stride, suffix='.dat')

    logger.info('taking intersection between science segments and idq segments')
    idqsegs = event.andsegments( [scisegs, idqsegs] )

    ### write segment file
    if opts.ignore_science_segments:
        idqseg_path = reed.idqsegascii(output_dir, '', gpsstart - lookback, lookback+stride)
    else:
        idqseg_path = reed.idqsegascii(output_dir, '_%s'%dq_name, gpsstart - lookback, lookback+stride)
    f = open(idqseg_path, 'w')
    for seg in idqsegs:
        print >> f, seg[0], seg[1]
    f.close()

    #===============================================================================================
    # find data
    #===============================================================================================
    ### find all *dat files, bin them according to classifier
    logger.info('finding all *dat files')
    datsD = defaultdict( list )
    for dat in reed.get_all_files_in_range(realtimedir, gpsstart-lookback, gpsstart+stride, pad=0, suffix='.dat' ):
        datsD[reed.extract_dat_name( dat )].append( dat )

    ### throw away any un-needed files
    for key in datsD.keys():
        if key not in classifiers:
            datsD.pop(key)
        else: ### throw out files that don't contain any science time
            datsD[key] = [ dat for dat in datsD[key] if event.livetime(event.andsegments([idqsegs, [reed.extract_start_stop(dat, suffix='.dat')]])) ]

    #===============================================================================================
    # build plots
    #===============================================================================================
    fignames = {'roc':{}, 'hst':{}, 'kde':{}, 'trg':{}, 'bwh':{}}
    ### roc -> Receiver Operating Characteristic curves
    ### hst -> raw histograms of distributions
    ### kde -> kde smoothed histograms of distributions
    ### trg -> histograms of gch features
    ### bwh -> bit-word histograms

    ### axes for overlay plots
    roc_figax = None
    hst_figax = None
    kde_figax = None

    ### bit-word storage
    gch_bitword = dict( (FAP, {}) for FAP in FAPthrs )
    cln_bitword = dict( (FAP, {}) for FAP in FAPthrs )

    #====================
    # build plots for classifiers separately and add to overlays
    # only plots that depend on this stride alone!
    #====================
    for classifier in classifiers:
        color = colors[classifier]

        ### write list of dats to cache file
        cache = reed.cache(output_dir, classifier, "_datcache%s"%usertag)
        logger.info('writing list of dat files to %s'%cache)
        f = open(cache, 'w')
        for dat in datsD[classifier]:
            print >>f, dat
        f.close()

        logger.info('building plots for %s'%classifier)

        ### load dat files
        output = reed.slim_load_datfiles(datsD[classifier], skip_lines=0, columns=columns)

        ### filter times by scisegs -> keep only the ones within scisegs
        out = np.array(event.include( [ [ float(output['GPS'][i]) ] + [output[column][i] for column in columns] for i in xrange(len(output['GPS'])) ], idqsegs, tcent=0 ))
        for ind, column in enumerate(columns):
            if column == 'GPS':
                output[column] = out[:,1+ind].astype(float)
            elif column == 'i':
                output['i'] = out[:,1+ind].astype(int)
            elif column == 'rank':
                output['rank'] = out[:,1+ind].astype(float)
            else:
                output[column] = out[:,1+ind]

        ### compute rcg from output
        r, c, g = reed.dat_to_rcg( output )
        dc, dg = reed.rcg_to_diff( c, g ) ### get the numbers at each rank

        ### dump into roc file
        roc = reed.roc(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writting %s'%roc)
        reed.rcg_to_file(roc, r, c, g)

        ### generate ROC plot
        rocfig = rsp.rocfig(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['roc'][classifier] = rocfig ### store for reference
        logger.info('  plotting %s'%rocfig)

        fig, ax = rsp.rcg_to_rocFig(c, g, color=color, label=classifier)
        ax.plot(faircoin, faircoin, 'k--')
        ax.legend(loc='best')

        fig.savefig(rocfig)
        rsp.close(fig)

        roc_figax = rsp.rcg_to_rocFig(c, g, color=color, label=classifier, figax=roc_figax) ### for overlay plot

        ### generate histogram, cdf 
        histfig = rsp.histfig(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['hst'][classifier] = histfig ### store for reference
        logger.info('  plotting %s'%histfig)

        if np.any(dc):
            fig, axh, axc = rsp.stacked_hist(r, dc, color=color, linestyle=cln_linestyle, label="%s cln"%classifier)
            hst_figax = rsp.stacked_hist(r, dc, color=color, linestyle=cln_linestyle, label="%s cln"%classifier, figax=hst_figax) ### for overlay plots
        if np.any(dg):
            fig, axh, axc = rsp.stacked_hist(r, dg, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=(fig, axh, axc))
            hst_figax = rsp.stacked_hist(r, dg, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=hst_figax)
        if np.any(dc) or np.any(dg):
            axh.set_xlim(xmin=0, xmax=1)
            axc.set_xlim(xmin=0, xmax=1)
            axc.legend(loc='best')
            fig.savefig(histfig)
            rsp.close(fig)

        ### compute kde estimates
        logger.info('  computing kde_pwg estimates')
        kde_cln = rsp.kde_pwg( kde, r, dc )
        kde_gch = rsp.kde_pwg( kde, r, dg )

        ### write kde points to file
        kde_cln_name = rsp.kdename(output_dir, classifier, ifo, "_cln%s"%usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writing %s'%kde_cln_name)
        np.save(event.gzopen(kde_cln_name, "w"), (kde, kde_cln))

        kde_gch_name = rsp.kdename(output_dir, classifier, ifo, "_gch%s"%usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writing %s'%kde_gch_name)
        np.save(event.gzopen(kde_gch_name, "w"), (kde, kde_gch))

        ### generate kde pdf, cdf (above)
        kdefig = rsp.kdefig(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['kde'][classifier] = kdefig ### store for reference
        logger.info(' plotting %s'%kdefig)

        fig, axh, axc = rsp.stacked_kde(kde, kde_cln, color=color, linestyle=cln_linestyle, label="%s cln"%classifier, figax=None)
        fig, axh, axc = rsp.stacked_kde(kde, kde_gch, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=(fig, axh, axc))
        axc.legend(loc='best')

        fig.savefig(kdefig)
        rsp.close(fig)

        kde_figax = rsp.stacked_kde(kde, kde_cln, color=color, linestyle=cln_linestyle, label="%s cln"%classifier, figax=kde_figax) ### for overlay plot
        kde_figax = rsp.stacked_kde(kde, kde_gch, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=kde_figax)

        ### figure out which gch and cln are removed below FAP thresholds
        rankthrs = []
        for FAP in FAPthrs:
            rankthr = min(r[c<=FAP*c[-1]]) ### smallest rank below FAP
            rankthrs.append( rankthr )
            cln_bitword[FAP][classifier], gch_bitword[FAP][classifier] = rsp.bin_by_rankthr(rankthr, output, columns=dat_columns)
                
        ### generate histograms of glitch parameters
        for column in dat_columns:
            paramhist = rsp.histfig(output_dir, classifier, ifo, "%s_%s"%(usertag, column), gpsstart-lookback, lookback+stride)
            fignames['trg'][(classifier, column)] = paramhist
            logger.info('  plotting %s'%paramhist)

            figax = None
            for rankthr, FAP in zip(rankthrs, FAPthrs):
                thr = min(r[c<=FAP*c[-1]]) ### smallest rank below FAP
                x = [float(output[column][i]) for i in xrange(len(output[column])) if output['rank'][i] >= thr]
                if x: ### will break if x is empty
                    figax = rsp.stacked_hist( x, np.ones_like(x), color=None, label="$FAP\leq%E$"%FAP, bmin=min(x), bmax=max(x))
            if figax: ### only save if we plotted something...
                fig, axh, axc = figax
                axh.set_xlabel(column)
                axc.set_xlim(axh.get_xlim())
                ax.legend(loc='best')
                fig.savefig(paramhist)
                rsp.close(fig)

    #====================
    # save overlays for this stride
    #====================

    ### roc overlay
    rocfig = rsp.rocfig(output_dir, "overlay", ifo, usertag, gpsstart-lookback, lookback+stride)
    fignames['roc']["overlay"] = rocfig ### store for reference
    logger.info('  plotting %s'%rocfig)
    fig, ax = roc_figax
    ax.plot(faircoin, faircoin, 'k--')
    ax.legend(loc='best')
    fig.savefig(rocfig)
    rsp.close(fig)

    ### histogram overlay
    if hst_figax: ### breaks if we have zero samples...
        histfig = rsp.histfig(output_dir, "overlay", ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['hst']["overlay"] = histfig ### store for reference
        logger.info('  plotting %s'%histfig)
        fig, axh, axc = hst_figax
        axh.set_xlim(xmin=0, xmax=1)
        axc.set_xlim(xmin=0, xmax=1)
        axc.legend(loc='best')
        fig.savefig(histfig)
        rsp.close(fig)

    ### kde overlay
    kdefig = rsp.kdefig(output_dir, "overlay", ifo, usertag, gpsstart-lookback, lookback+stride)
    fignames['kde']["overlay"] = kdefig ### store for reference
    logger.info(' plotting %s'%kdefig)
    fig, axh, axc = kde_figax
    axc.legend(loc='best')
    fig.savefig(kdefig)
    rsp.close(fig)

    #====================
    # bitword histograms
    #====================

    ### gch bitword histogram
    bitfig = rsp.bitfig(output_dir, ifo, "%s_gch"%usertag, gpsstart-lookback, lookback+stride)
    fignames['bwh']['gch'] = bitfig
    logger.info('  plotting %s'%bitfig)

    figax=None
    for FAP in FAPthrs:
        figax = rsp.bitword( gch_bitword[FAP], classifiers, label="FAP$\leq$%.3e"%FAP, figax=figax )
    if figax:
        fig, ax = figax
        ax.legend(loc='upper center')

        fig.savefig(bitfig)
        rsp.close(fig)

    ### cln bitword histogram
    bitfig = rsp.bitfig(output_dir, ifo, "%s_cln"%usertag, gpsstart-lookback, lookback+stride)
    fignames['bwh']['cln'] = bitfig
    logger.info('  plotting %s'%bitfig)

    figax=None
    for FAP in FAPthrs:
        figax = rsp.bitword( cln_bitword[FAP], classifiers, label="FAP$\leq$%.3e"%FAP, figax=figax )
    if figax:
        fig, ax = figax
        ax.legend(loc='upper center')

        fig.savefig(bitfig)
        rsp.close(fig)


    #=============================================
    # trending plots!
    # data depends on previous strides as well
    #=============================================

    logger.warning("WARNING: write trending plots!")
    """
    figure showing FAP calibration: 
        stated FAP vs deadtime <- pick up fap*npy.gz files for this!
        do this for both point estimates and 90% UL

TRENDING:

    NEED A GOOD WAY OF FINDING HISTORICAL DATA/FIGURES

    reference to historical plots that rely on only a single stride?
        overlay those plots on the same axis? (need to save plotting data into pickles for easy reference)
        just reproduce the old figures (point to them on disk) with an option to collapse a tab and hide these?

    first and second moment of residuals between historical ROC curves (trending)

    livetime trending plot: both instrument livetime (scisegs) and idq livetime (idqsegs)
    clean rate trending
    glitch rate trending
    eff at fixed FAP trending

    channel rank trending
    channel eff trending
    channel fap trending
    repeat for configurations?
    """

    #===============================================================================================
    # write html page
    #===============================================================================================

    """
REPORT:
time-stamp of page creation
command line options used to generate this page

science segments used: link to ascii file and xml file
idq segments used: link to ascii file and xml file

GWchannels, frequency band, signif_thr, etc
auxchannels, frequency bands, signif_thrs, etc
unsafe channels (not used, but report them for reference)

vetolist/configuration lists -> make human readable
segment lists -> a form to request segments?

table showing number of cln, gch for each classifier

all plots, figures built above
try to make them in expandable sections
    different trending ranges
    different types of plots
    else?
    """

    #=============================================
    # format index.html for summary_dir
    #=============================================

    """
    make this formatted so that it doesn't take forever to load.
    mirroring the current implementation is ok, but we can probably make the pages better using some standard python library?
    """

    #===============================================================================================

    gpsstart += stride

#===================================================================================================

import sys
sys.exit(0)











































def generate_html(
    path,
    gpsstart=False,
    gpsstop=False,
    gwchannel=False,
    tot_livetime=False,
    tot_glitches=False,
    roc_url=False,
    vetolist_url=False,
    ovl_segments=False,
    ovl_fdt=False,
    eff_url=False,
    FAP=False,
    roc_urls=False,
    kde_urls=False,
    chanlist_urls=False,
    stat_trends_urls=False,
    eff_trends_urls=False,
    stat_file=False,
    chan_trends_urls=False,
    ):
    """ writes the standard html summary page given the available information"""

    tzname = time.tzname[0]
    f = open(path, 'w')

    print >> f, '''<body>
<h1>iDQ Summary Page</h1>
<hr />'''

    if gpsstart and gpsstop:
        gpsstart_date = time.localtime(gpsstart + idq.gpsref
                + time.timezone)  # RE is not sure this is the correct way to go from gps time to local time
        gpsstop_date = time.localtime(gpsstop + idq.gpsref
                + time.timezone)  #
        print >> f, \
            '<p>data from : %2d:%2d:%2d %s %2d/%2d/%4d (%d GPSseconds) to %2d:%2d:%2d %s %2d/%2d/%4d (%d GPS seconds)</p>' \
            % (
            gpsstart_date.tm_hour,
            gpsstart_date.tm_min,
            gpsstart_date.tm_sec,
            tzname,
            gpsstart_date.tm_mon,
            gpsstart_date.tm_mday,
            gpsstart_date.tm_year,
            gpsstart,
            gpsstop_date.tm_hour,
            gpsstop_date.tm_min,
            gpsstop_date.tm_sec,
            tzname,
            gpsstop_date.tm_mon,
            gpsstop_date.tm_mday,
            gpsstop_date.tm_year,
            gpsstop,
            )
    if gwchannel:
        print >> f, '<p>Target channel : %s</p>' % gwchannel
    if tot_livetime:
        print >> f, '<p>total Livetime = %f seconds</p>' % tot_livetime
    if tot_glitches:
        print >> f, '<p>total No.glitches = %d events</p>' \
            % int(tot_glitches)
    if roc_url:
        print >> f, '<h2>Receiver Operating Characteristic Curve</h2>'
        print >> f, \
            '<img src="%s" alt="ROC curve" title="ROC curve" />' \
            % roc_url
    if stat_file:
        print >> f, \
            '<p>Basic statistic for each classifier can be found here: '
        print >> f, '<a href="%s">%s</a>, ' % (stat_file,
                stat_file.split('/')[-1])
        print >> f, '</p>'
    if roc_urls:
        print >> f, \
            '<p>Individual-classifier ROC curves can be found here: '
        for (url, classifier) in roc_urls:
            if url:
                print >> f, '<a href="%s">%s</a>, ' % (url, classifier)
        print >> f, '</p>'
    if chanlist_urls:
        print >> f, \
            '<p>Auxiliary channel performance for each classifier can be found here: '
        for (url, classifier) in chanlist_urls:
            if url:
                print >> f, '<a href="%s">%s<a/>, ' % (url, classifier)
        print >> f, '</p>'
    print >> f, '<hr />'
    if vetolist_url:
        print >> f, \
            '<p> Most recent ovl-vetolist can be found <a href="%s">here</a></p>' \
            % vetolist_url
    if ovl_segments and ovl_fdt:
        print >> f, \
            '<p> Most recent ovl-segments at %f FAP can be found <a href="%s">here</a></p>' \
            % (ovl_fdt, ovl_segments)
    if eff_url and FAP:
        print >> f, \
            '<h2>Efficiency at False-Alarm-Probability = %f</h2>' % FAP
        print >> f, \
            '<img src="%s" alt="EFF at FAP=%f" title="EFF at FAP=%f />' \
            % (eff_url, FAP, FAP)
    if kde_urls:
        print >> f, \
            '<h2>Estimates of Classifier Probability Density Functions</h2>'
        print >> f, \
            '<p>Individual-classifier estimates can be found here: '
        for (url, classifier) in kde_urls:
            if url:
                print >> f, '<a href="%s">%s</a>, ' % (url, classifier)
        print >> f, '</p>'

    if eff_trends_urls:
        print >> f, '<h2>Efficiency trends</h2>'
        for (url, descriptor) in eff_trends_urls:
            if url:
                print >> f, \
                    '<img src="%s" alt="efficiency trend %s" title="efficiency trend %s" />' \
                    % (url, descriptor, descriptor)
    if stat_trends_urls:
        print >> f, '<h2>Other statistical trends</h2>'
        print >> f, \
            '<p> Plots of livetime, glitch/clean rates can be found here: '
        for (url, descriptor) in stat_trends_urls:
            print >> f, '<a href="%s">%s</a>, ' % (url, descriptor)
        print >> f, '</p>'

    if chan_trends_urls:
        print >> f, '<h2>Channel performance trends</h2>'
        for (url, descriptor) in chan_trends_urls:
            if url:
                print >> f, \
                    '<img src="%s" alt="channel performance trend %s" title="channel performance trend %s" />' \
                    % (url, descriptor, descriptor)

    print >> f, '<hr />'

    # define current time

    c_time = time.localtime()
    print >> f, '<p>last updated %d:%d:%d %s %2d/%2d/%4d </p>' % (
        c_time.tm_hour,
        c_time.tm_min,
        c_time.tm_sec,
        tzname,
        c_time.tm_mon,
        c_time.tm_mday,
        c_time.tm_year,
        )
    print >> f, '</body>'

    f.close()
    return path

###
def path_to_url(path, base, remove):
    """ converts a path name to a url by removing all directory reference in "remove" and appending the remainder of "path" to "base" """

    path_list = path.split('/')
    for p in path_list[:-1]:
        if p not in remove and p != '':
            base += p + '/'
    return base + path_list[-1]


#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('-c', '--config', default='idq.ini', type='string', help='configuration file')
parser.add_option('-s', '--gps-start', default=False, type='int',
                  help='a GPS start time for the analysis. If default, gpsstart is calculated from the current time.')
parser.add_option('-e', '--gps-stop', default=False, type='int',
                  help='a GPS stop time for the analysis. If default, gpsstop is calculated from the current time.')
parser.add_option('-l', '--log-file', default='idq_summary.log', type='string', help='log file')

parser.add_option("", "--ignore-science-segments", default=False, action="store_true")

parser.add_option('-f','--force',default=False, action='store_true', help="forces *uroc cache file to be updated, even if we have no data. Use with caution.")

parser.add_option("", "--no-robot-cert", default=False, action="store_true")

(opts, args) = parser.parse_args()

#===================================================================================================
### setup logger to record process
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

#=================================================
### generate a dictionary for idq_summary specific options

myconf = dict(config.items('summary'))

stride = int(myconf['stride']) ### summary stride
delay = int(myconf['delay']) ### delay for summary job
                             ### provides a buffer for other jobs to finish

lookback = int(myconf['lookback']) ### how many summary strides we include in these plots

cluster_win = float(myconf['cluster_win']) ### used to cluster glitches

FAP = float(myconf['fap']) ### FAP at which we report segments, efficiency trends, etc

min_num_cln = float(myconf['min_num_cln'])
min_num_gch = float(myconf['min_num_gch'])

### THESE ARE NOT USED. CONSIDER REMOVING?
#symlink_path = myconf['symlink']
#url_base = myconf['url_base']
#url_remove = myconf['url_remove'].split()

gw_thr = float(myconf['gw_thr']) 

classifiers_switch_snr_signif = myconf['switch_snr_signif'].split() ### old "feature" with incorrectly labeled columns

#=================================================
### pull other parameters from config
sumdir = config.get('general', 'summarydir')
realtimedir = config.get('general', 'realtimedir')
traindir = config.get('general', 'traindir')

usertag = config.get('general', 'usertag')

gwchannel = config.get('general', 'gwchannel')

classifiers = config.get('general', 'classifiers').split()

vetolist_cache = config.get('general', 'ovl_train_cache')

unsafe_win = float(config.get('realtime', 'clean_window'))

columns = config.get('realtime', 'dat_columns').split() + ['rank']

kwtrgdir = config.get('general', 'kwtrgdir')

kde_num_samples = int(config.get('summary', 'kde_num_samples'))

### kleineWelle config
kwconfig = idq.loadkwconfig(config.get('general', 'kwconfig'))
kwbasename = kwconfig['basename']

#=================================================
### set classifier colors and labels
classifier_colors = [idq_s_p.classifier_colors(classifier) for classifier in classifiers]
classifier_labels = [idq_s_p.classifier_labels(classifier) for classifier in classifiers]

#=================================================
### set up ROBOT certificates
### IF ligolw_segement_query FAILS, THIS IS A LIKELY CAUSE
if opts.no_robot_cert:
    logger.info("Warning: running without a robot certificate. Your personal certificate may expire and this job may fail")
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
### current time and boundaries

t = int(idq.nowgps())
if not opts.gps_stop: ### stop time of this analysis
    print 'computing gpsstop from current time'
    gpsstop = (t - delay) / stride * stride  # require boundaries to be integer multiples of stride
else:
    gpsstop = opts.gps_stop / stride * stride
#print 'gpsstop = %d' % gpsstop

if not opts.gps_start:
    print 'computing gpsstart from gpsstop'
    gpsstart = gpsstop - stride
else:
    gpsstart = opts.gps_start / stride * stride # require boundaries to be integer mutliples of stride
#print 'gpsstart = %d' % gpsstart

#===================================================================================================
#
# MAIN
#
#===================================================================================================

# loop over all data ranges

while gpsstart < gpsstop:
    logger.info('-----------------------------------------------------------------'
                )
    logger.info('summarizing data from %d to %d' % (gpsstart, gpsstart + stride))

    ### output directory for this data
    this_sumdir = "%s/%d_%d"%(sumdir,gpsstart,gpsstart+stride)
    if not os.path.exists(this_sumdir):
        os.makedirs(this_sumdir)

    ### attempt to write an index.html
    index_html = "%s/index.html"%(this_sumdir)
    index_html_obj = open(index_html, "w")
    print >> index_html_obj, "attempted to build this page at:"
    c_time = time.localtime()
    tzname = time.tzname[0]
    print >> index_html_obj, '%d:%d:%d %s %2d/%2d/%4d' % (
        c_time.tm_hour,
        c_time.tm_min,
        c_time.tm_sec,
        tzname,
        c_time.tm_mon,
        c_time.tm_mday,
        c_time.tm_year,
        )
    print >> index_html_obj, "Target Channel: %s"%gwchannel
    index_html_obj.close()

    #=============================================
    # sciseg query
    #=============================================
    if opts.ignore_science_segments:
        logger.info("ignoring science segments")
        scisegs = [[gpsstart, gpsstart+stride]]

    else:
        logger.info("generating science segments")
        try:
            seg_xml_file = idq.segment_query(config, gpsstart, gpsstart+stride, url=config.get("get_science_segments","segdb"))

            lsctables.use_in(ligolw.LIGOLWContentHandler)
            xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

            seg_file = "%s/science_segements-%d-%d.xml.gz"%(this_sumdir, int(gpsstart), int(stride))
            logger.info("writting science segments to file : %s"%seg_file)
            ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))

            (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, config.get('get_science_segments', 'include'))

        except Exception as e:
            traceback.print_exc()
            logger.info("ERROR: segment generation failed. Skipping this summary period.")

            gpsstart += stride
            continue

    #=============================================
    # generating summary datfiles filtered by segments
    #=============================================
    ### get all relevant datfiles
    datfiles = idq.get_all_files_in_range(realtimedir, gpsstart, gpsstart+stride, suffix=".dat")

    for classifier in classifiers:
        ### follow standard datfile naming conventions 
        summary_dat = "%s/%s_%s_0-0_%sSummary-%d-%d.dat"%(this_sumdir, kwbasename, classifier, usertag, gpsstart, stride)

        logger.info("generating summary dat file for %s : %s"%(classifier, summary_dat))

        columns = ['GPS', 'i', 'rank', "signif", "SNR"]
        if classifier == "ovl":
            columns += ["vchan", "vthr", "vwin"]

        ### write columns to summary dat
        file_obj = open(summary_dat, "w")
        print >> file_obj, " ".join(columns)

        for datfile in datfiles:
            if classifier not in datfile: ### downselect based on classifier
                continue

            ### load only the required columns from datfile 
            output = idq.slim_load_datfile(datfile, skip_lines=0, columns=columns)

            ### convert to list of triggers
            triggers = [[output[c][i] for c in columns] for i in xrange(len(output["GPS"]))]
            triggers = [[float(l[0])]+l[1:] for l in triggers]

            ### window triggers with scisegs
            triggers = event.include( triggers, scisegs, tcent=0 ) ### tcent hard coded according to columns

            ### write surviving triggers to datfile
            for trigger in triggers:
                trigger = [str(trigger[0])] + trigger[1:]
                print >> file_obj, " ".join(trigger)

        file_obj.close()
        logger.info("Done.")
   
    #=============================================
    # generate *roc files
    #=============================================
    ### collect *dat files and merge into *roc files
    logger.info('generating *.roc files')
    roc_paths = idq.datfiles_to_roc(
        gpsstart,
        gpsstart + stride,
        columns=columns,
        classifiers=classifiers,
        basename=False,
        source_dir=this_sumdir,
#        source_dir=realtimedir,
        output_dir=this_sumdir,
        cluster_win=cluster_win,
        unsafe_win=unsafe_win,
        gw_thr=gw_thr,
        switch_snr_signif=classifiers_switch_snr_signif,
        )

    ### generate uniformly sampled ROC files (sampled 100 times)
    ### these are used for rank->FAP maps in the realtime job
    logger.info('generating *.uroc files')
    uroc_paths = []
    for (roc_path, classifier) in roc_paths:
        (uniform_ranks, uniform_ccln, uniform_cgch, tcln, tgch) = \
            idq_s_p.ROC_to_uniformROC(roc_path, num_samples=100)

        uniformROCfilename = roc_path[:-4] + '.uroc'  # suffix is for uniformly-sampled roc file

        uroc_paths.append((idq_s_p.rcg_to_ROC(
            uniformROCfilename,
            uniform_ranks,
            uniform_ccln,
            uniform_cgch,
            tcln,
            tgch,
            ), classifier, (tcln>=min_num_cln) and (tgch>=min_num_gch)))

    ### update uroc_cachefiles used in realtime job
    logger.info('updating *_uroc.cache files')
    for (uroc_path, classifier, include) in uroc_paths:
        if include or opts.force:
            logger.info('updating *_uroc.cache file for %s'%classifier)

            uroc_cachefilename = sumdir + '/' + classifier + '_uroc.cache'
            file = open(uroc_cachefilename, 'a')
            print >> file, uroc_path
            file.close()
        else:
           logger.info('not enough glitches or cleans were found to justify updating uroc.cache for %s'%classifier)

    #=============================================
    # generat ROC figures
    #=============================================
    logger.info('generating ROC figures')

    ### generate an individual ROC plot for each classifier
    roc_fig_paths = []
    for (ind, roc_path) in enumerate(roc_paths):
        (path, classifier) = roc_path
        try:
            figname = path[:-4] + '_roc'
            idq_s_p.ROC_to_ROC_plot([path],
                                    labels=[idq_s_p.classifier_labels(classifier)],
                                    colors=[idq_s_p.classifier_colors(classifier)],
                                    figure_name=figname, write=True)
            roc_fig_paths.append((figname + idq_s_p.fig_type,
                                 classifier))
        except:
            traceback.print_exc()
            roc_fig_paths.append(('', classifier))
            logger.info('WARNING: FAILED to build roc figure for '
                        + classifier)

    ### generate a combined ROC plot showing all classifiers
    try:
        figname = "%s/all-%d-%d_roc"%(this_sumdir,gpsstart,stride)
        roc_fig_path = figname + idq_s_p.fig_type
        idq_s_p.ROC_to_ROC_plot([l[0] for l in roc_paths],
                                labels=classifier_labels,
                                colors=classifier_colors,
                                figure_name=figname, write=True)
    except:
        traceback.print_exc()
        roc_fig_path = False
        logger.info('WARNING: FAILED to generate combined roc figure')
    logger.info('Done.')

    #=============================================
    # basic lists of channels, segments, etc
    #=============================================
    ### location of most recent vetolist:
    logger.info('finding pointers to OVL vetolists')
    vetolist_link = False
    try:
        vetolist_path = open(vetolist_cache, 'r').readlines()[-1].strip('\n')
        vetolist_link = this_sumdir + '/vetolist.eval'
        if os.path.lexists(vetolist_link):
            os.remove(vetolist_link)
        os.symlink(vetolist_path, vetolist_link)
    except:
        traceback.print_exc()
        logger.info('WARNING: FAILED to find most recent OVL vetolist')
    logger.info('Done')

    ### generate segments
    logger.info('generating OVL vetolist segments')
    vetolist_seg_path = False
    ovl_fdt = False
    try:
        vetolist_seg_path = this_sumdir + '/' + vetolist_path.split('/'
                )[-1] + '_FAP_<=' + str(FAP) + '.seg'
        (vetolist_segs, ovl_fdt) = idq.ovl.vetolist_to_segments(
            vetolist_path,
            gpsstart,
            gpsstart + stride,
            FAP,
            trg_dir=kwtrgdir,
            scisegs=False,
            output_filename=vetolist_seg_path,
            )
    except:
        traceback.print_exc()
        logger.info('WARNING: FAILED to generate OVL segments')
    logger.info('Done')

    #=============================================
    # algorithmic trending 
    #=============================================
    ### trending plots

    logger.info('generating trending plots')

    ### compute how far in time to look back
    lookbacktime = gpsstart - lookback * stride

    ### get stat summary files produced by earlier summary jobs
    stat_summary_files = idq.get_all_files_in_range(sumdir,
            lookbacktime, gpsstart + stride, pad=0, suffix='.stat')

    ### get the most recent stat file
    stat_summary_files.sort()
    latest_stat_summary_file = stat_summary_files[-1]

    ### generate trending plots for livetime, glitch and clean samples rates
    logger.info('generating summary statistic trending plots')
    try:
        stat_trends_plot_paths = idq_s_p.stat_to_trends_plot(
            stat_summary_files,
            classifiers,
            labels=classifier_labels,
            colors=classifier_colors,
            output_dir=this_sumdir,
            figure_basename='stat_trending',
            write=True,
            )
    except:
        traceback.print_exc()
        stat_trends_plot_paths = []
        logger.info('WARNING: FAILED to generate summary statistic trending plots.'
                    )
    logger.info('Done')

    ### get roc files
    roc_files_for_trending = idq.get_all_files_in_range(sumdir,
            lookbacktime, gpsstart + stride, pad=0, suffix='.roc')

    # define dictionary to hold them
    roc_files_dict = {}
    for classifier in classifiers:
        roc_files_dict[classifier] = []

    for file in roc_files_for_trending:
        classifier = file.split('/')[-1].split('-')[0]
        roc_files_dict[classifier].append(file)

    ### generate efficiency trending plot....
    logger.info('generating effciency trending plot')
    try:
        eff_trends_plot_paths = idq_s_p.ROC_to_eff_trends_plot(
            roc_files_dict,
            FAP,
            classifiers=classifiers,
            labels=classifier_labels,
            colors=classifier_colors,
            output_dir=this_sumdir,
            figure_basename='_trending_',
            write=True,
            )
    except:
        traceback.print_exc()
        eff_trends_plot_paths = []
        logger.info('WARNING: FAILED to generate efficiency trending plot.'
                    )
    logger.info('Done')

    #=============================================
    # channel/configuration specific performance
    #=============================================
    ### channel statistics
    logger.info('generating channel statistics')
    chanlist = []
    for classifier in classifiers:
        if classifier not in ['ovl']:
            logger.info('skipping ' + classifier)
            continue
        try:
            chanlist += idq.datfiles_to_chanlist(
                gpsstart,
                gpsstart + stride,
                columns=columns + ['vchan'],
                classifiers=[classifier],
                basename=False,
                source_dir=this_sumdir,
#                source_dir=realtimedir,
                output_dir=this_sumdir,
                cluster_win=cluster_win,
                unsafe_win=unsafe_win,
                gw_thr=gw_thr,
                switch_snr_signif=classifier
                    in classifiers_switch_snr_signif,
                )
        except:
            traceback.print_exc()
            chanlist.append(('', classifier))
            logger.info('WARNING: FAILED to generate chanlist  for '
                        + classifier)
    logger.info('Done')

    ### generate channel trending plots
    logger.info('generating channel performance trending plot')
    chanlist_trend = []
    for yvalue in ["rank", "eff", "fap"]:
        for classifier in classifiers:
            logger.info('%s %s channel performance trending'%(classifier, yvalue))

            if classifier not in ['ovl']:
                logger.info('skipping ' + classifier)
                continue

            try:
                figure_name = this_sumdir + '/' + classifier \
                    + '-%d-%d_channel_performance_%s_trends.png' \
                    % (lookbacktime, gpsstart + stride, yvalue)
                chan_perform_png = idq_s_p.chanlist_trending(
                    lookbacktime,
                    gpsstart + stride,
                    sumdir,
                    classifier=classifier,
                    figure_name=figure_name,
                    annotated=False,
                    yvalue=yvalue
                    )
                chanlist_trend.append((chan_perform_png, classifier))
            except:
                traceback.print_exc()
                chanlist_trend.append(('', classifier))
                logger.info('WARNING: FAILED to generated channel trending plot for '
                             + classifier)

    logger.info('Done')

    ### generate config statistics
    logger.info('generating config statistics')
    configlist = []
    for classifier in classifiers:
        if classifier != 'ovl':
            logger.info('skipping ' + classifier)
            continue
        try:
            configlist += idq.datfiles_to_configlist(
                gpsstart,
                gpsstart + stride,
                columns=columns + ['vchan', 'vthr', 'vwin'],
                classifiers=[classifier],
                basename=False,
                source_dir=this_sumdir,
#                source_dir=realtimedir,
                output_dir=this_sumdir,
                cluster_win=cluster_win,
                unsafe_win=unsafe_win,
                gw_thr=gw_thr,
                switch_snr_signif=classifier
                    in classifiers_switch_snr_signif,
                )
        except:
            traceback.print_exc()
            configlist.append((False, classifier, False))
            logger.info('WARNING: FAILED to generate config statistics for '
                         + classifier)
    logger.info('Done')

    #=============================================
    # mappings from rank to other statistics
    #=============================================

    # generate rank-->GPS maps
    logger.info('generating rank-->GPS maps')
    ranklist = []
    for classifier in classifiers:
        try:
            ranklist += idq.datfiles_to_ranklist(
                gpsstart,
                gpsstart + stride,
                columns=columns,
                classifiers=[classifier],
                basename=False,
                source_dir=this_sumdir,
#                source_dir=realtimedir,
                output_dir=this_sumdir,
                cluster_win=cluster_win,
                unsafe_win=unsafe_win,
                gw_thr=gw_thr,
                switch_snr_signif=classifier
                    in classifiers_switch_snr_signif,
                )
        except:
            traceback.print_exc()
            ranklist.append((False, classifier, False))
            logger.info('WARNING: FAILED to generate rank-->GPS map for '
                         + classifier)
    logger.info('Done')

    ### generate kde overlays (map rank -> p(r|g), p(r|c))
    logger.info('generating kde overlays')
    kde_fig_paths = []
    for (ind, roc_path) in enumerate(roc_paths):
        (path, classifier) = roc_path
        try:
            figname = path[:-4] + '_pwg_kde' + idq_s_p.fig_type
            fig = idq_s_p.ROC_to_pwg_kde_plot(path, write=False,
                    num_samples=kde_num_samples)
            idq_s_p.plt.savefig(figname)
            idq_s_p.plt.close(fig)
            kde_fig_paths.append((figname, classifier))
        except:
            traceback.print_exc()
            kde_fig_paths.append(('', classifier))
            logger.info('WARNING: FAILED to generate kdw plot for '
                        + classifier)
    logger.info('Done')


    ### glitch rates, glitch parameter histograms, pointers to the craziest glitches?

    ### plot residuals between old ROC curve and new ROC curve?

    #===============================================================================================
    # generate html summary page
    #===============================================================================================

    ### generate url for roc figure with all classifiers
    roc_url = False
    if roc_fig_path:
        roc_url = roc_fig_path.split(this_sumdir + '/')[-1]

    ### generate url for vetolist stuff
    vetolist_url = False
    if vetolist_link:
        vetolist_url = vetolist_link.split(this_sumdir + '/')[-1]
    if vetolist_seg_path:
        vetolist_seg_url = vetolist_seg_path.split(this_sumdir + '/')[-1]

    ### generate urls for figures for each classifier separately
    roc_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                for path in roc_fig_paths if path[0]]
    kde_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                for path in kde_fig_paths if path[0]]

    ### generate url for important channel lists
    chanlist_urls = [(path[-1].split(this_sumdir + '/')[-1], path[1])
                     for path in chanlist if path[0]]

    ### construct urls for trending plots
    stat_trends_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                        for path in stat_trends_plot_paths if path[0]]
    eff_trends_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                       for path in eff_trends_plot_paths if path[0]]

    ### construct url for stat summary file
    stat_summary_url = latest_stat_summary_file.split(this_sumdir + '/')[-1]

    ### construct url for channel performance file
    chanlist_trend_urls = [(path[0].split('/')[-1], path[1])
                           for path in chanlist_trend if path[0]]

    #=============================================
    # actuall build the html page
    #=============================================
    logger.info('generating html summary page')
    html_path = this_sumdir + '/' + str(gpsstart) + '_' + str(gpsstart
            + stride) + '-summary.html'
    generate_html(
        html_path,
        gpsstart=gpsstart,
        gpsstop=gpsstart + stride,
        gwchannel=gwchannel,
        roc_url=roc_url,
        vetolist_url=vetolist_url,
        roc_urls=roc_urls,
        kde_urls=kde_urls,
        FAP=FAP,
        ovl_segments=vetolist_seg_url,
        ovl_fdt=ovl_fdt,
        chanlist_urls=chanlist_urls,
        stat_trends_urls=stat_trends_urls,
        eff_trends_urls=eff_trends_urls,
        stat_file=stat_summary_url,
        chan_trends_urls=chanlist_trend_urls,
        )

    logger.info('Done')

    ### soft link html page to index.html
    index_html = "%s/index.html"%this_sumdir
    logger.info("soft linking %s -> %s"%(html_path, index_html))
    os.system("rm %s"%(index_html))
    os.system("ln -s %s %s"%(html_path, index_html))

    # update symbolic link
# ....iif os.path.lexists(symlink_path):
# ........os.remove(symlink_path)
# ....os.symlink(html_path, symlink_path)

    ### continue onto the next stride
    gpsstart += stride

#===================================================================================================
### (re)write index.html for sumdir once all pages are built
logger.info("writting top-level index.html page to point to individual directories")

index_html = "%s/index.html"%(sumdir)

file_obj = open(index_html, "w")
print >> file_obj, "<body>"

print >> file_obj, "<h1>iDQ summary pages Directory</h1>"

sumdirs = [l.strip("index.html") for l in sorted(glob.glob("%s/*_*/index.html"%sumdir))]
sumdirs.reverse()

for this_sumdir in sumdirs:
    this_sumdir = this_sumdir.strip("/").split("/")[-1]
    print >> file_obj, "<p>information about iDQ between %s and %s : "%tuple(this_sumdir.split("_"))
    print >> file_obj, "<a href=\"%s/\">here</a>"%(this_sumdir)
    print >> file_obj, "</p>"

print >> file_obj, "<hr />"

### print timestamp
c_time = time.localtime()
tzname = time.tzname[0]
print >> file_obj, '<p>last updated %d:%d:%d %s %2d/%2d/%4d </p>' % (
    c_time.tm_hour,
    c_time.tm_min,
    c_time.tm_sec,
    tzname,
    c_time.tm_mon,
    c_time.tm_mday,
    c_time.tm_year,
    )

print >> file_obj, "</body>"

file_obj.close()
