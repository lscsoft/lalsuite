# Copyright (C) 2015 Lindy Blackburn, Reed Essick
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

from lal import gpstime

from laldetchar.idq import idq
from laldetchar.idq import idq_summary_plots as isp
#from laldetchar.idq import reed as idq
#from laldetchar.idq import reed_summary_plots as isp 
from laldetchar.idq import event
from laldetchar.idq import calibration

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table

from glue import markup ### an easy way to write html pages. THIS MAY BE REPLACED???

from laldetchar import git_version

#===================================================================================================

__author__ = 'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = """This program generates summary html pages from iDQ pipeline output. The summmary pages provide variety of diagnostic and interpretational plots and data."""

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

parser.add_option('-t', '--trending', default=[], action='append', type='string', help='Number of seconds to look back for trending plots')

parser.add_option('-l', '--log-file', default='idq_summary.log', type='string', help='log file')

parser.add_option("", "--ignore-science-segments", default=False, action="store_true")

parser.add_option("", "--dont-cluster", default=False, action="store_true")

parser.add_option('-f','--force',default=False, action='store_true', help="forces *uroc cache file to be updated, even if we have no data. Use with caution.")

parser.add_option("", "--no-robot-cert", default=False, action="store_true")

parser.add_option("", "--no-html", default=False, action="store_true")

parser.add_option('', '--FAPthr', default=[], action="append", type='float', help='check calibration at this FAP value. This argument can be supplied multiple times to check multiple values.')

(opts, args) = parser.parse_args()

if opts.lookback != "infinity":
    lookback = int(opts.lookback)

if not opts.trending:
    opts.trending = ['infinity']

for ind, trend in enumerate(opts.trending):
    if trend != "infinity":
        opts.trending[ind] = int(trend)











### these should be set somewhere else!
html_width=500
html_height=500

thumbnail_width=100
thumbnail_height=100

num_trend_bins = 20

cln_linestyle='dashed' ### FIXME: pull from config file?
gch_linestyle='solid' 














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

mainidqdir = config.get('general', 'idqdir') ### get the main directory where idq pipeline is going to be running.

usertag = config.get('general', 'usertag')
if usertag:
    usertag = "_%s"%usertag

ifo = config.get('general', 'ifo')

channels = config.get('general','selected-channels')
unsafechannels = config.get('general', 'unsafe-channels')

gwchannel = config.get('general', 'gwchannel')
gch_kwsignif_thr = config.getfloat('general', 'gw_kwsignif_thr')
cln_kwsignif_thr = config.getfloat('realtime', 'clean_threshold')

#========================
# which classifiers
#========================
### ensure we have a section for each classifier and fill out dictionary of options
classifiersD, mla, ovl = idq.config_to_classifiersD( config )

### get combiners information and add these to classifiersD
combinersD, referenced_classifiers = idq.config_to_combinersD( config )
for combiner, value in combinersD.items():
    classifiersD[combiner] = value

classifiers = sorted(classifiersD.keys())

#if mla:
#    ### reading parameters from config file needed for mla
#    auxmvc_coinc_window = config.getfloat('build_auxmvc_vectors','time-window')
#    auxmc_gw_signif_thr = config.getfloat('build_auxmvc_vectors','signif-threshold')
#    auxmvc_selected_channels = config.get('general','selected-channels')
#    auxmvc_unsafe_channels = config.get('general','unsafe-channels')

colors = dict( zip( classifiers, 'r b g c m k y'.split() ) )
if len(colors) < len(classifiers):
    raise ValueError("ran out of colors for classifiers!")

#========================
# realtime 
#========================
realtimedir = config.get('general', 'realtimedir')### output directory for realtime predictions

dat_columns = config.get('realtime', 'dat_columns').split()
columns = dat_columns[:]
for c in "GPS i rank".split():
    if c not in columns:
        columns.append( c )
#columns = list(set( dat_columns + "GPS i rank".split() ) )

#========================
# summary jobs
#========================
summarydir = config.get('general', 'summarydir')
trenddir = "%s/trending"%(summarydir) ### used to store trending information

stride = config.getint('summary', 'stride')
delay = config.getint('summary', 'delay')

emaillist = config.get('warnings', 'summary') 
errorthr = config.getfloat('warnings', 'summary_errorthr')

kde_nsamples = config.getint('summary', 'kde_num_samples')
kde = np.linspace(0, 1, kde_nsamples) ### uniformly spaced ranks used to sample kde

faircoin = kde[:] ### used for ROC plots

cluster_key = config.get('summary', 'cluster_key')
if cluster_key not in columns:
    columns.append( cluster_key )
cluster_win = config.getfloat('summary', 'cluster_win')

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

    wait = gpsstart + stride + delay - idq.nowgps()
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

    ### dictionary to hold onto filenames we want to keep
    files = {'segs':{}, 'dat':{}, 'kde':{}, 'roc':{}, 'chn':{}}

    ### dictionary to hold data
    data = {}

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
            seg_xml_file = idq.segment_query(config, gpsstart - lookback , gpsstart + stride, url=segdb_url)

            ### write seg_xml_file to disk
            lsctables.use_in(ligolw.LIGOLWContentHandler)
            xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

            ### science segments xml filename
            seg_file = idq.segxml(output_dir, "_%s"%dq_name, gpsstart - lookback , lookback+stride)
            logger.info('writing science segments to file : '+seg_file)
            ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))
            files['segs']['sci'] = seg_file

            (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name) ### read in segments from xml file

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
    idqsegs = idq.get_idq_segments(realtimedir, gpsstart-lookback, gpsstart+stride, suffix='.dat')

    logger.info('taking intersection between science segments and idq segments')
    idqsegs = event.andsegments( [scisegs, idqsegs] )

    ### write segment file
    if opts.ignore_science_segments:
        idqseg_path = idq.idqsegascii(output_dir, '', gpsstart - lookback, lookback+stride)
    else:
        idqseg_path = idq.idqsegascii(output_dir, '_%s'%dq_name, gpsstart - lookback, lookback+stride)
    f = open(idqseg_path, 'w')
    for seg in idqsegs:
        print >> f, seg[0], seg[1]
    f.close()

    files['segs']['idq'] = idqseg_path

    #===============================================================================================
    # find data
    #===============================================================================================
    ### find all *dat files, bin them according to classifier
    logger.info('finding all *dat files')
    datsD = defaultdict( list )
    for dat in idq.get_all_files_in_range(realtimedir, gpsstart-lookback, gpsstart+stride, pad=0, suffix='.dat' ):
        datsD[idq.extract_dat_name( dat )].append( dat )

    ### throw away any un-needed files
    for key in datsD.keys():
        if key not in classifiers:
            datsD.pop(key)
        else: ### throw out files that don't contain any science time
            datsD[key] = [ dat for dat in datsD[key] if event.livetime(event.andsegments([idqsegs, [idq.extract_start_stop(dat, suffix='.dat')]])) ]

    #===============================================================================================
    # build plots
    #===============================================================================================
    fignames = {'linroc':{}, 'roc':{}, 'hst':{}, 'kde':{}, 'L':{}, 'trg':{}, 'bwh':{}, 'fap':{}, 'fapUL':{}}
    ### roc -> Receiver Operating Characteristic curves
    ### hst -> raw histograms of distributions
    ### kde -> kde smoothed histograms of distributions
    ### trg -> histograms of gch features
    ### bwh -> bit-word histograms

    ### axes for overlay plots
    roc_figax = None
    hst_figax = None
    kde_figax = None
    L_figax = None

    ### bit-word storage
    gch_bitword = dict( (FAP, {}) for FAP in opts.FAPthr )
    cln_bitword = dict( (FAP, {}) for FAP in opts.FAPthr )

    #====================
    # build plots for classifiers separately and add to overlays
    # only plots that depend on this stride alone!
    #====================
    for classifier in classifiers:
        color = colors[classifier]

        ### write list of dats to cache file
        cache = idq.cache(output_dir, classifier, "_datcache%s"%usertag)
        logger.info('writing list of dat files to %s'%cache)
        f = open(cache, 'w')
        for dat in datsD[classifier]:
            print >>f, dat
        f.close()

        files['dat'][classifier] = cache

        logger.info('building plots for %s'%classifier)

        ### load dat files
        output = idq.slim_load_datfiles(datsD[classifier], skip_lines=0, columns=columns)

        ### filter times by scisegs -> keep only the ones within scisegs
        output = idq.filter_datfile_output( output, idqsegs )

        if not opts.dont_cluster:
            output = idq.cluster_datfile_output( output, cluster_key=cluster_key, cluster_win=cluster_win)
            cluster_dat = idq.dat(output_dir, classifier, ifo, "clustered", usertag, gpsstart-lookback, lookback+stride)
            logger.info('  writing %s'%cluster_dat)
            idq.output_to_datfile( output, cluster_dat )
        else:
            cluster_dat = idq.dat(output_dir, classifier, ifo, "unclustered", usertag, gpsstart-lookback, lookback+stride)
            logger.info('  writing %s'%cluster_dat)
            idq.output_to_datfile( output, cluster_dat )


        ### compute rcg from output
        r, c, g = idq.dat_to_rcg( output )
        dc, dg = idq.rcg_to_diff( c, g ) ### get the numbers at each rank

        data[classifier] = {'num_cln':c[-1], 'num_gch':g[-1]}

        ### dump into roc file
        roc = idq.roc(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writting %s'%roc)
        idq.rcg_to_file(roc, r, c, g)

        files['roc'][classifier] = roc

        ### generate ROC plot
        rocfig = isp.rocfig(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['roc'][classifier] = rocfig ### store for reference
        logger.info('  plotting %s'%rocfig)

        fig, ax = isp.rcg_to_rocFig(c, g, color=color, label=classifier)
        ax.plot(faircoin, faircoin, 'k--')
        ax.legend(loc='best')

        fig.savefig(rocfig)
        isp.close(fig)

        rocfig = isp.rocfig(output_dir, classifier, ifo, "%s_lin"%usertag, gpsstart-lookback, lookback+stride)
        fignames['linroc'][classifier] = rocfig
        logger.info('  plotting %s'%rocfig)
        ax.set_xscale('linear')
        ax.set_xlim(xmin=0, xmax=1)
        fig.savefig(rocfig)
        isp.close(fig)

        roc_figax = isp.rcg_to_rocFig(c, g, color=color, label=classifier, figax=roc_figax) ### for overlay plot

        ### generate histogram, cdf 
        histfig = isp.histfig(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['hst'][classifier] = histfig ### store for reference
        logger.info('  plotting %s'%histfig)

        if np.any(dc):
            fig, axh, axc = isp.stacked_hist(r, dc, color=color, linestyle=cln_linestyle, label="%s cln"%classifier)
            hst_figax = isp.stacked_hist(r, dc, color=color, linestyle=cln_linestyle, label="%s cln"%classifier, figax=hst_figax) ### for overlay plots
        if np.any(dg):
            fig, axh, axc = isp.stacked_hist(r, dg, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=(fig, axh, axc))
            hst_figax = isp.stacked_hist(r, dg, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=hst_figax)
        if np.any(dc) or np.any(dg):
            axh.set_xlim(xmin=0, xmax=1)
            axc.set_xlim(xmin=0, xmax=1)
            axc.legend(loc='best')
            fig.savefig(histfig)
            isp.close(fig)

        ### compute kde estimates
        logger.info('  computing kde_pwg estimates')
        kde_cln = idq.kde_pwg( kde, r, dc )
        kde_gch = idq.kde_pwg( kde, r, dg )

        ### write kde points to file
        kde_cln_name = idq.kdename(output_dir, classifier, ifo, "_cln%s"%usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writing %s'%kde_cln_name)
        np.save(event.gzopen(kde_cln_name, "w"), (kde, kde_cln))

        kde_gch_name = idq.kdename(output_dir, classifier, ifo, "_gch%s"%usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writing %s'%kde_gch_name)
        np.save(event.gzopen(kde_gch_name, "w"), (kde, kde_gch))

        files['kde'][classifier] = {'cln':kde_cln_name, 'gch':kde_gch_name}

        ### generate kde pdf, cdf (above)
        kdefig = isp.kdefig(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['kde'][classifier] = kdefig ### store for reference
        logger.info('  plotting %s'%kdefig)

        fig, axh, axc = isp.stacked_kde(kde, kde_cln, color=color, linestyle=cln_linestyle, label="%s cln"%classifier, figax=None)
        fig, axh, axc = isp.stacked_kde(kde, kde_gch, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=(fig, axh, axc))
        axc.legend(loc='best')

        fig.savefig(kdefig)
        isp.close(fig)

        kde_figax = isp.stacked_kde(kde, kde_cln, color=color, linestyle=cln_linestyle, label="%s cln"%classifier, figax=kde_figax) ### for overlay plot
        kde_figax = isp.stacked_kde(kde, kde_gch, color=color, linestyle=gch_linestyle, label="%s gch"%classifier, figax=kde_figax)

        ### likelihood ratio
        Lfig = isp.Lfig(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        fignames['L'][classifier] = Lfig
        logger.info('  plotting %s'%Lfig)

        fig, axh, axc = isp.stacked_L(kde, kde_cln, kde_gch, color=color, linestyle='solid', label=classifier, figax=None)
        axc.legend(loc='best')

        fig.savefig(Lfig)
        isp.close(fig)

        L_figax = isp.stacked_L(kde, kde_cln, kde_gch, color=color, linestyle='solid', label=classifier, figax=L_figax) ### for overlay plot

        ### figure out which gch and cln are removed below FAP thresholds
        rankthrs = []
        for FAP in opts.FAPthr:
            rankthr = min(r[c<=FAP*c[-1]]) ### smallest rank below FAP
            rankthrs.append( rankthr )
            cln_bitword[FAP][classifier], gch_bitword[FAP][classifier] = idq.bin_by_rankthr(rankthr, output, columns=dat_columns)
                
        ### generate histograms of glitch parameters
        for column in dat_columns:
            paramhist = isp.histfig(output_dir, classifier, ifo, "%s_%s"%(usertag, column), gpsstart-lookback, lookback+stride)
            fignames['trg'][(classifier, column)] = paramhist
            logger.info('  plotting %s'%paramhist)

            figax = None
            for rankthr, FAP in zip(rankthrs, opts.FAPthr):
                thr = min(r[c<=FAP*c[-1]]) ### smallest rank below FAP
                x = [float(output[column][i]) for i in xrange(len(output[column])) if output['rank'][i] >= thr]
                if x: ### will break if x is empty
                    figax = isp.stacked_hist( x, np.ones_like(x), color=None, label="$FAP\leq%E$"%FAP, bmin=min(x), bmax=max(x), figax=figax)
            if figax: ### only save if we plotted something...
                fig, axh, axc = figax
                axh.set_xlabel(column)
                axc.set_xlim(axh.get_xlim())
                axc.legend(loc='best')
                fig.savefig(paramhist)
                isp.close(fig)

    #====================
    # save overlays for this stride
    #====================

    ### roc overlay
    rocfig = isp.rocfig(output_dir, "overlay", ifo, usertag, gpsstart-lookback, lookback+stride)
    fignames['roc']["overlay"] = rocfig ### store for reference
    logger.info('  plotting %s'%rocfig)
    fig, ax = roc_figax
    ax.plot(faircoin, faircoin, 'k--')
    ax.legend(loc='best')
    fig.savefig(rocfig)
    isp.close(fig)

    rocfig = isp.rocfig(output_dir, "overlay", ifo, "%s_lin"%usertag, gpsstart-lookback, lookback+stride)
    fignames['linroc']['overlay'] = rocfig
    logger.info('  plotting %s'%rocfig)
    ax.set_xscale('linear')
    ax.set_xlim(xmin=0, xmax=1)
    fig.savefig(rocfig)
    isp.close(fig)

    ### histogram overlay
    histfig = isp.histfig(output_dir, "overlay", ifo, usertag, gpsstart-lookback, lookback+stride)
    fignames['hst']["overlay"] = histfig ### store for reference
    if hst_figax:
        logger.info('  plotting %s'%histfig)
        fig, axh, axc = hst_figax
        axh.set_xlim(xmin=0, xmax=1)
        axc.set_xlim(xmin=0, xmax=1)
        axc.legend(loc='best')
        fig.savefig(histfig)
        isp.close(fig)

    ### kde overlay
    kdefig = isp.kdefig(output_dir, "overlay", ifo, usertag, gpsstart-lookback, lookback+stride)
    fignames['kde']["overlay"] = kdefig ### store for reference
    logger.info('  plotting %s'%kdefig)
    fig, axh, axc = kde_figax
    axc.legend(loc='best')
    fig.savefig(kdefig)
    isp.close(fig)

    ### L overlay
    Lfig = isp.Lfig(output_dir, "overlay", ifo, usertag, gpsstart-lookback, lookback+stride)
    fignames['L']["overlay"] = Lfig
    logger.info('  plotting %s'%Lfig)
    fig, axh, axc = L_figax
    axc.legend(loc='best')
    fig.savefig(Lfig)
    isp.close(fig)

    #====================
    # bitword histograms
    #====================

    ### gch bitword histogram
    bitfig = isp.bitfig(output_dir, ifo, "%s_gch"%usertag, gpsstart-lookback, lookback+stride)
    fignames['bwh']['gch'] = bitfig
    logger.info('  plotting %s'%bitfig)

    figax=None
    for FAP in opts.FAPthr:
        figax = isp.bitword( gch_bitword[FAP], classifiers, label="FAP$\leq$%.3e"%FAP, figax=figax )
    if figax:
        fig, ax = figax
        ax.legend(loc='upper center')

        fig.savefig(bitfig)
        isp.close(fig)

    ### cln bitword histogram
    bitfig = isp.bitfig(output_dir, ifo, "%s_cln"%usertag, gpsstart-lookback, lookback+stride)
    fignames['bwh']['cln'] = bitfig
    logger.info('  plotting %s'%bitfig)

    figax=None
    for FAP in opts.FAPthr:
        figax = isp.bitword( cln_bitword[FAP], classifiers, label="FAP$\leq$%.3e"%FAP, figax=figax )
    if figax:
        fig, ax = figax
        ax.legend(loc='upper center')

        fig.savefig(bitfig)
        isp.close(fig)

    #====================
    # fap calibration
    #====================
    logger.info('checking FAP calibration')
    ### find all *dat files, bin them according to classifier
    logger.info('  finding all *fap*npy.gz files')
    fapsD = defaultdict( list )
    for npygz in idq.get_all_files_in_range(realtimedir, gpsstart-lookback, gpsstart+stride, pad=0, suffix='.npy.gz' ):
        if "fap" in npygz:
            fapsD[idq.extract_fap_name( npygz )].append( npygz )

    ### throw away any un-needed files
    for key in fapsD.keys():
        if key not in classifiers:
            fapsD.pop(key)
        else: ### throw out files that don't contain any science time
            fapsD[key] = [ fap for fap in fapsD[key] if event.livetime(event.andsegments([idqsegs, [idq.extract_start_stop(fap, suffix='.npy.gz')]])) ]

    ### check calibration and build figures
    alerts = {}
    fap_figax = None
    fapUL_figax = None
    for classifier in classifiers:

        _times, timeseries = idq.combine_ts(fapsD[classifier], n=2) ### read in time-series

        times = []
        faps = []
        fapsUL = []
        for t, ts in zip(_times, timeseries):
            _t, _ts = idq.timeseries_in_segments(t, ts, idqsegs)
            if len(_ts):
                times.append( _t )
                faps.append( _ts[0] )
                fapsUL.append( _ts[1] )

        ### check point estimate calibration
        _, deadtimes, statedFAPs, errs = calibration.check_calibration(idqsegs, times, faps, opts.FAPthr)
#        errs = np.array([ d/F - 1.0 for d, F in zip(deadtimes, statedFAPs) if F ])

        ### check UL estimate calibration
        _, deadtimesUL, statedFAPsUL, errsUL = calibration.check_calibration(idqsegs, times, fapsUL, opts.FAPthr)
#        errsUL = np.array([ d/F - 1.0 for d, F in zip(deadtimesUL, statedFAPs) if F ])

        ### write output file
        calib_check = idq.calib_check(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writing %s'%calib_check)

        file_obj = open(calib_check, "w")
        print >> file_obj, "livetime = %.3f"%event.livetime(idqsegs)
        for FAPthr, deadtime, statedFAP, err, deadtimeUL, statedFAPUL, errUL in zip(opts.FAPthr, deadtimes, statedFAPs, errs, deadtimesUL, statedFAPsUL, errsUL):
                print >> file_obj, calibration.report_str%(FAPthr, statedFAP, deadtime, err , statedFAPUL, deadtimeUL, errUL )
        file_obj.close()

        if np.any(np.abs(errs) > errorthr) or np.any(np.abs(errsUL) > errorthr):
            alerts[classifier] = calib_check

        ### build pt. estimate figures
        color = colors[classifier]
        fig, ax = isp.calibration_scatter( deadtimes, statedFAPs, color=color, label=classifier, figax=None)
        figname = isp.calibfig(output_dir, ifo, classifier, usertag, gpsstart-lookback, lookback+stride)
        logger.info('  plotting %s'%figname)
        fignames['fap'][classifier] = figname
        ax.plot(faircoin, faircoin, 'k--')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc='best')
        fig.savefig(figname)
        isp.close(fig)

        ### for overlay
        fap_figax = isp.calibration_scatter( deadtimes, statedFAPs, color=color, label=classifier, figax=fap_figax)

        ### build upper-limit figures
        fig, ax = isp.calibration_scatter( deadtimesUL, statedFAPsUL, color=color, label=classifier, figax=None)
        figname = isp.calibfig(output_dir, ifo, classifier, "%s_UL"%usertag, gpsstart-lookback, lookback+stride)
        logger.info('  plotting %s'%figname)
        fignames['fapUL'][classifier] = figname
        ax.plot(faircoin, faircoin, 'k--')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc='best')
        fig.savefig(figname)
        isp.close(fig)

        fapUL_figax = isp.calibration_scatter( deadtimesUL, statedFAPsUL, color=color, label=classifier, figax=fapUL_figax)

    ### save overlays
    fig, ax = fap_figax
    figname = isp.calibfig(output_dir, ifo, 'overlay', usertag, gpsstart-lookback, lookback+stride)
    logger.info('  plotting %s'%figname)
    fignames['fap']['overlay'] = figname
    ax.plot(faircoin, faircoin, 'k--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    fig.savefig(figname)
    isp.close(fig)

    fig, ax = fapUL_figax
    figname = isp.calibfig(output_dir, ifo, 'overlay', "%s_UL"%usertag, gpsstart-lookback, lookback+stride)
    logger.info('  plotting %s'%figname)
    fignames['fapUL']['overlay'] = figname
    ax.plot(faircoin, faircoin, 'k--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='best')
    fig.savefig(figname)
    isp.close(fig)

    ### send alerts if needed
    if alerts: ### there are some interesting files
        alerts_keys = sorted(alerts.keys())
        alerts_keys_str = " ".join(alerts_keys)
        logger.warning('WARNING: found suspicous historical calibrations for : %s'%alerts_keys_str )
        if emaillist:
            email_cmd = "echo \"calibration check summary files are attached for: %s\" | mailx -s \"%s idq%s calibration warning in laldetchar-idq-summary\" %s \"%s\""%(alerts_keys_str, ifo, usertag, " ".join("-a \"%s\""%alerts[key] for key in alerts_keys), emaillist)
            logger.warning("  %s"%email_cmd)
            exit_code = os.system( email_cmd )
            if exit_code:
                logger.warning("WARNING: failed to send email!")

    #=============================================
    # compute channel performance for flavor=="ovl"
    #=============================================
    for classifier in classifiers:
        if classifiersD[classifier]['flavor'] != "ovl": ### channel information only available from flavor=='ovl'
            continue

        ### load dat files
        output = idq.slim_load_datfiles(datsD[classifier], skip_lines=0, columns='GPS i vchan'.split()+[cluster_key]) ### specific columns known to be in ovl dat files

        ### filter times by scisegs -> keep only the ones within scisegs
        output = idq.filter_datfile_output( output, idqsegs )
        if not opts.dont_cluster:
            output = idq.cluster_datfile_output( output, cluster_key=cluster_key, cluster_win=cluster_win)
        
        ### don't re-write clustered dat files here...

        ### extract channel performances
        performance, num_gch, num_cln = idq.dat_to_perf( output )

        ### write to file!
        chanperf = idq.chan_perf(output_dir, classifier, ifo, usertag, gpsstart-lookback, lookback+stride)
        logger.info('  writing : %s'%chanperf)
        idq.perf_to_file( performance, num_gch, num_cln, chanperf )

        files['chn'][classifier] = chanperf

    #===============================================================================================
    # trending plots!
    # data depends on previous strides as well
    # we search for and plot data from these periods binned automatically into shorter-scale averages
    #===============================================================================================
    logger.info('Begin: Trending')
    if not os.path.exists(trenddir):
        os.makedirs(trenddir)

    trending = {'sci':[], 
                'idq':[],
                'seg':[],
                'trg_rate':dict( [('overlay',[])] + [(classifier, []) for classifier in classifiers] ), 
                'effatfap':dict( [('overlay',[])] + [(classifier, []) for classifier in classifiers] ), 
                'rank':dict( (classifier, []) for classifier in classifiers ),
                'eff':dict( (classifier, []) for classifier in classifiers ),
                'fap':dict( (classifier, []) for classifier in classifiers )
                }

    for trend in opts.trending:
        ### compute the trending range
        if trend == "infinity":
            trend_start = global_start
        else:
            trend_start = gpsstart - trend

        logger.info('  Begin: %d - %d'%(trend_start, gpsstart+stride) )

        bindatas = [] ### store the data for plotting

        ### compute bins for this trend stride
        bins = np.linspace(trend_start, gpsstart+stride, num_trend_bins+1).astype(int)
        for s, e in zip(bins[:-1], bins[1:]):
            logger.info('    Begin: %d - %d'%(s , e) )

            bindata = {}

            bindata['range'] = (s, e)

            output_trenddir = "%s/%d_%d/"%(trenddir, s, e)

            if not os.path.exists(output_trenddir): ### we need to make this trending information
                logger.info('    making %s and all required data'%output_trenddir)
                os.makedirs(output_trenddir)

                #=====================================
                # science segments
                # we query the segdb right now, although that latency may be an issue...
                #=====================================
                if opts.ignore_science_segments:
                    logger.info('    analyzing data regardless of science segements')
                    scisegs = [[gpsstart-lookback, gpsstart+stride]] ### set segs to be this stride range
                    coveredsegs = [[gpsstart-lookback, gpsstart+stride]] ### set segs to be this stride range

                else:
                    logger.info('    querrying science segments')

                    try:
                        ### this returns a string
                        seg_xml_file = idq.segment_query(config, s, e, url=segdb_url)

                        ### write seg_xml_file to disk
                        lsctables.use_in(ligolw.LIGOLWContentHandler)
                        xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

                        ### science segments xml filename
                        seg_file = idq.segxml(output_trenddir, "_%s"%dq_name, s , e-s)
                        logger.info('    writing science segments to file : '+seg_file)
                        ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))

                        (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name) ### read in segments from xml file
 
                        bindata['sci'] = scisegs
 
                    except Exception as er:
                        traceback.print_exc()
                        logger.info('ERROR: segment generation failed. Skipping this trending period.')

                        if opts.force: ### we are require successful training or else we want errors
                            logger.info(traceback.print_exc())
                            raise er
                        else: ### we don't care if any particular training job fails
                            continue

                logger.info('    finding idq segments')
                idqsegs = idq.get_idq_segments(realtimedir, s, e, suffix='.dat')

                logger.info('    taking intersection between science segments and idq segments')
                idqsegs = event.andsegments( [scisegs, idqsegs] )

                bindata['idq'] = idqsegs

                ### write segment file
                if opts.ignore_science_segments:
                    idqseg_path = idq.idqsegascii(output_trenddir, '', s, e-s)
                else:
                    idqseg_path = idq.idqsegascii(output_trenddir, '_%s'%dq_name, s, e-s)
                f = open(idqseg_path, 'w')
                for seg in idqsegs:
                    print >> f, seg[0], seg[1]
                f.close()

                #=====================================
                # find data
                #=====================================
                ### find all *dat files, bin them according to classifier
                logger.info('    finding all *dat files')
                datsD = defaultdict( list )
                for dat in idq.get_all_files_in_range(realtimedir, s, e, pad=0, suffix='.dat' ):
                    datsD[idq.extract_dat_name( dat )].append( dat )

                ### throw away any un-needed files
                for key in datsD.keys():
                    if key not in classifiers:
                        datsD.pop(key)
                    else: ### throw out files that don't contain any science time
                        datsD[key] = [ dat for dat in datsD[key] if event.livetime(event.andsegments([idqsegs, [idq.extract_start_stop(dat, suffix='.dat')]])) ]

                bindata['rcg'] = {}
                bindata['chn'] = {}
                for classifier in classifiers:
                    ### write list of dats to cache file
                    cache = idq.cache(output_trenddir, classifier, "_datcache%s"%usertag)
                    logger.info('  writing list of dat files to %s'%cache)
                    f = open(cache, 'w')
                    for dat in datsD[classifier]:
                        print >>f, dat
                    f.close()

                    output = idq.slim_load_datfiles(datsD[classifier], skip_lines=0, columns=columns)
                    output = idq.filter_datfile_output( output, idqsegs )

                    if not opts.dont_cluster:
                        output = idq.cluster_datfile_output( output, cluster_key=cluster_key, cluster_win=cluster_win)
                        cluster_dat = idq.dat(output_trenddir, classifier, ifo, "clustered", usertag, s, e-s)
                        logger.info('  writing %s'%cluster_dat)
                        idq.output_to_datfile( output, cluster_dat )
                    else:
                        cluster_dat = idq.dat(output_trenddir, classifier, ifo, "unclustered", usertag, gpsstart-lookback, lookback+stride)
                        logger.info('  writing %s'%cluster_dat)
                        idq.output_to_datfile( output, cluster_dat )

                    r, c, g = idq.dat_to_rcg( output )
                    bindata['rcg'][classifier] = (r, c, g, c[-1], g[-1])

                    ### save roc file
                    roc = idq.roc(output_trenddir, classifier, ifo, usertag, s, e-s)
                    logger.info('  writing %s'%roc)
                    idq.rcg_to_file(roc, r, c, g)

                    ### build/save chan-perf files
                    for classifier in classifiers:
                        if classifiersD[classifier]['flavor'] != 'ovl':
                            continue

                        ### load dat files
                        output = idq.slim_load_datfiles(datsD[classifier], skip_lines=0, columns='GPS i vchan'.split()+[cluster_key]) ### specific columns known to be in ovl dat files

                        ### filter times by scisegs -> keep only the ones within scisegs
                        output = idq.filter_datfile_output( output, idqsegs )
                        if not opts.dont_cluster:
                            output = idq.cluster_datfile_output( output, cluster_key=cluster_key, cluster_win=cluster_win)

                        ### don't re-write clustered dat here...

                        ### extract channel performances!
                        performance, num_gch, num_cln = idq.dat_to_perf( output )

                        ### write to file!
                        chanperf = idq.chan_perf(output_trenddir, classifier, ifo, usertag, s, e-s)
                        logger.info('  writing : %s'%chanperf)
                        idq.perf_to_file( performance, num_gch, num_cln, chanperf )
                        bindata['chn'][classifier] = (performance, num_gch, num_cln)

            else: ### information already exists, just read it in
                logger.info('    %s exists. Reading in existing data'%output_trenddir)

                if not opts.ignore_science_segments:
                    seg_file = seg_file = idq.segxml(output_trenddir, "_%s"%dq_name, s , e-s)
                    logger.info('    reading %s'%seg_file)
                    bindata['sci'] = idq.extract_dq_segments(seg_file, dq_name)[0]

                    idq_file = idq.idqsegascii(output_trenddir, '_%s'%dq_name, s, e-s)
                else:
                    idq_file = idq.idqsegascii(output_trenddir, '', s, e-s)
               
                logger.info('    reading %s'%idq_file)
                try:
                    bindata['idq'] = np.loadtxt( idq_file )
                except IOError: ### thrown if file was empty
                    bindata['idq'] = []
 
                bindata['rcg'] = {}
                bindata['chn'] = {}
                for classifier in classifiers:
                    roc = idq.roc(output_trenddir, classifier, ifo, usertag, s, e-s)
                    logger.info('    reading %s'%roc)
                    bindata['rcg'][classifier] = idq.file_to_rcg( roc )
 
                    if classifiersD[classifier]['flavor'] == 'ovl':
                        chanperf = idq.chan_perf(output_trenddir, classifier, ifo, usertag, s, e-s)
                        logger.info('    reading %s'%chanperf)
                        bindata['chn'][classifier] = idq.file_to_perf( chanperf )

            ### store data for plotting!
            bindatas.append( bindata )

        #=========================================
        # build figures!
        #=========================================

        ### idq_seg rates
        ranges = []
        values = []
        for bindata in bindatas:
            s, e = bindata['range']
            segs = bindata['idq']
            uptime = event.livetime( segs )
            ranges.append( (s, e) )
            values.append( 1.0*uptime/(e-s) )

        fig, ax = isp.rates( ranges, values, color='b', label='idq-segs', figax=None )
        ax.set_ylim(ymin=0, ymax=1)
        ax.set_ylabel('duty cycle')
        ax.legend(loc='lower left')

        figname = isp.ratefig(output_dir, ifo, "%s_idq-segments"%(usertag), trend_start, gpsstart+stride-trend_start)
        logger.info('  plotting %s'%figname)
        fig.savefig(figname)
        isp.close(fig)

        trending['idq'].append( figname )

        ### plot for overlay
        uptime_overlay_figax = isp.rates( ranges, values, color='b', label='idq-segs', figax=None )

        if not opts.ignore_science_segments:
            ranges = []
            values = []
            for bindata in bindatas:
                s, e = bindata['range']
                segs = bindata['sci']
                uptime = event.livetime( segs )
                ranges.append( (s, e) )
                values.append( 1.0*uptime/(e-s) )

            fig, ax = isp.rates( ranges, values, color='g', label=dq_name.replace("_","\_"), figax=None)
            ax.set_ylim(ymin=0, ymax=1)
            ax.set_ylabel('duty cycle')
            ax.legend(loc='lower left')

            figname = isp.ratefig(output_dir, ifo, "%s_%s"%(usertag, dq_name), trend_start, gpsstart+stride-trend_start)
            logger.info('  plotting %s'%figname)
            fig.savefig(figname)
            isp.close(fig)

            trending['sci'].append( figname )

            ### plot overlay
            uptime_overlay_figax = isp.rates( ranges, values, color='g', label=dq_name.replace("_","\_"), figax=uptime_overlay_figax )

        fig, ax = uptime_overlay_figax
        ax.set_ylim(ymin=0, ymax=1)
        ax.set_ylabel('duty cycle')
        ax.legend(loc='lower left')

        figname = isp.ratefig(output_dir, ifo, "%s_segments"%(usertag), trend_start, gpsstart+stride-trend_start)
        logger.info('  plotting %s'%figname)
        fig.savefig(figname)
        isp.close(fig)

        trending['seg'].append( figname )

        ### glitch/clean rates and effatfap
        trg_rate_overlay_figax = None
        eff_fap_overlay_figax = [None]*len(opts.FAPthr)
        for classifier in classifiers:
            color = colors[classifier]

            ### extract data
            ranges = []
            gvalues = []
            cvalues = []
            effs = []
            for bindata in bindatas:
                s, e = bindata['range']
                r, c, g, _, _ = bindata['rcg'][classifier]
                ranges.append( (s, e) )
                gvalues.append( 1.0*g[-1]/(e-s) )
                cvalues.append( 1.0*c[-1]/(e-s) )

                if c[-1]:
                    fap = np.array(c, dtype='float')/c[-1]
                else:
                    fap = [0.0]
                if g[-1]:
                    eff = np.array(g, dtype='float')/g[-1]
                else:
                    eff = [0.0]
        
                ### iterate over opts.FAPthr -> efficiencies at these FAPs -> plot!
                _effs = []
                for fapthr in opts.FAPthr:
                    if np.any(fap<=fapthr):
                        _effs.append( np.max(eff[fap <= fapthr]) )
                    else:
                        _effs.append( 0.0 )
                effs.append( _effs )

            effs = np.array(effs )

            figname = isp.ratefig(output_dir, ifo, "_%s%s_trg"%(classifier, usertag), trend_start, gpsstart+stride-trend_start)
            trending['trg_rate'][classifier].append( figname )
            logger.info('  plotting %s'%figname)
            
            figax = isp.rates( ranges, cvalues, linestyle=cln_linestyle , color=color, label='%s cln'%classifier, figax=None)
            fig, ax = isp.rates( ranges, gvalues, linestyle=gch_linestyle , color=color, label='%s gch'%classifier, figax=figax)

            ax.set_ylabel('event rate [Hz]')
            ax.legend(loc='best')

            fig.savefig(figname)
            isp.close(fig)

            trg_rate_overlay_figax = isp.rates( ranges, cvalues, linestyle=cln_linestyle , color=color, label='%s cln'%classifier, figax=trg_rate_overlay_figax)
            trg_rate_overlay_figax = isp.rates( ranges, gvalues, linestyle=gch_linestyle , color=color, label='%s gch'%classifier, figax=trg_rate_overlay_figax)

            ### effatfap
            figax = None
            c = 'b g r m c y k'.split()
            cind = 0
            for ind, fap in enumerate(opts.FAPthr):
                figax = isp.rates( ranges, effs[:,ind], color=c[cind], label='%s fap<=%.3e'%(classifier, fap), figax=figax)
                eff_fap_overlay_figax[ind] = isp.rates( ranges, effs[:,ind], color=color, label=classifier, figax=eff_fap_overlay_figax[ind])
                cind = (cind+1)%len(c)
            if figax:
                fig, ax = figax
                ax.set_ylim(ymin=0, ymax=1)
                ax.set_ylabel('Glitch Detection Efficiency')
                ax.legend(loc='lower left')

            figname = isp.ratefig(output_dir, ifo, "_%s%s_eff"%(classifier, usertag), trend_start, gpsstart+stride-trend_start)
            trending['effatfap'][classifier].append( figname )
            logger.info('  plotting %s'%figname)
            if figax:
                fig.savefig(figname)
                isp.close(fig)

        ### clean/glitch rates overlay
        fig, ax = trg_rate_overlay_figax
        ax.set_ylabel('event rate [Hz]')
        ax.legend(loc='best')

        figname = isp.ratefig(output_dir, ifo, "_overlay%s_trg"%(usertag), trend_start, gpsstart+stride-trend_start)
        trending['trg_rate']['overlay'].append( figname )
        logger.info('  plotting %s'%figname)
        fig.savefig( figname )
        isp.close(fig)

        eff_fap_overlays = []
        ### effatfap overlay
        for (fig, ax), fap in zip(eff_fap_overlay_figax, opts.FAPthr):
            ax.set_ylim(ymin=0, ymax=1)
            ax.set_ylabel('Glitch Detection Efficiency')
            ax.legend(loc='lower left')

            figname = isp.ratefig(output_dir, ifo, "_overlay%s_eff-fap<%.3e"%(usertag, fap), trend_start, gpsstart+stride-trend_start)
            logger.info('  plotting %s'%figname)
            eff_fap_overlays.append( figname )
            fig.savefig( figname )
            isp.close( fig )

        trending['effatfap']['overlay'].append( eff_fap_overlays )

        ### channel performance
        for classifier in classifiers:
            if classifiersD[classifier]['flavor'] != 'ovl':
                continue

            ranges = []
            performances = []
            num_gchs = []
            num_clns = []
            vchans = set()
            for bindata in bindatas:
                ranges.append( bindata['range'] )
                performance, gch, cln = bindata['chn'][classifier]
                for vchan in performance.keys():
                    vchans.add( vchan )
                performances.append( performance )
                num_gchs.append( gch )
                num_clns.append( cln )
            vchans = sorted([vchan for vchan in vchans if vchan != "none"])

            (figr, axr), (fige, axe), (figf, axf) = isp.channel_performance( vchans, ranges, performances, num_gchs, num_clns, figax=None, cmap=None)

            rnkname = isp.chanfig(output_dir, ifo, classifier, "eff-fap", usertag, trend_start, gpsstart+stride-trend_start)
            effname = isp.chanfig(output_dir, ifo, classifier, "eff", usertag, trend_start, gpsstart+stride-trend_start)
            fapname = isp.chanfig(output_dir, ifo, classifier, "fap", usertag, trend_start, gpsstart+stride-trend_start)

            logger.info('  plotting %s'%rnkname)
            trending['rank'][classifier].append( rnkname )
            figr.savefig(rnkname)
            isp.close(figr)

            logger.info('  plotting %s'%effname)
            trending['eff'][classifier].append( effname )
            fige.savefig(effname)
            isp.close(fige)

            logger.info('  plotting %s'%fapname)
            trending['fap'][classifier].append( fapname )
            figf.savefig(fapname)
            isp.close(figf)












    """
TRENDING:

    first and second moment of residuals between historical ROC curves (trending)

    """




















    #===============================================================================================
    # write html page
    #===============================================================================================
    if not opts.no_html:












        """
REPORT:

vetolist/configuration lists -> make human readable
segment lists -> a form to request segments?

        """
















        logger.info('writing html for %s'%output_dir)

        title = "idq%s : %d-%d"%(usertag, gpsstart, stride)
        header = ""

        creation_time = idq.nowgps()

        process_params = __file__ ### the script currently executing
        for key, value in vars(opts).items():
            key = key.replace("_","-")

            if isinstance(value, str) and os.path.exists(value):
                value = os.path.abspath( value )
            
            if isinstance(value, list):
                for v in value:
                    process_params += " --%s %s"%(key, v)
            elif isinstance(value, bool):
                if value:
                    process_params += " --%s"%(key)
            else:
                process_params += " --%s %s"%(key, value)

        footer = "page created on %s (%d) via:<br/>%s"%(gpstime.tconvert(creation_time), creation_time, process_params)

        page = markup.page()
        page.init( title=title, header=header, footer=footer)

        lensumdir = len(output_dir)

        page.h1( "idq%s : %s (%d) - %s (%d)"%(usertag, gpstime.tconvert( gpsstart-lookback), gpsstart-lookback, gpstime.tconvert(gpsstart+stride), gpsstart+stride ) )

        page.hr( )

        #=========================================
        # info about parameters
        #=========================================

        ### copy channels, unsafe_channels into output_dir
        os.system("cp %s %s"%(channels, output_dir) )
        os.system("cp %s %s"%(unsafechannels, output_dir) )

        minipage = markup.page()
        minipage.a( 'AUX channels', href="./%s"%(channels.split('/')[-1]) )
        minipage.br( )
        minipage.a( 'unsafe AUX channels', href='./%s'%(unsafechannels.split("/")[-1]) )
        minipage.br( )
        minipage.p( ('GW channel : %s'%gwchannel, 'gch_thr : %.1f'%gch_kwsignif_thr, 'cln_thr : %.1f'%cln_kwsignif_thr) )

        page.p( str(minipage) )

        #=========================================
        # pointers to data
        #=========================================

        ### pointers to segments!
        page.hr( )

        minipage = markup.page()
        if files['segs'].has_key('sci'):
            minipage.a( dq_name, href="./%s"%files['segs']['sci'][lensumdir:] )
            minipage.br( )
        minipage.a( 'idq_segments', href="./%s"%files['segs']['idq'][lensumdir:] )

        page.p( (str(minipage)) )

        ### pointers to data storage
        page.hr( )

        minipage = markup.page()
        for classifier in classifiers[:-1]:
            minipage.a( classifier, href="./%s"%files['dat'][classifier][lensumdir:] )
            minipage.addcontent( ", " )
        minipage.a( classifiers[-1], href="./%s"%files['dat'][classifiers[-1]][lensumdir:] )

        page.p( "datfile caches: %s"%str(minipage) )

        minipage = markup.page()
        for classifier in classifiers[:-1]:
            minipage.a( classifier, href="./%s"%files['roc'][classifier][lensumdir:] )
            minipage.addcontent( ", " )
        minipage.a( classifiers[-1], href="./%s"%files['roc'][classifiers[-1]][lensumdir:] )

        page.p( "rocfiles: %s"%str(minipage) )

        minipage = markup.page()
        for classifier in classifiers[:-1]:
            minipage.a( classifier, href="./%s"%files['kde'][classifier]['cln'][lensumdir:] )
            minipage.addcontent( ", ")
        minipage.a( classifiers[-1], href="./%s"%files['kde'][classifiers[-1]]['cln'][lensumdir:] )

        page.p( "cln kde files: %s"%str(minipage) )

        minipage = markup.page()
        for classifier in classifiers[:-1]:
            minipage.a( classifier, href="./%s"%files['kde'][classifier]['gch'][lensumdir:] )
            minipage.addcontent( ", ")
        minipage.a( classifiers[-1], href="./%s"%files['kde'][classifiers[-1]]['gch'][lensumdir:] )

        page.p( "gch kde files: %s"%str(minipage) )

        #=========================================
        # number of samples for this stride
        #=========================================

        minipage = markup.page()

        micropage = markup.page()
        micropage.td( 'classifier' )
        micropage.td( 'No. cln' )
        micropage.td( 'No. gch' )

        minipage.tr( str(micropage) )

        for classifier in classifiers:

             micropage = markup.page()

             micropage.td( classifier )
             micropage.td( "%d"%data[classifier]['num_cln'] )
             micropage.td( "%d"%data[classifier]['num_gch'] )

             minipage.tr( str(micropage) )

        page.table( str(minipage) )

        #=========================================
        # images for this stride 
        #=========================================
        page.hr( )

        for figflavor in 'linroc roc hst kde L fap fapUL'.split():
            page.img( width=html_width, height=html_height, alt="%s overlay"%figflavor, src="./%s"%fignames[figflavor]['overlay'][lensumdir:] )
            page.br( )

            minipage = markup.page()
            for classifier in classifiers[:-1]:
                minipage.a( classifier, href="./%s"%fignames[figflavor][classifier][lensumdir:] )
                minipage.addcontent( ", " )
            minipage.a( classifiers[-1], href="./%s"%fignames[figflavor][classifiers[-1]][lensumdir:] )

            page.p( "%s figures: %s"%(figflavor, str(minipage)) )

        page.img( width=html_width, height=html_height, alt='bitword cln', src="./%s"%fignames['bwh']['cln'][lensumdir:] )
        page.img( width=html_width, height=html_height, alt='bitword gch', src="./%s"%fignames['bwh']['gch'][lensumdir:] )
        page.br( )
    
        ### print trg histograms, or at least links!
        for column in dat_columns:
            minipage = markup.page()
            for classifier in classifiers:
                minipage.a( classifier, href="./%s"%fignames['trg'][(classifier, column)][lensumdir:] )
                minipage.addcontent( ", " )
            minipage.a( classifiers[-1], href="./%s"%fignames['trg'][(classifiers[-1], column)][lensumdir:] )
            
            page.p( "%s : %s"%(column,str(minipage)) )

        page.hr( )

        ### print channel performance files
        chn_classifiers = sorted(files['chn'].keys())
        if chn_classifiers:
            minipage = markup.page()
            for classifier in chn_classifiers[:-1]:
                minipage.a( classifier, href="./%s"%files['chn'][classifier][lensumdir:] )
                minipage.addcontent( ", ")
            minipage.a( chn_classifiers[-1], href="./%s"%files['chn'][chn_classifiers[-1]][lensumdir:] )

            page.p( "channel performance : %s"%str(minipage) )

        page.hr( )

        #=========================================
        # trending images
        #=========================================

        ### seg trending plots
        for figname in trending['seg']:
            page.img( width=html_width, height=html_height, alt="duty cycle trending", src="./%s"%figname[lensumdir:] )
#        page.br( )

        for segtype in ['idq', 'sci']:
            if trending[segtype]:
                minipage = markup.page()
                for figname in trending[segtype][:-1]:
                    label = "%d_%d"%tuple(idq.extract_start_stop( figname, suffix=".png" ))
                    minipage.a( label, href="./%s"%figname[lensumdir:] )    
                    minipage.addcontent(", ")
                figname = trending[segtype][-1]
                label = "%d_%d"%tuple(idq.extract_start_stop( figname, suffix=".png" ))
                minipage.a( label, href="./%s"%figname[lensumdir:] )

            page.p( "%s : %s"%(segtype, str(minipage)) )

        page.hr( )

        ### trigger rates
        minipages = dict( (classifier, markup.page()) for classifier in classifiers )
        N = len(opts.trending)
        for ind in xrange(N-1):
            figname = trending['trg_rate']['overlay'][ind]
            page.img( width=html_width, height=html_height, alt="trigger rates", src="./%s"%figname[lensumdir:] )
            for classifier in classifiers:
                figname = trending['trg_rate'][classifier][ind]
                label = "%d_%d"%tuple(idq.extract_start_stop( figname, suffix=".png" ))
                minipages[classifier].a( label, href="./%s"%figname[lensumdir:] )
                minipages[classifier].addcontent(", ")
        figname = trending['trg_rate']['overlay'][-1]
        page.img( width=html_width, height=html_height, alt="trigger rates", src="./%s"%figname[lensumdir:] )
#        page.br( )       

        for classifier in classifiers:
            figname = trending['trg_rate'][classifier][-1]
            label = "%d_%d"%tuple(idq.extract_start_stop( figname, suffix=".png" ))
            minipages[classifier].a( label, href="./%s"%figname[lensumdir:] )
            page.p( "%s : %s"%(classifier, str(minipages[classifier])) )

        page.hr( )

        ### effatfap
        for ind, fap in enumerate(opts.FAPthr):
            minipage = markup.page()
            minipage.br( )
            for trendind in xrange(N):
                figname = trending['effatfap']['overlay'][trendind][ind]
                minipage.img( width=html_width, height=html_height, alt="efficiency at fixed fap", src="./%s"%figname[lensumdir:] )
            page.p( "fap <= %.3e : %s"%(fap, str(minipage) ) )

        minipages = dict( (classifier, markup.page()) for classifier in classifiers )
        for ind in xrange(N-1):
            for classifier in classifiers:
                figname = trending['effatfap'][classifier][ind]
                label = "%d_%d"%tuple(idq.extract_start_stop( figname, suffix=".png" ))
                minipages[classifier].a( label, href="./%s"%figname[lensumdir:] )
                minipages[classifier].addcontent(", ")
        for classifier in classifiers:
            figname = trending['effatfap'][classifier][-1]
            label = "%d_%d"%tuple(idq.extract_start_stop( figname, suffix=".png" ))
            minipages[classifier].a( label, href="./%s"%figname[lensumdir:] )
            page.p( "%s : %s"%(classifier, str(minipages[classifier])) )

        page.hr( )

        ### channel performance plots
        for classifier in classifiers:
            if classifiersD[classifier]['flavor'] != 'ovl':
                continue

            page.h1( "%s channel performance"%classifier )

            ### rank
            page.p( "Rank" )
            for figname in trending['rank'][classifier]:
                page.img( width=html_width, height=html_height, alt="%s channel rank"%classifier, src="./%s"%figname[lensumdir:] )
            page.br()

            ### eff
            page.p( "Glitch Detection Efficiency" )
            for figname in trending['eff'][classifier]:
                page.img( width=html_width, height=html_height, alt="%s channel eff"%classifier, src="./%s"%figname[lensumdir:] )
            page.br( )

            ### fap
            page.p( "False Alarm Probability" )
            for figname in trending['fap'][classifier]:
                page.img( width=html_width, height=html_height, alt="%s channel fap"%classifier, src="./%s"%figname[lensumdir:] )
#            page.br( )

            page.hr( )

        #=========================================
        ### write page to disk
        html = idq.sumhtml( output_dir, usertag, gpsstart-lookback, lookback+stride)
        logger.info('writing %s'%html)
        file_obj = open(html, "w")
        print >> file_obj, page
        file_obj.close()

        #===========================================================================================
        # format index.html for summary_dir
        #===========================================================================================

        logger.info('updating index.html for %s'%summarydir)

        title = "idq%s"%usertag

        index = markup.page()
        index.init( title=title, header=header, footer=footer )

        ### iterate over existing subdirectories
        for subdir in sorted(glob.glob("%s/*_*/"%summarydir), reverse=True):
            sub = subdir.split('/')[-2]
            subindex = glob.glob("%s/summary-index*html"%subdir)
            if subindex:
                index.a( sub, href="./%s"%("/".join(subindex[0].split("/")[-2:])) )
            else:
                index.a( sub, href="./%s"%sub )
            roc = glob.glob("%s/*overlay*_ROC*png"%subdir)

            if roc:
                index.img( width=thumbnail_width, height=thumbnail_height, alt="%s roc overlay"%sub, src="./%s"%("/".join(roc[0].split('/')[-2:]) ) )
            index.br( )

        index.hr( )

        ### write page to disk
        file_obj = open("%s/index.html"%summarydir, "w")
        print >> file_obj, index
        file_obj.close()

    #===============================================================================================

    gpsstart += stride

#===================================================================================================
if opts.lockfile:
    idq.release( opts.lockfile ) ### unlock lockfile
