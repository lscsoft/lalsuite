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

import ConfigParser
from optparse import *
import sys
import os
import time
from laldetchar.idq import idq
import numpy
from laldetchar.idq import idq_summary_plots as idq_s_p
import traceback
import logging
from laldetchar import git_version

__author__ = \
    'Lindy Blackburn (<lindy.blackburn@ligo.org>), Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date


def generate_html(
    path,
    gpsstart=False,
    gpsstop=False,
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
    if tot_livetime:
        print >> f, '<p>total Livetime = %f seconds</p>' % tot_livetime
    if tot_glitches:
        print >> f, '<p>total No.glitches = %d events</p>' \
            % int(tot_glitches)
    if roc_url:
        print >> f, '<h2>Reciever Operating Characteristic Curve</h2>'
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


def path_to_url(path, base, remove):
    """ converts a path name to a url by removing all directory reference in "remove" and appending the remainder of "path" to "base" """

    path_list = path.split('/')
    for p in path_list[:-1]:
        if p not in remove and p != '':
            base += p + '/'
    return base + path_list[-1]


################################################################

description = \
    """This program generates summary html pages from iDQ pipeline output. The summmary pages provide variety of diagnostic and interpretational plots and data."""

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
parser.add_option('-l', '--log-file', default='idq_summary.log',
                  type='string', help='log file')
(opts, args) = parser.parse_args()

################################################################

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

config = ConfigParser.SafeConfigParser()
config.read(opts.config)

# generate a dictionary for idq_summary specific options

myconf = dict(config.items('idq_summary'))
stride = int(myconf['stride'])
delay = int(myconf['delay'])
lookback = int(myconf['lookback'])
cluster_win = float(myconf['cluster_win'])
FAP = float(myconf['fap'])

symlink_path = myconf['symlink']
url_base = myconf['url_base']
url_remove = myconf['url_remove'].split()

gw_thr = float(myconf['gw_thr'])
classifiers_switch_snr_signif = myconf['switch_snr_signif'].split()

# other parameters from config

sumdir = config.get('general', 'summarydir')
realtimedir = config.get('general', 'realtimedir')
dailydir = config.get('general', 'dailydir')
classifiers = config.get('general', 'classifiers').split()
vetolist_cache = config.get('general', 'ovl_train_cache')
unsafe_win = float(config.get('idq_realtime', 'clean_window'))
columns = config.get('idq_realtime', 'dat_columns').split() + ['rank']
kwtrgdir = config.get('general', 'kwtrgdir')

kde_num_samples = int(config.get('idq_summary', 'kde_num_samples'))

# set classifier colors and labels

classifier_colors = [idq_s_p.classifier_colors(classifier)
                     for classifier in classifiers]
classifier_labels = [idq_s_p.classifier_labels(classifier)
                     for classifier in classifiers]

# current time and boundaries

t = int(idq.nowgps())
if not opts.gps_stop:
    print 'computing gpsstop from current time'
    gpsstop = (t - delay) / stride * stride  # gpsstop time for this analysis
else:
    gpsstop = opts.gps_stop / stride * stride
print 'gpsstop = %d' % gpsstop
if not opts.gps_start:
    print 'computing gpsstart from gpsstop'
    gpsstart = gpsstop - stride
else:
    gpsstart = opts.gps_start / stride * stride
print 'gpsstart = %d' % gpsstart

######
# MAIN
######

# loop over all data ranges

while __name__ == '__main__' and gpsstart < gpsstop:
    logger.info('-----------------------------------------------------------------'
                )
    logger.info('summarizing data from %d to %d' % (gpsstart, gpsstart
                + stride))

    # generate output directory for this data

    this_sumdir = sumdir + '/' + str(gpsstart) + '_' + str(gpsstart
            + stride)
    if not os.path.exists(this_sumdir):
        os.makedirs(this_sumdir)

    # collect *dat files and merge into *roc files

    logger.info('generating *.roc files')
    roc_paths = idq.datfiles_to_roc(
        gpsstart,
        gpsstart + stride,
        columns=columns,
        classifiers=classifiers,
        basename=False,
        source_dir=realtimedir,
        output_dir=this_sumdir,
        cluster_win=cluster_win,
        unsafe_win=unsafe_win,
        gw_thr=gw_thr,
        switch_snr_signif=classifiers_switch_snr_signif,
        )

    # generate uniformly sampled ROC files (sampled 100 times)

    logger.info('generating *.uroc files')
    uroc_paths = []
    for (roc_path, classifier) in roc_paths:
        (uniform_ranks, uniform_ccln, uniform_cgch, tcln, tgch) = \
            idq_s_p.ROC_to_uniformROC(roc_path, num_samples=100)
        uniformROCfilename = roc_path[:-4] + '.uroc'  # suffix is for uniformly-sampled roc file
        uroc_paths.append((idq_s_p.rcg_to_ROC(
            roc_path[:-4] + '.uroc',
            uniform_ranks,
            uniform_ccln,
            uniform_cgch,
            tcln,
            tgch,
            ), classifier))

    # update uroc_cachefiles!

    logger.info('updating *_uroc.cache files')
    for (uroc_path, classifier) in uroc_paths:
        uroc_cachefilename = sumdir + '/' + classifier + '_uroc.cache'
        file = open(uroc_cachefilename, 'a')
        print >> file, uroc_path
        file.close()

    logger.info('generating ROC figures')

    # generate an individual ROC plot for each classifier

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

    # generate a combined ROC plot showing all classifiers

    try:
        figname = this_sumdir + '/all-' + str(gpsstart) + '-' \
            + str(stride) + '_roc'
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

    # location of most recent vetolist:

    logger.info('finding pointers to OVL vetolists')
    vetolist_link = False
    try:
        vetolist_path = open(vetolist_cache, 'r'
                             ).readlines()[-1].strip('\n')
        vetolist_link = this_sumdir + '/vetolist.eval'
        if os.path.lexists(vetolist_link):
            os.remove(vetolist_link)
        os.symlink(vetolist_path, vetolist_link)
    except:
        traceback.print_exc()
        logger.info('WARNING: FAILED to find most recent OVL vetolist')
    logger.info('Done')

    # generate segments

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

    # trending plots

    logger.info('generating trending plots')

    # compute how far in time to look back

    lookbacktime = gpsstart - lookback * stride

    # get stat summary files

    stat_summary_files = idq.get_all_files_in_range(sumdir,
            lookbacktime, gpsstart + stride, pad=0, suffix='.stat')

    # get the latest stat file

    stat_summary_files.sort()
    latest_stat_summary_file = stat_summary_files[-1]

    # generate trending plots for livetime, glitch and clean samples rates

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

    # get roc files

    roc_files_for_trending = idq.get_all_files_in_range(sumdir,
            lookbacktime, gpsstart + stride, pad=0, suffix='.roc')

    # define dictionary to hold them

    roc_files_dict = {}
    for classifier in classifiers:
        roc_files_dict[classifier] = []

    for file in roc_files_for_trending:
        classifier = file.split('/')[-1].split('-')[0]
        roc_files_dict[classifier].append(file)

    # generate efficiency trending plot....

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

    # channel statistics

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
                source_dir=realtimedir,
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

    # generate channel trending plot

    logger.info('generating channel performance trending plot')
    chanlist_trend = []
    for classifier in classifiers:
        if classifier not in ['ovl']:
            logger.info('skipping ' + classifier)
            continue
        try:
            figure_name = this_sumdir + '/' + classifier \
                + '-%d-%d_channel_performance_trends.png' \
                % (lookbacktime, gpsstart + stride)
            chan_perform_png = idq_s_p.chanlist_trending(
                lookbacktime,
                gpsstart + stride,
                sumdir,
                classifier=classifier,
                figure_name=figure_name,
                annotated=True,
                )
            chanlist_trend.append((chan_perform_png, classifier))
        except:
            traceback.print_exc()
            chanlist_trend.append(('', classifier))
            logger.info('WARNING: FAILED to generated channel trending plot for '
                         + classifier)

    logger.info('Done')

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
                source_dir=realtimedir,
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

    # generate config statistics

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
                source_dir=realtimedir,
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

    # generate kde overlays

    logger.info('generating kde overlays')
    kde_fig_paths = []
    for (ind, roc_path) in enumerate(roc_paths):
        (path, classifier) = roc_path
        try:
            figname = path[:-4] + '_pwg_kde' + idq_s_p.fig_type
            fig = idq_s_p.ROC_to_pwg_kde_plot(path, write=False,
                    num_samples=kde_num_samples)

# ............print figname
            # ################## FIND AND PLOT OLD pwg_kde FROM TRAINING FILES ##############################################
            # train_roc_path = dailydir +"/"+_trange+"/"+classifier+"/*roc"  <------ glob pattern
            # eval, kde_c, kde_g = idq_s_p.ROC_to_pwg_kde(train_roc_path)
            # idq_s_p.plt.plot(eval, kde_c, color=idq.__kde_c_color, linewidth=2, alpha=0.5)
            # idq_s_p.plt.plot(eval, kde_g, color=idq.__kde_g_color, linewidth=2, alpha=0.5)
            # ###############################################################################################################

            idq_s_p.plt.savefig(figname)
            idq_s_p.plt.close(fig)
            kde_fig_paths.append((figname, classifier))
        except:
            traceback.print_exc()
            kde_fig_paths.append(('', classifier))
            logger.info('WARNING: FAILED to generate kdw plot for '
                        + classifier)

    logger.info('Done')

    # glitch rates, glitch parameter histograms, pointers to the craziest glitches

    # important channel trending plots, figure of merit for important channels

    # plot residuals between old ROC curve and new ROC curve

    # channel appearence and importance plots (surface plots for all channels)

    # generate html summary page

    roc_url = False
    vetolist_url = False
    if roc_fig_path:
        roc_url = roc_fig_path.split(this_sumdir + '/')[-1]
    if vetolist_link:
        vetolist_url = vetolist_link.split(this_sumdir + '/')[-1]
    if vetolist_seg_path:
        vetolist_seg_url = vetolist_seg_path.split(this_sumdir + '/'
                )[-1]

    roc_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                for path in roc_fig_paths if path[0]]
    kde_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                for path in kde_fig_paths if path[0]]
    chanlist_urls = [(path[-1].split(this_sumdir + '/')[-1], path[1])
                     for path in chanlist if path[0]]

    # construct urls for trending plots

    stat_trends_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                        for path in stat_trends_plot_paths if path[0]]
    eff_trends_urls = [(path[0].split(this_sumdir + '/')[-1], path[1])
                       for path in eff_trends_plot_paths if path[0]]

    # construct url for stat summary file

    stat_summary_url = latest_stat_summary_file.split(this_sumdir + '/'
            )[-1]

    # construct url for channel performance file

    chanlist_trend_urls = [(path[0].split('/')[-1], path[1])
                           for path in chanlist_trend if path[0]]

    logger.info('generating html summary page')
    html_path = this_sumdir + '/' + str(gpsstart) + '_' + str(gpsstart
            + stride) + '-summary.html'
    generate_html(
        html_path,
        gpsstart=gpsstart,
        gpsstop=gpsstart + stride,
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

    # update symbolic link
# ....iif os.path.lexists(symlink_path):
# ........os.remove(symlink_path)
# ....os.symlink(html_path, symlink_path)

    gpsstart += stride

