# Copyright (C) 2013 Reed Essick
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
import numpy as np
import re as re
import json

from ligo.gracedb.rest import GraceDb

from laldetchar.idq import idq
from laldetchar.idq import idq_summary_plots as isp

from laldetchar.idq import idq_gdb_utils as igu

from ConfigParser import SafeConfigParser

from optparse import OptionParser

from laldetchar import git_version

#===================================================================================================

__author__ = 'Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id__date__ = git_version.date

description = \
    """ Program generates a summary of iDQ performance around an event."""

#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'% git_version.verbose_msg,
        usage='%prog [options]',
        description=description)

parser.add_option('-v',
        '--verbose',
        default=False,
        action='store_true')

parser.add_option(
        '-c', '--config',
        default='idq.ini',
        type='string',
        )

parser.add_option(
        '-s',
        '--gps-start',
        dest='start',
        default=0,
        type='float',
        help='the gps start time of the time range of interest')

parser.add_option(
        '-e',
        '--gps-end',
        dest='end',
        default=0,
        type='float',
        help='the gps end time of the time range of interest')

parser.add_option('-g',
        '--gracedb-id',
        default=None,
        type='string',
        help='GraceDB ID')

parser.add_option('',
    '--skip-gracedb-upload',
    default=False,
    action='store_true',
    help='skip steps involving communication with GraceDB. Automatically set to True if --gracedb-id==None')

parser.add_option('-C',
        '--classifier',
        default='ovl',
        type='string',
        help='the classifier used to generate the timeseries data. Default="ovl"')

(opts, args) = parser.parse_args()

opts.skip_gracedb_upload = (opts.gracedb_id==None) or opts.skip_gracedb_upload

#=================================================
# read relevant stuff from config file
#=================================================
config = SafeConfigParser()
config.read( opts.config )

ifo = config.get('general','ifo')
tag = config.get('general','usertag')

if tag and opts.gracedb_id:
    filetag = "_%s_%s"%(tag, opts.gracedb_id)
elif opts.gracedb_id:
    filetag = "_%s"%(opts.gracedb_id)
elif tag:
    filetag = "_%s"%(tag)
else:
    filetag = ""

realtimedir = config.get('general','realtimedir')
gdbdir = config.get('gdb general','main_gdb_dir')

if not opts.skip_gracedb_upload:
    if config.has_option('gdb general', 'gdb_url'):
        gracedb = GraceDb(config.get('gdb general', 'gdb_url'))
    else:
        gracedb = GraceDb()

if config.has_option(opts.classifier, 'plotting_label'):
    plotting_label = config.get(opts.classifier, 'plotting_label')
else:
    plotting_label = opts.classifier

#===================================================================================================

rank_channame  = idq.channame(ifo, opts.classifier, "%s_rank"%tag)
fap_channame   = idq.channame(ifo, opts.classifier, "%s_fap"%tag)
fapUL_channame = idq.channame(ifo, opts.classifier, "%s_fapUL"%tag)

#===================================================================================================

'''
Need to generate:
  ROC curves 
    use both dat and gwf files to estimate this. Difference should primarily be in FAP estiamtes
    annotate with the min-FAP associated with this event
  calibration checks 
    use both dat and gwf files to estimate this.
  likelihood ratios
  active channels (over a wider timescale than timeseries?)
  information on when the last training was performed, how much data was available
  when the last calibration was performed, how much data was available
  iDQ glitch/clean rates
  else?

Essentially, a faster summary job without the html wrapper and uploaded to GraceDB.
  Should we just write an html document and post a link? (I tend to think no...)
  Should everything be tagged data_quality, or will that crowd the pull-down section too much?
'''


raise ValueError("this executable still needs to be written...")

















































# get all *.gwf files in range

if opts.verbose:
    print "Finding relevant *.gwf files"
rank_filenames = []
fap_filenames = []
all_files = idq.get_all_files_in_range(realtimedir, opts.plotting_gps_start, opts.plotting_gps_end, pad=0, suffix='.gwf')
for filename in all_files:
    if opts.classifier == idq.extract_fap_name(filename): # and ifo in filename: ### this last bit not needed?
        if 'rank' in filename:
            rank_filenames.append(filename)
        if 'fap' in filename:
            fap_filenames.append(filename)

rank_filenames.sort()
fap_filenames.sort()

if (not rank_filenames) or (not fap_filenames): # we couldn't find either rank or fap files
    # exit gracefully
    if opts.verbose:
        print "no iDQ timeseries for %s at %s"%(opts.classifier, ifo)
    if not opts.skip_gracedb_upload:
        gracedb.writeLog(opts.gracedb_id, message="No iDQ timeseries for %s at %s"%(opts.classifier, ifo))
    sys.exit(0)

#===================================================================================================
### process FAP files

if opts.verbose:
    print "reading fap timeseries from:"
    for filename in fap_filenames:
        print '\t' + filename

#fig = isp.plt.figure(figsize=isp.rank_timeseries_figsize)
#rnk_ax = fig.add_axes(isp.rank_timeseries_axpos)
#fap_ax = rnk_ax.twinx()
fig = isp.plt.figure(figsize=isp.rank_splittimeseries_figsize)
rnk_ax = fig.add_axes(isp.rank_splittimeseries_axpos)
fap_ax = fig.add_axes(isp.fap_splittimeseries_axpos)

fap_figax = (fig, fap_ax)
rnk_figax = (fig, rnk_ax)

if opts.gps:
    to = opts.gps
#    rnk_ax.plot( np.zeros((2, )), rnk_ax.get_ylim(), ':k', linewidth=2, alpha=0.5)
#    rnk_ax.text( 0.0, 0.8, "%.3f"%opts.gps, ha='center', va='center')
else:
    to = opts.plotting_gps_start

### check calibration type -> include UpperLimits when possible
if config.get('calibration','mode')=='dat': 
    f_times, f_timeseries = idq.combine_gwf(fap_filenames, [fap_channame, fapUL_channame])
    fUL_timeseries = [f[1] for f in f_timeseries]
    f_timeseries = [f[0] for f in f_timeseries]

    ### write combined data to disk
    if opts.verbose:
        print "writing combined fap frames to disk"
    for t, ts, tS in zip(f_times, f_timeseries, fUL_timeseries):
        truth = (opts.plotting_gps_start <= t)*(t <= opts.plotting_gps_end)
        t = t[truth]
        ts = ts[truth]
        tS = tS[truth]

        start = int(t[0])
        dt = t[1]-t[0]
        dur = int(len(t)*dt)
        fapfr = idq.gdb_timeseriesgwf( gdbdir , opts.classifier, ifo, "_fap%s"%filetag, start, dur)                      
        if opts.verbose:
            print "    %s"%fapfr
        idq.timeseries2frame( fapfr, {fap_channame:ts, fapUL_channame:tS}, t[0], dt )
        if not opts.skip_gracedb_upload:
            message = "iDQ fap timeseries for %s at %s within [%d, %d] :"%(opts.classifier, ifo, start, start+dur)
            if opts.verbose:
                print "    %s"%message
            gracedb.writeLog( opts.gracedb_id, message=message, filename=fapfr ) #, tagname=['data_quality'] )

    ### post min-fap value
    if opts.verbose:
        print "finding minimum FAP observed within [%.3f, %.3f]"%(opts.start, opts.end)
    min_fap = 1.0
    for (t, ts) in zip(f_times, f_timeseries):
        # ensure time series only fall within desired range
        truth = (opts.start <= t) * (t <= opts.end)
        ts = ts[truth]
        t = t[truth]

        if len(t): # some surviving data
            # generate and write summary statistics
            (f_min, f_max, f_mean, f_stdv) = idq.stats_ts(ts)

            # update min_fap
            if min_fap > f_min:
                min_fap = f_min

    # upload minimum fap observed within opts.start, opts.end
    jsonfilename = idq.gdb_minFap_json(gdbdir, opts.classifier, ifo, "_minFAP%s"%filetag, int(opts.start), int(opts.end-opts.start))
    file_obj = open(jsonfilename, "w")
    file_obj.write( json.dumps( {fap_channame:{"min":min_fap, "start":opts.start, "end":opts.end}} ) )
    file_obj.close()

    if not opts.skip_gracedb_upload:
        message = "minimum glitch-FAP for %s at %s within [%.3f, %.3f] is %.3e"%(opts.classifier, ifo, opts.start, opts.end, min_fap)
        if opts.verbose:
            print "    %s"%message
        gracedb.writeLog( opts.gracedb_id, message=message, filename=jsonfilename, tagname=['data_quality'] )

#    ### compute statistics
#    if opts.verbose:
#        print "computing fap timeseries statistics"
#    raise StandardError("WRITE ME")
#
#    ### write statistics to disk
#    if opts.verbose:
#        print "writing fap timeseries statistics to disk" ### ascii and json
#    raise StandardError("WRITE ME")

    ### plot
    if opts.verbose:
        print "plotting fap timeseries"
    fap_figax = isp.rank_timeseries_filled( f_times, fUL_timeseries, baseline=1, start=opts.plotting_gps_start, end=opts.plotting_gps_end, to=opts.gps, color='c', linestyle='dashed', figax=fap_figax, shade_alpha=0.0, zorder=9, alpha=0.5 )
#    fap_figax = isp.rank_timeseries( f_times, fUL_timeseries, start=opts.plotting_gps_start, end=opts.plotting_gps_end, to=opts.gps, color='c', linestyle='solid', figax=fap_figax, shade_alpha=0.0, zorder=9, alpha=1.0, linewidth=1.0 )
    fap_figax = isp.rank_timeseries( f_times, f_timeseries, start=opts.plotting_gps_start, end=opts.plotting_gps_end, to=opts.gps, color='b', linestyle='solid', figax=fap_figax, shade_alpha=0.1, zorder=10, alpha=0.75 )

else:
    f_times, f_timeseries = idq.combine_gwf(fap_filenames, [fap_channame])

    ### write combined data to disk
    if opts.verbose:
        print "writing combined fap frames to disk"
    for t, ts in zip(f_times, f_timeseries):
        truth = (opts.plotting_gps_start <= t)*(t <= opts.plotting_gps_end)
        t = t[truth]
        ts = ts[truth]

        start = int(t[0])
        dt = t[1]-t[0]
        dur = int(len(t)*dt)
        fapfr = idq.gdb_timeseriesgwf( gdbdir , opts.classifier, ifo, "_fap%s"%filetag, start, dur)
        if opts.verbose:
            print "    %s"%fapfr
        idq.timeseries2frame( fapfr, {fap_channame:ts}, t[0], dt )
        if not opts.skip_gracedb_upload:
            message = "iDQ fap frame for %s at %s within [%d, %d] :"%(opts.classifier, ifo, start, start+dur)
            if opts.verbose:
                print "    %s"%message
            gracedb.writeLog( opts.gracedb_id, message=message, filename=fapfr ) #, tagname=['data_quality'] )

    ### post min-fap value
    if opts.verbose:
        print "finding minimum FAP observed within [%.3f, %.3f]"%(opts.start, opts.end)
    min_fap = 1.0
    for (t, ts) in zip(f_times, f_timeseries):
        # ensure time series only fall within desired range
        truth = (opts.start <= t) * (t <= opts.end)
        ts = ts[truth]
        t = t[truth]

        if len(t): # some surviving data
            # generate and write summary statistics
            (f_min, f_max, f_mean, f_stdv) = idq.stats_ts(ts)

            # update min_fap
            if min_fap > f_min:
                min_fap = f_min

    # upload minimum fap observed within opts.start, opts.end
    jsonfilename = idq.gdb_minFap_json(gdbdir, opts.classifier, ifo, "_minFAP%s"%filetag, int(opts.start), int(opts.end-opts.start))
    file_obj = open(jsonfilename, "w")
    file_obj.write( json.dumps( {fap_channame:{"min":min_fap, "start":opts.start, "end":opts.end}} ) )
    file_obj.close()

    if not opts.skip_gracedb_upload:
        message = "minimum glitch-FAP for %s at %s within [%.3f, %.3f] is %.3e"%(opts.classifier, ifo, opts.start, opts.end, min_fap)
        if opts.verbose:
            print "    %s"%message
        gracedb.writeLog( opts.gracedb_id, message=message, filename=jsonfilename, tagname=['data_quality'] )

#    ### compute statistics
#    if opts.verbose:
#        print "computing fap timeseries statistics"
#    raise StandardError("WRITE ME")
#
#    ### write statistics to disk
#    if opts.verbose:
#        print "writing fap timeseries statistics to disk" ### ascii and json
#    raise StandardError("WRITE ME")

    ### plot
    if opts.verbose:
        print "plotting fap timeseries"
    fap_figax = isp.rank_timeseries( f_times, f_timeseries, start=opts.plotting_gps_start, end=opts.plotting_gps_end, to=opts.gps, color='b', linestyle='solid', figax=fap_figax, shade_alpha=0.1, zorder=10, alpha=0.75 )

#===================================================================================================
### process rank files

if opts.verbose:
    print "reading rank timeseries from:"
    for filename in rank_filenames:
        print '\t' + filename
(r_times, r_timeseries) = idq.combine_gwf(rank_filenames, [rank_channame])

### write combined data to disk
if opts.verbose:
    print "writing combined rank frames to disk:"
for t, ts in zip(r_times, r_timeseries):
    truth = (opts.plotting_gps_start <= t)*(t <= opts.plotting_gps_end)
    t = t[truth]
    ts = ts[truth]

    start = int(t[0])
    dt = t[1]-t[0]
    dur = int(len(t)*dt)
    rnkfr = idq.gdb_timeseriesgwf( gdbdir , opts.classifier, ifo, "_rank%s"%filetag, start, dur)
    if opts.verbose:
        print "    %s"%rnkfr
    idq.timeseries2frame( rnkfr, {rank_channame:ts}, t[0], dt )
    if not opts.skip_gracedb_upload:
        message = "iDQ glitch-rank frame for %s at %s within [%d, %d] :"%(opts.classifier, ifo, start, start+dur)
        if opts.verbose:
            print "    %s"%message
        gracedb.writeLog( opts.gracedb_id, message=message, filename=rnkfr ) #, tagname=['data_quality'] )

### post max-rank value
if opts.verbose:
    print "finding max rank observed within [%.3f, %.3f]"%(opts.start, opts.end)
    max_rnk = 0.0
    for (t, ts) in zip(r_times, r_timeseries):
        # ensure time series only fall within desired range
        truth = (opts.start <= t) * (t <= opts.end)
        ts = ts[truth]
        t = t[truth]

        if len(t): # some surviving data
            # generate and write summary statistics
            (r_min, r_max, r_mean, r_stdv) = idq.stats_ts(ts)

            # update min_fap
            if max_rnk < r_max:
                max_rnk = r_max

    # upload minimum fap observed within opts.start, opts.end
    jsonfilename = idq.gdb_minFap_json(gdbdir, opts.classifier, ifo, "_maxRANK%s"%filetag, int(opts.start), int(opts.end-opts.start))
    file_obj = open(jsonfilename, "w")
    file_obj.write( json.dumps( {fap_channame:{"min":min_fap, "start":opts.start, "end":opts.end}} ) )
    file_obj.close()

    if not opts.skip_gracedb_upload:
        message = "maximum glitch-RANK for %s at %s within [%.3f, %.3f] is %.3e"%(opts.classifier, ifo, opts.start, opts.end, min_fap)
        if opts.verbose:
            print "    %s"%message
        gracedb.writeLog( opts.gracedb_id, message=message, filename=jsonfilename, tagname=['data_quality'] )

#### compute statistics
#if opts.verbose:
#    print "computing rank timeseries statistics"
#raise StandardError("WRITE ME")
#
#### write statistics to disk
#if opts.verbose:
#    print "writing fap timeseries statistics to disk" ### ascii and json
#raise StandardError("WRITE ME")

### add to plot
if opts.verbose:
    print "plotting rank timeseries"
rnk_figax = isp.rank_timeseries( r_times, r_timeseries, start=opts.plotting_gps_start, end=opts.plotting_gps_end, to=opts.gps, color='r', linestyle='solid', figax=rnk_figax, shade_alpha=0.1, zorder=1, alpha=0.75 )

### look up flavor-specific stuff
if config.has_section(opts.classifier) and config.has_option(opts.classifier, 'flavor'): ### only allow classifiers, not combiners
    flavor = config.get(opts.classifier, 'flavor')
    if flavor == "ovl":
        if opts.verbose:
            print "  looking up extra data for flavor=ovl"

        time_ordered_vconfigs, channel_ordered_vconfigs, vconfig_keys = igu.extract_ovl_vconfigs( rank_filenames, rank_channame, config.get('general','traindir'), opts.plotting_gps_start, opts.plotting_gps_end, metric=config.get(opts.classifier, 'metric') )

        ### make plot
#        scfig, ax = isp.ovl_vconfig_stripchart( time_ordered_vconfigs, to=opts.gps, mapIt=vconfig_keys )
        scfig, ax = isp.ovl_vconfig_stripchart( channel_ordered_vconfigs, to=opts.gps, mapIt=vconfig_keys )

        ax.set_xlim(xmin=opts.plotting_gps_start-to, xmax=opts.plotting_gps_end-to)
        ax.set_title('iDQ (possible) active channels for %s at %s'%(plotting_label, ifo))

        figname = isp.ovlstripchart(gdbdir, ifo, opts.classifier, filetag, opts.plotting_gps_start, opts.plotting_gps_end-opts.plotting_gps_start, figtype="png")
        if opts.verbose:
            print "  saving : %s"%figname
        scfig.savefig(figname)
        isp.plt.close( scfig )

        if not opts.skip_gracedb_upload:
            message = "iDQ channel strip chart for %s at %s between [%.3f, %.3f]"%(opts.classifier, ifo, opts.plotting_gps_start, opts.plotting_gps_end)
            if opts.verbose:
                print "  %s"%message
            gracedb.writeLog(opts.gracedb_id, message=message, filename=figname, tagname='data_quality')
            
        ### json file containing segments
        json_cov = []
        for key, value in channel_ordered_vconfigs.items():
            key = [ dict((k,K[vconfig_keys[k]]) for k in vconfig_keys) for K in key]
            json_cov.append( {'vconfigs':key, 'segments':value} )

        ### perhaps we should make this a fancy xml table of some sort?
        jsonfilename = idq.gdb_ovlstripchart_json(gdbdir, opts.classifier, ifo, filetag, int(opts.plotting_gps_start), int(opts.plotting_gps_end-opts.plotting_gps_start))
        if opts.verbose:
            print "  writing : %s"%jsonfilename
        file_obj = open(jsonfilename, "w")
        file_obj.write( json.dumps( json_cov ) )
        file_obj.close()
        if not opts.skip_gracedb_upload:
            message = "iDQ (possible) active channels for %s at %s between [%.3f, %.3f]"%(opts.classifier, ifo, opts.plotting_gps_start, opts.plotting_gps_end)
            if opts.verbose:
                print "  %s"%message
            gracedb.writeLog(opts.gracedb_id, message=message, filename=jsonfilename, tagname='data_quality')

    elif opts.verbose:
        print "flavor=%s does not require any special extra data"%flavor

#===================================================================================================

### add annotations to the plot

if opts.gch_xml:
    ### annotate glitches
    gps_times = [ gps + gps_ns * 1e-9 for (gps, gps_ns) in igu.get_glitch_times(opts.gch_xml)]
    gps_times = np.asarray(gps_times)
    # channel annotation looks awfull on the static plot, must find a better way to visualize it
    #rnk_ax.text(gps - opts.plotting_gps_start, 0.1, auxchannel, ha="left", va="bottom", rotation=45)
    rnk_ax.plot(gps_times - to, 0.1*np.ones(len(gps_times)), marker="x", markerfacecolor="g", \
        markeredgecolor = "g", linestyle="none")

### shade region outside of opts.start, opts.end
if opts.start != opts.plotting_gps_start:
    for ax in [fap_ax, rnk_ax]:
        ax.fill_between( [opts.plotting_gps_start - to, opts.start-to],
            np.ones((2, ))*1e-5,
            np.ones((2, )),
            color='k',
            edgecolor='none',
            alpha=0.25,
            )

if opts.end != opts.plotting_gps_end:
    for ax in [fap_ax, rnk_ax]:
        ax.fill_between( [opts.end-to, opts.plotting_gps_end-to],
            np.ones((2, ))*1e-5,
            np.ones((2, )),
            color='k',
            edgecolor='none',
            alpha=0.25,
            )

#fig.suptitle('%s at %s'%(opts.classifier, ifo))
fap_ax.set_title('iDQ glitch False Alarm Probability for %s at %s'%(plotting_label, ifo))

rnk_ax.set_xlim(xmin=opts.plotting_gps_start-to, xmax=opts.plotting_gps_end-to)
fap_ax.set_xlim(rnk_ax.get_xlim())

isp.plt.setp(fap_ax.get_xticklabels(), visible=False)

rnk_ax.set_ylim(ymin=0-2e-2, ymax=1+2e-2)
fap_ax.set_ylim(ymin=1e-5*10**(-5*2e-2), ymax=10**(5*2e-2))

rnk_ax.set_xlabel('Time since %.3f'%to, visible=True)
fap_ax.set_xlabel('Time since %.3f'%to, visible=False)

rnk_ax.set_ylabel('%s rank' % plotting_label, color='r')
fap_ax.set_ylabel('%s FAP' % plotting_label, color='b')

rnk_ax.set_yscale('linear')
fap_ax.set_yscale('log')

rnk_ax.yaxis.set_label_position('left')
rnk_ax.yaxis.tick_left()

fap_ax.yaxis.set_label_position('right')
fap_ax.yaxis.tick_right()

### save figure
figname = isp.timeseriesfig(gdbdir, ifo, opts.classifier, filetag, int(opts.plotting_gps_start), int(opts.plotting_gps_end-opts.plotting_gps_start))
if opts.verbose:
    print '  saving : %s'%figname
fig.savefig(figname)
isp.plt.close(fig)

if not opts.skip_gracedb_upload:
    message = "iDQ fap and glitch-rank timeseries plot for " + opts.classifier + " at "+ifo+":"
    if opts.verbose:
        print "  %s"%message
    gracedb.writeLog(opts.gracedb_id, message=message, filename=figname, tagname='data_quality')

if opts.verbose:
    print "Done"

