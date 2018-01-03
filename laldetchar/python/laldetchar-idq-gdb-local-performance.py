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
import glob
import numpy as np
import re as re
import json

from collections import defaultdict

from ligo.gracedb.rest import GraceDb

from laldetchar.idq import event
from laldetchar.idq import idq
from laldetchar.idq import idq_summary_plots as isp

from laldetchar.idq import idq_gdb_utils as igu

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
#from glue.ligolw import table

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

parser.add_option('-F', 
        '--FAPthr',
        default=[],
        action='append',
        type='float')

parser.add_option('-S',
        '--KWsignifThr',
        default=[],
        action='append',
        type='float')

parser.add_option('',
        '--ignore-science-segments',
        default=False, 
        action='store_true')

(opts, args) = parser.parse_args()

opts.skip_gracedb_upload = (opts.gracedb_id==None) or opts.skip_gracedb_upload

if not opts.FAPthr:
    opts.FAPthr.append( 0.1 )

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
traindir = config.get('general','traindir')
calibrationdir = config.get('general','calibrationdir')

gdbdir = config.get('gdb general','main_gdb_dir')

gwchannel = config.get('general', 'gwchannel')

if not opts.KWsignifThr:
    opts.KWsignifThr.append( config.getfloat('general', 'gw_kwsignif_thr') )

GWkwconfig = idq.loadkwconfig(config.get('data_discovery', 'GWkwconfig'))
GWkwbasename = GWkwconfig['basename']
GWgdsdir = config.get('data_discovery', 'GWgdsdir')

kwstride = int(float(GWkwconfig['stride']))

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

### science segments
if opts.ignore_science_segments:
    if opts.verbose:
        print 'analyzing data regardless of science segements'
    scisegs = [[opts.start, opts.end]] ### set segs to be this stride range
    coveredsegs = [[opts.start, opts.end]] ### set segs to be this stride range

else:
    ### load settings for accessing dmt segment files
    dq_name = config.get('get_science_segments', 'include')
    segdb_url = config.get('get_science_segments', 'segdb')

    if opts.verbose:
        print 'querrying science segments'

    ### this returns a string
    seg_xml_file = idq.segment_query(config, opts.start , opts.end, url=segdb_url)

    ### write seg_xml_file to disk
    lsctables.use_in(ligolw.LIGOLWContentHandler)
    xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

    ### science segments xml filename
    seg_file = idq.segxml(gdbdir, "%s_%s"%(filetag, dq_name), opts.start , opts.end-opts.start )
    if opts.verbose:
        print '  writing science segments to file : '+seg_file
    ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))

    (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name) ### read in segments from xml file

if opts.verbose:
    print 'finding idq segments'
idqsegs = idq.get_idq_segments(realtimedir, opts.start, opts.end, suffix='.dat')

if opts.verbose:
    print 'taking intersection between science segments and idq segments'
idqsegs = event.andsegments( [scisegs, idqsegs] )

### write segment file
if opts.ignore_science_segments:
    idqseg_path = idq.idqsegascii(gdbdir, filetag, opts.start, opts.end-opts.start)
else:
    idqseg_path = idq.idqsegascii(gdbdir, '%s_%s'%(filetag, dq_name), opts.start, opts.end-opts.start)
if opts.verbose:
    print "  writing : "+idqseg_path
f = open(idqseg_path, 'w')
for seg in idqsegs:
    print >> f, seg[0], seg[1]
f.close()

#=================================================

rank_channame  = idq.channame(ifo, opts.classifier, "%s_rank"%tag)
fap_channame   = idq.channame(ifo, opts.classifier, "%s_fap"%tag)
fapUL_channame = idq.channame(ifo, opts.classifier, "%s_fapUL"%tag)

flavor = config.get(opts.classifier, 'flavor')
if config.has_option(opts.classifier, 'plotting_label'):
    plotting_label = config.get(opts.classifier, 'plotting_label')
else:
    plotting_label = opts.classifier

#===================================================================================================

### Find all FAP files
if opts.verbose:
    print "finding all fap*gwf files"
faps = [fap for fap in idq.get_all_files_in_range( realtimedir, opts.start, opts.end, pad=0, suffix='.gwf') if ('fap' in fap) and (opts.classifier==idq.extract_fap_name( fap )) and event.livetime(event.andsegments([[idq.extract_start_stop(fap, suffix=".gwf")], idqsegs])) ]

### compute total time covered
#T = event.livetime( [idq.extract_start_stop(fap, suffix='.gwf') for fap in faps] )*1.0
T = event.livetime( idqsegs )*1.0

### combine timeseries and generate segments
if opts.verbose:
    print "generating segments from %d fap files"%(len(faps))
segs = dict( (fapThr, [[], 1.0]) for fapThr in opts.FAPthr )
t, ts = idq.combine_gwf(faps, [fap_channame])
for t, ts in zip(t, ts):

    t, ts = idq.timeseries_in_segments( t, ts, idqsegs )

    for fapThr in opts.FAPthr:
        s, minFAP = idq.timeseries_to_segments(t, -ts, -fapThr) # we want FAP <= FAPthr <--> -FAP >= -FAPthr
        s = event.andsegments( [s, idqsegs] ) ### necessary because of how timeseries_to_segments may interact with timeseries_in_segments

        segs[fapThr][0] += s
        if minFAP!=None:
            segs[fapThr][1] = min(segs[fapThr][1], -minFAP)
if opts.verbose:
    print "computing associated deadtimes"
dt = [event.livetime(segs[fapThr][0])/T for fapThr in opts.FAPthr]
maxFAP = [segs[fapThr][1] for fapThr in opts.FAPthr]

### write json for calibration check
jsonfilename = idq.gdb_calib_json( gdbdir, ifo, opts.classifier, filetag, opts.start, opts.end-opts.start )
if opts.verbose:
    print "  %s"%jsonfilename
file_obj = open(jsonfilename, "w")
file_obj.write( json.dumps( {opts.classifier:{'nominal FAP':opts.FAPthr, 'maximum reported FAP':maxFAP, 'observed deadtime':dt, 'duration':T} } ) )
file_obj.close()

if not opts.skip_gracedb_upload:
    message = "iDQ calibration sanity check for %s at %s within [%.3f, %.3f]"%(opts.classifier, ifo, opts.start, opts.end)
    if opts.verbose:
        print "    "+message
    gracedb.writeLog( opts.gracedb_id, message=message, filename=jsonfilename )

### plot calibration check
fig = isp.plt.figure()
ax = fig.add_axes( isp.default_axpos )

if np.any(np.array(maxFAP)>0) and np.any(np.array(dt)>0): 
    ax.loglog(maxFAP, dt, marker='o', linestyle='none')
else:
    ax.plot(maxFAP, dt, marker='o', linestyle='none')

ax.set_xlabel('Nominal FAP')
ax.set_ylabel('Observed deadtime')

ax.grid(True)
xmin = 1e-4
xmax = 1.0
ax.set_xlim(xmin=xmin, xmax=xmax)
ax.set_ylim(ax.get_xlim())

faircoin = np.linspace( xmin, xmax, 1001 )
ax.plot( faircoin, faircoin, color='k' )

ax.set_title('iDQ local FAP calibration for %s at %s'%(plotting_label, ifo))

figname = isp.calibfig( gdbdir, ifo, opts.classifier, filetag, opts.start, opts.end-opts.start )
if opts.verbose:
    print "  %s"%figname
fig.savefig(figname)
isp.plt.close(fig)

if not opts.skip_gracedb_upload:
    message = "iDQ calibration sanity check figure for %s at %s within [%.3f, %.3f]"%(opts.classifier, ifo, opts.start, opts.end)
    if opts.verbose:
        print "    "+message
    gracedb.writeLog( opts.gracedb_id, message=message, filename=figname, tagname=['data_quality'] )

### discover KW triggers, compute efficiencies, plot ROC curves
if opts.verbose:
    print "retrieving gw triggers"
gwtrgs = event.include( idq.retrieve_kwtrigs(GWgdsdir, GWkwbasename, opts.start, opts.end-opts.start, kwstride, sleep=0, ntrials=1, verbose=False)[gwchannel], idqsegs, tcent=event.col_kw['tcent'] )

fig = isp.plt.figure()
ax = fig.add_axes( isp.default_axpos )

jsonD = {}
for signif in opts.KWsignifThr:
    trgs = [trg for trg in gwtrgs if trg[event.col_kw['signif']] >= signif]
    N = len(trgs)*1.0
    if N:
        eff = []
        for fapThr in opts.FAPthr:
            s = segs[fapThr][0]
            eff.append( len(event.include( trgs, s, tcent=event.col_kw['tcent'] ))/N )
    else:
        eff = [0.0 for fapThr in opts.FAPthr]

#    color = ax.step( dt, eff, label='%d events with KWsignif $\geq %.1f$'%(N, signif), where='post' )[0].get_color()
    color = ax.plot( dt, eff, label='%d events with KWsignif $\geq %.1f$'%(N, signif) )[0].get_color()

    jsonD[signif] = {'observed deadtime':dt, 'observed efficiency':eff, 'number of glitches':N, 'duration':T}

    l, h = [], []
    for e in eff:
        cr = idq.binomialCR( e*N, N, conf=0.68 )
        l.append(cr[0])
        h.append(cr[1])
    ax.fill_between( dt, l, h, color=color, alpha=0.25, edgecolor='none' )    

### add curve from dat files
if opts.verbose:
    print "finding all *dat files"
dats = [dat for dat in idq.get_all_files_in_range( realtimedir, opts.start, opts.end, pad=0, suffix='.dat') if (opts.classifier==idq.extract_dat_name( dat )) and event.livetime(event.andsegments([[idq.extract_start_stop(dat, suffix=".dat")], idqsegs]))]
if opts.verbose:
    print "reading samples from %d dat files"%len(dats)
output = idq.slim_load_datfiles( dats, skip_lines=0, columns=['GPS', 'i', 'rank'])
output = idq.filter_datfile_output( output, idqsegs )
#output['GPS'] = [float(l) for l in output['GPS']] ### not necessary because values are cast appropriately within idq.filter_datfile_output
#output['i'] = [float(l) for l in output['i']]
#output['rank'] = [float(l) for l in output['rank']]

r, c, g, = idq.dat_to_rcg( output )

if g[-1] and c[-1]:
    color = ax.step( 1.0*c/c[-1], 1.0*g/g[-1], label='datfiles: $N_c=%d$, $N_g=%d$'%(c[-1], g[-1]), linewidth=2, where='post' )[0].get_color()
    for G, C in zip(g, c):
        c_cr = idq.binomialCR( C, c[-1], conf=0.68 )
        g_cr = idq.binomialCR( G, g[-1], conf=0.68 )
        ax.fill_between( c_cr, [g_cr[0]]*2, [g_cr[1]]*2, color=color, alpha=0.25, edgecolor='none' )

jsonD['dat'] = {'rank':list(r), 'cumulative cleans':list(c), 'cumulative glitches':list(g)}

### finish plot
ax.set_xlabel('Deadtime (FAP)')
ax.set_ylabel('Efficiency')

ax.set_xscale('log')
ax.set_yscale('linear')

ax.grid(True)
xmin = 1e-4
xmax = 1.0
ax.set_xlim(xmin=xmin, xmax=xmax)
ax.set_ylim(ymin=0.0, ymax=1.0)

faircoin = np.linspace( xmin, xmax, 1001 )
ax.plot( faircoin, faircoin, color='k' )

ax.legend(loc='best')

ax.set_title('iDQ local ROC curves for %s at %s'%(plotting_label, ifo))

figname = isp.rocfig( gdbdir, opts.classifier, ifo, filetag, opts.start, opts.end-opts.start )
if opts.verbose:
    print "  %s"%figname
fig.savefig(figname)
isp.plt.close(fig)

if not opts.skip_gracedb_upload:
    message = "iDQ local ROC figure for %s at %s within [%.3f, %.3f]"%(ifo, opts.classifier, opts.start, opts.end)
    if opts.verbose:
        print "    "+message
    gracedb.writeLog( opts.gracedb_id, message=message, filename=figname, tagname=['data_quality'] )

jsonfilename = idq.gdb_roc_json(  gdbdir, opts.classifier, ifo, filetag, opts.start, opts.end-opts.start )
if opts.verbose:
    print "  %s"%jsonfilename
file_obj = open(jsonfilename, "w")
file_obj.write( json.dumps( {opts.classifier:jsonD} ) )
file_obj.close()

if not opts.skip_gracedb_upload:
    message = "iDQ local ROC curves for %s at %s within [%.3f, %.3f]"%(opts.classifier, ifo, opts.start, opts.end)
    if opts.verbose:
        print "    "+message
    gracedb.writeLog( opts.gracedb_id, message=message, filename=jsonfilename )

#=================================================

### determine vital statistics about training and calibration 

### extract trained, calib ranges and associate them with segments
trainD = defaultdict( list )
calibD = defaultdict( list )
for fap in faps:
    seg = idq.extract_start_stop( fap, suffix=".gwf" )
    trained, calib = idq.extract_timeseries_ranges( fap )
    trainD[tuple(trained)].append( seg )
    calibD[tuple(calib)].append( seg )

if opts.ignore_science_segments:
    dq_name = config.get('get_science_segments', 'include')

### fix segments, extract info
jsonD = {}
for (calib_start, calib_end), segs in calibD.items():
    segs = event.fixsegments( segs )
    thisD = {'used':segs}

    ### extract livetime
    sciseg_files = glob.glob("%s/*_%d/science_segments*%s-*-*.xml.gz"%(calibrationdir, calib_end, dq_name))
    if len(sciseg_files) > 1:
        raise ValueError("something odd with %s scisegs in : %s/*_%d/"%(dq_name, calibrationdir, calib_end))
    elif sciseg_files: ### len == 1
        thisD['%s livetime'%dq_name] = event.livetime( idq.extract_dq_segments( sciseg_files[0], dq_name)[0] )
        segstart, segend = idq.extract_start_stop( sciseg_files[0], suffix='.xml.gz')
        thisD['start'] = int(segstart)
        thisD['end'] = int(segend)
    else:
        thisD['%s livetime'%dq_name] = None
        thisD['start'] = None
        thisD['end'] = None

    ### extract Ng, Nc
    ### format is the same for all classifiers, independent of flavor, so we just write the command here instead of delegating...
    rocs = [ _ for _ in glob.glob("%s/*_%d/%s_%s-*-*.uroc"%(calibrationdir, calib_end, ifo, opts.classifier)) if (len(_.split(opts.classifier)[-1].split('-'))==3)]
    if len(rocs) != 1:
        raise ValueError("something odd with %s uroc files in : %s/*_%d/"%(opts.classifier, calibrationdir, calib_end))
    else:
        _, _, _, Nc, Ng = idq.file_to_rcg( rocs[0] )
        thisD['Num glitches'] = int(Ng)
        thisD['Num cleans'] = Nc

    ### add to dictionary
    jsonD["%d_%d"%(calib_start, calib_end)] = thisD

### write jsonD to file
jsonfilename = idq.useSummary_json( gdbdir, ifo, opts.classifier, "%s_calibStats"%filetag, opts.start, opts.end-opts.start)
if opts.verbose:
    print "writing : %s"%jsonfilename
file_obj = open(jsonfilename, "w")
file_obj.write( json.dumps( jsonD ) )
file_obj.close()
if not opts.skip_gracedb_upload:
    message = "iDQ local calibration vital statistics for %s at %s within [%.3f, %.3f]"%(opts.classifier, ifo, opts.start, opts.end)
    if opts.verbose:
        print "    "+message
    gracedb.writeLog( opts.gracedb_id, message=message, filename=jsonfilename, tagname=['data_quality'] )


### repeat for training
jsonD = {}
for (train_start, train_end), segs in trainD.items():
    segs = event.fixsegments( segs )
    thisD = {'used':segs}

    ### extract livetime
    sciseg_files = glob.glob("%s/%d_%d/science_segments*%s*.xml.gz"%(traindir, train_start, train_end, dq_name))
    if len(sciseg_files) > 1:
        raise ValueError("something odd with %s scisegs in : %s/%d_%d/"%(dq_name, traindir, train_start, train_end))
    elif sciseg_files: ### len == 1
        thisD['%s livetime'%dq_name] = event.livetime( idq.extract_dq_segments( sciseg_files[0], dq_name)[0] )
        segstart, segend = idq.extract_start_stop( sciseg_files[0], suffix='.xml.gz')
        thisD['start'] = int(segstart)
        thisD['end'] = int(segend)
    else:
        thisD['%s livetime'%dq_name] = None
        thisD['start'] = None
        thisD['end'] = None

    ### extract Ng, Nc
    if flavor == "ovl":
        vetolist = glob.glob("%s/%d_%d/%s/ovl/*vetolist.eval"%(traindir, train_start, train_end, opts.classifier))
        if len(vetolist) != 1:
            raise ValueError("something odd with %s vetolists : %s/%d_%d/%s/ovl/"%(opts.classifier, traindir, train_start, train_end))
        else:
            thisD['Num glitches'] = int( event.loadstringtable( vetolist[0] )[0][idq.ovl.vD['#gwtrg']] ) 
    elif flavor in idq.mla_flavors:
        pat = glob.glob("%s/%d_%d/%s_mla*.pat"%(traindir, train_start, train_end, ifo))
        if len(pat) != 1:
            raise ValueError("something odd with %s patfile : %s/%d_%d/"%(opts.classifier, traindir, train_start, train_end))
        else:
            i = [float(i) for i in idq.slim_load_datfiles( pat[0], columns=['i'] )]
            Ng = np.sum( i )
            Nc = len(i) - Ng
            thisD['Num glitches'] = int(Ng)
            thisD['Num cleans'] = int(Nc)
    else:
        raise ValueError("do not know how to extract training information for flavor=%s"%(flavor))

    ### add to dictionary
    jsonD["%d_%d"%(train_start, train_end)] = thisD

### write jsonD to file
jsonfilename = idq.useSummary_json( gdbdir, ifo, opts.classifier, "%s_trainStats"%filetag, opts.start, opts.end-opts.start)
if opts.verbose:
    print "writing : %s"%jsonfilename
file_obj = open(jsonfilename, "w")
file_obj.write( json.dumps( jsonD ) )
file_obj.close()
if not opts.skip_gracedb_upload:
    message = "iDQ local training vital statistics for %s at %s within [%.3f, %.3f]"%(opts.classifier, ifo, opts.start, opts.end)
    if opts.verbose:
        print "    "+message
    gracedb.writeLog( opts.gracedb_id, message=message, filename=jsonfilename, tagname=['data_quality'] )


#=================================================

'''
  likelihood ratios
  active channels (over a wider timescale than timeseries?)
  iDQ glitch/clean rates
  else?
'''

