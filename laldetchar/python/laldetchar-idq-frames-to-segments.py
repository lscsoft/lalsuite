# Copyright (C) 2015 Reed Essick
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
import logging
import ConfigParser

import numpy as np
from collections import defaultdict

from laldetchar.idq import idq
#from laldetchar.idq import reed as idq

from laldetchar.idq import event
from pylal import Fr

from laldetchar import git_version

from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table

from optparse import OptionParser

#===================================================================================================

__prog__ = 'laldetchar-idq-frames-to-segments'

__author__ = \
    'Reed Essick (<reed.essick@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

description = \
""" 
an executable that knows how to discover idq FAP frames based on an iDQ config file and converts them into segments.
outputs an xml document with ligolw segment tables and supporting metadata
"""

#===================================================================================================

parser = OptionParser(version='Name: %%prog\n%s'%git_version.verbose_msg,
                                usage='%prog [options] gpsstart gpsend',
                                description=description)
parser.add_option('-c', '--config-file', dest='config_file',
    help='configuration file', metavar='FILE', default='idq.ini' )

parser.add_option("", "--realtime-log", default=None, type="string",
    help="if supplied, waits for realtime-log to pass the end time before proceeding.")

parser.add_option("", "--max-wait", default=128, type="float", 
    help="the maximum amount of time we wait while parsing realtime-log")

parser.add_option('-l', '--log-file', default='idq_frames_to_segments.log',
                  type='string', help='log file')

parser.add_option('', '--ignore-science-segments',
    default=False, action="store_true",
    help='generate segments regardless of the science segment content.'    
    )

parser.add_option("", "--no-robot-cert",
    default=False, action="store_true",
    help="do not use robot cert with segment query in training jobs"
    )

parser.add_option('-C', '--classifier', default=[], action='append', type='string', \
    help='generate segments for this classifier. Can generate segments. \
    for more than one classifier by repeating this option'
    )

parser.add_option('-F', '--FAPthr', default=[], action='append', type='float', \
    help='a threshold used to generate segments from time-series. \
    Can generate multiple sets of segments corresponding to \
    different thresholds by repeating this option.' \
    )

parser.add_option('-R', '--right-padding', default=0, type='float', \
    help="transform all segments from [start, end] -> [start, end+right_padding]"
    )
parser.add_option('-L', '--left-padding', default=0, type='float', \
    help="transform all segments from [start, end] -> [start-left_padding, end]"
    )
parser.add_option('-T', '--t-lag', default=0, type='float', \
    help="transform all segments from [start, end] -> [start+t_lag, end+t_lag]"
    )

parser.add_option('-t', '--tag', default="", type="string")
parser.add_option('-o', '--output-dir', default=".", type="string", \
    help="the directory into which the resulting segment files will be written" \
    )

(opts, args) = parser.parse_args()

### ensure output directory exists
if not os.path.exists(opts.output_dir):
    os.makedirs(opts.output_dir)

if opts.tag:
    opts.tag="_"+opts.tag

### parser gps times from arguments
if len(args) != 2:
    raise ValueError("please supply exactly two arguments: gpsstart gpsend")
startgps, endgps = [int(l) for l in args]

#===================================================================================================
### setup logger to record processes
logger = idq.setup_logger('idq_logger', opts.log_file, sys.stdout, format='%(asctime)s %(message)s')

sys.stdout = idq.LogFile(logger)
sys.stderr = idq.LogFile(logger)

#===================================================================================================
### read global configuration file

config = ConfigParser.SafeConfigParser()
config.read(opts.config_file)

mainidqdir = config.get('general', 'idqdir') ### get the main directory where idq pipeline is going to be running.

usertag = config.get('general', 'usertag')

ifo = config.get('general', 'ifo')

#========================
# which classifiers
#========================
### ensure we have a section for each classifier and fill out dictionary of options
classifiersD, mla, ovl = idq.config_to_classifiersD( config )
combinersD, _ = idq.config_to_combinersD( config )
classifiers = sorted(classifiersD.keys() + combinersD.keys())

### check that all requested classifiers exist...
for classifier in opts.classifier:
    if classifier not in classifiers:
        logger.info('WARNING: %s not present in config file : %s. skipping...'%(classifier, opts.config_file))
        opts.classifier.remove(classifier)

opts.classifier = sorted(list(set(opts.classifier)))

### check that FAPthrs make sense
opts.FAPthr = [FAPthr for FAPthr in opts.FAPthr if (FAPthr >=0) and (FAPthr <= 1)]

opts.FAPthr = sorted(list(set(opts.FAPthr)))

#========================
# check that there is anything to do...
#========================

if not opts.classifier:
    logger.info('WARNING: no valid classifiers specified. Nothing to do...')
    sys.exit(0)

if not opts.FAPthr:
    logger.info('WARNING: no valid FAPthr specified. Nothing to do...')
    sys.exit(0)

#========================
# realtime job
#========================
ts_fs = config.getfloat('realtime','sampling_rate') ### sampling frequency for time-series files
ts_dt = 1.0/ts_fs

realtimedir = config.get('general', 'realtimedir')### output directory for realtime predictions

#========================
# science segments 
#========================
if not opts.ignore_science_segments:
    ### load settings for accessing dmt segment files
    dq_name = config.get('get_science_segments', 'include')
    segdb_url = config.get('get_science_segments', 'segdb')

#========================
### wait for realtime.log to pass endgps

if opts.realtime_log:
    logger.info("parsing %s to extract idq-realtime state"%opts.realtime_log)

    realtime_log = open(realtime_log, "r") ### open realtime log for reading
    realtime_log.seek(0, 2) ### go to end of file

    ### wait until realtime has passed gps_end+delay
    past, dead, timed_out = idq.block_until(endgps+opts.left_padding-opts.t_lag, realtime_log, max_wait=opts.max_wait, timeout=2*opts.max_wait)

    if past:
        logger.info("found realtime stride starting after t=%.3f+%.3f-%.3f"%(endgps,opts.left_padding,opts.t_lag))
    elif timed_out:
        logger.info("WARNING: could not find a recent enough stride in %s after searching for %.1f seconds. Realtime process may be behind"%(opts.realtime_log, 2*opts.max_wait))
    else: # dead
        logger.info("WARNING: no new iDQ information was reported to %s after waiting %.1f seconds. Realtime process may be dead."%(opts.realtime_log, opts.max_wait))

#===================================================================================================
### begin the analysis

stride = endgps - startgps

lookup_startgps = startgps - opts.right_padding - opts.t_lag
lookup_endgps = endgps + opts.left_padding - opts.t_lag

#========================
# find science segments
#=======================
if opts.ignore_science_segments:
    logger.info('analyzing data regardless of science segements')
    scisegs = [[startgps, endgps]] ### set segs to be this stride range
    coveredsegs = [[startgps, endgps]] ### set segs to be this stride range

else:
    logger.info('Begin: querrying science segments')

    ### this returns a string
    seg_xml_file = idq.segment_query(config, startgps, endgps, url=segdb_url)

    ### write seg_xml_file to disk
    lsctables.use_in(ligolw.LIGOLWContentHandler)
    xmldoc = ligolw_utils.load_fileobj(seg_xml_file, contenthandler=ligolw.LIGOLWContentHandler)[0]

    ### science segments xml filename
    seg_file = idq.segxml(opts.output_dir, "_%s"%dq_name, startgps, stride)
    logger.info('writing science segments to file : '+seg_file)
    ligolw_utils.write_filename(xmldoc, seg_file, gz=seg_file.endswith(".gz"))

    (scisegs, coveredseg) = idq.extract_dq_segments(seg_file, dq_name) ### read in segments from xml file

### modify scisegs to account for shifts
logger.info('modifying scisegs to account for shifts')
modified_scisegs = [(s-opts.right_padding-opts.t_lag, e+opts.left_padding-opts.t_lag) for s, e in scisegs]

### find idq segments
logger.info('finding idq segments')
idqsegs = idq.get_idq_segments(realtimedir, lookup_startgps, lookup_endgps, suffix='.dat')

logger.info('taking intersection between modified science segments and idq segments')
idqsegs = event.andsegments( [modified_scisegs, idqsegs] )

### write segment file
if opts.ignore_science_segments:
    idqseg_path = idq.idqsegascii(opts.output_dir, '', startgps, stride)
else:
    idqseg_path = idq.idqsegascii(opts.output_dir, '_%s'%dq_name, startgps , stride)
f = open(idqseg_path, 'w')
for seg in idqsegs:
    print >> f, seg[0], seg[1]
f.close()

#========================
# go findeth the frame data
#========================
logger.info('  finding all *fap*.gwf files')
fapsD = defaultdict( list )
for fap in [fap for fap in  idq.get_all_files_in_range(realtimedir, lookup_startgps, lookup_endgps, pad=0, suffix='.gwf') if "fap" in fap]:
    fapsD[idq.extract_fap_name( fap )].append( fap )

### throw away files we will never need
for key in fapsD.keys():
    if key not in opts.classifier: ### throw away unwanted files
        fapsD.pop(key)
    else: ### keep only files that overlap with scisegs
        fapsD[key] = [ fap for fap in fapsD[key] if event.livetime(event.andsegments([idqsegs, [idq.extract_start_stop(fap, suffix='.gwf')]])) ]

#========================
# iterate through classifiers -> generate segments
#========================

### set up xml document
from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables

from glue.ligolw.utils import process

xmldoc = ligolw.Document()
xml_element = ligolw.LIGO_LW()
xmldoc.appendChild( xml_element )

proc_id = process.register_to_xmldoc( xmldoc, __prog__, opts.__dict__, version=__version__ ).process_id

segdef = lsctables.New( lsctables.SegmentDefTable, columns=["process_id", "segment_def_id", "ifos", "name", "version", "comment"] )
segsum = lsctables.New( lsctables.SegmentSumTable, columns=["process_id", "segment_sum_id", "start_time", "start_time_ns", "end_time", "end_time_ns", "comment", "segment_def_id"] )
segtab = lsctables.New( lsctables.SegmentTable, columns=["process_id", "segment_def_id", "segment_id", "start_time", "start_time_ns", "end_time", "end_time_ns"] )

xml_element.appendChild( segdef )
xml_element.appendChild( segsum )
xml_element.appendChild( segtab )

### iterate through classifiers
for classifier in opts.classifier:
    logger.info('Begin: generating segments for %s'%classifier)

    faps = fapsD[classifier]
    logger.info('  found %d files'%(len(faps)))

    ### need to load in time-series from frames here!
    chan = idq.channame(ifo, classifier, "%s_fap"%usertag)
    t, ts = idq.combine_gwf(faps, [chan]) ### loads in the data from frames
   
    logger.info('  found %d continous segments'%(len(t)))
 
    ### set up segdef row
    fap2segdef_id = {}
    for FAPthr in opts.FAPthr:
        segdef_id = segdef.get_next_id()
        
        segdef_row = lsctables.SegmentDef()
        segdef_row.process_id = proc_id
        segdef_row.segment_def_id = segdef_id
        segdef_row.ifos = ifo
        segdef_row.name = classifier
        segdef_row.version = 1
        segdef_row.comment = 'FAPthr=%.9e'%FAPthr

        segdef.append( segdef_row )

        fap2segdef_id[FAPthr] = segdef_id

    ### iterate through contiguous segments:
    for T, TS in zip(t, ts):
        if not len(T):
            logger.info('len(T)=0, skipping...')
            continue

        dT = T[1]-T[0] ### assumes we have more than one sample...
        start = max(scisegs[0][0], T[0]+opts.t_lag)
        end = min(scisegs[-1][-1], T[-1]+dT+opts.t_lag)

        start_time = int(start)
        start_time_ns = (start-start_time)*1e9
        end_time = int(end)
        end_time_ns = (end-end_time)*1e9

        dur = end_time - start_time

        ### iterate through FAPthr
        for FAPthr in opts.FAPthr:
            print "    FAPthr : %.6e"%(FAPthr)

            segdef_id = fap2segdef_id[FAPthr]
            ### set up segment summary row
            segsum_row = lsctables.SegmentSum()

            segsum_row.segment_sum_id = segsum.get_next_id()

            segsum_row.comment = ''
            segsum_row.process_id = proc_id
            segsum_row.segment_def_id = segdef_id
            segsum_row.start_time    = start_time
            segsum_row.start_time_ns = start_time_ns
            segsum_row.end_time      = end_time
            segsum_row.end_time_ns   = end_time_ns

            segsum.append( segsum_row )

            ### generate segments for this threshold
            segs, min_TS = idq.timeseries_to_segments(T[:], -TS[:], -FAPthr)

            if opts.right_padding != 0: ### transform all segments from [start, end] -> [start, end+right_padding]
                logger.info('      moving right edge of segments: end->end+%.6f'%(opts.right_padding))
                for ind, (s, e) in enumerate(segs):
                    segs[ind] = [s, e + opts.right_padding]

            if opts.left_padding != 0: ### transform all segments from [start, end] -> [start-left_padding, end]
                logger.info('      moving left edge of segments: start->start+%.6f'%(opts.right_padding))
                for ind, (s, e) in enumerate(segs):
                    segs[ind] = [s - opts.left_padding, e]

            if (opts.right_padding != 0) or (opts.left_padding != 0): ### clean up segments
                logger.info('      ensuring segments still make sense after moving the edges')
                good_segs = []
                for (s, e) in segs:
                    if e > s: ### require the segments to still make sense after we fustz with the left/right edges
                        good_segs.append( [s, e] )

            if opts.t_lag != 0: ### transform all segments from [start, end] -> [start+t_lag, end+t_lag]
                logger.info('      shifting all segments: [s, e]->[s+%.6f, e+%.6f]'%(opts.t_lag, opts.t_lag))
                segs = np.array(segs) + opts.t_lag

            ### take intersection of segs and scisegs
            logger.info('taking intersection of vetosegs and scisegs')
            segs = event.andsegments( [scisegs, segs] )

            logger.info('    writing segments to SegmentTable')
            ### add segments to table 
            for s, e in segs:
                seg_row = lsctables.Segment()

                seg_row.segment_id = segtab.get_next_id()
                
                seg_row.process_id = proc_id
                seg_row.segment_def_id = segdef_id
                seg_row.start_time = int(s)
                seg_row.start_time_ns = (s-seg_row.start_time)*1e9
                seg_row.end_time = int(e)
                seg_row.end_time_ns = (e-seg_row.end_time)*1e9

                segtab.append( seg_row )

    logger.info('Done: generating segments for %s'%classifier)

### write to disk
filename = idq.frame2segmentxml( opts.output_dir, opts.tag, startgps, endgps-startgps )
logger.info('writing segments to : %s'%filename)
ligolw_utils.write_filename( xmldoc, filename, gz=filename.endswith('.gz') )
