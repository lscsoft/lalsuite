#!/usr/bin/python

# Copyright (C) 2013 Ruslan Vaulin
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

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

import os
import sys 
from optparse import *
import glob
from pylal import auxmvc_utils
#from pylal import git_version
import numpy
import random
import event
import idq
#import idq_daily

######################################################################################
#
#                          MAIN
#
######################################################################################

usage = """
        This code uses triggers from GW and auxiliary channels to build auxmvc feature vectors. The vectors are later passed as input to multivariate classifiers.         
		"""
parser = OptionParser(version="", usage=usage) 
parser.add_option("", "--trigger-files", help="glob for files with triggers")
parser.add_option("", "--trigger-dir", help="path to directory containing (possibly sub-directories with) trigger files.")
parser.add_option("", "--time-window",default=0.1,type="float",help="time window to applied on samples (in seconds).") 
parser.add_option("", "--signif-threshold",default=35.0,type="float",help="threshold on trigger significance in the main channel") 
parser.add_option("", "--channels", default=None, help="file with the list of auxiliary channels to be used in the analysis.\
                                                        If not given, all channels in the input data are used.")
parser.add_option("", "--unsafe-channels", default=None, help="file with the list of unsafe channels.")
parser.add_option("", "--main-channel", help="Name of the channel in which artifacts should be detected. \
                                              Triggers from this channel are used as indicators for building auxmvc feature vector.\
											  Typically, main channel is set to GW (DARM) channel of the detector.")
parser.add_option("", "--dq-segments", type="string", default=None, help="xml file with dq segments")
parser.add_option("", "--dq-segments-name", type="string", help="name of dq segments, --include-segments option to \
                                                                ligolw_segment_query that generated the xml segment file.")
parser.add_option("", "--max-glitch-samples", default=None, type="int",help="Upper limit on number of glitch samples in the training set.")
parser.add_option("", "--clean-samples-rate",default=0.1,type="float",help="Rate of clean samples in Hz")
parser.add_option("", "--max-clean-samples", default=None, type="int",help="Upper limit on number of clean samples in the training set.")
parser.add_option("", "--filter-out-unclean", action="store_true", default=False, help="get rid of clean samples that are too close to glitches.")
parser.add_option("", "--gps-start-time", default="", type="string",help="GPS start time of the period  to be processed")
parser.add_option("", "--gps-end-time", default="", type="string",help="GPS end time of the period  to be processed")
parser.add_option("", "--output-file",type="string", default="auxmvc_samples.pat", help="path and name of the output file to save auxmvc feature vectors")
parser.add_option("", "--verbose", action="store_true", default=False, help="run in verbose mode")
(opts,args) = parser.parse_args()


main_channel = opts.main_channel

if opts.gps_start_time and opts.gps_end_time:
	gps_start_time = int(opts.gps_start_time)
	gps_end_time = int(opts.gps_end_time)


# get list of trigger files
if opts.trigger_dir:
	trig_files = idq.get_all_files_in_range(opts.trigger_dir, gps_start_time, gps_end_time, pad=0, suffix='.trg')
else:
	trig_files = glob.glob(opts.trigger_files)

if not trig_files:
	print "Warning: Empty list of trigger files, exiting without doing anything."
	sys.exit(0)



if opts.verbose: 
	print "Loading triggers ..."
# load triggers from files
trigger_dict = event.loadkwm(trig_files) 

if opts.verbose:
	print "Done."	


if not trigger_dict:
	print "Warning: No triggers in the input files, exiting without doing anything."
	sys.exit(0)

if opts.dq_segments:
	# load dq segments
	(dq_segments, covered_segments) = idq.extract_dq_segments(open(opts.dq_segments, "r"), opts.dq_segments_name)
	# sort and merge segments
	dq_segments = event.fixsegments(dq_segments)
else:
	dq_segments = None
  

# construct auxmvc feature vectors
auxmvc_vectors = idq.build_auxmvc_vectors(trigger_dict, main_channel, opts.time_window, opts.signif_threshold, opts.output_file,\
  gps_start_time=gps_start_time, gps_end_time=gps_end_time, channels=opts.channels, unsafe_channels=opts.unsafe_channels,\
  science_segments = dq_segments, clean_samples_rate=opts.clean_samples_rate, filter_out_unclean=opts.filter_out_unclean,\
  max_clean_samples = opts.max_clean_samples, max_glitch_samples = opts.max_glitch_samples)

clean_samples = auxmvc_vectors[numpy.nonzero(auxmvc_vectors['i']==0)[0],:]
glitch_samples = auxmvc_vectors[numpy.nonzero(auxmvc_vectors['i']==1)[0],:]

if opts.verbose:
	print "total number of glitch samples in the set: " + str(len(glitch_samples))
	print "total number of clean samples in the set: " +str(len(clean_samples))






















