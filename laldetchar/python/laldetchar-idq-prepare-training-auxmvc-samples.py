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
from laldetchar.idq import auxmvc_utils
import numpy
from laldetchar.idq import event
from laldetchar.idq import idq

from laldetchar import git_version

__author__ = 'Ruslan Vaulin <ruslan.vaulin@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

######################################################################################
#
#                          MAIN
#
######################################################################################

description = \
    """This program finds *.pat files that were generated during the low-latency evaluation step and that fall  within the specified time segment. It merges them into a single file used for training of classifiers."""
parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('', '--source-directory',
                  help='path to where *.pat files of low-latency evaluation step are located'
                  )
parser.add_option('', '--dq-segments', type='string', default='',
                  help='xml file with dq segments')
parser.add_option('', '--dq-segments-name', type='string',
                  help='name of dq segments, --include-segments option given to ligolw_segment_query during generation of the xml dq segment file.'
                  )
parser.add_option('', '--gps-start-time', default='', type='string',
                  help='gps start time of the period  to be processed')
parser.add_option('', '--gps-end-time', default='', type='string',
                  help='gps end time of the period  to be processed')
parser.add_option('', '--max-clean-samples', default=None, type='int',
                  help='upper limit on number of clean samples in the training set'
                  )
parser.add_option('', '--max-glitch-samples', default=None, type='int',
                  help='upper limit on number of glitch samples in the training set.'
                  )
parser.add_option('', '--output-file', type='string',
                  default='training_samples.pat',
                  help='full path and name of the output file into which the training samples will be saved'
                  )
parser.add_option('', '--verbose', action='store_true', default=False,
                  help='run in verbose mode')

(opts, args) = parser.parse_args()

gps_start_time = int(opts.gps_start_time)
gps_end_time = int(opts.gps_end_time)

# get all *.pat files in the specififed range

patfiles = idq.get_all_files_in_range(opts.source_directory,
        gps_start_time, gps_end_time, suffix='.pat')

if len(patfiles) == 0:
    print 'No *.pat files found in the gps range ' \
        + str(gps_start_time) + ' - ' + str(gps_end_time)
    print 'Exiting with status 2.'
    sys.exit(2)

# load auxmvc vector samples

auxmvc_samples = auxmvc_utils.ReadMVSCTriggers(patfiles,
        Classified=False)

if opts.dq_segments:

    # load dq segments

    (dq_segments, covered_segments) = \
        idq.extract_dq_segments(open(opts.dq_segments, 'r'),
                                opts.dq_segments_name)

    # sort and merge segments

    dq_segments = event.fixsegments(dq_segments)

    # keep only samples that fall in the dq_segments

    auxmvc_samples = \
        auxmvc_utils.get_samples_in_segments(auxmvc_samples,
            dq_segments)

# get clean and glitch samples

random_samples = auxmvc_samples[numpy.nonzero(auxmvc_samples['i']
                                == 0)[0], :]
clean_samples = auxmvc_utils.get_clean_samples(random_samples)

# apply upper limit on the number of clean samples if given

if opts.max_clean_samples:
    if len(clean_samples) > opts.max_clean_samples:
        clean_samples = clean_samples[-opts.max_clean_samples:]

glitch_samples = auxmvc_samples[numpy.nonzero(auxmvc_samples['i']
                                == 1)[0], :]

# apply upper limit on the number of glitch samples if given

if opts.max_glitch_samples:
    if len(glitch_samples) > opts.max_glitch_samples:
        glitch_samples = glitch_samples[-opts.max_glitch_samples:]

if opts.verbose:
    print 'total number of glitch samples in training set: ' \
        + str(len(glitch_samples))
    print 'total number of clean samples in training set: ' \
        + str(len(clean_samples))

auxmvc_samples = numpy.concatenate((glitch_samples, clean_samples))

# construct name for MVSC training pat file....

output_file_name = opts.output_file

# save training samples into a file

auxmvc_utils.WriteMVSCTriggers(auxmvc_samples,
                               output_filename=output_file_name,
                               Classified=False)

