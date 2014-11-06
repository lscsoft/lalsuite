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

from optparse import *
from laldetchar.idq import auxmvc_utils
import numpy

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
    """
        This code takes *.pat and *.dat files (input and output of the MLA classifier evaluation) and combines them into a single output file adding the variables that were  excluded during classification. By default it overwrites the original *.dat file.
        """

parser = OptionParser(version='Name: %%prog\n%s'
                      % git_version.verbose_msg, usage='%prog [options]'
                      , description=description)
parser.add_option('', '--pat-file', type='string',
                  help='*.pat file, containing unclassified samples'
                  )
parser.add_option('', '--dat-file', type='string',
                  help='*.dat file, containing classified samples'
                  )
parser.add_option('', '--excluded-variables', type='string',
                  help='comma separated list of variables that were excluded during classification'
                  )
parser.add_option('', '--output-file', type='string', default=None,
                  help='output file')
(opts, args) = parser.parse_args()

# load samples from pat file

auxmvc_pat_samples = auxmvc_utils.ReadMVSCTriggers([opts.pat_file],
        Classified=False)

# load samples from dat file

auxmvc_dat_samples = auxmvc_utils.ReadMVSCTriggers([opts.dat_file],
        Classified=True)

excluded_vars = opts.excluded_variables.split(',')

vars = []

if 'GPS_s' and 'GPS_ms' in excluded_vars:
    vars.append('GPS')
    excluded_vars.remove('GPS_s')
    excluded_vars.remove('GPS_ms')

vars.extend(['i', 'w'])

if 'unclean' in excluded_vars:
    vars.append('unclean')
    excluded_vars.remove('unclean')

for var in excluded_vars:
    vars.append(var)

excluded_vars.append('unclean')

vars.append('rank')

for var in auxmvc_dat_samples.dtype.names:
    if not var in ['index', 'i', 'w', 'Bagger']:
        vars.append(var)

datformats = []

for var in vars:
    if var in ['i', 'unclean']:
        datformats.append('i')
    else:
        datformats.append('g8')

dat_samples = numpy.empty(auxmvc_dat_samples.shape,
                          dtype={'names': vars, 'formats': datformats})

dat_samples['GPS'] = auxmvc_pat_samples['GPS_s'] \
    + auxmvc_pat_samples['GPS_ms'] * 10 ** -3

vars.remove('GPS')

dat_samples['rank'] = auxmvc_dat_samples['Bagger']

vars.remove('rank')

for var in vars:
    if var in excluded_vars:
        dat_samples[var] = auxmvc_pat_samples[var]
    else:
        dat_samples[var] = auxmvc_dat_samples[var]

if not opts.output_file:
    opts.output_file = opts.dat_file

auxmvc_utils.WriteMVSCTriggers(dat_samples,
                               output_filename=opts.output_file,
                               Classified=True)

