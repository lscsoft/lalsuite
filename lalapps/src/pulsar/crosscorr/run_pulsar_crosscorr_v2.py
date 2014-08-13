#!/usr/bin/env python
# script to call lalapps_pulsar_crosscorr_v2 (for use with e.g. condor)

# Copyright (C) 2014 John Whelan

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with with program; see the file COPYING. If not, write to the
# Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

from __future__ import division
from argparse import ArgumentParser
from subprocess import check_call
from ConfigParser import SafeConfigParser

ap = ArgumentParser()
ap.add_argument("--configFile", action="store", required=True,
                    help='Configuration file containing arguments')
ap.add_argument("--numJobs", action="store", type=int, required=True,
                    help='Number of jobs into which run is split')
ap.add_argument("--jobNum", action="store", type=int, required=True,
                    help='Number of this job')

args = ap.parse_args()

cp = SafeConfigParser()
cp.optionxform = str
cp.read(args.configFile)

fMin = cp.getfloat('param-space','f_min')
fMax = cp.getfloat('param-space','f_max')
fFullBand = fMax-fMin
fBand = fFullBand / args.numJobs
fStart = fMin + args.jobNum * fBand

argstring = ''
argstring += '--fStart %f ' % fStart
argstring += '--fBand %f ' % fBand

# Pass along the rest of the arguments from the ini file
argstring += ' '.join(['--%s %s' % a for a in cp.items('raw-program-arguments')])

program = 'lalapps_pulsar_crosscorr_v2'
print program, argstring
