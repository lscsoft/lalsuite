#
# Copyright (C) 2017  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Average several sky maps according to the Burst group established by:
 * LIGO-T1700050-v7 (https://dcc.ligo.org/LIGO-T1700050-v7)
 * LIGO-T1700078-v2 (https://dcc.ligo.org/LIGO-T1700078-v2)

The sky maps are first upsampled to a common resolution and then averaged
pixel-by-pixel. Some FITS header values are copied from the first input
sky map.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Command line interface
import argparse
from lalinference.bayestar import command
parser = command.ArgumentParser()
parser.add_argument(
    'input', metavar='INPUT.fits[.gz]', type=argparse.FileType('rb'),
    nargs='+', help='Input FITS files')
parser.add_argument(
    '-o', '--output', metavar='OUTPUT.fits[.gz]',
    help='Ouptut FITS file [default: derived from input filenames]')
opts = parser.parse_args()


# Late imports
import os.path
import healpy as hp
import numpy as np
from lalinference.io import fits


# Read all input sky maps
probs, metas = zip(*(fits.read_sky_map(file.name, nest=True)
                     for file in opts.input))

# Determine output HEALPix resolution
npix = max(len(prob) for prob in probs)
nside = hp.npix2nside(npix)

# Average sky maps
prob = np.mean([hp.ud_grade(prob, nside,
                            order_in='NESTED', order_out='NESTED',
                            power=-2) for prob in probs], axis=0)

# Describe the processing of this file
history = ['Arithmetic mean of the following {} sky maps:'.format(len(probs))]
history += ['    ' + file.name for file in opts.input]

# Create metadata dictionary for FITS header
meta = dict(creator=parser.prog, nest=True, history=history)

# Copy some header values from the first sky map
copy_keys = ('objid', 'url', 'instruments', 'gps_time')
meta.update({key: metas[0][key] for key in copy_keys if key in metas[0]})

# Determine output filename
if opts.output is not None:
    output = opts.output
else:
    output = '_plus_'.join(os.path.split(file.name)[1].partition('.fits')[0]
                      for file in opts.input) + '.fits.gz'

# Write output filename to stdout
print(output)

# Write output file
fits.write_sky_map(output, prob, **meta)
