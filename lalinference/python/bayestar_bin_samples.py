#
# Copyright (C) 2011-2015  Leo Singer
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
Convert a list of posterior samples to a HEALPix FITS image using an adaptively
refined HEALPix tree, subdividing each node if it contains more than
--samples-per-bin posterior samples. By default, the resolution is set to the
deepest resolution of the tree.

If an input filename is specified, then the posterior samples are read from it.
If no input filename is provided, then input is read from stdin.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Command line interface

from lalinference.bayestar import command
parser = command.ArgumentParser()
parser.add_argument('--nside', '-n', type=int, default=-1,
    help='Write final sky map at this HEALPix lateral resolution [default: auto]')
parser.add_argument('--max-nside', '-m', type=int, default=-1,
    help='Stop subdivision at this maximum HEALPix lateral resolution [default: max 64-bit]')
parser.add_argument('--samples-per-bin', type=int, default=30,
    help='Samples per bin [default: %(default)s]')
parser.add_argument('--objid',
    help='Event ID to be stored in FITS header [default: %(default)s]')
parser.add_argument('input', metavar='INPUT.dat', default='-',
    help='name of input posterior samples file [default: stdin]')
parser.add_argument('-o', '--output', metavar='OUTPUT.fits[.gz]', required=True,
    help='name of output FITS file [required]')
opts = parser.parse_args()


# Late imports.
import numpy as np
from lalinference.io import fits
import lalinference.bayestar.postprocess

samples = np.recfromtxt(opts.input, names=True)
theta = 0.5*np.pi - samples['dec']
phi = samples['ra']

p = lalinference.bayestar.postprocess.adaptive_healpix_histogram(
    theta, phi, opts.samples_per_bin,
    nside=opts.nside, max_nside=opts.max_nside)

# Write output to FITS file.
fits.write_sky_map(opts.output, p,
    creator=parser.prog, objid=opts.objid,
    gps_time=samples['time'].mean())
