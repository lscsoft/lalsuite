#
# Copyright (C) 2011-2017  Leo Singer
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
from __future__ import division


# Command line interface.
from lalinference.bayestar import command
parser = command.ArgumentParser()
parser.add_argument(
    '--nside', '-n', type=int, default=-1,
    help='Write final sky map at this HEALPix lateral resolution '
    '[default: auto]')
parser.add_argument(
    '--max-nside', '-m', type=int, default=-1,
    help='Stop subdivision at this maximum HEALPix lateral resolution '
    '[default: max 64-bit]')
parser.add_argument(
    '--samples-per-bin', type=int, default=30,
    help='Samples per bin [default: %(default)s]')
parser.add_argument(
    '--two-step', action='store_true', default=False,
    help='Partition the samples into one sub-population for laying out bin '
    'boundaries and one sub-population for evaluating densities instead of '
    'using the whole population for both steps')
parser.add_argument(
    '--objid',
    help='Event ID to be stored in FITS header [default: %(default)s]')
parser.add_argument(
    'input', metavar='INPUT.dat', default='-',
    help='name of input posterior samples file [default: stdin]')
parser.add_argument(
    '-o', '--output', metavar='OUTPUT.fits[.gz]', required=True,
    help='name of output FITS file [required]')
opts = parser.parse_args()


# Late imports.
import healpy as hp
import numpy as np
from astropy.table import Table
from astropy.utils.misc import NumpyRNGContext
from lalinference.io import fits
from lalinference.io import hdf5
from lalinference.bayestar import distance
from lalinference.bayestar import moc
from lalinference.bayestar.sky_map import derasterize
from lalinference.healpix_tree import adaptive_healpix_histogram


def pixstats(samples, max_nside, nside, ipix):
    step = (max_nside // nside) ** 2
    i0 = np.searchsorted(samples['ipix'], step * ipix)
    i1 = np.searchsorted(samples['ipix'], step * (ipix + 1))
    n = i1 - i0
    if n == 0:
        if 'dist' in samples.colnames:
            return 0.0, np.inf, 0.0
        else:
            return 0.0
    else:
        probdensity = n / len(samples) / hp.nside2pixarea(nside)
        if 'dist' in samples.colnames:
            dist = samples['dist'][i0:i1]
            return probdensity, np.mean(dist), np.std(dist)
        else:
            return probdensity


# Read samples.
samples = hdf5.read_samples(opts.input)

# If two-step binning is requested, then randomly partition the samples.
if opts.two_step:
    with NumpyRNGContext(0):  # Fix random seed to make results reproducible
        ranking_samples, samples = np.array_split(
            np.random.permutation(samples), 2)
    ranking_samples = Table(ranking_samples)
    hist_samples = Table(samples)
else:
    ranking_samples = hist_samples = samples

# Place the histogram bins.
theta = 0.5*np.pi - hist_samples['dec']
phi = hist_samples['ra']
p = adaptive_healpix_histogram(
    theta, phi, opts.samples_per_bin,
    nside=opts.nside, max_nside=opts.max_nside, nest=True)

# Evaluate per-pixel density.
p = derasterize(Table([p], names=['PROB']))
order, ipix = moc.uniq2nest(p['UNIQ'])
nside = hp.order2nside(order.astype(int))
max_order = order.max().astype(int)
max_nside = hp.order2nside(max_order)
ranking_samples['ipix'] = hp.ang2pix(max_nside, theta, phi, nest=True)
ranking_samples.sort('ipix')
result = np.transpose(
    [pixstats(ranking_samples, max_nside, n, i) for n, i in zip(nside, ipix)])

# Add distance info if necessary.
if 'dist' in ranking_samples.colnames:
    p['PROBDENSITY'], distmean, diststd = result
    p['DISTMU'], p['DISTSIGMA'], p['DISTNORM'] = \
        distance.moments_to_parameters(distmean, diststd)

    # Add marginal distance moments
    p.meta['distmean'] = np.mean(samples['dist'])
    p.meta['diststd'] = np.std(samples['dist'])
else:
    p['PROBDENSITY'] = result

# Write output to FITS file.
fits.write_sky_map(
    opts.output, p, nest=True,
    creator=parser.prog, objid=opts.objid, gps_time=samples['time'].mean())
