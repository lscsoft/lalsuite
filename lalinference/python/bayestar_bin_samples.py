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
from lalinference.io import fits
from lalinference.bayestar import distance
from lalinference.bayestar import moc
from lalinference.bayestar.sky_map import derasterize
from lalinference.healpix_tree import adaptive_healpix_histogram

samples = Table.read(opts.input, format='ascii')
theta = 0.5*np.pi - samples['dec']
phi = samples['ra']
if 'distance' in samples.colnames:
    samples.rename_column('distance', 'dist')

p = adaptive_healpix_histogram(
    theta, phi, opts.samples_per_bin,
    nside=opts.nside, max_nside=opts.max_nside, nest=True)


def diststats(samples, max_nside, nside, ipix):
    step = (max_nside // nside) ** 2
    i0 = np.searchsorted(samples['ipix'], step * ipix)
    i1 = np.searchsorted(samples['ipix'], step * (ipix + 1))
    if i0 == i1:
        return np.inf, 0.0
    else:
        dist = samples['dist'][i0:i1]
        return np.mean(dist), np.std(dist)

if 'dist' in samples.colnames:
    p = derasterize(Table([p], names=['PROB']))
    order, ipix = moc.uniq2nest(p['UNIQ'])
    nside = hp.order2nside(order.astype(int))
    max_order = order.max().astype(int)
    max_nside = hp.order2nside(max_order)
    samples['ipix'] = hp.ang2pix(max_nside, theta, phi, nest=True)
    samples.sort('ipix')
    distmean, diststd = np.transpose(
        [diststats(samples, max_nside, n, i) for n, i in zip(nside, ipix)])

    p['DISTMU'], p['DISTSIGMA'], p['DISTNORM'] = \
        distance.moments_to_parameters(distmean, diststd)

    # Add marginal distance moments
    p.meta['distmean'] = np.mean(samples['dist'])
    p.meta['diststd'] = np.std(samples['dist'])

# Write output to FITS file.
fits.write_sky_map(
    opts.output, p, nest=True,
    creator=parser.prog, objid=opts.objid, gps_time=samples['time'].mean())
