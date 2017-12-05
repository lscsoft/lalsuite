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
Overplot contours for a large number of sky maps, in geographical coordinates.
This can reveal patterns due to the priors (i.e., the network antenna pattern).
"""


# Command line interface

from lalinference.bayestar import command
parser = command.ArgumentParser(parents=[command.figure_parser])
parser.add_argument(
    '--contour', metavar='PERCENT', type=float, default=90,
    help='plot contour enclosing this percentage of'
    ' probability mass [default: %(default)s]')
parser.add_argument(
    '--alpha', metavar='ALPHA', type=float, default=0.1,
    help='alpha blending for each sky map [default: %(default)s]')
parser.add_argument(
    'fitsfilenames', metavar='GLOB.fits[.gz]', nargs='+', action='glob',
    help='Input FITS filenames and/or globs')
parser.set_defaults(colormap=None)
opts = parser.parse_args()

# Late imports

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import healpy as hp
import lal
from lalinference.io import fits
from lalinference import plot
from glue.text_progress_bar import ProgressBar


fig = plt.figure(frameon=False)
ax = plt.axes(projection='mollweide')
ax.grid()

progress = ProgressBar()

progress.max = len(opts.fitsfilenames)

matplotlib.rc('path', simplify=True, simplify_threshold=1)

if opts.colormap is None:
    colors = ['k'] * len(opts.fitsfilenames)
else:
    colors = matplotlib.cm.get_cmap(opts.colormap)
    colors = colors(np.linspace(0, 1, len(opts.fitsfilenames)))
for count_records, (color, fitsfilename) in enumerate(
        zip(colors, opts.fitsfilenames)):
    progress.update(count_records, fitsfilename)
    skymap, metadata = fits.read_sky_map(fitsfilename, nest=None)
    nside = hp.npix2nside(len(skymap))
    gmst = lal.GreenwichMeanSiderealTime(metadata['gps_time']) % (2*np.pi)

    indices = np.argsort(-skymap)
    region = np.empty(skymap.shape)
    region[indices] = 100 * np.cumsum(skymap[indices])
    plot.healpix_contour(
        region, nest=metadata['nest'], dlon=-gmst,
        colors=[color], linewidths=0.5, levels=[opts.contour],
        alpha=opts.alpha)

progress.update(-1, 'saving figure')

# Add a white outline to all text to make it stand out from the background.
plot.outline_text(ax)

opts.output()
