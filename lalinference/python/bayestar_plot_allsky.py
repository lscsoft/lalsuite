#
# Copyright (C) 2011  Leo Singer
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
__doc__ = """
Plot an all-sky map on a Mollweide projection.
By default, plot in celestial coordinates (RA, Dec).

To plot in geographic coordinates (longitude, latitude) with major
coastlines overlaid, provide the --geo flag.

Public-domain cartographic data is courtesy of Natural Earth
(http://www.naturalearthdata.com) and processed with MapShaper
(http://www.mapshaper.org).
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Command line interface

from optparse import Option, OptionParser
from lalinference.bayestar import command
from matplotlib import cm
colormap_choices = sorted(cm.cmap_d.keys())
parser = OptionParser(
    description = __doc__,
    usage = "%prog [options] [INPUT]",
    option_list = [
        Option("-o", "--output", metavar="FILE.{pdf,png}",
            help="name of output file [default: plot to screen]"),
        Option("--colormap", default="jet", choices=colormap_choices,
            metavar='|'.join(colormap_choices),
            help="name of matplotlib colormap [default: %default]"),
        Option("--figure-width", metavar="INCHES", type=float, default=8.,
            help="width of figure in inches [default: %default]"),
        Option("--figure-height", metavar="INCHES", type=float, default=6.,
            help="height of figure in inches [default: %default]"),
        Option("--dpi", metavar="PIXELS", type=int, default=300,
            help="resolution of figure in dots per inch [default: %default]"),
        Option("--contour", metavar="PERCENT", type=float, action="append",
            default=[], help="plot contour enclosing this percentage of"
            + " probability mass [may be specified multiple times, default: none]"),
        Option("--colorbar", default=False, action="store_true",
            help="Show colorbar [default: %default]"),
        Option("--radec", nargs=2, metavar="RA Dec", type=float, action="append",
            default=[], help="right ascension (deg) and declination (deg) to mark"
            + " [may be specified multiple times, default: none]"),
        Option("--geo", action="store_true", default=False,
            help="Plot in geographic coordinates, (lat, lon) instead of (RA, Dec) [default: %default]")
    ]
)
opts, args = parser.parse_args()
infilename = command.get_input_filename(parser, args)

# Late imports

# Choose a matplotlib backend that is suitable for headless
# rendering if output to file is requested
import matplotlib
if opts.output is not None:
    matplotlib.use('agg')

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import healpy as hp
import lal
from lalinference import fits
from lalinference import plot

fig = plt.figure(figsize=(opts.figure_width, opts.figure_height), frameon=False)
ax = plt.subplot(111,
    projection='mollweide' if opts.geo else 'astro mollweide')
ax.cla()
ax.grid()

skymap, metadata = fits.read_sky_map(infilename)
nside = hp.npix2nside(len(skymap))

if opts.geo:
    dlon = -lal.GreenwichMeanSiderealTime(lal.LIGOTimeGPS(metadata['gps_time'])) % (2*np.pi)
else:
    dlon = 0

# Convert sky map from probability to probability per square degree.
probperdeg2 = skymap / hp.nside2pixarea(nside, degrees=True)

# Plot sky map.
vmax = probperdeg2.max()
plot.healpix_heatmap(probperdeg2, dlon=dlon,
    vmin=0., vmax=vmax, cmap=plt.get_cmap(opts.colormap))

if opts.colorbar:
    # Plot colorbar.
    cb = plot.colorbar(vmax)

    # Set colorbar label.
    cb.set_label(r'prob. per deg$^2$')

# Add contours.
if opts.contour:
    indices = np.argsort(-skymap)
    region = np.empty(skymap.shape)
    region[indices] = 100 * np.cumsum(skymap[indices])
    cs = plot.healpix_contour(region, dlon=dlon,
        colors='k', linewidths=0.5, levels=opts.contour)
    fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
    plt.clabel(cs, fmt=fmt, fontsize=6, inline=True)

if opts.geo:
    geojson_filename = os.path.join(os.path.dirname(plot.__file__),
        'ne_simplified_coastline.json')
    with open(geojson_filename, 'r') as geojson_file:
        geojson = json.load(geojson_file)
    for shape in geojson['geometries']:
        verts = np.deg2rad(shape['coordinates'])
        plt.plot(verts[:, 0], verts[:, 1], color='black', linewidth=1, alpha=0.25)

# Add markers (e.g., for injections or external triggers).
for ra, dec in np.deg2rad(opts.radec):
    # Convert the right ascension to either a reference angle from -pi to pi
    # or a wrapped angle from 0 to 2 pi, depending on the version of Matplotlib.
    # FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
    if opts.geo or plot.mpl_version < '1.2.0':
        ra = plot.reference_angle(ra + dlon)
    else:
        ra = plot.wrapped_angle(ra + dlon)
    ax.plot(ra, dec, '*', markerfacecolor='white', markeredgecolor='black', markersize=10)

# If we are using a new enough version of matplotlib, then
# add a white outline to all text to make it stand out from the background.
plot.outline_text(ax)

fig.patch.set_alpha(0.)
ax.patch.set_alpha(0.)
ax.set_alpha(0.)

if opts.output is None:
    plt.show()
else:
    plt.savefig(opts.output, dpi=opts.dpi)
