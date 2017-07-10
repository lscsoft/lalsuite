#
# Copyright (C) 2012-2017  Leo Singer
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
Plotting tools for drawing HEALPix images
"""
from __future__ import division

import functools
import matplotlib
from distutils.version import LooseVersion
from matplotlib import text
from matplotlib import ticker
from matplotlib import patheffects
import numpy as np
import healpy as hp

__all__ = ('heatmap', 'contour', 'contourf', 'healpix_heatmap',
           'healpix_contour', 'healpix_contourf', 'colorbar', 'outline_text')


def heatmap(func, *args, **kwargs):
    "Plot a function on the sphere using the current geographic projection."""
    from matplotlib import pyplot as plt

    # Get current axis.
    ax = plt.gca()

    # Set up a regular grid tiling the bounding box of the axes.
    x = np.arange(ax.bbox.x0, ax.bbox.x1 + 0.5, 0.5)
    y = np.arange(ax.bbox.y0, ax.bbox.y1 + 0.5, 0.5)
    xx, yy = np.meshgrid(x, y)

    # Get axis data limits.
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Retrieve the inverse transform of the current axes (which converts
    # display coodinates to data coordinates).
    itrans = ax.transData.inverted()

    # Get the longitude and latitude of every point in the bounding box.
    lons, lats = itrans.transform(np.column_stack((xx.ravel(), yy.ravel()))).T

    # Create a mask that selects only the pixels that fall inside the map
    # boundary.
    mask = \
        np.isfinite(lons) & np.isfinite(lats) & (lons >= xmin) & (lons <= xmax)
    zz = np.ma.array(np.empty(lons.shape), mask=~mask)

    # Evaluate the function everywhere that the mask is set.
    zz[mask] = func(lons[mask], lats[mask])

    # Plot bitmap using imshow.
    if LooseVersion(matplotlib.__version__) < LooseVersion('2.0'):
        # FIXME: workaround for old behavior of imshow().
        # Remove this once we require matplotlib >= 2.0.
        # See also:
        #   * https://bugs.ligo.org/redmine/issues/5152
        #   * https://github.com/matplotlib/matplotlib/issues/7903
        aximg = plt.imshow(
            zz.reshape(xx.shape), aspect=ax.get_aspect(),
            interpolation='nearest',
            origin='upper', extent=(xmin, xmax, ymax, ymin), *args, **kwargs)
    else:
        aximg = plt.imshow(
            zz.reshape(xx.shape), aspect=ax.get_aspect(),
            interpolation='nearest',
            origin='upper', extent=(0, 1, 1, 0), transform=ax.transAxes,
            *args, **kwargs)

    # Hide masked-out values by displaying them in transparent white.
    aximg.cmap.set_bad('w', alpha=0.)

    # Done.
    return aximg


def contour(func, *args, **kwargs):
    "Plot a function on the sphere using the current geographic projection."""
    from matplotlib import pyplot as plt

    # Get current axis.
    ax = plt.gca()

    # Set up a regular grid tiling in right ascension and declination
    x = np.linspace(*ax.get_xlim(), num=500)
    y = np.linspace(*ax.get_ylim(), num=500)
    xx, yy = np.meshgrid(x, y)

    # Evaluate the function everywhere.
    zz = func(xx, yy)

    # Add contour plot
    ax = plt.contour(xx, yy, zz, *args, **kwargs)

    # Done.
    return ax


def contourf(func, *args, **kwargs):
    "Plot a function on the sphere using the current geographic projection."""
    from matplotlib import pyplot as plt

    # Get current axis.
    ax = plt.gca()

    # Set up a regular grid tiling in right ascension and declination
    x = np.linspace(*ax.get_xlim(), num=500)
    y = np.linspace(*ax.get_ylim(), num=500)
    xx, yy = np.meshgrid(x, y)

    # Evaluate the function everywhere.
    zz = func(xx, yy)

    # Add contour plot
    ax = plt.contourf(xx, yy, zz, *args, **kwargs)

    # Done.
    return ax


def _healpix_lookup(map, lon, lat, nest=False, dlon=0):
    """Look up the value of a HEALPix map in the pixel containing the point
    with the specified longitude and latitude."""
    nside = hp.npix2nside(len(map))
    return map[hp.ang2pix(nside, 0.5 * np.pi - lat, lon - dlon, nest=nest)]


def healpix_heatmap(map, *args, **kwargs):
    """Produce a heatmap from a HEALPix map."""
    mpl_kwargs = dict(kwargs)
    dlon = mpl_kwargs.pop('dlon', 0)
    nest = mpl_kwargs.pop('nest', False)
    return heatmap(
        functools.partial(_healpix_lookup, map, nest=nest, dlon=dlon),
        *args, **mpl_kwargs)


def healpix_contour(map, *args, **kwargs):
    """Produce a contour plot from a HEALPix map."""
    mpl_kwargs = dict(kwargs)
    dlon = mpl_kwargs.pop('dlon', 0)
    nest = mpl_kwargs.pop('nest', False)
    return contour(
        functools.partial(_healpix_lookup, map, nest=nest, dlon=dlon),
        *args, **mpl_kwargs)


def healpix_contourf(map, *args, **kwargs):
    """Produce a contour plot from a HEALPix map."""
    mpl_kwargs = dict(kwargs)
    dlon = mpl_kwargs.pop('dlon', 0)
    nest = mpl_kwargs.pop('nest', False)
    return contourf(
        functools.partial(_healpix_lookup, map, nest=nest, dlon=dlon),
        *args, **mpl_kwargs)


def colorbar():
    from matplotlib import pyplot as plt

    usetex = matplotlib.rcParams['text.usetex']
    locator = ticker.AutoLocator()
    formatter = ticker.ScalarFormatter(useMathText=not usetex)
    formatter.set_scientific(True)
    formatter.set_powerlimits((1e-1, 100))

    # Plot colorbar
    cb = plt.colorbar(
        orientation='horizontal', shrink=0.4,
        ticks=locator, format=formatter)

    if cb.orientation == 'vertical':
        axis = cb.ax.yaxis
    else:
        axis = cb.ax.xaxis

    # Move order of magnitude text into last label.
    ticklabels = [label.get_text() for label in axis.get_ticklabels()]
    # Avoid putting two '$' next to each other if we are in tex mode.
    if usetex:
        fmt = '{{{0}}}{{{1}}}'
    else:
        fmt = u'{0}{1}'
    ticklabels[-1] = fmt.format(ticklabels[-1], formatter.get_offset())
    axis.set_ticklabels(ticklabels)
    last_ticklabel = axis.get_ticklabels()[-1]
    last_ticklabel.set_horizontalalignment('left')

    # Draw edges in colorbar bands to correct thin white bands that
    # appear in buggy PDF viewers. See:
    # https://github.com/matplotlib/matplotlib/pull/1301
    cb.solids.set_edgecolor("face")

    # Done.
    return cb


def outline_text(ax):
    """Add a white outline to all text to make it stand out from the
    background."""
    effects = [patheffects.withStroke(linewidth=2, foreground='w')]
    for artist in ax.findobj(text.Text):
        artist.set_path_effects(effects)
