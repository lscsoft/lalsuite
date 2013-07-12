#
# Copyright (C) 2012  Leo Singer
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
Plotting tools for drawing skymaps
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"
__all__ = ("AstroMollweideAxes", "reference_angle", "make_rect_poly", "heatmap")


import warnings
import functools

# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
import matplotlib
import distutils.version
mpl_version = distutils.version.LooseVersion(matplotlib.__version__)

from matplotlib.axes import Axes
from matplotlib import text
from matplotlib.ticker import Formatter, FixedFormatter, FixedLocator
from matplotlib.projections import projection_registry
from matplotlib.transforms import Transform, Affine2D
from matplotlib.projections.geo import MollweideAxes
from mpl_toolkits.basemap import _geoslib as geos
from matplotlib import pyplot as plt
import numpy as np
import healpy as hp


# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
if mpl_version >= '1.3.0':
    FixedMollweideAxes = MollweideAxes
elif mpl_version < '1.2.0':
    class FixedMollweideAxes(MollweideAxes):
        """Patched version of matplotlib's Mollweide projection that implements a
        correct inverse transform."""

        name = 'fixed mollweide'

        class FixedMollweideTransform(MollweideAxes.MollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.InvertedFixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

        class InvertedFixedMollweideTransform(MollweideAxes.InvertedMollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.FixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

            def transform(self, xy):
                x = xy[:, 0:1]
                y = xy[:, 1:2]

                sqrt2 = np.sqrt(2)
                sintheta = y / sqrt2
                with np.errstate(invalid='ignore'):
                    costheta = np.sqrt(1. - 0.5 * y * y)
                longitude = 0.25 * sqrt2 * np.pi * x / costheta
                latitude = np.arcsin(2 / np.pi * (np.arcsin(sintheta) + sintheta * costheta))
                return np.concatenate((longitude, latitude), 1)
            transform.__doc__ = Transform.transform.__doc__

        def _get_core_transform(self, resolution):
            return self.FixedMollweideTransform(resolution)
else:
    class FixedMollweideAxes(MollweideAxes):
        """Patched version of matplotlib's Mollweide projection that implements a
        correct inverse transform."""

        name = 'fixed mollweide'

        class FixedMollweideTransform(MollweideAxes.MollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.InvertedFixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

        class InvertedFixedMollweideTransform(MollweideAxes.InvertedMollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.FixedMollweideTransform(self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

            def transform_non_affine(self, xy):
                x = xy[:, 0:1]
                y = xy[:, 1:2]

                sqrt2 = np.sqrt(2)
                sintheta = y / sqrt2
                with np.errstate(invalid='ignore'):
                    costheta = np.sqrt(1. - 0.5 * y * y)
                longitude = 0.25 * sqrt2 * np.pi * x / costheta
                latitude = np.arcsin(2 / np.pi * (np.arcsin(sintheta) + sintheta * costheta))
                return np.concatenate((longitude, latitude), 1)
            transform_non_affine.__doc__ = Transform.transform_non_affine.__doc__

        def _get_core_transform(self, resolution):
            return self.FixedMollweideTransform(resolution)


# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
if mpl_version < '1.2.0':
    class AstroMollweideAxes(FixedMollweideAxes):
        """Mollweide axes with phi axis flipped and in hours instead of degrees."""

        name = 'astro mollweide'

        class RaFormatter(Formatter):
            # Copied from matplotlib.geo.GeoAxes.ThetaFormatter and modified
            def __init__(self, round_to=1.0):
                self._round_to = round_to

            def __call__(self, x, pos=None):
                hours = (x / np.pi) * 12.
                hours = round(15 * hours / self._round_to) * self._round_to / 15
                return r"%0.0f$^\mathrm{h}$" % hours

        def set_longitude_grid(self, degrees):
            # Copied from matplotlib.geo.GeoAxes.set_longitude_grid and modified
            number = (360.0 / degrees) + 1
            self.xaxis.set_major_locator(
                FixedLocator(
                    np.linspace(-np.pi, np.pi, number, True)[1:-1]))
            self._longitude_degrees = degrees
            self.xaxis.set_major_formatter(self.RaFormatter(degrees))

        def _set_lim_and_transforms(self):
            # Copied from matplotlib.geo.GeoAxes._set_lim_and_transforms and modified
            FixedMollweideAxes._set_lim_and_transforms(self)

            # This is the transform for latitude ticks.
            yaxis_stretch = Affine2D().scale(np.pi * 2.0, 1.0).translate(-np.pi, 0.0)
            yaxis_space = Affine2D().scale(-1.0, 1.1)
            self._yaxis_transform = \
                yaxis_stretch + \
                self.transData
            yaxis_text_base = \
                yaxis_stretch + \
                self.transProjection + \
                (yaxis_space + \
                 self.transAffine + \
                 self.transAxes)
            self._yaxis_text1_transform = \
                yaxis_text_base + \
                Affine2D().translate(-8.0, 0.0)
            self._yaxis_text2_transform = \
                yaxis_text_base + \
                Affine2D().translate(8.0, 0.0)

        def _get_affine_transform(self):
            transform = self._get_core_transform(1)
            xscale, _ = transform.transform_point((-np.pi, 0))
            _, yscale = transform.transform_point((0, np.pi / 2.0))
            return Affine2D() \
                .scale(0.5 / xscale, 0.5 / yscale) \
                .translate(0.5, 0.5)
else:
    class AstroMollweideAxes(FixedMollweideAxes):
        """Mollweide axes with phi axis flipped and in hours from 24 to 0 instead of
        in degrees from -180 to 180."""

        name = 'astro mollweide'

        def cla(self):
            super(AstroMollweideAxes, self).cla()
            self.set_xlim(0, 2*np.pi)

        def set_xlim(self, *args, **kwargs):
            Axes.set_xlim(self, 0., 2*np.pi)
            Axes.set_ylim(self, -np.pi / 2.0, np.pi / 2.0)

        def _get_core_transform(self, resolution):
            return Affine2D().translate(-np.pi, 0.) + super(AstroMollweideAxes, self)._get_core_transform(resolution)

        class RaFormatter(Formatter):
            # Copied from matplotlib.geo.GeoAxes.ThetaFormatter and modified
            def __init__(self, round_to=1.0):
                self._round_to = round_to

            def __call__(self, x, pos=None):
                hours = (x / np.pi) * 12.
                hours = round(15 * hours / self._round_to) * self._round_to / 15
                return r"%0.0f$^\mathrm{h}$" % hours

        def set_longitude_grid(self, degrees):
            # Copied from matplotlib.geo.GeoAxes.set_longitude_grid and modified
            number = (360.0 / degrees) + 1
            self.xaxis.set_major_locator(
                FixedLocator(
                    np.linspace(0, 2*np.pi, number, True)[1:-1]))
            self._longitude_degrees = degrees
            self.xaxis.set_major_formatter(self.RaFormatter(degrees))

        def _set_lim_and_transforms(self):
            # Copied from matplotlib.geo.GeoAxes._set_lim_and_transforms and modified
            super(AstroMollweideAxes, self)._set_lim_and_transforms()

            # This is the transform for latitude ticks.
            yaxis_stretch = Affine2D().scale(np.pi * 2.0, 1.0)
            yaxis_space = Affine2D().scale(-1.0, 1.1)
            self._yaxis_transform = \
                yaxis_stretch + \
                self.transData
            yaxis_text_base = \
                yaxis_stretch + \
                self.transProjection + \
                (yaxis_space + \
                 self.transAffine + \
                 self.transAxes)
            self._yaxis_text1_transform = \
                yaxis_text_base + \
                Affine2D().translate(-8.0, 0.0)
            self._yaxis_text2_transform = \
                yaxis_text_base + \
                Affine2D().translate(8.0, 0.0)

        def _get_affine_transform(self):
            transform = self._get_core_transform(1)
            xscale, _ = transform.transform_point((0, 0))
            _, yscale = transform.transform_point((0, np.pi / 2.0))
            return Affine2D() \
                .scale(0.5 / xscale, 0.5 / yscale) \
                .translate(0.5, 0.5)


projection_registry.register(AstroMollweideAxes)


def wrapped_angle(a):
    """Convert an angle to a reference angle between 0 and 2*pi."""
    return np.mod(a, 2 * np.pi)


def reference_angle(a):
    """Convert an angle to a reference angle between -pi and pi."""
    a = np.mod(a, 2 * np.pi)
    return np.where(a <= np.pi, a, a - 2 * np.pi)


def reference_angle_deg(a):
    """Convert an angle to a reference angle between -180 and 180 degrees."""
    a = np.mod(a, 360)
    return np.where(a <= 180, a, a - 360)


def subdivide_vertices(vertices, subdivisions):
    """Subdivide a list of vertices by inserting subdivisions additional vertices
    between each original pair of vertices using linear interpolation."""
    subvertices = np.empty((subdivisions * len(vertices), vertices.shape[1]))
    frac = np.atleast_2d(np.arange(subdivisions + 1, dtype=float) / subdivisions).T.repeat(vertices.shape[1], 1)
    for i in range(len(vertices)):
        subvertices[i*subdivisions:(i+1)*subdivisions] = frac[:0:-1, :] * np.expand_dims(vertices[i-1, :], 0).repeat(subdivisions, 0)  + frac[:-1, :] * np.expand_dims(vertices[i, :], 0).repeat(subdivisions, 0)
    return subvertices


# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
def cut_dateline(vertices):
    """Cut a polygon across the dateline, possibly splitting it into multiple
    polygons.  Vertices consist of (longitude, latitude) pairs where longitude
    is always given in terms of a reference angle (between -pi and pi).

    This routine is not meant to cover all possible cases; it will only work for
    convex polygons that extend over less than a hemisphere."""

    out_vertices = []

    # Ensure that the list of vertices does not contain a repeated endpoint.
    if (vertices[0, :] == vertices[-1, :]).all():
        vertices = vertices[:-1, :]

    def count_dateline_crossings(phis):
        n = 0
        for i in range(len(phis)):
            if crosses_dateline(phis[i - 1], phis[i]):
                n += 1
        return n

    def crosses_dateline(phi0, phi1):
        """Test if the segment consisting of v0 and v1 croses the meridian."""
        phi0, phi1 = sorted((phi0, phi1))
        return phi1 - phi0 > np.pi

    dateline_crossings = count_dateline_crossings(vertices[:, 0])
    if dateline_crossings % 2:
        # Determine index of the (unique) line segment that crosses the dateline.
        for i in range(len(vertices)):
            v0 = vertices[i - 1, :]
            v1 = vertices[i, :]
            if crosses_dateline(v0[0], v1[0]):
                delta_lat = abs(reference_angle(v1[0] - v0[0]))
                lat = (np.pi - abs(v0[0])) / delta_lat * v0[1] + (np.pi - abs(v1[0])) / delta_lat * v1[1]
                out_vertices += [np.vstack((vertices[:i, :], [
                    [np.sign(v0[0]) * np.pi, lat],
                    [np.sign(v0[0]) * np.pi, np.sign(lat) * np.pi / 2],
                    [-np.sign(v0[0]) * np.pi, np.sign(lat) * np.pi / 2],
                    [-np.sign(v0[0]) * np.pi, lat],
                ], vertices[i:, :]))]
                break
    elif dateline_crossings:
        frame_poly = geos.Polygon(np.array([[-np.pi, np.pi/2], [-np.pi, -np.pi/2], [np.pi, -np.pi/2], [np.pi, np.pi/2]]))
        poly = geos.Polygon(np.vstack((vertices[:, 0] % (2 * np.pi), vertices[:, 1])).T)
        if poly.intersects(frame_poly):
            out_vertices += [p.get_coords() for p in poly.intersection(frame_poly)]
        poly = geos.Polygon(np.vstack((vertices[:, 0] % (-2 * np.pi), vertices[:, 1])).T)
        if poly.intersects(frame_poly):
            out_vertices += [p.get_coords() for p in poly.intersection(frame_poly)]
    else:
        out_vertices += [vertices]

    return out_vertices


def cut_prime_meridian(vertices):
    """Cut a polygon across the prime meridian, possibly splitting it into multiple
    polygons.  Vertices consist of (longitude, latitude) pairs where longitude
    is always given in terms of a wrapped angle (between 0 and 2*pi).

    This routine is not meant to cover all possible cases; it will only work for
    convex polygons that extend over less than a hemisphere."""

    out_vertices = []

    # Ensure that the list of vertices does not contain a repeated endpoint.
    if (vertices[0, :] == vertices[-1, :]).all():
        vertices = vertices[:-1, :]

    # Ensure that the longitudes are wrapped from 0 to 2*pi.
    vertices = np.vstack((wrapped_angle(vertices[:, 0]), vertices[:, 1])).T

    def count_meridian_crossings(phis):
        n = 0
        for i in range(len(phis)):
            if crosses_meridian(phis[i - 1], phis[i]):
                n += 1
        return n

    def crosses_meridian(phi0, phi1):
        """Test if the segment consisting of v0 and v1 croses the meridian."""
        # If the two angles are in [0, 2pi), then the shortest arc connecting
        # them crosses the meridian if the difference of the angles is greater
        # than pi.
        phi0, phi1 = sorted((phi0, phi1))
        return phi1 - phi0 > np.pi

    # Count the number of times that the polygon crosses the meridian.
    meridian_crossings = count_meridian_crossings(vertices[:, 0])

    if meridian_crossings % 2:
        # If there are an odd number of meridian crossings, then the polygon
        # encloses the pole. Any meridian-crossing edge has to be extended
        # into a curve following the nearest polar edge of the map.
        for i in range(len(vertices)):
            v0 = vertices[i - 1, :]
            v1 = vertices[i, :]
            # Loop through the edges until we find one that crosses the meridian.
            if crosses_meridian(v0[0], v1[0]):
                # If this segment crosses the meridian, then fill it to
                # the edge of the map by inserting new line segments.

                # Find the latitude at which the meridian crossing occurs by
                # linear interpolation.
                delta_lon = abs(reference_angle(v1[0] - v0[0]))
                lat = abs(reference_angle(v0[0])) / delta_lon * v0[1] + abs(reference_angle(v1[0])) / delta_lon * v1[1]

                # Find the closer of the left or the right map boundary for
                # each vertex in the line segment.
                lon_0 = 0. if v0[0] < np.pi else 2*np.pi
                lon_1 = 0. if v1[0] < np.pi else 2*np.pi

                # Set the output vertices to the polar cap plus the original
                # vertices.
                out_vertices += [np.vstack((vertices[:i, :], [
                    [lon_0, lat],
                    [lon_0, np.sign(lat) * np.pi / 2],
                    [lon_1, np.sign(lat) * np.pi / 2],
                    [lon_1, lat],
                ], vertices[i:, :]))]

                # Since the polygon is assumed to be convex, the only possible
                # odd number of meridian crossings is 1, so we are now done.
                break
    elif meridian_crossings:
        # Since the polygon is assumed to be convex, if there is an even number
        # of meridian crossings, we know that the polygon does not enclose
        # either pole. Then we can use ordinary Euclidean polygon intersection
        # algorithms.

        # Construct polygon representing map boundaries in longitude and latitude.
        frame_poly = geos.Polygon(np.array([[0., np.pi/2], [0., -np.pi/2], [2*np.pi, -np.pi/2], [2*np.pi, np.pi/2]]))

        # Intersect with polygon re-wrapped to lie in [pi, 3*pi).
        poly = geos.Polygon(np.vstack((reference_angle(vertices[:, 0]) + 2 * np.pi, vertices[:, 1])).T)
        if poly.intersects(frame_poly):
            out_vertices += [p.get_coords() for p in poly.intersection(frame_poly)]

        # Intersect with polygon re-wrapped to lie in [-pi, pi).
        poly = geos.Polygon(np.vstack((reference_angle(vertices[:, 0]), vertices[:, 1])).T)
        if poly.intersects(frame_poly):
            out_vertices += [p.get_coords() for p in poly.intersection(frame_poly)]
    else:
        # Otherwise, there were zero meridian crossings, so we can use the
        # original vertices as is.
        out_vertices += [vertices]

    # Done!
    return out_vertices


def make_rect_poly(width, height, theta, phi, subdivisions=10):
    """Create a Polygon patch representing a rectangle with half-angles width
    and height rotated from the north pole to (theta, phi)."""

    # Convert width and height to radians, then to Cartesian coordinates.
    w = np.sin(np.deg2rad(width))
    h = np.sin(np.deg2rad(height))

    # Generate vertices of rectangle.
    v = np.array([[-w, -h], [w, -h], [w, h], [-w, h]])

    # Subdivide.
    v = subdivide_vertices(v, subdivisions)

    # Project onto sphere by calculating z-coord from normalization condition.
    v = np.hstack((v, np.sqrt(1. - np.expand_dims((v * v).sum(1), 1))))

    # Transform vertices.
    v = np.dot(v, hp.rotator.euler_matrix_new(phi, theta, 0, Y=True))

    # Convert to spherical polar coordinates.
    thetas, phis = hp.vec2ang(v)

    # FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
    if mpl_version < '1.2.0':
        # Return list of vertices as longitude, latitude pairs.
        return np.vstack((reference_angle(phis), 0.5 * np.pi - thetas)).T
    else:
        # Return list of vertices as longitude, latitude pairs.
        return np.vstack((wrapped_angle(phis), 0.5 * np.pi - thetas)).T


def heatmap(func, *args, **kwargs):
    "Plot a function on the sphere using the current geographic projection."""

    # Get current axis.
    ax = plt.gca()

    # Set up a regular grid tiling the bounding box of the axes.
    x = np.arange(ax.bbox.x0, ax.bbox.x1 + 0.5, 0.5)
    y = np.arange(ax.bbox.y0, ax.bbox.y1 + 0.5, 0.5)
    xx, yy = np.meshgrid(x, y)

    # Get axis data limits.
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Retrieve the inverse transform of the current axes (which converts display
    # coodinates to data coordinates).
    itrans = ax.transData.inverted()

    # Get the longitude and latitude of every point in the bounding box.
    lons, lats = itrans.transform(np.vstack((xx.flatten(), yy.flatten())).T).T

    # Create a mask that selects only the pixels that fall inside the map boundary.
    mask = np.isfinite(lons) & np.isfinite(lats) & (lons >= xmin) & (lons <= xmax)
    zz = np.ma.array(np.empty(lons.shape), mask=~mask)

    # Evaluate the function everywhere that the mask is set.
    zz[mask] = func(lons[mask], lats[mask])

    # Plot bitmap using imshow.
    aximg = plt.imshow(zz.reshape(xx.shape), aspect=ax.get_aspect(),
        extent=(xmin, xmax, ymax, ymin), *args, **kwargs)

    # Hide masked-out values by displaying them in transparent white.
    aximg.cmap.set_bad('w', alpha=0.)

    # Done.
    return aximg


def contour(func, *args, **kwargs):
    "Plot a function on the sphere using the current geographic projection."""

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


def _healpix_lookup(map, lon, lat):
    """Look up the value of a HEALPix map in the pixel containing the point
    with the specified longitude and latitude."""
    nside = hp.npix2nside(len(map))
    return map[hp.ang2pix(nside, 0.5 * np.pi - lat, lon)]


def healpix_heatmap(map, *args, **kwargs):
    """Produce a heatmap from a HEALPix map."""
    return heatmap(functools.partial(_healpix_lookup, map),
        *args, **kwargs)


def healpix_contour(map, *args, **kwargs):
    """Produce a contour plot from a HEALPix map."""
    return contour(functools.partial(_healpix_lookup, map),
        *args, **kwargs)


def colorbar(vmax):
    # Work out a good tick spacing for colorbar.  Why is this so complicated?
    base = int(np.floor(np.log10(vmax)))
    dtick = 10. ** base
    if vmax / dtick < 2:
        dtick *= 0.25
    elif vmax / dtick < 5:
        dtick *= 0.5
    if vmax % dtick == 0:
        ticks = np.arange(0, vmax + 0.5 * dtick, dtick)
    else:
        ticks = np.arange(0, vmax, dtick)
    ticklabels = ['$%g$' % (tick / 10.**base) for tick in ticks]
    if '.' in ticklabels[-1]: ticklabels[-1] = r'$\;\;\;\;$' + ticklabels[-1]
    else: ticklabels[-1] = r'$\;\;\;\,\,$' + ticklabels[-1]
    ticklabels[-1] += r'$\times 10^{%d}$' % base
    formatter = FixedFormatter(ticklabels)

    # Plot colorbar
    cb = plt.colorbar(orientation='horizontal', ticks=ticks, format=formatter,
        shrink=0.4)

    # Adjust appearance of colorbar tick labels
    for tick, ticklabel in zip(cb.ax.get_xticks(), cb.ax.get_xticklabels()):
        ticklabel.set_verticalalignment('baseline')
        ticklabel.set_y(-1.5)

    # Done.
    return cb


def outline_text(ax):
    """If we are using a new enough version of matplotlib, then
    add a white outline to all text to make it stand out from the background."""
    try:
        # Try to import matplotlib.patheffects (requires matplotlib 1.0+).
        from matplotlib import patheffects
    except ImportError:
        # If import failed, print a warning and do nothing.
        warnings.warn("This version of matplotlib does not support path effects.")
    else:
        # Otherwise, add the path effects.
        effects = [patheffects.withStroke(linewidth=2, foreground='w')]
        for artist in ax.findobj(text.Text):
            artist.set_path_effects(effects)
