#
# Copyright (C) 2012-2018  Leo Singer
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
Axes subclasses for all-sky maps
"""
from __future__ import division

import warnings

from matplotlib.axes import Axes
from matplotlib.ticker import Formatter, FixedLocator
from matplotlib.projections import projection_registry
from matplotlib.transforms import Transform, Affine2D
from matplotlib.projections.geo import LambertAxes, MollweideAxes
import numpy as np

# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
import matplotlib
import distutils.version
mpl_version = distutils.version.LooseVersion(matplotlib.__version__)

__all__ = ('AstroDegreesMollweideAxes', 'AstroHoursMollweideAxes',
           'AstroMollweideAxes', 'AstroLambertAxes')


# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
if mpl_version >= '1.3.0':
    FixedMollweideAxes = MollweideAxes
elif mpl_version < '1.2.0':
    raise NotImplemented(
        'This module requires matplotlib >= 1.2.0. '
        'You have matplotlib {}.'.format(mpl_version))
else:
    class FixedMollweideAxes(MollweideAxes):
        """Patched version of matplotlib's Mollweide projection that implements a
        correct inverse transform."""

        class FixedMollweideTransform(MollweideAxes.MollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.InvertedFixedMollweideTransform(
                    self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

        class InvertedFixedMollweideTransform(
                MollweideAxes.InvertedMollweideTransform):

            def inverted(self):
                return FixedMollweideAxes.FixedMollweideTransform(
                    self._resolution)
            inverted.__doc__ = Transform.inverted.__doc__

            def transform_non_affine(self, xy):
                x = xy[:, 0:1]
                y = xy[:, 1:2]

                sqrt2 = np.sqrt(2)
                sintheta = y / sqrt2
                with np.errstate(invalid='ignore'):
                    costheta = np.sqrt(1. - 0.5 * y * y)
                longitude = 0.25 * sqrt2 * np.pi * x / costheta
                latitude = np.arcsin(
                    2 / np.pi * (np.arcsin(sintheta) + sintheta * costheta))
                return np.concatenate((longitude, latitude), 1)
            transform_non_affine.__doc__ = \
                Transform.transform_non_affine.__doc__

        def _get_core_transform(self, resolution):
            return self.FixedMollweideTransform(resolution)

    projection_registry.register(FixedMollweideAxes)


class AstroDegreesMollweideAxes(FixedMollweideAxes):
    """Mollweide axes with phi axis flipped and in degrees from 360 to 0
    instead of in degrees from -180 to 180."""

    name = 'astro degrees mollweide'

    def cla(self):
        super(AstroDegreesMollweideAxes, self).cla()
        self.set_xlim(0, 2*np.pi)

    def set_xlim(self, *args, **kwargs):
        Axes.set_xlim(self, 0., 2*np.pi)
        Axes.set_ylim(self, -np.pi / 2.0, np.pi / 2.0)

    def _get_core_transform(self, resolution):
        return Affine2D().translate(-np.pi, 0.) + \
            super(AstroDegreesMollweideAxes, self)._get_core_transform(
                resolution)

    def set_longitude_grid(self, degrees):
        # Copied from matplotlib.geo.GeoAxes.set_longitude_grid and modified
        super(AstroDegreesMollweideAxes, self).set_longitude_grid(degrees)
        number = (360.0 / degrees) + 1
        self.xaxis.set_major_locator(
            FixedLocator(
                np.linspace(0, 2*np.pi, number, True)[1:-1]))

    def _set_lim_and_transforms(self):
        # Copied from matplotlib.geo.GeoAxes._set_lim_and_transforms
        # and modified
        super(AstroDegreesMollweideAxes, self)._set_lim_and_transforms()

        # This is the transform for latitude ticks.
        yaxis_stretch = Affine2D().scale(np.pi * 2.0, 1.0)
        yaxis_space = Affine2D().scale(-1.0, 1.1)
        self._yaxis_transform = \
            yaxis_stretch + \
            self.transData
        yaxis_text_base = \
            yaxis_stretch + \
            self.transProjection + \
            (yaxis_space +
             self.transAffine +
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

projection_registry.register(AstroDegreesMollweideAxes)


class AstroHoursMollweideAxes(AstroDegreesMollweideAxes):
    """Mollweide axes with phi axis flipped and in hours from 24 to 0 instead of
    in degrees from -180 to 180."""

    name = 'astro hours mollweide'

    class RaFormatter(Formatter):
        # Copied from matplotlib.geo.GeoAxes.ThetaFormatter and modified
        def __init__(self, round_to=1.0):
            self._round_to = round_to

        def __call__(self, x, pos=None):
            hours = (x / np.pi) * 12.
            hours = round(15 * hours / self._round_to) * self._round_to / 15
            return r"%0.0f$^\mathrm{h}$" % hours

    def set_longitude_grid(self, degrees):
        super(AstroHoursMollweideAxes, self).set_longitude_grid(degrees)
        self.xaxis.set_major_formatter(self.RaFormatter(degrees))

projection_registry.register(AstroHoursMollweideAxes)


# For backward compatibility
class AstroMollweideAxes(AstroHoursMollweideAxes):

    name = 'astro mollweide'

    def __init__(self, *args, **kwargs):
        warnings.warn("The AstroMollweideAxes ('astro mollweide') class has "
                      "been deprecated. Please use AstroHoursMollweideAxes "
                      "('astro hours mollweide') instead.", stacklevel=2)
        super(AstroMollweideAxes, self).__init__(*args, **kwargs)

projection_registry.register(AstroMollweideAxes)


class AstroLambertAxes(LambertAxes):
    name = 'astro lambert'

    def cla(self):
        super(AstroLambertAxes, self).cla()
        self.set_xlim(0, 2*np.pi)

    def set_xlim(self, *args, **kwargs):
        Axes.set_xlim(self, 0., 2*np.pi)
        Axes.set_ylim(self, -np.pi / 2.0, np.pi / 2.0)

    def _get_core_transform(self, resolution):
        return Affine2D().translate(-np.pi, 0.).scale(-1, 1) + super(
            AstroLambertAxes, self)._get_core_transform(resolution)

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
        self.xaxis.set_major_formatter(self.RaFormatter(degrees))

projection_registry.register(AstroLambertAxes)

from lalinference.bayestar.deprecation import warn
warn('ligo.skymap.plot.allsky')
