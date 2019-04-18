#
# Copyright (C) 2013-2017  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
Convenience function to produce a sky map from LIGO-LW rows.
"""
import numpy as np
from astropy.table import Column, Table
from astropy import units as u
from .. import moc
from .. import healpix_tree

from lalinference.bayestar.deprecation import warn
warn('ligo.skymap.bayestar')


def rasterize(skymap):
    skymap = Table(moc.rasterize(skymap), meta=skymap.meta)
    skymap.rename_column('PROBDENSITY', 'PROB')
    skymap['PROB'] *= 4 * np.pi / len(skymap)
    skymap['PROB'].unit = u.pixel ** -1
    return skymap


def derasterize(skymap):
    skymap.rename_column('PROB', 'PROBDENSITY')
    skymap['PROBDENSITY'] *= len(skymap) / (4 * np.pi)
    skymap['PROBDENSITY'].unit = u.steradian ** -1
    nside, _, ipix, _, _, value = zip(
        *healpix_tree.reconstruct_nested(skymap))
    nside = np.asarray(nside)
    ipix = np.asarray(ipix)
    # FIXME: replace with np.stack() when Numpy 1.10.0 is on all
    # of the LIGO Data Grid clusters
    value = np.hstack(value)
    uniq = (4 * np.square(nside) + ipix).astype(np.uint64)
    old_units = [column.unit for column in skymap.columns.values()]
    skymap = Table(value, meta=skymap.meta)
    for old_unit, column in zip(old_units, skymap.columns.values()):
        column.unit = old_unit
    skymap.add_column(Column(uniq, name='UNIQ'), 0)
    return skymap
