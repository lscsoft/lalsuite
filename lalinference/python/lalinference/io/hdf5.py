# Copyright (C) 2016  Leo Singer
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
Reading HDF5 posterior sample chain HDF5 files.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"
__all__ = ("read_samples",)


import h5py
import lalinference
from astropy.table import Column, Table


_colname_map = (('rightascension', 'ra'),
                ('declination', 'dec'),
                ('distance', 'dist'))
_mcmc_path = 'lalinference/lalinference_mcmc/posterior_samples'
_nest_path = 'lalinference/lalinference_nest/nested_samples'


def _remap_colnames(table):
    for old_name, new_name in _colname_map:
        if old_name in table.colnames:
            table.rename_column(old_name, new_name)


def read_samples(filename, path=None):
    """Read an HDF5 sample chain file.

    DEPRECATED! This function will eventually be deleted once text formatted
    sample chains are retired.

    Parameters
    ----------
    filename : str
        The path of the HDF5 file on the filesystem.
    path : str, optional
        The path of the dataset within the HDF5 file. By default, try the
        conventional path for lalinference_mcmc and lalinference_nest.

    Returns
    -------
    table : `astropy.table.Table`
        The sample chain as an Astropy table.
    """
    with h5py.File(filename, 'r') as f:
        if path is not None:
            table = f[path]
        else:
            try:
                try:
                    table = f[_mcmc_path]
                except KeyError:
                    table = f[_nest_path]
            except KeyError:
                raise KeyError('The HDF5 file did not contain a dataset called '
                               '`{0}` or `{1}`. Try providing the dataset path '
                               'explicitly.'.format(_mcmc_path, _nest_path))
        table = Table.read(table)

    # Restore vary types.
    for column, vary in zip(table.columns.values(), table.meta.pop('vary')):
        column.meta['vary'] = vary

    # Restore fixed columns from table attributes.
    for key, value in table.meta.items():
        table.add_column(Column([value] * len(table), name=key,
                         meta={'vary': lalinference.LALINFERENCE_PARAM_FIXED}))

    # Delete table attributes.
    for key in table.meta:
        del table.meta[key]

    # Normalize column names.
    _remap_colnames(table)

    # Done!
    return table
