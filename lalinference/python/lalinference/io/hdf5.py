# Copyright (C) 2016  Leo Singer, John Veitch
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
__all__ = ('read_samples', 'write_samples')


import numpy as np
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


def read_samples(filename, path=None, tablename=None):
    """Read an HDF5 sample chain file.

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
	# Look for a given path
        if path is not None:
            table = f[path]
        # Look for a given table name
        elif tablename is not None:
	    table = f.visititems(lambda name,val: val if name.rsplit('/')[-1]==tablename else None)
	# Look for some common paths
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


def write_samples(table, filename, path, metadata=None):
    """Write an HDF5 sample chain file.

    Parameters
    ----------
    table : `astropy.table.Table`
        The sample chain as an Astropy table.
    filename : str
        The path of the HDF5 file on the filesystem.
    path : str
        The path of the dataset within the HDF5 file.
    metadata: dict (Optional)
        Dictionary of (path, value) pairs of metadata attributes
        to add to the output file

    Example... first some imports:
    >>> from lalinference import LALINFERENCE_PARAM_LINEAR as LINEAR
    >>> from lalinference import LALINFERENCE_PARAM_CIRCULAR as CIRCULAR
    >>> from lalinference import LALINFERENCE_PARAM_FIXED as FIXED
    >>> from lalinference import LALINFERENCE_PARAM_OUTPUT as OUTPUT

    Check that we catch columns that are supposed to be FIXED but are not:
    >>> table = Table([
    ...     Column(np.arange(10), name='foo', meta={'vary': FIXED})
    ... ])
    >>> write_samples(table, 'bar.hdf5', 'bat/baz')
    Traceback (most recent call last):
        ...
    AssertionError: 
    Arrays are not equal
    Column {0} is a `fixed` column, but its values are not identical
    (mismatch 100.0%)
     x: Column([1, 2, 3, 4, 5, 6, 7, 8, 9])
     y: array(0)

    And now try writing an arbitrary example to a temporary file:
    >>> import os.path
    >>> from lalinference.bayestar.command import TemporaryDirectory
    >>> table = Table([
    ...     Column(np.ones(10), name='foo', meta={'vary': FIXED}),
    ...     Column(np.arange(10), name='bar', meta={'vary': LINEAR}),
    ...     Column(np.arange(10) * np.pi, name='bat', meta={'vary': CIRCULAR}),
    ...     Column(np.arange(10), name='baz', meta={'vary': OUTPUT})
    ... ])
    >>> with TemporaryDirectory() as dir:
    ...     write_samples(table, os.path.join(dir, 'test.hdf5'), 'bat/baz')
    """
    # Copy the table so that we do not modify the original.
    table = table.copy()

    # Reconstruct table attributes.
    vary = []
    for colname, column in table.columns.items():
        if column.meta['vary'] == lalinference.LALINFERENCE_PARAM_FIXED:
            np.testing.assert_array_equal(column[1:], column[0],
                                          'Column {0} is a `fixed` column, '
                                          'but its values are not identical')
            table.meta[colname] = column[0]
            del table[colname]
        else:
            vary.insert(0, column.meta.pop('vary'))
    table.meta['vary'] = np.asarray(vary)
    table.write(filename, format='hdf5', path=path)
    if metadata:
      with h5py.File(filename) as hdf:
        for internal_path, attributes in metadata.items():
	  for key, value in attributes.items():
	    try:
	      hdf[internal_path].attrs[key] = value
            except:
	      print('Unable to set metadata %s[%s] = %s'%(internal_path,key,str(value)))
