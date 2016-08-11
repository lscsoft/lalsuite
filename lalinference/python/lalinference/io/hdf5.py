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
from astropy.table import Column, Table
from lalinference import LALInferenceHDF5PosteriorSamplesDatasetName \
    as POSTERIOR_SAMPLES
from lalinference import LALINFERENCE_PARAM_LINEAR as LINEAR
from lalinference import LALINFERENCE_PARAM_CIRCULAR as CIRCULAR
from lalinference import LALINFERENCE_PARAM_FIXED as FIXED
from lalinference import LALINFERENCE_PARAM_OUTPUT as OUTPUT


_colname_map = (('rightascension', 'ra'),
                ('declination', 'dec'),
                ('distance', 'dist'),
                ('polarisation','psi'),
                ('mchirp', 'mc'),
                ('a_spin1', 'a1'),
                ('a_spin2','a2'),
                ('tilt_spin1', 'tilt1'),
                ('tilt_spin2', 'tilt2'))


def _remap_colnames(table):
    for old_name, new_name in _colname_map:
        if old_name in table.colnames:
            table.rename_column(old_name, new_name)


def _find_table(group, tablename):
    """Recursively search an HDF5 group or file for a dataset by name.

    Parameters
    ----------
    group : `h5py.File` or `h5py.Group`
        The file or group to search
    tablename : str
        The name of the table to search for

    Returns
    -------
    table : `h5py.Dataset`
        The dataset whose name is `tablename`

    Raises
    ------
    KeyError
        If the table is not found or if multiple matching tables are found

    Check that we can find a file by name:
    >>> import os.path
    >>> from lalinference.bayestar.command import TemporaryDirectory
    >>> table = Table(np.eye(3), names=['a', 'b', 'c'])
    >>> with TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     table.write(filename, path='foo/bar', append=True)
    ...     table.write(filename, path='foo/bat', append=True)
    ...     table.write(filename, path='foo/xyzzy/bat', append=True)
    ...     with h5py.File(filename, 'r') as f:
    ...         _find_table(f, 'bar')
    <HDF5 dataset "bar": shape (3,), type "|V24">

    Check that an exception is raised if the table is not found:
    >>> with TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     table.write(filename, path='foo/bar', append=True)
    ...     table.write(filename, path='foo/bat', append=True)
    ...     table.write(filename, path='foo/xyzzy/bat', append=True)
    ...     with h5py.File(filename, 'r') as f:
    ...         _find_table(f, 'plugh')
    Traceback (most recent call last):
        ...
    KeyError: 'Table not found: plugh'

    Check that an exception is raised if multiple tables are found:
    >>> with TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     table.write(filename, path='foo/bar', append=True)
    ...     table.write(filename, path='foo/bat', append=True)
    ...     table.write(filename, path='foo/xyzzy/bat', append=True)
    ...     with h5py.File(filename, 'r') as f:
    ...         _find_table(f, 'bat')
    Traceback (most recent call last):
        ...
    KeyError: 'Multiple tables called bat exist: foo/xyzzy/bat, foo/bat'
    """
    results = {}

    def visitor(key, value):
        _, _, name = key.rpartition('/')
        if name == tablename:
            results[key] = value

    group.visititems(visitor)

    if len(results) == 0:
        raise KeyError('Table not found: {0}'.format(tablename))

    if len(results) > 1:
        raise KeyError('Multiple tables called {0} exist: {1}'.format(
            tablename, ', '.join(results.keys())))

    table, = results.values()
    return table


def read_samples(filename, path=None, tablename=POSTERIOR_SAMPLES):
    """Read an HDF5 sample chain file.

    Parameters
    ----------
    filename : str
        The path of the HDF5 file on the filesystem.
    path : str, optional
        The path of the dataset within the HDF5 file.
    tablename : str, optional
        The name of table to search for recursively within the HDF5 file.
        By default, search for 'posterior_samples'.

    Returns
    -------
    table : `astropy.table.Table`
        The sample chain as an Astropy table.

    >>> import os.path
    >>> from lalinference.bayestar.command import TemporaryDirectory
    >>> table = Table([
    ...     Column(np.ones(10), name='foo', meta={'vary': FIXED}),
    ...     Column(np.arange(10), name='bar', meta={'vary': LINEAR}),
    ...     Column(np.arange(10) * np.pi, name='bat', meta={'vary': CIRCULAR}),
    ...     Column(np.arange(10), name='baz', meta={'vary': OUTPUT})
    ... ])
    >>> with TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     write_samples(table, filename, 'foo/bar/posterior_samples')
    ...     len(read_samples(filename))
    10
    """
    with h5py.File(filename, 'r') as f:
        if path is not None: # Look for a given path
            table = f[path]
        else: # Look for a given table name
            table = _find_table(f, tablename)
        table = Table.read(table)

    # Restore vary types.
    for column, vary in zip(table.columns.values(), table.meta.pop('vary')):
        column.meta['vary'] = vary

    # Restore fixed columns from table attributes.
    for key, value in table.meta.items():
        table.add_column(Column([value] * len(table), name=key,
                         meta={'vary': FIXED}))

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
    metadata: dict (optional)
        Dictionary of (path, value) pairs of metadata attributes
        to add to the output file

    Check that we catch columns that are supposed to be FIXED but are not:
    >>> table = Table([
    ...     Column(np.arange(10), name='foo', meta={'vary': FIXED})
    ... ])
    >>> write_samples(table, 'bar.hdf5', 'bat/baz') # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    AssertionError: 
    Arrays are not equal
    Column foo is a fixed column, but its values are not identical
    ...

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
        if column.meta['vary'] == FIXED:
            np.testing.assert_array_equal(column[1:], column[0],
                                          'Column {0} is a fixed column, but '
                                          'its values are not identical'
                                          .format(column.name))
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
                    except KeyError:
                        raise KeyError(
                            'Unable to set metadata {0}[{1}] = {2}'.format(
                            internal_path, key, value))
