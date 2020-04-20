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

import numpy as np
import h5py
from astropy.table import Column, Table
from lalinference import LALInferenceHDF5PosteriorSamplesDatasetName \
    as POSTERIOR_SAMPLES
from lalinference import LALINFERENCE_PARAM_FIXED as FIXED
from lalinference import LALINFERENCE_PARAM_OUTPUT as OUTPUT

__all__ = ('read_samples', 'write_samples', 'extract_metadata')


def _identity(x):
    return x


_colname_map = (('rightascension', 'ra', _identity),
                ('declination', 'dec', _identity),
                ('logdistance', 'dist', np.exp),
                ('distance', 'dist', _identity),
                ('polarisation', 'psi', _identity),
                ('chirpmass', 'mc', _identity),
                ('a_spin1', 'a1', _identity),
                ('a_spin2', 'a2', _identity),
                ('tilt_spin1', 'tilt1', _identity),
                ('tilt_spin2', 'tilt2', _identity))


def _remap_colnames(table):
    for old_name, new_name, func in _colname_map:
        if old_name in table.colnames:
            table[new_name] = func(table.columns.pop(old_name))


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
    >>> from tempfile import TemporaryDirectory
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
    KeyError: 'Multiple tables called bat exist: foo/bat, foo/xyzzy/bat'
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
            tablename, ', '.join(sorted(results.keys()))))

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

    Test reading a file written using the Python API:
    >>> import os.path
    >>> from tempfile import TemporaryDirectory
    >>> table = Table([
    ...     Column(np.ones(10), name='foo', meta={'vary': FIXED}),
    ...     Column(np.arange(10), name='bar', meta={'vary': LINEAR}),
    ...     Column(np.arange(10) * np.pi, name='bat', meta={'vary': CIRCULAR}),
    ...     Column(np.arange(10), name='baz', meta={'vary': OUTPUT})
    ... ])
    >>> with TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     write_samples(table, filename, path='foo/bar/posterior_samples')
    ...     len(read_samples(filename))
    10

    Test reading a file that was written using the LAL HDF5 C API:
    >>> table = read_samples('test.hdf5')
    >>> table.colnames
    ['uvw', 'opq', 'lmn', 'ijk', 'def', 'abc', 'rst', 'ghi']
    """
    with h5py.File(filename, 'r') as f:
        if path is not None:  # Look for a given path
            table = f[path]
        else:  # Look for a given table name
            table = _find_table(f, tablename)
        table = Table.read(table)

    # Restore vary types.
    for i, column in enumerate(table.columns.values()):
        column.meta['vary'] = table.meta['FIELD_{0}_VARY'.format(i)]

    # Restore fixed columns from table attributes.
    for key, value in table.meta.items():
        # Skip attributes from H5TB interface
        # (https://www.hdfgroup.org/HDF5/doc/HL/H5TB_Spec.html).
        if key == 'CLASS' or key == 'VERSION' or key == 'TITLE' or key.startswith('FIELD_'):
            continue
        if key in table.colnames:
            # This is handled separately as rename_duplicate can trigger a bug in astropy < 2.0.16
            table.add_column(Column([value] * len(table), name=key,
                         meta={'vary': FIXED}), rename_duplicate=True)
        else:
            table.add_column(Column([value] * len(table), name=key,
                         meta={'vary': FIXED}))

    # Delete remaining table attributes.
    table.meta.clear()

    # Normalize column names.
    _remap_colnames(table)

    # Done!
    return table


def write_samples(table, filename, metadata=None, **kwargs):
    """Write an HDF5 sample chain file.

    Parameters
    ----------
    table : `astropy.table.Table`
        The sample chain as an Astropy table.
    filename : str
        The path of the HDF5 file on the filesystem.
    metadata: dict (optional)
        Dictionary of (path, value) pairs of metadata attributes
        to add to the output file
    kwargs: dict
        Any keyword arguments for `astropy.table.Table.write`.

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

    And now try writing an arbitrary example to a temporary file.
    >>> import os.path
    >>> from tempfile import TemporaryDirectory
    >>> table = Table([
    ...     Column(np.ones(10), name='foo', meta={'vary': FIXED}),
    ...     Column(np.arange(10), name='bar', meta={'vary': LINEAR}),
    ...     Column(np.arange(10) * np.pi, name='bat', meta={'vary': CIRCULAR}),
    ...     Column(np.arange(10), name='baz', meta={'vary': OUTPUT})
    ... ])
    >>> with TemporaryDirectory() as dir:
    ...     write_samples(
    ...         table, os.path.join(dir, 'test.hdf5'), path='bat/baz')
    """
    # Copy the table so that we do not modify the original.
    table = table.copy()

    # Make sure that all tables have a 'vary' type.
    for column in table.columns.values():
        if 'vary' not in column.meta:
            if np.all(column[0] == column[1:]):
                column.meta['vary'] = FIXED
            else:
                column.meta['vary'] = OUTPUT
    # Reconstruct table attributes.
    for colname, column in tuple(table.columns.items()):
        if column.meta['vary'] == FIXED:
            np.testing.assert_array_equal(column[1:], column[0],
                                          'Column {0} is a fixed column, but '
                                          'its values are not identical'
                                          .format(column.name))
            table.meta[colname] = column[0]
            del table[colname]
    for i, column in enumerate(table.columns.values()):
        table.meta['FIELD_{0}_VARY'.format(i)] = column.meta['vary']
    table.write(filename, format='hdf5', **kwargs)
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

def update_metadata(metadata, level, attrs, strict_versions, collision='raise'):
        """Updates the sub-dictionary 'key' of 'metadata' with the values from
        'attrs', while enforcing that existing values are equal to those with
        which the dict is updated.
        """
        if level not in metadata:
            metadata[level] = {}
        for key in attrs:
            if key in metadata[level]:
                    if collision == 'raise':
                        if attrs[key]!=metadata[level][key]:
                            if key == 'version' and not strict_versions:
                                continue
                            else:
                                raise ValueError(
                                    'Metadata mismatch on level %r for key %r:\n\t%r != %r'
                                    % (level, key, attrs[key], metadata[level][key]))
                    elif collision == 'append':
                        if isinstance(metadata[level][key], list):
                            metadata[level][key].append(attrs[key])
                        else:
                            metadata[level][key] = [metadata[level][key], attrs[key]]
                    elif collision == 'ignore':
                        pass
                    else:
                        raise ValueError('Invalid value for collision: %r' % collision)
            else:
                metadata[level][key] = attrs[key]
        return

def extract_metadata(filename, metadata, log_noise_evidences=[], log_max_likelihoods=[], nlive=[], dset_name=None, nest=False, strict_versions=True):
        """
        Extract metadata from HDF5 sample chain file

        Parameters
        ----------
        filename : str
            The path of the HDF5 file on the filesystem.
        metadata : dict
            Dict into which to place metadata
        log_noise_evidences : array (optional)
            Array into which to place log noise evidences (if nest = True)
        log_max_likelihoods : array (optional)
            Array into which to place log max likelihoods (if nest = True)
        nlive : array (optional)
            Array into which to place number of live points (if nest = True)
        return_run_identifier : Boolean (optional : default False)
            Whether to return the run identifier
        nest : Boolean (optional : default False)
            Whether to output quantities that only exist for nest runs

        Returns
        -------
        run_identifier : str
            The run identifier
        """
        with h5py.File(filename, 'r') as hdf:
            # walk down the groups until the actual data is reached, storing
            # metadata for each step.
            current_level = '/lalinference'
            group = hdf[current_level]
            update_metadata(metadata, current_level, group.attrs, strict_versions)

            if len(hdf[current_level].keys()) != 1:
                raise KeyError('Multiple run-identifiers found: %r'
                               % list(hdf[current_level].keys()))
            # we ensured above that there is only one identifier in the group.
            run_identifier = list(hdf[current_level].keys())[0]

            current_level = '/lalinference/' + run_identifier
            group = hdf[current_level]
            update_metadata(metadata, current_level, group.attrs, strict_versions, collision='append')

            if nest:
                # store the noise evidence and max likelihood seperately for later use
                log_noise_evidences.append(group.attrs['log_noise_evidence'])
                log_max_likelihoods.append(group.attrs['log_max_likelihood'])
                nlive.append(group.attrs['number_live_points'])

            # storing the metadata under the posterior_group name simplifies
            # writing it into the output hdf file.
            if dset_name == None:
                dset_name = POSTERIOR_SAMPLES
            current_level = '/lalinference/' + run_identifier + '/' + dset_name
            current_level_posterior = '/lalinference/' + run_identifier + '/' + POSTERIOR_SAMPLES
            group = hdf[current_level]
            update_metadata(metadata, current_level_posterior, group.attrs, strict_versions, collision='ignore')

            return run_identifier
