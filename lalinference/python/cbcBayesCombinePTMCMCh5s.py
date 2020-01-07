# -*- coding: utf-8 -*-
#
#       cbcBayesCombinePTMCMCh5s.py
#
#       Copyright 2016
#       Ben Farr <benjamin.farr@ligo.org>
#
#       Combine multiple HDF5s from the same parallel-tempered run into a common hdf5
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

from optparse import OptionParser

import h5py

from lalinference import git_version
__author__="Ben Farr <benjamin.farr@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

USAGE='''%prog [options] PTMCMC_datafile.hdf5 [PTMCMC_datafile2.hdf5 ...]
Combine chains from a parallel-tempered MCMC run spread across several HDF5 files.
'''

if __name__ == '__main__':
    parser = OptionParser(USAGE)
    parser.add_option(
        '-o', '--outfile', type='string', default=None,
        help='Output file for posterior samples. If None, file containing T=1 chain will be used', metavar='combined_chains.hdf5')
    opts, args = parser.parse_args()

    datafiles = args

    group_id = '/lalinference/lalinference_mcmc'

    outfile = opts.outfile

    # find the file containing the coldest chain (i.e. a posterior_samples group)
    for datafile in datafiles:
        rootfile = datafile
        possible_root = h5py.File(datafile, 'a')
        if group_id + '/posterior_samples' in possible_root:
            possible_root.close()
            break
        else:
            possible_root.close()

    if outfile:
        if outfile != rootfile:
            try:
                outp = h5py.File(outfile, 'w-')

                # Copy contents of the root hdf5 to preserve metadata
                with h5py.File(rootfile, 'r') as inp:
                    for group in inp.keys():
                       inp.copy(group, outp['/'])

            except IOError:
                assert outfile in datafiles, \
                    "Trying to write to an existing file that isn't being explicitly combined"

                outp = h5py.File(outfile, 'a')

                # Copy root over now for consistent behavior
                with h5py.File(rootfile, 'r') as inp:
                    for chain in inp[group_id].keys():
                        chain_id = group_id + '/' + chain
                        inp.copy(chain_id, outp[group_id])

            # make sure the target group exists
            try:
                outp.create_group(group_id)
            except ValueError:
                pass
    else:
        outfile = rootfile

    for datafile in datafiles:
        if datafile == outfile:
            continue
        with h5py.File(datafile, 'r') as inp:
            for chain in inp[group_id].keys():
                chain_id = group_id + '/' + chain
                if chain_id not in outp:
                   inp.copy(chain_id, outp[group_id])
    outp.close()
