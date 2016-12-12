# Copyright (C) 2016 Ian Harry
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# the knowledge that it will not be of any use whatsoever. It is distributed
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
The code reads in a set of compressed template banks and combines them into a
single template bank.
"""

import argparse
import numpy
import h5py
import logging
__author__  = "Ian Harry <ian.harry@ligo.org>"
__program__ = "lalapps_cbc_sbank_hdf5_bankcombiner"

parser = argparse.ArgumentParser(description=__doc__[1:])
parser.add_argument("--output-file", type=str,
                    help="Output hdf bank file.")
parser.add_argument("--input-filenames", nargs='*', default=None,
                    action="store",
                    help="List of input hdf bank files.")
parser.add_argument("--verbose", action="store_true", default=False)

args = parser.parse_args()

attrs_dict = None
items_dict = None
approx_map_dict = {}
approx_map_dict['counter'] = 1

for file_name in args.input_filenames:
    hdf_fp = h5py.File(file_name, 'r')
    if 'empty_file' in hdf_fp.attrs:
        continue
    if attrs_dict is None:
        attrs_dict = {}
        for key, item in hdf_fp.attrs.items():
            attrs_dict[key] = item
    if items_dict is None:
        items_dict = {}
        for item, entries in hdf_fp.items():
            items_dict[item] = entries[:]
    else:
        curr_items = set(items_dict.keys())
        new_items = set(hdf_fp.keys())
        # This does the XOR check of the two sets of keys.
        # Basically we demand that the two files must have the same items.
        if set(curr_items).symmetric_difference(new_items):
            err_msg = "All input files must contain the same data structures. "
            err_msg += "File {} ".format(file_name)
            err_msg += "contains fields {} ".format(new_items)
            err_msg += "other files contain {}.".format(curr_items)
            raise ValueError(err_msg)
        for item, entries in hdf_fp.items():
            items_dict[item] = numpy.append(items_dict[item], entries[:])
    hdf_fp.close()

out_fp = h5py.File(args.output_file, 'w')
if attrs_dict is None:
    out_fp.attrs['empty_file'] = True
else:
    for item, value in items_dict.items():
        out_fp[item] = value
    for item, value in attrs_dict.items():
        out_fp.attrs[item] = value
out_fp.close()
