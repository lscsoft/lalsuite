# Copyright (C) 2011 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation,.

## \defgroup laldetchar_py_triggers Triggers
## \ingroup laldetchar_python
"""A collection of tools for reading, manipulating, and writing triggers
in a standard format.

Functions are provided to read events from file, with internal methods
for reading from `ASCII`, `ROOT`, and `LIGO_LW` XML files.

For each of those event trigger generators (ETG) that do not (or did not)
write natively to `LIGO_LW` XML, an I/O module is provided to import
events from their native format, either `ASCII` or `ROOT`, into standard
`LIGO_LW` tables, imported from the `glue.ligolw` package [1].
"""
# [1]: <a href="https://www.lsc-group.phys.uwm.edu/daswg/projects/glue/doc/glue.ligolw-module.html" ref="external">glue.ligolw library</a>
# \author Duncan Macleod (<duncan.macleod@ligo.org>)
#
# ### Synopsis ###
#
# ~~~
# from laldetchar import triggers
# ~~~
#
# ### Example ###
#
# \code
# from laldetchar import triggers
# trigs = triggers.load_triggers("mytrigfile.root", "omicron")
# \endcode

import sys
import re as _re
import warnings as _warnings
_warnings.filterwarnings("ignore", "column name", UserWarning)

from . import utils

from glue import (segments, lal as _cache)

from laldetchar import git_version as version
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = version.id
__date__ = version.date

_re_xml = _re.compile("(xml|xml.gz)\Z")
_re_root = _re.compile("root\Z")

## \addtogroup laldetchar_py_triggers
#@{

def from_file(filename, etg, columns=None, start=None, end=None, **kwargs):
    """Reads the triggers for the stated trigger generator from the
    given file.

    @param filename
        path to file containing triggers
    @param etg
        the name of the parent trigger generator
    @param columns
        a list of valid LIGO_LW column names for the new table
        (defaults to all)
    @param start
        minimum GPS time for returned triggers
    @param end
        maximum GPS time for returned triggers
    @param kwargs
        other keyword arguments passed to the relevant trigger reader
        for this ETG

    @returns a LIGO_LW table containing the triggers
    """
    if _re_xml.search(filename):
        return utils.from_ligolw(filename, etg, columns=columns,
                                 start=start, end=end, **kwargs)
    elif _re_root.search(filename):
        return utils.from_root(filename, etg, columns=columns,
                               start=start, end=end, **kwargs)
    else:
        return utils.from_ascii(filename, etg, columns=columns,
                                start=start, end=end, **kwargs)


def from_files(filelist, etg, columns=None, start=None, end=None,
               verbose=False, **kwargs):
    """Read the triggers for the stated trigger generator from each of
    the files in filelist into a single LIGO_LW table.

    @param filelist
        list of file names, or `glue.lal.Cache` containing file entries
    @param etg
        the name of the parent trigger generator
    @param columns
        a list of valid LIGO_LW column names for the new table
        (defaults to all)
    @param start
        minimum GPS time for returned triggers
    @param end
        maximum GPS time for returned triggers
    @param verbose UNDOCUMENTED
    @param kwargs
        other keyword arguments passed to the relevant trigger reader
        for this ETG

    @returns a LIGO_LW table containing the triggers
    """
    if verbose:
        N = float(len(filelist))
        def print_verbose(i, final=False):
            v = "Reading triggers... %d/%d (%.1f%%)" % (i, N, i/N * 100)
            if final:
                sys.stdout.write("%s\n" % v)
            else:
                sys.stdout.write("%s\r" % v)
            sys.stdout.flush()
    if isinstance(filelist, _cache.Cache):
        span = segments.segment(start is not None and start or
                                segments.NegInfinity,
                                end is not None and end or segments.PosInfinity)
        filelist = filelist.sieve(segment=span).pfnlist()
    if len(filelist) == 0:
        if (kwargs.has_key('channel') and
            (not isinstance(kwargs['channel'], basestring) and
             kwargs['channel'] is not None)):
            channels = kwargs['channel']
            return dict((c, utils.new_ligolw_table(etg, columns=columns)) for
                        c in channels)
        else:
            return utils.new_ligolw_table(etg, columns=columns)
    if verbose:
        print_verbose(0)
    out = from_file(filelist[0], etg, columns=columns, start=start, end=end,
                    **kwargs)
    if isinstance(out, dict):
        def extend(in_):
            for key,tab in in_.iteritems():
                out[key].extend(tab)
    else:
        extend = out.extend
    if verbose:
        print_verbose(1)
    for i,fp in enumerate(filelist[1:]):
        extend(from_file(fp, etg, columns=columns, start=start, end=end,
                         **kwargs))
        if verbose:
            print_verbose(i+2)
    if verbose:
        print_verbose(N, final=True)
    return out

# close doxygen
##
#	\defgroup laldetchar_py_triggers_cwb		Coherent WaveBurst
#	\defgroup laldetchar_py_triggers_excesspower	ExcessPower
#	\defgroup laldetchar_py_triggers_kleinewelle	KleineWelle
#	\defgroup laldetchar_py_triggers_omega		Omega
#	\defgroup laldetchar_py_triggers_omicron	Omicron
#	\defgroup laldetchar_py_triggers_utils		Utilities
#@}
