# Copyright (C) 2013 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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

"""This module provides a method to return trigger files in a similar
manner to gw_data_find
"""

import glob
import os
import numpy
import re
import sys

from glue import segments
from glue.lal import (Cache, CacheEntry)

from laldetchar import git_version
__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

TRIGFIND_BASE_PATH = "/home/detchar/triggers"

def find_trigger_urls(channel, etg, gpsstart, gpsend, verbose=False, **kwargs):
    """Find the paths of trigger files that represent the given
    observatory, channel, and ETG (event trigger generator) for a given
    GPS [start, end) segment.
    """
    # special case for KW
    if etg.lower() in ['kw', 'kleinewelle']:
        from .kw import find_dmt_cache
        ifo = channel.split(':')[0]
        kwargs.setdefault('extension', 'xml')
        kwargs.setdefault('check_files', True)
        return find_dmt_cache(gpsstart, gpsend, ifo, **kwargs)
    elif etg.lower() == 'omega':
        from .omega import find_dmt_cache
        ifo = channel.split(':')[0]
        kwargs.setdefault('check_files', True)
        return find_dmt_cache(gpsstart, gpsend, ifo, **kwargs)
    elif etg.lower() == 'omicron':
        etg = '?micron'

    # construct search
    span = segments.segment(gpsstart, gpsend)
    ifo, channel = channel.split(':', 1)
    trigtype = "%s_%s" % (channel, etg.lower())
    epoch = '*'
    searchbase = os.path.join(TRIGFIND_BASE_PATH, epoch, ifo, trigtype)
    gpsdirs = numpy.arange(int(str(gpsstart)[:5]), int(str(gpsend)[:5])+1)
    trigform = ('%s-%s_%s-%s-*.xml*'
                % (ifo, re.sub('-', '_', channel), etg.lower(), '[0-9]'*10))

    # perform and cache results
    out = Cache()
    for gpsdir in gpsdirs:
        gpssearchpath = os.path.join(searchbase, str(gpsdir), trigform)
        if verbose:
            sys.stdout.write("Searching %s..."
                             % os.path.split(gpssearchpath)[0])
            sys.stdout.flush()
        gpscache = Cache(map(CacheEntry.from_T050017,
                             glob.glob(os.path.join(searchbase, str(gpsdir),
                                                    trigform))))
        out.extend(gpscache.sieve(segment=span))
        if verbose:
            sys.stdout.write(" %d found\n" % len(gpscache.sieve(segment=span)))
    out.sort(key=lambda e: e.path)

    return out
