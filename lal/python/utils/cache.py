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

"""Modules extending the Cache file functionality from LAL
"""

import tempfile
import os

from .. import git_version
from ..lal import CacheImport

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

__all__ = ['lalcache_from_gluecache']

def lalcache_from_gluecache(cache):
    """Convert a glue.lal.Cache object to a lal.Cache object.
    Writes cache to temporary file and reads to Cache.

    @param cache
        LAL cache object from GLUE to convert
        type cache glue.lal.Cache

    @returns a lal.Cache object representing the same data
    """
    with tempfile.NamedTemporaryFile(delete=False) as t:
        cache = cache
        for e in cache:
            e.segment = type(e.segment)(int(e.segment[0]), int(e.segment[1]))
        cache.tofile(t)
        frcache = CacheImport(t.name)
    os.remove(t.name)
    return frcache
