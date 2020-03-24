#!/usr/bin/env python
#
# Copyright (C) 2019 Duncan Macleod
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

"""Tests for lal.utils.cache

See utils_cache_verify.py for more tests of the same module
"""

import os
import sys

import pytest

from glue.lal import Cache as GlueCache

from lal import (
    LIGOTimeGPS,
    utils as lal_utils,
)

GlueCache.entry_class = lal_utils.CacheEntry


def test_lalcache_from_gluecache():
    files = [
        "X-TEST-0-1.gwf",
        "X-TEST-1-1.gwf",
    ]
    gcache = GlueCache.from_urls(files, coltype=LIGOTimeGPS)
    try:
        lcache = lal_utils.lalcache_from_gluecache(gcache)
    finally:
        for fp in files:
            if os.path.isfile(fp):
                os.remove(fp)
    assert lcache.length == len(gcache)
    assert lcache.list.url == (
        "file://localhost{}".format(os.path.abspath(files[0]))
    )


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-utils-cache.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
