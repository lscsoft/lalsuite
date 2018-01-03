# Copyright (C) 2012 Duncan M. Macleod
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

## \addtogroup laldetchar_py_triggers_excesspower
"""Utitilies to find and manipulate triggers from (GSTLAL) ExcessPower
"""
#
# ### Synopsis ###
#
# ~~~
# from laldetchar.triggers import excesspower
# ~~~
# \author Duncan Macleod (<duncan.macleod@ligo.org>)

import os

from datetime import datetime
from glob import glob

import lal

from glue.lal import (Cache, CacheEntry)
from glue.segments import segment as Segment

from laldetchar import git_version

__author__ = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

ER3_RUN_DIRECTORY = "/home/detchar/excesspower/ER3"

## \addtogroup laldetchar_py_triggers_excesspower
#@{

def find_online_cache(start, end, channel, **kwargs):
    """Find ExcessPower files from the online GSTLAL analysis
    for the given span

    @param start
        GPS start time for search
    @param end
        GPS end time for search
    @param channel UNDOCUMENTED
    @param kwargs UNDOCUMENTED
    'ifo' observatory for search
    'clustering'
        tag for clustering stage to search, default: unclustered
    'check_files'
        check that the returned files can be read on disk, default False
    """
    out = Cache()

    # set base directory
    directory = kwargs.pop("directory", ER3_RUN_DIRECTORY)
    ifo,channel = channel.split(":", 1)
    channel_dir = os.path.join(directory, ifo, "%s_excesspower" % channel)

    glob_query = "%s-%s_excesspower-*.xml" % (ifo, channel.replace("-", "_"))

    span = Segment(start, end)

    # optimise
    append = out.append
    splitext = os.path.splitext
    isfile = os.path.isfile
    pjoin = os.path.join
    intersects = span.intersects
    from_T050017 = CacheEntry.from_T050017


    # loop over days gathering files
    t = start // 1e4 * 1e4
    while t < end:
        gps_dir = os.path.join(channel_dir, "%.6s" % t)
        if os.path.isdir(gps_dir):
            file_list = glob(os.path.join(gps_dir, glob_query))
            for f in file_list:
                e = from_T050017(f)
                if intersects(e.segment):
                   append(e)
        t += 1e4

    out.sort(key=lambda e: e.segment[0])
    return out

##@}
