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

"""Utitilies to find and manipulate triggers from Daily ihope
"""

import os

from datetime import datetime
from glob import glob

import lal

from glue.lal import (Cache, CacheEntry)
from glue.segments import segment as Segment

from laldetchar import (git_version, triggers)

__author__ = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


def find_daily_cache(start, end, ifo, clustering=None, check_files=False,
                     **kwargs):
    """Find Daily ihope files from the daily runs for the given span

    @param start
        GPS start time for search
    @param end
        GPS end time for search
    @param ifo
        observatory for search
    @param clustering
        tag for clustering stage to search, default: unclustered
    @param check_files
        check that the returned files can be read on disk, default False
    @param kwargs UNDOCUMENTED
    """
    out = Cache()

    # set clustering tag
    if clustering==None or clustering.upper()=='UNCLUSTERED':
        file_tag='INSPIRAL_UNCLUSTERED'
    elif clustering.upper() in ["100MS", "100MILLISEC"]:
        file_tag='INSPIRAL_100MILLISEC_CLUSTERED'
    elif clustering.upper() in ["30MS", "30MILLISEC"]:
        file_tag='INSPIRAL_30MILLISEC_CLUSTERED'
    elif clustering.upper() in ["16S", "16SECOND"]:
        file_tag='INSPIRAL_16SEC_CLUSTERED'

    # set base directory
    directory = kwargs.pop("directory", os.path.expanduser("~cbc/ihope_daily"))

    # work out days
    span = Segment(start, end)
    start = int(start)
    start_d = lal.UTCToGPS(datetime(*lal.GPSToUTC(start)[:6]).replace(
                                       hour=0, minute=0, second=0).timetuple())
    days = []
    day = start_d
    while day <= end:
        days.append(day)
        day+=86400

    # optimise
    append = out.append
    splitext = os.path.splitext
    isfile = os.path.isfile
    pjoin = os.path.join
    intersects = span.intersects
    from_T050017 = CacheEntry.from_T050017

    # loop over days gathering files
    for day in days:
        utc = datetime(*lal.GPSToUTC(day)[:6])
        day_path = pjoin(directory, utc.strftime("%Y%m"),
                         utc.strftime("%Y%m%d"))
        day_cache = os.path.join(day_path, "%s-%s.cache" % (ifo, file_tag))
        if isfile(day_cache):
            with open(day_cache, "r") as f:
                filenames = Cache.fromfile(f).pfnlist()
        else:
            filenames = glob(os.path.join(day_path,
                                               ("%s-%s-*.xml.gz"
                                                % (ifo, file_tag))))
        for filename in filenames:
            e = from_T050017(filename)
            if intersects(e.segment):
                append(e)

    out.sort(key=lambda e: e.path)
    return out

def from_file(fileobj, start=None, end=None, columns=None):
    """@returns a SnglInspiralTable of events from the given fileobj

    @param fileobj UNDOCUMENTED
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param columns
        set of valid SnglInspiral columns to read from data
    """
    return triggers.load_table_from_fileobj(fileobj, "daily_ihope",
                                            start=start, end=end,
                                            columns=columns)


def from_lal_cache(cache, start=None, end=None, columns=None, verbose=False):
    """Read a SnglInspiralTable from a Cache of daily ihope XML files

    @param cache
        glue.lal.Cache of filepaths from which to read the data
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param columns
        set of valid SnglInspiral columns to read from data
    @param verbose
        print verbose progress, default False

    @returns a glue.lsctables.SnglBurstTable object representing the data
    """
    return triggers.load_table_from_lal_cache(cache, etg="daily_ihope",
                                              columns=columns,
                                              verbose=verbose,
                                              start=start, end=end)


def from_files(filelist, start=None, end=None, columns=None, verbose=False):
    """Read a SnglInspiralTable from a list of daily ihope XML files

    @param filelist
        glue.lal.Cache of filepaths from which to read the data
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param columns
        set of valid SnglInspiral columns to read from data
    @param verbose
        print verbose progress, default False

    @returns a glue.lsctables.SnglBurstTable object representing the data
    """
    return triggers.load_table_from_files(filelist, etg="daily_ihope",
                                          columns=columns, verbose=verbose,
                                          start=start, end=end)
