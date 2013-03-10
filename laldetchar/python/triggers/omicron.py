# Copyright (C) 2012 Duncan Macleod
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

"""Read and manipulate Omicron events from ROOT files
"""

from __future__ import division

import os
import re

from ROOT import TChain

from glue.ligolw import ilwd,lsctables,table,utils
from glue.segments import segment as Segment

from lal import LIGOTimeGPS
from laldetchar import git_version

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

OMICRON_COLUMNS = ["process_id", "search", "channel", "ifo",
                   "peak_time", "peak_time_ns", "start_time", "start_time_ns",
                   "stop_time", "stop_time_ns", "duration",
                   "central_freq", "flow", "fhigh", "bandwidth",
                   "snr", "amplitude", "confidence"]


def get_sngl_burst(root_event, columns=OMICRON_COLUMNS):
    """@returns a LIGO_LW SnglBurst event with attributes seeded from
    the given Omicron ROOT tree event
    """
    sb = lsctables.SnglBurst()
    if "search" in columns:
        sb.search = u"omicron"

    flow = root_event.fstart
    fhigh = root_event.fend
    if "flow" in columns:
        sb.flow = flow
    if "fhigh" in columns:
        sb.fhigh = fhigh
    if "bandwidth" in columns:
        sb.bandwidth = fhigh-flow
    if "central_freq" in columns:
        sb.central_freq = root_event.frequency

    peak_time = LIGOTimeGPS(root_event.time)
    if "time" in columns or "peak_time" in columns:
        sb.peak_time = peak_time.gpsSeconds
    if "time" in columns or "peak_time_ns" in columns:
        sb.peak_time_ns = peak_time.gpsNanoSeconds
    start_time = LIGOTimeGPS(root_event.tstart)
    if "start_time" in columns:
        sb.start_time = start_time.gpsSeconds
    if "start_time_ns" in columns:
        sb.start_time_ns = start_time.gpsNanoSeconds
    stop_time = LIGOTimeGPS(root_event.tend)
    if "stop_time" in columns:
        sb.stop_time = stop_time.gpsSeconds
    if "stop_time_ns" in columns:
        sb.stop_time_ns = stop_time.gpsNanoSeconds
    if "duration" in columns:
        sb.duration = float(stop_time-start_time)

    if "snr" in columns:
        sb.snr = root_event.snr
    if "amplitude" in columns:
        sb.amplitude = root_event.snr**2 / 2.
    if "confidence" in columns:
        sb.confidence = root_event.snr

    return sb


def from_root_file(filename, start=None, end=None, ifo=None, channel=None,
                   columns=OMICRON_COLUMNS):
    """Read a SnglBurstTable from an Omicron ROOT file

    @param filename
        file path from which to read the data
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param ifo
        observatory that produced the given data
    @param channel
        source channel that produced the given data
    @param columns
        set of valid SnglBurst columns to read from data

    @returns a Multi or Sngl BurstTable object representing the data
    """
    if channel and not ifo and re.match("[A-Z]\d:", channel):
       ifo = channel[:2]

    # set columns
    if columns is None:
        columns = lsctables.SnglBurst.__slots__
        usercolumns = False
    else:
        usercolumns = True
    columns = set(columns)

    if start or end:
        if start is None:
            start = -numpy.inf
        if end is None:
            end = numpy.inf
        span = Segment(start, end)
        columns.update(["peak_time", "peak_time_ns"])
        check_time = True
    else:
        check_time = False

    # read file and generate triggers
    root_tree = TChain("triggers")
    root_tree.Add(filename)

    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    append = out.append

    # generate table
    if usercolumns:
        for c in out.columnnames:
            if c.lower() not in columns:
                idx = out.columnnames.index(c)
                out.columnnames.pop(idx)
                out.columntypes.pop(idx)

    # read table
    nevents = root_tree.GetEntries()
    for i in range(nevents):
        root_tree.GetEntry(i)
        t = get_sngl_burst(root_tree, columns=columns)
        if check_time and (float(t.get_peak()) not in span):
            continue
        if ifo:
            t.ifo = ifo
        if channel:
            t.channel = channel
        append(t)

    return out


def from_file(filename, start=None, end=None, ifo=None, channel=None,
              columns=OMICRON_COLUMNS):
    """Read a SnglBurstTable from an Omicron ROOT or ASCII file

    @param filename
        file path from which to read the data
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param ifo
        observatory that produced the given data
    @param channel
        source channel that produced the given data
    @param columns
        set of valid SnglBurst columns to read from data

    @returns a Multi or Sngl BurstTable object representing the data
    """
    if os.path.splitext(filename)[1] == ".root":
         return from_root_file(filename, start=start, end=end, ifo=ifo,
                               channel=channel, columns=columns)
    else:
         raise RuntimeError("Unable to parse Omicron-format ASCII files, "
                            "please code up the 'from_ascii_file' for "
                            "Omicron and patch this module...")


def from_files(filelist, start=None, end=None, ifo=None, channel=None,
               columns=OMICRON_COLUMNS, verbose=False):
    """Read a BurstTable from a list of Omicron-format ROOT or ASCII files

    @param filelist
        list of filepaths from which to read the data
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param ifo
        observatory that produced the given data
    @param channel
        source channel that produced the given data
    @param columns
        set of valid SnglBurst columns to read from data
    @param verbose
        print verbose progress, default False

    @returns a glue.lsctables.SnglBurstTable object representing the data
    """
    if verbose:
        sys.stdout.write("Extracting KW triggers from %d files...     \r"
                         % len(filelist))
        sys.stdout.flush()
        num = len(filelist)

    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    extend = out.extend

    # load results
    for i,fp in enumerate(filelist):
        extend(from_file(fp, start=start, end=end, columns=columns,
                         ifo=ifo, channel=channel))
        if verbose:
            progress = int((i+1)/num*100)
            sys.stdout.write("Extracting KW triggers from %d files... "
                             "%.2d%%\r" % (num, progress))
            sys.stdout.flush()
    if verbose:
        sys.stdout.write("Extracting KW triggers from %d files... "
                         "100%%\n" % (num))
        sys.stdout.flush()

    return out


def from_lal_cache(cache, start=None, end=None, ifo=None, channel=None,
                   columns=OMICRON_COLUMNS, verbose=False):
    """Read a SnglBurstTable from a Cache of Omicron ROOT files

    @param cache
        glue.lal.Cache of filepaths from which to read the data
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param ifo
        observatory that produced the given data
    @param channel
        source channel that produced the given data
    @param columns
        set of valid SnglBurst columns to read from data
    @param verbose
        print verbose progress, default False

    @returns a Multi or Sngl BurstTable object representing the data
    """
    return from_files(cache.pfnlist(), start=start, end=end, ifo=ifo,
                      channel=channel, columns=columns, verbose=verbose)
