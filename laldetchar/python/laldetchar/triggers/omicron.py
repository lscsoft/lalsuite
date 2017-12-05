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

## \addtogroup laldetchar_py_triggers_omicron
"""Read Omicron events from ROOT files
"""
#
# ### Synopsis ###
#
# ~~~
# from laldetchar.triggers import omicron
# ~~~
# \author Duncan Macleod (<duncan.macleod@ligo.org>)

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

OMICRON_COLUMNS = ["search", "peak_time", "peak_time_ns", "start_time",
                   "start_time_ns",
                   "stop_time", "stop_time_ns", "duration",
                   "central_freq", "peak_frequency",
                   "flow", "fhigh", "bandwidth",
                   "snr", "amplitude", "confidence"]

_re_comment = re.compile("[#%]")
_re_delim = re.compile("[\t\,\s]+")

## \addtogroup laldetchar_py_triggers_omicron
#@{


def ascii_trigger(line, columns=OMICRON_COLUMNS):
    """Parse a line of `ASCII` text into a `SnglBurst` object

    @param line
        single line of `ASCII` text from an Omicron file
    @param columns
        a list of valid `LIGO_LW` column names to load (defaults to all)

    @returns a `SnglBurst` built from the `ASCII` data
    """
    if isinstance(line, str):
        dat = map(float, _re_delim.split(line.rstrip()))
    else:
        dat = map(float, line)

    if len(dat) == 5:
        (peak, freq, duration, band, snr) = dat
        peak = LIGOTimeGPS(peak)
        start = peak - duration/2.
        stop = peak + duration/2.
    else:
        raise ValueError("Wrong number of columns in ASCII line. "
                         "Cannot read.")

    t = lsctables.SnglBurst()
    t.search = u"omicron"

    if 'start_time' in columns:
        t.start_time = start.gpsSeconds
    if 'start_time_ns' in columns:
        t.start_time_ns = start.gpsNanoSeconds
    if 'time' in columns or 'peak_time' in columns:
        t.peak_time = peak.gpsSeconds
    if 'time' in columns or 'peak_time_ns' in columns:
        t.peak_time_ns = peak.gpsNanoSeconds
    if 'stop_time' in columns:
        t.stop_time = stop.gpsSeconds
    if 'stop_time_ns' in columns:
        t.stop_time_ns  = stop.gpsNanoSeconds
    if 'duration' in columns:
        t.duration = duration

    if 'central_freq' in columns:
        t.central_freq = freq
    if 'peak_frequency' in columns:
        t.peak_frequency = freq
    if 'bandwidth' in columns:
        t.bandwidth = band
    if 'flow' in columns:
        t.flow = freq - band/2.
    if 'fhigh' in columns:
        t.fhigh = freq + band/2.

    if 'snr' in columns:
        t.snr = snr
    if 'amplitude' in columns:
        t.amplitude = snr ** 2 / 2.
    if 'confidence' in columns:
        t.confidence = snr
    return t


def root_trigger(root_event, columns=OMICRON_COLUMNS):
    """Parse a `ROOT` tree entry into a `SnglBurst` object.

    @param root_event
        `ROOT` `TChain` object
    @param columns
        a list of valid `LIGO_LW` column names to load (defaults to all)

    @returns a `SnglBurst` built from the `ROOT` data
    """

    t = lsctables.SnglBurst()
    t.search = u"omicron"

    flow = root_event.fstart
    fhigh = root_event.fend
    if 'flow' in columns:
        t.flow = flow
    if 'fhigh' in columns:
        t.fhigh = fhigh
    if 'bandwidth' in columns:
        t.bandwidth = fhigh-flow
    if 'central_freq' in columns:
        t.central_freq = root_event.frequency
    if 'peak_frequency' in columns:
        t.peak_frequency = root_event.frequency

    peak_time = LIGOTimeGPS(root_event.time)
    if 'time' in columns or 'peak_time' in columns:
        t.peak_time = peak_time.gpsSeconds
    if 'time' in columns or 'peak_time_ns' in columns:
        t.peak_time_ns = peak_time.gpsNanoSeconds
    start_time = LIGOTimeGPS(root_event.tstart)
    if 'start_time' in columns:
        t.start_time = start_time.gpsSeconds
    if 'start_time_ns' in columns:
        t.start_time_ns = start_time.gpsNanoSeconds
    stop_time = LIGOTimeGPS(root_event.tend)
    if 'stop_time' in columns:
        t.stop_time = stop_time.gpsSeconds
    if 'stop_time_ns' in columns:
        t.stop_time_ns = stop_time.gpsNanoSeconds
    if 'duration' in columns:
        t.duration = float(stop_time-start_time)

    if 'snr' in columns:
        t.snr = root_event.snr
    if 'amplitude' in columns:
        t.amplitude = root_event.snr**2 / 2.
    if 'confidence' in columns:
        t.confidence = root_event.snr

    return t


def from_ascii(filename, columns=None, start=None, end=None,
               channel=None):
    """Read Omicron triggers from an ASCII file

    Lines in the file are parsed one-by-one, excluding obvious comments
    with each converted to an `OmicronTrigger` (a sub-class of
    `SnglBurst`).

    @param filename
        path to the `ASCII` file
    @param columns
        a list of valid `LIGO_LW` column names to load (defaults to all)
    @param start
        minimum GPS time for returned triggers
    @param end
        maximum GPS time for returned triggers
    @param channel
        name of the source data channel for these events.

    @returns a `LIGO_LW` table containing the triggers
    """
    if channel and re.match("[A-Z]\d:", channel):
        ifo = channel[:2]
    else:
        ifo = None

    if not columns:
        columns = OMICRON_COLUMNS
    if columns:
        columns = set(columns)
        if 'peak_time' not in columns and (start or end):
            columns.add("peak_time")
            columns.add("peak_time_ns")
        if 'snr' in columns:
            columns.add('amplitude')
        if channel:
            columns.update(['ifo', 'channel'])
    if (start or end):
        start = start or segments.NegInfinity
        end = end or segments.PosInfinity
        span = Segment(start, end)
        check_time = True
    else:
        check_time = False

    # generate table
    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    columns = out.columnnames

    # read file and generate triggers
    append = out.append
    next_id = out.get_next_id
    with open(filename, 'r') as f:
        for line in f:
            if _comment.match(line):
               continue
            t = ascii_trigger(line, columns=columns)
            t.event_id = next_id()
            if channel:
                t.ifo = ifo
                t.channel = channel
            if not check_time or (float(t.get_peak()) in span):
                append(t)

    return out

def from_root(filename, columns=None, start=None, end=None,
              channel=None):
    """Read Omicron triggers from a `ROOT` file

    @param filename
        file path from which to read the data
    @param start
        GPS start time after which to restrict returned events
    @param end
        GPS end time before which to restrict returned
    @param channel
        source channel that produced the given data
    @param columns
        set of valid SnglBurst columns to read from data

    @returns a `SnglBurstTable` representing the data
    """
    if channel and re.match("[A-Z]\d:", channel):
        ifo = channel[:2]
    else:
        ifo = None

    if not columns:
        columns = OMICRON_COLUMNS
    if columns:
        columns = set(columns)
        if 'peak_time' not in columns and (start or end):
            columns.add("peak_time")
            columns.add("peak_time_ns")
        if 'snr' in columns:
            columns.add('amplitude')
        if channel:
            columns.update(['ifo', 'channel'])
    if (start or end):
        start = start or segments.NegInfinity
        end = end or segments.PosInfinity
        span = Segment(start, end)
        check_time = True
    else:
        check_time = False

    # generate table
    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    columns = out.columnnames
    append = out.append
    next_id = out.get_next_id

    # read file and generate triggers
    root_tree = TChain("triggers")
    root_tree.Add(filename)
    nevents = root_tree.GetEntries()
    for i in range(nevents):
        root_tree.GetEntry(i)
        t = root_trigger(root_tree, columns=columns)
        t.event_id = next_id()
        if channel:
            t.ifo = ifo
            t.channel = channel
        if not check_time or (float(t.get_peak()) in span):
            append(t)
    return out

##@}
