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
# with this program; if not, write to the Free Software Foundation,

## \addtogroup laldetchar_py_triggers_kleinewelle
"""Read KleineWelle events from `ASCII` files

The module provides a method to import triggers from `ASCII` into a
`SnglBurstTable`. Eventually, writing to `ASCII` will be supported also.
"""
#
# ### Synopsis ###
#
# ~~~
# from laldetchar.triggers import kw
# ~~~
# \author Duncan Macleod (<duncan.macleod@ligo.org>)

from __future__ import division


import os
import sys
import numpy
import re
import bisect

from socket import getfqdn
from lal.lal import LIGOTimeGPS

from glue import segments
from glue.lal import (Cache, CacheEntry)
from glue.ligolw import lsctables

from laldetchar import git_version

__author__  = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# set global variables
_comment = re.compile('[#%]')
_delim   = re.compile('[\t\,\s]+')

KLEINEWELLE_COLUMNS = ["search", "peak_time", "peak_time_ns", "start_time",
                       "start_time_ns", "stop_time", "stop_time_ns",
                       "duration", "central_freq", "peak_frequency",
                       "bandwidth", "snr", "amplitude", "confidence"]


## \addtogroup laldetchar_py_triggers_kleinewelle
#@{


def ascii_trigger(line, columns=KLEINEWELLE_COLUMNS):
    """Parse a line of `ASCII` text into a `SnglBurst` object

    @param line
        single line of `ASCII` text from a KleineWelle file
    @param columns
        a list of valid `LIGO_LW` column names to load (defaults to all)

    @returns a `SnglBurst` built from the `ASCII` data
    """
    t = lsctables.SnglBurst()
    t.search = u"kleinewelle"

    dat = line.rstrip().split()
    if len(dat) == 9:
        channel = re.sub("_", ":", dat.pop(-1), 1)
        if re.search("_\d+_\d+\Z", channel):
            channel = channel.rsplit("_", 2)[0]
        if 'channel' in columns:
            t.channel = channel
        if 'ifo' in columns and re.match("[A-Z]\d:", channel):
            t.ifo = channel[:2]
        elif 'ifo' in columns:
            t.ifo = None
    if len(dat) == 8:
        (start, stop, peak, freq, energy,
         amplitude, n_pix, sig) = list(map(float, dat))
        start = LIGOTimeGPS(start)
        stop = LIGOTimeGPS(stop)
        peak = LIGOTimeGPS(peak)
        duration = float(stop-start)
        snr = (amplitude-n_pix)**(1/2)
        band = 0
    else:
        raise ValueError("Wrong number of columns in ASCII line. "
                         "Cannot read.")

    # set times
    if 'start_time' in columns:
        t.start_time = start.gpsSeconds
    if 'start_time_ns' in columns:
        t.start_time_ns = start.gpsNanoSeconds
    if 'peak_time' in columns:
        t.peak_time = peak.gpsSeconds
    if 'peak_time_ns' in columns:
        t.peak_time_ns = peak.gpsNanoSeconds
    if 'stop_time' in columns:
        t.stop_time = stop.gpsSeconds
    if 'stop_time_ns' in columns:
        t.stop_time_ns  = stop.gpsNanoSeconds

    # set ms times
    if 'ms_start_time' in columns:
        t.ms_start_time = start.gpsSeconds
    if 'ms_start_time_ns' in columns:
        t.ms_start_time_ns = start.gpsNanoSeconds
    if 'ms_stop_time' in columns:
        t.ms_stop_time = stop.gpsSeconds
    if 'ms_stop_time_ns' in columns:
        t.ms_stop_time_ns = stop.gpsNanoSeconds

    # set other times
    if 'duration' in columns:
        t.duration = duration
    if 'ms_duration' in columns:
        t.ms_duration = duration

    # set frequencies
    if 'central_freq' in columns:
        t.central_freq = freq
    if 'peak_frequency' in columns:
        t.peak_frequency = freq
    if 'peak_frequency_eror' in columns:
        t.peak_frequency_error = 0
    if 'bandwidth' in columns:
        t.bandwidth = band
    if 'flow' in columns:
        t.flow = freq - band/2.
    if 'fhigh' in columns:
        t.fhigh = freq + band/2.

    # set ms frequencies
    if 'ms_bandwidth' in columns:
        t.ms_bandwidth = band
    if 'ms_flow' in columns:
        t.ms_flow = freq - band/2.
    if 'ms_fhigh' in columns:
        t.ms_fhigh = freq + band/2.

    # set amplitude and snr
    if 'snr' in columns:
        t.snr = snr
    if 'ms_snr' in columns:
        t.ms_snr = snr
    if 'amplitude' in columns:
        t.amplitude = amplitude
    if "confidence" in columns:
        t.confidence = sig
    if "ms_confidence" in columns:
        t.ms_confidence = sig
    if "hrss" in columns:
        t.hrss = 0
    if "ms_hrss" in columns:
        t.ms_hrss = 0

    # set other params
    if 'n_pix' in columns or 'param_one_value' in columns:
        t.param_one_name = "n_pix"
        t.param_one_value = n_pix
    if 'signifiance' in columns or 'param_two_value' in columns:
        t.param_two_name = "significance"
        t.param_two_value = sig

    return t


def from_ascii(filename, columns=None, start=None, end=None,
               channel=None):
    """Read KleineWelle triggers from an `ASCII` file

    Lines in the file are parsed one-by-one, excluding obvious comments
    with each converted to an `KleineWelleTrigger` (a sub-class of
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
    if isinstance(channel, basestring):
        channels = [channel]
    elif channel is None:
        channels = [None]
    else:
        channels = channel

    if not columns:
        columns = KLEINEWELLE_COLUMNS
    if columns:
        if channel:
            columns.append('channel')
            columns.append('ifo')
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
        span = segments.segment(start, end)
        check_time = True
    else:
        check_time = False

    out = dict()

    for ch in channels:
        # generate table
        out[ch] = lsctables.New(lsctables.SnglBurstTable, columns=columns)
        columns = out[ch].columnnames

    # read file and generate triggers
    append = dict((ch, out[ch].append) for ch in channels)
    next_id = out[ch].get_next_id
    with open(filename, 'r') as f:
        for line in f:
            if _comment.match(line):
               continue
            t = ascii_trigger(line, columns=columns)
            if channel and t.channel not in channels:
                continue
            t.event_id = next_id()
            if not check_time or (float(t.get_peak()) in span):
                append[t.channel](t)

    if isinstance(channel, basestring) or channel is None:
        return out[channel]
    else:
        return out


# TODO: remove following lines when a good trigfind method is in place
##@}

#
# =============================================================================
#
# Event finding
#
# =============================================================================
#

# find KW files from S6-style online analysis
def find_online_cache(start, end, ifo, mask='DOWNSELECT',
                      check_files=False, **kwargs):
    """Find KW Online files for the given GPS period.

    @param start
        GPS start time for search
    @param end
        GPS end time for search
    @param ifo
        observatory for search
    @param mask
        description tag of KW ASCII to search
    @param check_files
        check that the returned files can be read on disk, default False
    @param kwargs UNDOCUMENTED
    """
    out = Cache()

    # verify host
    host = { 'G1':'atlas', 'H1':'ligo-wa', 'H2':'ligo-wa', 'L1':'ligo-la'}
    if (not kwargs.has_key('directory') and not
            re.search(host[ifo],getfqdn())):
        sys.stderr.write("WARNING: KW online files are not available for "
                         "IFO=%s on this host." % ifo)
        sys.stderr.flush()
        return out

    span = segments.segment(start,end)
    # get parameters
    dt = kwargs.pop("duration", 64)
    overlap = kwargs.pop("overlap", 8)
    if ifo == "G1":
        directory = kwargs.pop("directory", "/home/omega/online/G1/segments")
        epoch = kwargs.pop("epoch", 983669456)
    else:
        directory = kwargs.pop("directory",\
                               "/home/omega/online/%s/archive/S6/segments"
                               % ifo)

    # optimise
    append = out.append
    splitext = os.path.splitext
    isfile = os.path.isfile
    intersects = span.intersects
    segment = segments.segment
    from_T050017 = CacheEntry.from_T050017

    # get times
    start_time = int(start-numpy.mod(start-epoch, dt-overlap))
    t = start_time

    if ifo == "G1":
        def _omega_file(gps, ifo):
            return ("%s/%.5d/%.10d-%.10d/%s-KW_TRIGGERS_%s-%.10d-%d.txt"
                    % (basedir, gps/100000, gps, gps+dt, ifo, mask, gps, dt))
    else:
        def _omega_file(gps, ifo):
            return ("%s/%s-%s/%s-KW_TRIGGERS_%s-%s-%s.txt"
                    % (basedir, gps, gps+dt, ifo, mask, gps, dt))

    # loop over time segments constructing file paths
    while t<end:
        fp = _omega_file(t, ifo)
        if (intersects(segment(t, t+dt)) and
               (not check_files or isfile(fp))):
            append(from_T050017(fp))
        t += dt - overlap
    out.sort(key=lambda e: e.path)

    return out


# find files from aLIGO DMTKW
def find_dmt_cache(start, end, ifo, extension="xml", check_files=False,
                   **kwargs):
    """Find DMT KW files for the given GPS period.

    @param start
        GPS start time for search
    @param end
        GPS end time for search
    @param ifo
        observatory for search
    @param extension UNDOCUMENTED
    @param check_files
        check that the returned files can be read on disk, default False
    @param kwargs UNDOCUMENTED
    """
    out = Cache()

    # verify host
    host = { 'G1':'atlas', 'H1':'ligo-wa', 'H2':'ligo-wa', 'L1':'ligo-la'}
    if (not kwargs.has_key('directory') and not
            re.search(host[ifo],getfqdn())):
        sys.stderr.write("WARNING: KW online files are not available for "
                         "IFO=%s on this host." % ifo)
        sys.stderr.flush()
        return out

    span = segments.segment(start,end)

    # set known epochs
    known_epochs = [1026263104]

    # get parameters
    dt = int(kwargs.pop("duration", 64))
    epoch = kwargs.pop("epoch", known_epochs)
    filetag = kwargs.pop("filetag", "KW_TRIGGERS")
    dirtag = filetag.endswith("_TRENDS") and filetag[:-7] or filetag
    try:
        iter(epoch)
    except TypeError:
        epoch = [int(epoch)]
    overlap = int(kwargs.pop("overlap", 0))
    directory = kwargs.pop("duration",
                           "/gds-%s/dmt/triggers/%s-%s"
                           % (ifo.lower(), ifo[0].upper(), dirtag))

    # optimise
    append = out.append
    splitext = os.path.splitext
    isfile = os.path.isfile
    intersects = span.intersects
    segment = segments.segment
    from_T050017 = CacheEntry.from_T050017

    # get times
    epoch_idx = bisect.bisect_right(epoch, start)-1
    next_epoch = len(epoch) >= epoch_idx+2  and epoch[epoch_idx+1] or 0
    start_time = int(start-numpy.mod(start-epoch[epoch_idx], dt-overlap))
    t = start_time

    def _kw_file(gps, ifo):
        return ("%s/%s-%s-%.5s/"
                "%s-%s-%.10d-%d.%s"
                % (directory, ifo.upper()[0], dirtag, gps,
                   ifo.upper()[0], filetag, gps, dt, extension))

    # loop over time segments constructing file paths
    while t<end:
        fp = _kw_file(t, ifo)
        if (intersects(segment(t, t+dt)) and
               (not check_files or isfile(fp))):
            append(from_T050017(fp))
        t += dt - overlap
        if next_epoch and t > next_epoch:
            t = next_epoch
            epoch_idx += 1
            next_epoch = len(epoch) >= epoch_idx+2  and epoch[epoch_idx+1] or 0
    out.sort(key=lambda e: e.path)

    return out

#@}
