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

## \addtogroup laldetchar_py_triggers_omega
"""Read and write `ASCII` files written by the Omega-pipeline.
"""
#
# ### Synopsis ###
#
# ~~~
# from laldetchar.triggers import omega
# ~~~
# \author Duncan Macleod (<duncan.macleod@ligo.org>)

from __future__ import division

import os
import sys
import numpy
import re
import bisect

from socket import getfqdn
from lal import (LIGOTimeGPS, lalStrainUnit)

from glue import segments
from glue.lal import (Cache, CacheEntry)
from glue.ligolw import lsctables

from laldetchar import git_version

__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# set global variables
_comment = re.compile('[#%]')
_delim   = re.compile('[\t\,\s]+')

## \addtogroup laldetchar_py_triggers_omega
#@{

OMEGA_COLUMNS = ['process_id', 'event_id', 'start_time', 'start_time_ns',
                 'peak_time', 'peak_time_ns', 'stop_time', 'stop_time_ns',
                 'duration', 'central_freq', 'peak_frequency',
                 'bandwidth', 'snr', 'amplitude']


def ascii_trigger(line, columns=OMEGA_COLUMNS):
    """Parse a line of ASCII text into a `SnglBurst` trigger

    @param line
        single line of `ASCII` text
    @param columns
        a list of `LIGO_LW` columns to laod, defaults to
        `OMEGA_COLUMNS`

    @returns a `SnglBurst` built from the `ASCII` data
    """
    if isinstance(line, str):
        dat = map(float, _delim.split(line.rstrip()))
    else:
        dat = map(float, line)

    # map to known formats
    if len(dat)==11:
        (peak, freq, duration, band, amplitude,
         cls, cle, cln, av_freq, av_band, err_freq) =  dat
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        peak     = LIGOTimeGPS(peak)
        snr      = (2*amplitude)**(1/2)
        clusters = True
    elif len(dat)==8:
        (peak, freq, duration, band, amplitude, cls, cle, cln) = dat
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        peak     = LIGOTimeGPS(peak)
        av_freq  = freq
        av_band  = band
        err_freq = 0
        snr      = (2*amplitude)**(1/2)
        clusters = True
    elif len(dat)==5:
        (peak, freq, duration, band, amplitude) = dat
        peak = LIGOTimeGPS(peak)
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        av_freq  = freq
        av_band  = band
        err_freq = 0
        snr      = (2*amplitude)**(1/2)
        clusters = False
    else:
        raise ValueError("Wrong number of columns in ASCII line. "
                         "Cannot read.")

    t = lsctables.SnglBurst()
    t.search = u"omega"
    t.event_id = lsctables.SnglBurstTable.get_next_id()

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
        t.peak_frequency = av_freq
    if 'peak_frequency_eror' in columns:
        t.peak_frequency_error = err_freq
    if 'bandwidth' in columns:
        t.bandwidth = av_band
    if 'flow' in columns:
        t.flow = freq-band/2
    if 'fhigh' in columns:
        t.fhigh = freq+band/2

    # set ms frequencies
    if 'ms_bandwidth' in columns:
        t.ms_bandwidth = band
    if 'ms_flow' in columns:
        t.ms_flow = freq-band/2
    if 'ms_fhigh' in columns:
        t.ms_fhigh = freq+band/2

    # set amplitude and snr
    if 'snr' in columns:
        t.snr = snr
    if 'ms_snr' in columns:
        t.ms_snr = snr
    if 'amplitude' in columns:
        t.amplitude = amplitude

    return t


def from_ascii(filename, columns=None, start=None, end=None,
               channel=None):
    """Read Omega triggers from an ASCII file

    Lines in the file are parsed one-by-one, excluding obvious comments
    with each converted to an `OmegaTrigger` (a sub-class of
    `SnglBurst`).

    @param filename
        path to the ASCII file
    @param columns
        a list of valid LIGO_LW column names to load (defaults to all)
    @param start
        minimum GPS time for returned triggers
    @param end
        maximum GPS time for returned triggers
    @param channel
        name of the source data channel for these events

    @returns a LIGO_LW table containing the triggers
    """
    if channel and re.match("[A-Z]\d:", channel):
        ifo = channel[:2]
    else:
        ifo = None

    if not columns:
        columns = OMEGA_COLUMNS
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
        span = segments.segment(start, end)
        check_time = True
    else:
        check_time = False

    # record amplitude if recording SNR

    # generate table
    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    columns = out.columnnames

    # read file and generate triggers
    append = out.append
    with open(filename, 'r') as f:
        for line in f:
            if _comment.match(line):
               continue
            t = ascii_trigger(line, columns=columns)
            if channel:
                t.ifo = ifo
                t.channel = channel
            if not check_time or (float(t.get_peak()) in span):
                append(t)

    return out

##@}

# TODO: remove the functions below if a good trigfind solution is
# implemented

#
# =============================================================================
#
# Event finding
#
# =============================================================================
#

# find Omega files from S6-style online analysis
def find_online_cache(start, end, ifo, mask='DOWNSELECT',
                      check_files=False, **kwargs):
    """Find Omega Online files for the given GPS period.

    @param start
        GPS start time for search
    @param end
        GPS end time for search
    @param ifo
        observatory for search
    @param mask
        description tag of Omega ASCII to search
    @param check_files
        check that the returned files can be read on disk, default False
    @param kwargs UNDOCUMENTED
    """
    out = Cache()

    # verify host
    host = { 'G1':'atlas', 'H1':'ligo-wa', 'H2':'ligo-wa', 'L1':'ligo-la'}
    if (not kwargs.has_key('directory') and not
            re.search(host[ifo],getfqdn())):
        sys.stderr.write("WARNING: Omega online files are not available for "
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
            return ("%s/%.5d/%.10d-%.10d/%s-OMEGA_TRIGGERS_%s-%.10d-%d.txt"
                    % (basedir, gps/100000, gps, gps+dt, ifo, mask, gps, dt))
    else:
        def _omega_file(gps, ifo):
            return ("%s/%s-%s/%s-OMEGA_TRIGGERS_%s-%s-%s.txt"
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


# find files from aLIGO DMTOmega
def find_dmt_cache(start, end, ifo, check_files=False, **kwargs):
    """Find DMTOmega files for the given GPS period.

    @param start
        GPS start time for search
    @param end
        GPS end time for search
    @param ifo
        observatory for search
    @param check_files
        check that the returned files can be read on disk, default False
    @param kwargs UNDOCUMENTED
    """
    out = Cache()

    # verify host
    host = { 'G1':'atlas', 'H1':'ligo-wa', 'H2':'ligo-wa', 'L1':'ligo-la'}
    if (not kwargs.has_key('directory') and not
            re.search(host[ifo],getfqdn())):
        sys.stderr.write("WARNING: Omega online files are not available for "
                         "IFO=%s on this host." % ifo)
        sys.stderr.flush()
        return out

    span = segments.segment(start,end)

    # set known epochs
    known_epochs = {1031340854:55, 1041657635:55, 1041669472:55,
                    1041682187:55, 1044093810:38, 1044111232:38, 1044111282:38,
                    1044112180:38, 1057700030:38, 1057722672:38}

    # get parameters
    epoch = kwargs.pop("epoch", sorted(known_epochs.keys()))
    dt = kwargs.pop("duration", 55)
    try:
        iter(epoch)
    except TypeError:
        epoch = [epoch]
    overlap = kwargs.pop("overlap", 0)
    directory = kwargs.pop("duration",
                           "/gds-%s/dmt/triggers/%s-Omega_Triggers"
                           % (ifo.lower(), ifo[0].upper()))

    # optimise
    append = out.append
    splitext = os.path.splitext
    isfile = os.path.isfile
    intersects = span.intersects
    segment = segments.segment
    from_T050017 = CacheEntry.from_T050017

    # get times
    epoch_idx = bisect.bisect_right(epoch, start)-1
    print epoch_idx
    try:
        dt = known_epochs[epoch[epoch_idx]]
    except KeyError:
        dt = 38
    next_epoch = len(epoch) >= epoch_idx+2  and epoch[epoch_idx+1] or 0
    start_time = int(start-numpy.mod(start-epoch[epoch_idx], dt-overlap))
    t = start_time

    def _omega_file(gps, ifo, deltaT):
        return ("%s/%s-OMEGA_TRIGGERS_CLUSTER-%.5s/"
                "%s-OMEGA_TRIGGERS_CLUSTER-%.10d-%d.xml"
                % (directory, ifo.upper(), gps, ifo.upper(), gps, deltaT))

    # loop over time segments constructing file paths
    while t<end:
        fp = _omega_file(t, ifo, dt)
        if (intersects(segment(t, t+dt)) and
               (not check_files or isfile(fp))):
            append(from_T050017(fp))
        t += dt - overlap
        if next_epoch and t > next_epoch:
            try:
                dt = known_epochs[next_epoch]
            except KeyError:
                dt = 55
            t = next_epoch
            epoch_idx += 1
            next_epoch = len(epoch) >= epoch_idx+2  and epoch[epoch_idx+1] or 0
    out.sort(key=lambda e: e.path)

    return out
