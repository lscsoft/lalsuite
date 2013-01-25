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

"""Utilites for manipulating the output of the Omega-pipeline
"""

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


#
# =============================================================================
#
# Event reading
#
# =============================================================================
#

# convert ascii line to omega trigger
def trigger(line, columns=lsctables.SnglBurst.__slots__, virgo=False):
    """Read a SnglBurst event from a line of Omega-format ASCII

    @param line
        string of ASCII to process, or list of attributes in ASCII order
    @param columns
        list of valid columns to load, defaults to ALL for SnglBurst
    @param virgo
        identifies ASCII as in Virgo-Omega format, default False

    @returns a glue.lsctables.SnglBurst object representing this event
    """
    if isinstance(line, str):
        dat = map(float, _delim.split(line.rstrip()))
    else:
        dat = map(float, line)

    # map to known formats
    if virgo:
        (start, stop, peak, freq, band, cln, cle, snr) = dat
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        peak     = LIGOTimeGPS(peak)
        av_freq  = freq
        av_band  = band
        err_freq = 0
        clusters = False
    elif len(dat)==11:
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
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        av_freq  = freq
        av_band  = band
        err_freq = 0
        snr      = (2*amplitude)**(1/2)
        clusters = False
    else:
        raise ValueError("Wrong number of columns in ASCII line. Cannot read.")
    
    # set object
    t = lsctables.SnglBurst()

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
     
    # set other params
    if 'cluster_size' in columns or 'param_one_value' in columns:
        t.param_one_name = 'cluster_size'
        if clusters:
            t.param_one_value = cls
        else:
            t.param_one_value = numpy.NaN
    if 'cluster_norm_energy' in columns or 'param_two_value' in columns:
        t.param_two_name = 'cluster_norm_energy'
        if clusters:
            t.param_two_value = cle
        else:
            t.param_two_value = numpy.NaN
    if 'cluster_number' in columns or 'param_three_value' in columns:
        t.param_three_name = 'cluster_number'
        if clusters:
            t.param_three_value = cln
        else:
            t.param_three_value = numpy.NaN

    return t


# read triggers from file
def from_file(fobj, start=None, end=None, ifo=None, channel=None,
              columns=None, virgo=False):
    """Read a SnglBurstTable from Omega-format ASCII file

    @param fobj
        file object from which to read the data
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
    @param virgo
        identifies ASCII as in Virgo-Omega format, default False

    @returns a glue.lsctables.SnglBurstTable object representing the data
    """
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
        span = segments.segment(start, end)
        columns.update(["peak_time", "peak_time_ns"])
        check_time = True
    else:
        check_time = False

    # record amplitude if recording SNR
    if 'snr' in columns:
        columns.add("amplitude")

    # generate table
    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    append = out.append

    # remove unused names and types
    if usercolumns:
        for c in out.columnnames:
            if c.lower() not in columns:
                idx = out.columnnames.index(c)
                out.columnnames.pop(idx)
                out.columntypes.pop(idx)

    # read file and generate triggers
    for i,line in enumerate(fobj):
        if _comment.match(line):
            continue
        t = trigger(line, columns=columns, virgo=virgo)
        if not check_time or (check_time and float(t.get_peak()) in span):
            append(t)

    return out


# read triggers from a list of files 
def from_files(filelist, start=None, end=None, ifo=None, channel=None,
               columns=None, verbose=False, virgo=False):
    """Read a SnglBurstTable from a list of Omega-format ASCII files

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
    @param virgo
        identifies ASCII as in Virgo-Omega format, default False

    @returns a glue.lsctables.SnglBurstTable object representing the data
    """
    if verbose:
        sys.stdout.write("Extracting Omega triggers from %d files...     \r"
                         % len(filelist))
        sys.stdout.flush()
        num = len(filelist)/100
    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    extend = out.extend
    for i,fp in enumerate(filelist):
        with open(fp, "r") as f:
            extend(from_file(f, start=start, end=end, columns=columns,\
                             ifo=ifo, channel=channel, virgo=virgo))
        if verbose:
            progress = int((i+1)/num)
            sys.stdout.write("Extracting Omega triggers from %d files... "
                             "%.2d%%\r" % (num, progress))
            sys.stdout.flush()
    if verbose:
        sys.stdout.write("Extracting Omega triggers from %d files... "
                         "100%%\r" % (num))
        sys.stdout.flush()
    return out

def from_lal_cache(cache, start=None, end=None, ifo=None, channel=None,
                   columns=None, verbose=False, virgo=False):
    """Read a SnglBurstTable from a Cache of Omega-format ASCII files
    
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
    @param virgo
        identifies ASCII as in Virgo-Omega format, default False

    @returns a glue.lsctables.SnglBurstTable object representing the data
    """
    return from_files(cache.pfnlist(), start=start, end=end, ifo=ifo,
                      channel=channel, columns=columns, verbose=verbose,
                      virgo=virgo)

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
    known_epochs = [1031340854, 1041657635, 1041669472, 1041682187]

    # get parameters
    dt = kwargs.pop("duration", 55)
    epoch = kwargs.pop("epoch", known_epochs)
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
    next_epoch = len(epoch) >= epoch_idx+2  and epoch[epoch_idx+1] or 0
    start_time = int(start-numpy.mod(start-epoch[epoch_idx], dt-overlap))
    t = start_time

    def _omega_file(gps, ifo):
        return ("%s/%s-OMEGA_TRIGGERS_CLUSTER-%.5s/"
                "%s-OMEGA_TRIGGERS_CLUSTER-%.10d-%d.xml"
                % (directory, ifo.upper(), gps, ifo.upper(), gps, dt))

    # loop over time segments constructing file paths
    while t<end:
        fp = _omega_file(t, ifo)
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
