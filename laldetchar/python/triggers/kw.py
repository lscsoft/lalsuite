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

"""Utilites for manipulating the output of the KleineWelle algorithm
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
_xml = re.compile("(xml|xml.gz)\Z")


#
# =============================================================================
#
# Event reading
#
# =============================================================================
#

# convert ascii line to omega trigger
def trigger(line, columns=lsctables.SnglBurst.__slots__, channel=None):
    """Read a SnglBurst event from a line of KW-format ASCII

    @param line
        string of ASCII to process, or list of attributes in ASCII order
    @param columns
        list of valid columns to load, defaults to ALL for SnglBurst

    @returns a glue.lsctables.SnglBurst object representing this event
    """
    if isinstance(line, str):
        dat = _delim.split(line.rstrip())

    if len(dat) == 9:
        c = re.sub("_", ":", dat.pop(-1), 1)
        if re.search("_\d+_\d+\Z", c):
            c = c.rsplit("_", 2)[0]
        if channel and channel not in c:
            return
        elif not channel:
            channel = c
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
        raise ValueError("Wrong number of columns in ASCII line. Cannot read.")
    
    # set object
    t = lsctables.SnglBurst()
    t.channel = channel

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


# read triggers from file
def from_file(fobj, start=None, end=None, ifo=None, channel=None,
              columns=None):
    """Read a SnglBurstTable from KW-format ASCII file

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
    if not channel or isinstance(channel, str):
        channels = [channel]
        return_dict = False
    else:
        return_dict = True
        channels = channel

    out = dict((c, lsctables.New(lsctables.SnglBurstTable, columns=columns)) for
               c in channels)

    # remove unused names and types
    if usercolumns:
        for key,tab in out.iteritems():
            for c in out[key].columnnames:
                if c.lower() not in columns:
                    idx = out[key].columnnames.index(c)
                    out[key].columnnames.pop(idx)
                    out[key].columntypes.pop(idx)

    # read file and generate triggers
    if _xml.search(fobj.name):
        pass
    else:
        for i,line in enumerate(fobj):
            if _comment.match(line):
                continue
            t = trigger(line, columns=columns)
            t.ifo = ifo
            t.search = "KleineWelle"
            if (t and t.channel in channels and (not check_time or
                       (check_time and float(t.get_peak()) in span))):
                if not channel:
                    out[channel].append(t)
                else:
                    out[t.channel].append(t)

    if return_dict:
        return out
    else:
        return out[channels[0]]    


# read triggers from a list of files 
def from_files(filelist, start=None, end=None, ifo=None, channel=None,
               columns=None, verbose=False):
    """Read a SnglBurstTable from a list of KW-format ASCII files

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

    # generate table for multi channel read
    if not channel or isinstance(channel, str):
        channels = [channel]
        return_dict = False
    else:
        channels = channel
        return_dict = True
    out = dict((c, lsctables.New(lsctables.SnglBurstTable, columns=columns)) for
               c in channels)

    # load results
    for i,fp in enumerate(filelist):
        with open(fp, "r") as f:
            new_trigs = from_file(f, start=start, end=end, columns=columns,\
                                  ifo=ifo, channel=channels)
            for channel in channels:
                out[channel].extend(new_trigs[channel])
        if verbose:
            progress = int((i+1)/num*100)
            sys.stdout.write("Extracting KW triggers from %d files... "
                             "%.2d%%\r" % (num, progress))
            sys.stdout.flush()
    if verbose:
        sys.stdout.write("Extracting KW triggers from %d files... "
                         "100%%\n" % (num))
        sys.stdout.flush()

    if return_dict:
        return out
    else:
        return out[channel]

def from_lal_cache(cache, start=None, end=None, ifo=None, channel=None,
                   columns=None, verbose=False):
    """Read a SnglBurstTable from a Cache of KW-format ASCII files
    
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

    @returns a glue.lsctables.SnglBurstTable object representing the data
    """
    return from_files(cache.pfnlist(), start=start, end=end, ifo=ifo,
                      channel=channel, columns=columns, verbose=verbose)

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
def find_dmt_cache(start, end, ifo, format="xml", check_files=False, **kwargs):
    """Find DMT KW files for the given GPS period.

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
        sys.stderr.write("WARNING: KW online files are not available for "
                         "IFO=%s on this host." % ifo)
        sys.stderr.flush()
        return out

    span = segments.segment(start,end)

    # set known epochs
    known_epochs = [1026263104]

    # get parameters
    dt = kwargs.pop("duration", 64)
    epoch = kwargs.pop("epoch", known_epochs)
    try:
        iter(epoch)
    except TypeError:
        epoch = [epoch]
    overlap = kwargs.pop("overlap", 0)
    directory = kwargs.pop("duration",
                           "/gds-%s/dmt/triggers/%s-KW_TRIGGERS"
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
        return ("%s/%s-KW_TRIGGERS-%.5s/"
                "%s-KW_TRIGGERS-%.10d-%d.%s"
                % (directory, ifo.upper()[0], gps, ifo.upper()[0], gps, dt,
                   format))

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
