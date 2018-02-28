#!/usr/bin/env python

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

"""
Auxiliary functions for running and postprocessing Omega pipeline code.
"""

from __future__ import division

import os
import sys
import numpy
import re

from socket import getfqdn
from lal import LIGOTimeGPS, lalStrainUnit

from pylal import git_version
from glue.ligolw import lsctables
from glue import segments
from glue.lal import Cache,CacheEntry

__author__  = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = git_version.id
__date__    = git_version.date

# =============================================================================
# convert ascii line to omega trigger
# =============================================================================

_comment = re.compile('[#%]')
_delim   = re.compile('[\t\,\s]+')

def trigger(line, columns=lsctables.SnglBurst.__slots__, virgo=False, ifo=None, channel=None):
    """
    Convert a line from an Omega-format ASCII file into a SnglBurst object.
    """
 
    if isinstance(line, str):
        dat = map(float, _delim.split(line.rstrip()))
    else:
        dat = map(float, line)

    # map to known formats
    if virgo:
        start, stop, peak, freq, band, cln, cle, snr =  dat
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        peak     = LIGOTimeGPS(peak)
        av_freq  = freq
        av_band  = band
        err_freq = 0
        clusters = False
    elif len(dat)==11:
        peak, freq, duration, band, amplitude, cls, cle, cln, av_freq,\
            av_band, err_freq =  dat
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        peak     = LIGOTimeGPS(peak)
        snr      = (2*amplitude)**(1/2)
        clusters = True
    elif len(dat)==8:
        peak, freq, duration, band, amplitude, cls, cle, cln =  dat
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        peak     = LIGOTimeGPS(peak)
        av_freq  = freq
        av_band  = band
        err_freq = 0
        snr      = (2*amplitude)**(1/2)
        clusters = True
    elif len(dat)==5:
        peak, freq, duration, band, amplitude =  dat
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        av_freq  = freq
        av_band  = band
        err_freq = 0
        snr      = (2*amplitude)**(1/2)
        clusters = False
    elif len(dat)==4:
        peak, av_freq, amplitude, chisq = dat
        snr      = (amplitude)**(1/2)
        peak     = LIGOTimeGPS(peak)
    else:
        raise ValueError("Wrong number of columns in ASCII line. Cannot read.")
    
    # set object
    t = lsctables.SnglBurst()

    # set columns that are same for all triggers
    if 'ifo' in columns: t.ifo = ifo
    if 'channel' in columns: t.channel = channel
    if 'search' in columns: t.search = 'omega'

    # set times
    if 'start_time' in columns:     t.start_time    = start.gpsSeconds
    if 'start_time_ns' in columns:  t.start_time_ns = start.gpsNanoSeconds
    if 'peak_time' in columns:      t.peak_time     = peak.gpsSeconds
    if 'peak_time_ns' in columns:   t.peak_time_ns  = peak.gpsNanoSeconds
    if 'stop_time' in columns:      t.stop_time     = stop.gpsSeconds
    if 'stop_time_ns' in columns:   t.stop_time_ns  = stop.gpsNanoSeconds
 
    # set ms times
    if 'ms_start_time' in columns:     t.ms_start_time    = start.gpsSeconds
    if 'ms_start_time_ns' in columns:  t.ms_start_time_ns = start.gpsNanoSeconds
    if 'ms_stop_time' in columns:      t.ms_stop_time     = stop.gpsSeconds
    if 'ms_stop_time_ns' in columns:   t.ms_stop_time_ns  = stop.gpsNanoSeconds

    # set other times
    if 'duration' in columns:     t.duration    = duration
    if 'ms_duration' in columns:  t.ms_duration = duration

    # set frequencies
    if 'central_freq' in columns:         t.central_freq         = freq 
    if 'peak_frequency' in columns:       t.peak_frequency       = av_freq 
    if 'peak_frequency_eror' in columns:  t.peak_frequency_error = err_freq
    if 'bandwidth' in columns:            t.bandwidth            = av_band
    if 'flow' in columns:                 t.flow                 = freq-band/2
    if 'fhigh' in columns:                t.fhigh                = freq+band/2

    # set ms frequencies
    if 'ms_bandwidth' in columns: t.ms_bandwidth = band
    if 'ms_flow' in columns:      t.ms_flow      = freq-band/2
    if 'ms_fhigh' in columns:     t.ms_fhigh     = freq+band/2

    # set amplitude and snr
    if 'snr' in columns:        t.snr        = snr
    if 'ms_snr' in columns:     t.ms_snr     = snr
    if 'amplitude' in columns:  t.amplitude  = amplitude
     
    # set other params
    if 'cluster_size' in columns or 'param_one_value' in columns:
        t.param_one_name      = 'cluster_size'
        if clusters:
            t.param_one_value = cls
        else:
            t.param_one_value = numpy.NaN
    if 'cluster_norm_energy' in columns or 'param_two_value' in columns:
        t.param_two_name      = 'cluster_norm_energy'
        if clusters:
            t.param_two_value = cle
        else:
            t.param_two_value = numpy.NaN
    if 'cluster_number' in columns or 'param_three_value' in columns:
        t.param_three_name      = 'cluster_number'
        if clusters:
            t.param_three_value = cln
        else:
            t.param_three_value = numpy.NaN

    return t

# =============================================================================
# read triggers from file
# =============================================================================

def fromfile(fobj, start=None, end=None, ifo=None, channel=None,\
                                    columns=None, virgo=False):

    """
    Load triggers from an Omega format text file into a SnglBurstTable object.
    Use start and end to restrict the returned triggers, and give ifo and
    channel to fill those columns in the table.

    If columns is given as a list, only those columns in the table will be
    filled. This is advisable to speed up future operations on this table.

    Arguments :

        fname : file or str
            file object or filename path to read with numpy.loadtext

    Keyword arguments :

        start : float
            minimum peak time for returned triggers
        end : float
            maximum peak time for returned triggers
        ifo : str
            name of IFO to fill in table
        channel : str
            name of channel to fill in table
        columns : iterable
            list of columnnames to populate in table
        virgo : [ True | False ]
            fobj written in Virgo OmegaOnline format 
    """

    # set columns
    if columns==None:
        columns = lsctables.SnglBurst.__slots__
        usercolumns = False
    else:
        usercolumns = True
    if start or end:
        if not start:
            start = -numpy.inf
        if not end:
            end     = numpy.inf
        span = segments.segment(start, end)
        if 'peak_time' not in columns: columns.append('peak_time')
        if 'peak_time_ns' not in columns: columns.append('peak_time_ns')
        check_time = True
    else:
        check_time = False

    # record amplitude if recording SNR
    if 'snr' in columns and not 'amplitude' in columns:
        columns.append('amplitude')

    #
    # generate table
    #

    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    append = out.append

    # remove unused names and types
    if usercolumns:
        # map types
        if isinstance(columns[0], str):
            columns = map(str.lower, columns)
        if isinstance(columns[0], unicode):
            columns = map(unicode.lower, columns)
        # remove unused
        for c in out.columnnames:
            if c.lower() not in columns:
                idx = out.columnnames.index(c)
                out.columnnames.pop(idx)
                out.columntypes.pop(idx)

    #
    # read file
    #

    # force filename not file object
    if hasattr(fobj, 'readline'):
        fh = fobj
    else:
        fh = open(fobj, 'r')

    # read file and generate triggers
    for i,line in enumerate(fh):
        if _comment.match(line): continue
        t = trigger(line, columns=columns, virgo=virgo, channel=channel, ifo=ifo)
        if not check_time or (check_time and float(t.get_peak()) in span):
            append(t)

    # close file if we opened it
    if not hasattr(fobj, 'readline'):
        fh.close()

    return out

# =============================================================================
# Return frequency series from omega triggers
# =============================================================================

def tofrequencyseries(bursttable, fcol='peak_frequency', pcol=None,\
                      name="", epoch=LIGOTimeGPS(), deltaF=0, f0=0,\
                      unit=lalStrainUnit):
    """
    Returns a numpy.array and REAL8FrequencySeries built from these
    OmegaSpectrum triggers. The array holds the discrete frequencies at
    which the sectrum was calculated and the series holds the data and
    associated metadata.

    If pcol is not given, the series data is the square of the SNR of each
    'trigger'.
    """

    freq = bursttable.getColumnByName('peak_frequency') 
    if pcol:
        data = bursttable.getColumnByName('pcol')
    else:
        data = bursttable.getColumnByName('snr')**2
    freq,data = list(map(numpy.asarray, zip(*sorted(zip(freq, data)))))

    if int(epoch) == 0 and len(bursttable) != 0:
        epoch    = LIGOTimeGPS(float(bursttable[0].get_time()))
    if deltaF == 0 and len(bursttable) > 1:
        deltaF = freq[1]-freq[0] 
    if f0 == 0 and len(bursttable) != 0:
        f0 = freq[0]

    series = seriesutils.fromarray(data, name=name, epoch=epoch, deltaT=deltaF,\
                                   f0=f0, unit=unit, frequencyseries=True)
    return freq,series

# =============================================================================
# Get LAL cache of omega files from their expected location
# =============================================================================

def get_cache(start, end, ifo, channel, mask='DOWNSELECT', checkfilesexist=False,\
              **kwargs):
    """
    Returns a glue.lal.Cache contatining CacheEntires for all omega online
    trigger files between the given start and end time for the given ifo.
    """
    cache = Cache()
    
    # verify host
    host = { 'G1':'atlas', 'H1':'ligo-wa', 'H2':'ligo-wa', 'L1':'ligo-la'}
    if (not kwargs.has_key('directory') and not re.search(host[ifo],getfqdn())):
        sys.stderr.write("warning: Omega online files are not available for "+\
                         "IFO=%s on this host." % ifo)
        sys.stderr.flush()
        return cache

    span = segments.segment(start,end)
    if ifo == 'G1':
        if channel:
            kwargs.setdefault('directory', '/home/omega/online/%s/segments' % channel.replace(':','_'))
        else:
            kwargs.setdefault('directory', '/home/omega/online/G1/segments')
        kwargs.setdefault('epoch', 0)
    else:
        kwargs.setdefault('directory',\
                            '/home/omega/online/%s/archive/S6/segments' % ifo)
        kwargs.setdefault('epoch', 931211808)
    kwargs.setdefault('duration', 64)
    kwargs.setdefault('overlap', 8)

    # optimise
    append       = cache.append
    splitext     = os.path.splitext
    isfile   = os.path.isfile
    intersects   = span.intersects
    segment      = segments.segment
    from_T050017 = CacheEntry.from_T050017
    basedir      = kwargs['directory']
    basetime     = kwargs['epoch']
    triglength   = kwargs['duration']
    overlap      = kwargs['overlap']

    # get times
    start_time = int(start-numpy.mod(start-basetime,triglength-overlap))
    t = start_time

    # loop over time segments constructing file paths and appending to the cache
    while t<end:
        if ifo == 'G1':
            trigfile = '%s/%.5d/%.10d-%10.d/%s-OMEGA_TRIGGERS_%s-%.10d-%d.txt'\
                % (basedir, t/100000, t, t+triglength, ifo, mask, t, triglength)
        else:
            trigfile = '%s/%.10d-%10.d/%s-OMEGA_TRIGGERS_%s-%.10d-%d.txt'\
                % (basedir, t, t+triglength, ifo, mask, t, triglength)
        if intersects(segment(t, t+triglength))\
        and (not checkfilesexist or isfile(trigfile)):
            append(from_T050017(trigfile))
        t+=triglength-overlap

    cache.sort(key=lambda e: e.path)

    return cache

# =============================================================================
# Read files
# =============================================================================

def fromfiles(filelist, start=None, end=None, columns=None, verbose=False,\
              virgo=False, channel=None):
    """
    Read omega triggers from a list of ASCII filepaths.
    """
    # set up counter
    if verbose:
        sys.stdout.write("Extracting Omega triggers from %d files...     "\
                         % len(filelist))
        sys.stdout.flush()
        delete = '\b\b\b'
        num = len(filelist)/100

    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    extend = out.extend

    for i,fp in enumerate(filelist):
        with open(fp, "r") as f:
            extend(fromfile(f, start=start, end=end, columns=columns,\
                            virgo=virgo, channel=channel))
        if verbose:
            progress = int((i+1)/num)
            sys.stdout.write('%s%.2d%%' % (delete, progress))
            sys.stdout.flush()

    if verbose:
        sys.stdout.write("\n")

    return out

def fromlalcache(cache, start=None, end=None, columns=None, verbose=False,\
                 virgo=False,channel=None):
     """
     Read omega triggers from a glue.lal.Cache list of ASCII CacheEntries.
     """
     return fromfiles(cache.pfnlist(), start=start, end=end, columns=columns,\
                      verbose=verbose, virgo=virgo, channel=channel)
