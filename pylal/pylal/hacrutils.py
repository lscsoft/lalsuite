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
Auxiliary functions for running and postprocessing HACR pipeline code.
"""

from __future__ import division

import MySQLdb
import sys
import re
import datetime

from lal import LIGOTimeGPS, GPSToUTC

from pylal import git_version

from glue.ligolw import lsctables

__author__  = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = git_version.id
__date__    = git_version.date

# =============================================================================
# convert ascii line to HACR trigger
# =============================================================================

_comment = re.compile('[#%]')
_delim   = re.compile('[\t\,\s]+')

def trigger(line, columns=lsctables.SnglBurst.__slots__):
    """
    Convert a line from an Omega-format ASCII file into a SnglBurst object.
    """

    if isinstance(line, str):
        dat = map(float, _delim.split(line.rstrip()))
    else:
        dat = map(float, line)

    if len(dat)==8:
        peak, peak_offset, freq, band, duration, n_pix, snr, totPower = dat
        peak = LIGOTimeGPS(peak+peak_offset)
        start    = LIGOTimeGPS(peak-duration/2)
        stop     = LIGOTimeGPS(peak+duration/2)
        av_freq  = freq
        av_band  = band
        err_freq = 0
    else:
        raise ValueError("Wrong number of columns in ASCII line. Cannot read.")

    # set object
    t = lsctables.SnglBurst()

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
    if 'snr' in columns:       t.snr       = snr
    if 'ms_snr' in columns:    t.ms_snr    = snr
    #if 'amplitude' in columns: t.amplitude = totPower**(1/2)

    # set other params
    if 'numPixels' in columns or 'param_two_value' in columns:
        t.param_two_name   = 'numPixels'
        t.param_two_value  = n_pix
    if 'totPower' in columns or 'param_three_value' in columns:
        t.param_three_name = 'cluster_number'
        t.param_two_name   = totPower

    return t

# =============================================================================
# Find databases
# =============================================================================

def find_databases(host, user="reader", passwd="readonly", match=None):
    """
    Query for all databases on the given host
    """
    connection = MySQLdb.connect(host=host, user=user, passwd=passwd)
    cursor = connection.cursor()
    cursor.execute("SHOW DATABASES")
    out = cursor.fetchall()
    connection.close()
    databases = [db[0] for db in out]
    if isinstance(match, str):
        match = re.compile(match) 
    if match:
        databases = [db for db in databases if match.search(db)]
    return databases

def find_geo_databases(start_time, end_time, host, user="reader",\
                       passwd="readonly"):
    """
    Query for all databases named geoYYYYMM whose year and month match the given
    [start_time, end_time) interval.
    """
    # find all relevant databases
    match = "\Ageo\d\d\d\d\d\d\Z"
    alldbs = find_databases(host, user=user, passwd=passwd, match=match)

    # format dates
    start_date = datetime.date(*GPSToUTC(int(start_time))[:3])
    end_date   = datetime.date(*GPSToUTC(int(end_time))[:3])

    # pick out those databases whose dates look correct
    databases = []
    d  = start_date
    dt = datetime.timedelta(days=28)
    while d<=end_date:
        db = 'geo%s' % d.strftime("%Y%m")
        if db not in databases:
            databases.append(db)
        d += dt

    return databases

# =============================================================================
# Find databases
# =============================================================================

def find_channels(host, database=None, user="reader", passwd="readonly",\
                  match=None):
    """
    Query for all channels in the given HACR database on the given host
    """
    if not database:
        database = find_databases(host, user=user, passwd=passwd)[-1]
    connection = MySQLdb.connect(host=host, user=user, passwd=passwd,\
                                 db=database)
    cursor = connection.cursor()
    query = "select channel from job where monitorName = 'chacr'"
    cursor.execute(query)
    out = cursor.fetchall()
    connection.close()
    channels = [ch[0] for ch in out]
    if isinstance(match, str): 
        match = re.compile(match)
    if match:
        channels = [ch for ch in channels if match.search(ch)]
    return channels

# =============================================================================
# Find databases
# =============================================================================

def get_triggers(start_time, end_time, channel, host, columns=None, pid=None,\
                 user="reader", passwd="readonly", verbose=False):
    """
    Query the given HACR MySQL host for HACR triggers in a given
    [start_time, stop_time) semi-open interval. Returns a LIGOLw SnglBurstTable.
    """
    ifo = channel[0:2]
    if ifo != "G1":
        raise NotImplementedError("Access to non-GEO databases hasn't been "+\
                                  "implemented yet.")

    #
    # generate table
    #

    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    append = out.append

    # remove unused names and types
    if columns:
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
    else:
        columns = lsctables.SnglBurst.__slots__
    
    # get the relevant databases
    databases  = find_geo_databases(start_time, end_time, host, user=user,\
                                    passwd=passwd)

    #
    # connect
    # 

    connection = MySQLdb.connect(host=host, user=user, passwd=passwd)
    cursor     = connection.cursor()

    #
    # find triggers
    #

    for db in databases:
        cursor.execute("use %s" % db)

        # verify channel in database
        query = "select process_id, channel from job where monitorName = "+\
                "'chacr' and channel = '%s'" % channel
        numchan = int(cursor.execute(query))
        result  = cursor.fetchall()
        if numchan == 0:
            raise AttributeError("%s not found in database %s.\n"\
                                 % (channel, db))
        elif numchan > 1:
            sys.stderr.write("warning: Multiple process_ids found for "+\
                             "channel %s in monitor chacr. Triggers from all "\
                             "processes will be amalgamated.\n" % channel)

        # get process ids of relevance
        process_ids = [int(process[0]) for process in result\
                       if pid==None or int(process[0]) == pid]

        # loop over process ids
        for pid in process_ids:
            query = "select gps_start, gps_offset, freq_central, bandwidth, "\
                    "duration, num_pixels, snr, totPower from mhacr where "\
                    "process_id = %s AND gps_start >= %s AND gps_start < %s "\
                    "order by gps_start asc" % (pid, start_time, end_time)
            ntrigs = cursor.execute(query)
            result = cursor.fetchall()
            if ntrigs:
                for row in result:
                    append(trigger(row, columns=columns))

    connection.close()
    return out
