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

## \addtogroup laldetchar_py_triggers_cwb
"""Read Coherent WaveBurst events from `ROOT` files
"""
#
# ### Synopsis ###
#
# ~~~
# from laldetchar.triggers import cwb
# ~~~
# \author Duncan Macleod (<duncan.macleod@ligo.org>)

from __future__ import division

import re

from ROOT import TChain

from glue.ligolw import lsctables

from lal import LIGOTimeGPS
from laldetchar import git_version

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

CWB_DETECTOR_INDEX = ["L1", "H1", "H2", "G1", "T1", "V1", "A1"]
CWB_MULTI_COLUMNS = ['process_id', 'event_id', 'ifos',
                     'start_time', 'start_time_ns', 'duration',
                     'peak_time', 'peak_time_ns', 'central_freq',
                     'bandwidth', 'snr', 'confidence']
CWB_SNGL_COLUMNS = ['process_id', 'event_id', 'search', 'ifo', 'peak_time',
                    'peak_time_ns', 'start_time', 'start_time_ns',
                    'duration', 'time_lag', 'peak_frequency', 'bandwidth',
                    'central_freq', 'snr', 'confidence', 'hrss', 'tfvolume']

# open doxygen
## \addtogroup laldetchar_py_triggers_cwb
#@{


def root_multi_trigger(root_event, columns=CWB_MULTI_COLUMNS):
    """Parse a multi-detector Coherent WaveBurst `ROOT` tree entry
    into a `MultiBurst` object.

    @param root_event
        `ROOT` `TChain` object
    @param columns
        a list of valid `LIGO_LW` column names to load (defaults to all)

    @returns a `MultiBurst` built from the `ROOT` data
    """
    ifos = get_ifos(root_event)
    first_ifo_idx = CWB_DETECTOR_INDEX.index(list(ifos)[0])

    mb = lsctables.MultiBurst()
    if 'process_id' in columns:
        mb.process_id = lsctables.ProcessID(root_event.run)
    if 'event_id' in columns:
        mb.event_id = lsctables.MultiBurstTable.get_next_id()
    if 'ifos' in columns:
        mb.set_ifos(ifos)

    peak_time = LIGOTimeGPS(list(root_event.time)[first_ifo_idx])
    if 'peak_time' in columns:
        mb.peak_time = peak_time.gpsSeconds
    if 'peak_time_ns' in columns:
        mb.peak_time_ns = peak_time.gpsNanoSeconds
    start_time = LIGOTimeGPS(list(root_event.start)[first_ifo_idx])
    if 'start_time' in columns:
        mb.start_time = start_time.gpsSeconds
    if 'start_time_ns' in columns:
        mb.start_time_ns = start_time.gpsNanoSeconds
    if 'duration' in columns:
        mb.duration = float(list(root_event.duration)[first_ifo_idx])

    fmin = min(root_event.low)
    fmax = min(root_event.high)
    if 'central_freq' in columns:
        mb.central_freq = list(root_event.frequency)[0]
    if 'bandwidth' in columns:
        mb.bandwidth = fmax-fmin

    if 'snr' in columns:
        mb.snr = min(root_event.rho)
    if 'confidence' in columns:
        mb.confidence = root_event.likelihood

    return mb

def root_sngl_trigger(root_event, instrument, columns=CWB_SNGL_COLUMNS):
    """Parse a multi-detector Coherent WaveBurst `ROOT` tree entry
    into a `SnglBurst` object for the given instrument.

    @param root_event
        `ROOT` `TChain` object
    @param instrument UNDOCUMENTED
    @param columns
        a list of valid `LIGO_LW` column names to load (defaults to all)

    @returns a `SnglBurst` built from the `ROOT` data
    """
    ifo_idx = CWB_DETECTOR_INDEX.index(instrument)

    sb = lsctables.SnglBurst()
    if 'process_id' in columns:
        sb.process_id = lsctables.ProcessID(root_event.run)
    if 'event_id' in columns:
        sb.event_id = lsctables.SnglBurstTable.get_next_id()
    if 'search' in columns:
        sb.search = u'waveburst'

    if 'ifo' in columns:
        sb.ifo = instrument

    peak_time = LIGOTimeGPS(list(root_event.time)[ifo_idx])
    if 'peak_time' in columns:
        sb.peak_time = peak_time.gpsSeconds
    if 'peak_time_ns' in columns:
        sb.peak_time_ns = peak_time.gpsNanoSeconds
    start_time = LIGOTimeGPS(list(root_event.start)[ifo_idx])
    if 'start_time' in columns:
        sb.start_time = start_time.gpsSeconds
    if 'start_time_ns' in columns:
        sb.start_time_ns = start_time.gpsNanoSeconds
    stop_time = LIGOTimeGPS(list(root_event.stop)[ifo_idx])
    if 'stop_time' in columns:
        sb.stop_time = stop_time.gpsSeconds
    if 'stop_time_ns' in columns:
        sb.stop_time_ns = stop_time.gpsNanoSeconds
    if 'duration' in columns:
        sb.duration = float(stop_time-start_time)
    if 'time_lag' in columns:
        sb.time_lag = list(root_event.lag)[ifo_idx]

    flow = list(root_event.low)[ifo_idx]
    fhigh = list(root_event.high)[ifo_idx]
    if 'peak_frequency' in columns:
        sb.peak_frequency = list(root_event.frequency)[ifo_idx]
    if 'flow' in columns:
        sb.flow = flow
    if 'fhigh' in columns:
        sb.fhigh = fhigh
    if 'bandwidth' in columns:
        sb.bandwidth = list(root_event.bandwidth)[ifo_idx]
    if 'central_freq' in columns:
        sb.central_freq = flow + list(root_event.bandwidth)[ifo_idx]/2.0

    if 'snr' in columns:
        sb.snr = list(root_event.snr)[ifo_idx]**(1./2.)
    if 'confidence' in columns:
        sb.confidence = root_event.likelihood

    if 'hrss' in columns:
        sb.hrss = list(root_event.hrss)[ifo_idx]
    if 'tfvolume' in columns:
        sb.tfvolume = list(root_event.volume)[ifo_idx]
    return sb


def from_root(filename, columns=None, start=None, end=None, ifo=None,
              channel=None):
    """Read cWB triggers from a `ROOT` file

    If channel is not given, the `ROOT` file will be read into a
    MultiBurstTable containing the coherent information. If ifo is
    given, a SnglBurstTable will be returned, containing the
    single-detector information for that observatory.

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

    @returns a `MultiBurstTable` or `SnglBurstTable` representing the data
    """
    if channel and re.match("[A-Z]\d:", channel) and not ifo:
        ifo = channel[:2]

    if not columns and ifo:
        columns = CWB_SNGL_COLUMNS
    elif not columns:
        columns = CWB_MULTI_COLUMNS
    if columns:
        columns = set(columns)
        if 'peak_time' not in columns and (start or end):
            columns.add("peak_time")
            columns.add("peak_time_ns")
        if 'snr' in columns:
            columns.add('amplitude')
    if (start or end):
        start = start or segments.NegInfinity
        end = end or segments.PosInfinity
        span = segments.segment(start, end)
        check_time = True
    else:
        check_time = False

    # generate table
    if ifo:
        out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    else:
        out = lsctables.New(lsctables.MultiBurstTable, columns=columns)
    columns = out.columnnames
    append = out.append

    # read file and generate triggers
    root_tree = TChain("waveburst")
    root_tree.Add(filename)
    nevents = root_tree.GetEntries()
    for i in range(nevents):
        root_tree.GetEntry(i)
        if ifo:
            t = root_sngl_trigger(root_tree, ifo, columns=columns)
            t.channel = channel
        else:
            t = root_multi_trigger(root_tree, columns=columns)
        if not check_time or (float(t.get_peak()) in span):
            append(t)
    return out


def get_ifos(root_event):
    """Find the instrument set for this cWB `ROOT` tree event

    @param root_event
        a cWB `ROOT` tree event

    @returns a set of instruments
    """
    ifo_index = list(root_event.ifo)
    ndim = root_event.ndim
    return set([CWB_DETECTOR_INDEX[ifo_index[i]-1] for
                i in range(root_event.ndim)])


##@}
