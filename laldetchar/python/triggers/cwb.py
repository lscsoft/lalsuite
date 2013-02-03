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
 
"""Read and manipulate Coherent WaveBurst events from ROOT files
"""

from __future__ import division

import re

from ROOT import TChain

from glue.ligolw import ilwd,lsctables,table,utils

from lal import LIGOTimeGPS
from laldetchar import git_version

__author__ = "Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

CWB_DETECTOR_INDEX = ["L1", "H1", "H2", "G1", "T1", "V1", "A1"]


def get_ifos_from_index(cwb_event):
    """@returns the instrument set for this cWB ROOT tree event
    """
    ifo_index = list(cwb_event.ifo)
    ndim = cwb_event.ndim
    return [CWB_DETECTOR_INDEX[ifo_index[i]-1] for i in range(cwb_event.ndim)]


def get_multi_burst(cwb_event,
                    columns=lsctables.MultiBurstTable.validcolumns.keys()):
    """@returns a LIGOLw MultiBurst event with attributes seeded from
    the given cWB ROOT tree event
    """
    ndim = cwb_event.ndim
    mb = lsctables.MultiBurst()

    ifos = get_ifos_from_index(cwb_event)
    first_ifo_idx = min(sorted(range(len(ifos)), key=ifos.__getitem__))

    if "ifos" in columns:
        mb.set_ifos(ifos)

    if "process_id" in columns:
        mb.process_id = ilwd.get_ilwdchar_class("multi_burst", "process_id")(
                            cwb_event.run)

    if "snr" in columns:
        mb.snr = min(cwb_event.rho) 

    fmin = min(cwb_event.low)
    fmax = min(cwb_event.high)
    if "peak_frequency" in columns:
        mb.peak_frequency = list(cwb_event.frequency)[0]
    if "bandwidth" in columns:
        mb.bandwidth = fmax-fmin
    if "central_freq" in columns:
        mb.central_freq = fmin + mb.bandwidth/2.0

    peak_time = LIGOTimeGPS(list(cwb_event.time)[first_ifo_idx])
    if "peak_time" in columns:
        mb.peak_time = peak_time.gpsSeconds
    if "peak_time_ns" in columns:
        mb.peak_time_ns = peak_time.gpsNanoSeconds
    start_time = LIGOTimeGPS(list(cwb_event.start)[first_ifo_idx])
    if "start_time" in columns:
        mb.start_time = start_time.gpsSeconds
    if "start_time_ns" in columns:
        mb.start_time_ns = start_time.gpsNanoSeconds
    stop_time = LIGOTimeGPS(list(cwb_event.stop)[first_ifo_idx])
    #mb.stop_time = stop_time.gpsSeconds
    #mb.stop_time_ns = stop_time.gpsNanoSeconds
    if "duration" in columns:
        mb.duration = float(stop_time-start_time)

    return mb


def get_sngl_burst(cwb_event, ifo,
                   columns=lsctables.MultiBurstTable.validcolumns.keys()):
    """@returns a LIGOLw SnglBurst event with attributes seeded from
    the given cWB ROOT tree event
    """
    ndim = cwb_event.ndim
    sb = lsctables.SnglBurst()

    ifo_idx = CWB_DETECTOR_INDEX.index(ifo)
    if ifo in columns:
        sb.ifo = ifo

    if "process_id" in columns:
        sb.process_id = ilwd.get_ilwdchar_class("multi_burst", "process_id")(
                            cwb_event.run)
    if "search" in columns:
        sb.search = u"waveburst"

    if "snr" in columns:
        sb.snr = list(cwb_event.rho)[ifo_idx]
    if "confidence" in columns:
        sb.confidence = cwb_event.likelihood
    
    if "peak_frequency" in columns:
        sb.peak_frequency = list(cwb_event.frequency)[ifo_idx]
    if "flow" in columns:
        sb.flow = list(cwb_event.low)[ifo_idx]
    if "fhigh" in columns:
        sb.fhigh = list(cwb_event.high)[ifo_idx]
    if "bandwidth" in columns:
        sb.bandwidth = list(cwb_event.bandwidth)[ifo_idx]
    if "central_freq" in columns:
        sb.central_freq = sb.flow + sb.bandwidth/2.0

    peak_time = LIGOTimeGPS(list(cwb_event.time)[ifo_idx])
    if "peak_time" in columns:
        sb.peak_time = peak_time.gpsSeconds
    if "peak_time_ns" in columns:
        sb.peak_time_ns = peak_time.gpsNanoSeconds
    start_time = LIGOTimeGPS(list(cwb_event.start)[ifo_idx])
    if "start_time" in columns:
        sb.start_time = start_time.gpsSeconds
    if "start_time_ns" in columns:
        sb.start_time_ns = start_time.gpsNanoSeconds
    stop_time = LIGOTimeGPS(list(cwb_event.stop)[ifo_idx])
    if "stop_time" in columns:
        sb.stop_time = stop_time.gpsSeconds
    if "stop_time_ns" in columns:
        sb.stop_time_ns = stop_time.gpsNanoSeconds
    if "duration" in columns:
        sb.duration = float(stop_time-start_time)
    if "time_lag" in columns:
        sb.time_lag = list(cwb_event.lag)[ifo_idx]

    if "hrss" in columns:
        sb.hrss = list(cwb_event.hrss)[ifo_idx]

    return sb


def from_root_file(filename, start=None, end=None, ifo=None, channel=None,
                   columns=None):
    """Read a {Multi,Sngl}BurstTable from cWB ROOT file

    If ifo is not given, the ROOT file will be read into a MultiBurstTable
    containing the coherent information. If ifo is given, a SnglBurstTable
    will be returned, containing the single-detector information for that
    observatory.

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
        if ifo:
            columns = lsctables.SnglBurst.__slots__
        else:
            columns = lsctables.MultiBurst.__slots__
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

    # generate table
    if usercolumns:
        for c in out.columnnames:
            if c.lower() not in columns:
                idx = out.columnnames.index(c)
                out.columnnames.pop(idx)
                out.columntypes.pop(idx)

    # read file and generate triggers
    root_tree = TChain("waveburst")
    root_tree.Add(filename)
    
    if ifo:
        out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    else:
        out = lsctables.New(lsctables.MultiBurstTable, columns=columns)
    append = out.append

    nevents = root_tree.GetEntries()
    for i in range(nevents):
        root_tree.GetEntry(i)
        if ifo:
            t = get_sngl_burst(root_tree, ifo, columns=columns)
        else:
            t = get_multi_burst(root_tree, columns=columns)
        if check_time and (float(t.get_peak()) not in span):
            continue
        if channel:
            t.channel = channel
        append(t)

    return out


def from_file(filename, start=None, end=None, ifo=None, channel=None,
              columns=None):
    """Read a {Multi,Sngl}BurstTable from a cWB ROOT or ASCII file

    If ifo is not given, the file will be read into a MultiBurstTable
    containing the coherent information. If ifo is given, a SnglBurstTable
    will be returned, containing the single-detector information for that
    observatory.

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
         raise RuntimeError("Unable to parse cWB-format ASCII files, "
                            "please code up the 'from_ascii_file' for "
                            "cWB and patch this module...")


def from_files(filelist, start=None, end=None, ifo=None, channel=None,
               columns=None, verbose=False):
    """Read a BurstTable from a list of cWB-format ROOT or ASCII files

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

    if ifo:
        out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    else:
        out = lsctables.New(lsctables.MultiBurstTable, columns=columns)
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
                   columns=None, verbose=False):
    """Read a BurstTable from a Cache of cWB ROOT files
    
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
