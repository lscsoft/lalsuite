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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# ==============================================================================
# Preamble
# ==============================================================================

import os, sys, re, operator, math
from StringIO import StringIO
from glue import segments
from glue.ligolw import ligolw, lsctables, table, utils
from glue.ligolw.utils import segments as ligolw_segments
from glue.segmentdb import query_engine, segmentdb_utils
from glue.ligolw.utils import process as ligolw_process
from pylal.dq.dqTriggerUtils import def_get_time

from scipy.stats import poisson

LIGOTimeGPS = lsctables.LIGOTimeGPS

# Some boilerplate to make segmentlists picklable
import copy_reg
copy_reg.pickle(type(segments.segment(0, 1)), \
                lambda x:(segments.segment, (x[0], x[1])))
copy_reg.pickle(type(segments.segmentlist([])), 
                lambda x:(segments.segmentlist, ([y for y in x], )))

from glue import git_version

__author__  = "Andrew P Lundgren <andrew.lundgren@ligo.org>, Duncan Macleod <duncan.macleod@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides useful segment and veto tools for data quality investigations.
"""

# ==============================================================================
# Function to load segments from an xml file
# ==============================================================================

def fromsegmentxml(file, dict=False, id=None):

  """
    Read a glue.segments.segmentlist from the file object file containing an
    xml segment table.

    Arguments:

      file : file object
        file object for segment xml file

    Keyword Arguments:

      dict : [ True | False ]
        returns a glue.segments.segmentlistdict containing coalesced
        glue.segments.segmentlists keyed by seg_def.name for each entry in the
        contained segment_def_table. Default False
      id : int
        returns a glue.segments.segmentlist object containing only those
        segments matching the given segment_def_id integer
        
  """

  # load xmldocument and SegmentDefTable and SegmentTables
  xmldoc, digest = utils.load_fileobj(file, gz=file.name.endswith(".gz"), contenthandler = lsctables.use_in(ligolw.LIGOLWContentHandler))
  seg_def_table  = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
  seg_table      = table.get_table(xmldoc, lsctables.SegmentTable.tableName)

  if dict:
    segs = segments.segmentlistdict()
  else:
    segs = segments.segmentlist()

  seg_id = {}
  for seg_def in seg_def_table:
    seg_id[int(seg_def.segment_def_id)] = str(seg_def.name)
    if dict:
      segs[str(seg_def.name)] = segments.segmentlist()

  for seg in seg_table:
    if dict:
      segs[seg_id[int(seg.segment_def_id)]]\
          .append(segments.segment(seg.start_time, seg.end_time))
      continue
    if id and int(seg.segment_def_id)==id:
      segs.append(segments.segment(seg.start_time, seg.end_time))
      continue
    segs.append(segments.segment(seg.start_time, seg.end_time))

  if dict:
   for seg_name in seg_id.values():
     segs[seg_name] = segs[seg_name].coalesce()
  else:
    segs = segs.coalesce()

  xmldoc.unlink()

  return segs

# ==============================================================================
# Write to segment xml file
# ==============================================================================

def tosegmentxml(file, segs):

  """
    Write the glue.segments.segmentlist object segs to file object file in xml
    format with appropriate tables.
  """

  # generate empty document
  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessTable))
  xmldoc.childNodes[-1].appendChild(lsctables.New(lsctables.ProcessParamsTable))

  # append process to table
  process = ligolw_process.append_process(xmldoc,\
                                  program='pylal.dq.dqSegmentUtils',\
                                  version=__version__,\
                                  cvs_repository='lscsoft',\
                                  cvs_entry_time=__date__)

  gpssegs = segments.segmentlist()
  for seg in segs:
    gpssegs.append(segments.segment(LIGOTimeGPS(seg[0]), LIGOTimeGPS(seg[1])))

  # append segs and seg definer
  segments_tables = ligolw_segments.LigolwSegments(xmldoc)
  segments_tables.add(ligolw_segments.LigolwSegmentList(active=gpssegs))
  # finalise
  segments_tables.coalesce()
  segments_tables.optimize()
  segments_tables.finalize(process)
  ligolw_process.set_process_end_time(process)

  # write file
  with utils.SignalsTrap():
    utils.write_fileobj(xmldoc, file, gz=False)

# ==============================================================================
# Function to load segments from a csv file
# ==============================================================================

def fromsegmentcsv(csvfile):

  """
    Read a glue.segments.segmentlist object from the file object file containin
    a comma separated list of segments.
  """

  def CSVLineToSeg(line):
    tstart, tend = map(LIGOTimeGPS, line.split(','))
    return segments.segment(tstart, tend)

  segs = segments.segmentlist([CSVLineToSeg(line) for line in csvfile])

  return segs.coalesce()

# ==============================================================================
# Function to parse a segment list for CBC analysable segments
# ==============================================================================

def CBCAnalyzableSegs(seglist):

  """
    Remove any segments shorter than 2064 seconds from seglist because ihope
    won't analyze them.
  """

  return segments.segmentlist([seg for seg in seglist if abs(seg) >= 2064])

# ==============================================================================
# Function to pad a list of segments given start and end paddings
# ==============================================================================

def pad_segmentlist(seglist, start_pad, end_pad=None):

  """
    Given a veto segmentlist, start pad, and end pad, pads and coalesces the
    segments. Convention is to expand segments with positive padding, contract
    with negative. Any segments that are not big enough to be contracted
    appropriately are removed.
  """

  if not end_pad: end_pad = start_pad

  padded = lambda seg: segments.segment(seg[0]-start_pad, seg[1]+end_pad)

  seglist = segments.segmentlist([padded(seg) for seg in seglist if \
                                  abs(seg) > start_pad+end_pad])

  return seglist.coalesce()

# ==============================================================================
# Function to crop a list of segments
# ==============================================================================

def crop_segmentlist(seglist, end_chop=30):

  """
    Given a segmentlist and time to chop, removes time from the end of each
    segment (defaults to 30 seconds).
  """

  chopped = lambda seg: segments.segment(seg[0], max(seg[0], seg[1] - end_chop))
  seglist = segments.segmentlist([chopped(seg) for seg in seglist])
  return seglist.coalesce()

# =============================================================================
# Function to return segments in given gps time range
# =============================================================================

def grab_segments(start, end, flag,\
                  segment_url='https://segdb.ligo.caltech.edu',\
                  segment_summary=False):

  """
    Returns a segmentlist containing the segments during which the given flag
    was active in the given period.

    Arguments:

      start : int
        GPS start time
      end : int
        GPS end time
      flag : string
        'IFO:NAME:VERSION' format string

    Keyword arguments:

      segment_url : string
        url of segment database to query, default https://segdb.ligo.caltech.edu
      segment_summary : [ True | False ]
        also return the glue.segments.segmentlist defining the valid span of the
        returned segments
  """

  # set times
  start = int(math.floor(start))
  end   = int(math.ceil(end))

  # set query engine
  connection        = segmentdb_utils.setup_database(segment_url)
  engine            = query_engine.LdbdQueryEngine(connection)

  # format flag name
  if isinstance(flag, basestring):
    flags = flag.split(',')
  else:
    flags = flag

  segdefs = []
  for f in flags:
    spec = f.split(':')
    if len(spec) < 2 or len(spec) > 3:
      raise AttributeError, "Included segements must be of the form "+\
                            "ifo:name:version or ifo:name:*"

    ifo  = spec[0]
    name = spec[1]

    if len(spec) is 3 and spec[2] is not '*':
      version = int(spec[2])
      if version < 1:
        raise AttributeError, "Segment version numbers must be greater than zero"
    else:
      version = '*'

    # expand segment definer
    segdefs += segmentdb_utils.expand_version_number(engine, (ifo, name, version, \
                                                            start, end, 0, 0))

  # query database and return
  segs = segmentdb_utils.query_segments(engine, 'segment', segdefs)
  segs = [s.coalesce() for s in segs]
  if segment_summary:
    segsums = segmentdb_utils.query_segments(engine, 'segment_summary', segdefs)
    #segsums = reduce(operator.or_, segsums).coalesce()
    segsums = [s.coalesce() for s in segsums]
    segsummap = [segments.segmentlist() for f in flags]
    for segdef,segsum in zip(segdefs, segsums):
      try:
        fidx = flags.index(':'.join(map(str, segdef[:3])))
      except ValueError:
        fidx = flags.index(':'.join(segdef[:2]))
      segsummap[fidx].extend(segsum)
    if flag == flags[0]:
      return segs[0],segsummap[0]
    else:
      return segs,segsummap
  if flag == flags[0]:
    return segs[0]
  else:
    return segs

# ==============================================================================
# Dump flags from segment database
# ==============================================================================

def dump_flags(ifos=None, segment_url=None, match=None, unmatch=None,\
               latest=False):

  """
    Returns the list of all flags defined in the database.

    Keyword rguments:
      ifo : [ str | list ]
        list of ifos to query, or str for single ifo
      segment_url : str 
        url of segment database, defaults to contents of S6_SEGMENT_SERVER
        environment variable
      match : [ str | regular pattern ]
        regular expression to search against returned flag names, e.g, 'UPV'
      unmatch : [ str | regular pattern ]
        regular expression to negatively search against returned flag names
  """

  if isinstance(ifos, str):
    ifos = [ifos]

  # get url
  if not segment_url:
    segment_url = os.getenv('S6_SEGMENT_SERVER')

  # open connection to LDBD(W)Server
  myClient = segmentdb_utils.setup_database(segment_url)

  reply = StringIO(myClient.query(squery))
  xmldoc, digest = utils.load_fileobj(reply)
  seg_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)

  # sort table by ifo, name and version
  seg_def_table.sort(key=lambda flag: (flag.ifos[0], flag.name, \
                                       flag.version), reverse=True)

  flags = lsctables.New(type(seg_def_table))

  for row in seg_def_table:

    # test re match
    if match and not re.search(match, row.name):  continue

    # test re unmatch
    if unmatch and re.search(unmatch, row.name):  continue

    # only append latest versions of multiple flags
    flatest=True
    if latest:
      # get all flags with same ifo and name
      vflags = [f for f in flags if row.name==f.name and\
                row.get_ifos()==f.get_ifos()]
      # if later version in list, move on
      for f in vflags:
        if f.version>=row.version:
          flatest=False
          break
    if not flatest:
      continue

    # append those flags matching ifos requirement
    for ifo in ifos:
      if ifo in row.get_ifos():
        flags.append(row)
        break

  return flags

# ==============================================================================
# Calculate Poisson safety
# ==============================================================================

def poisson_safety(segs, injTable, livetime):

  """
    Return a tuple containing the number of vetoed injections, the number
    expected, and the Poisson safety probability based on the number of
    injections vetoed relative to random chance according to Poisson statistics.

    Arguments:

      segs : glue.segments.segmentlist
        list of segments to be tested
      injTable : glue.ligolw.table.Table
        table of injections
      livetime : [ float ]
        livetime of search
  """

  deadtime = segs.__abs__()

  get_time = def_get_time(injTable.tableName)
  injvetoed = len([ inj for inj in injTable if get_time(inj) in segs ])

  injexp = len(injTable)*float(deadtime) / float(livetime)

  prob = 1 - poisson.cdf(injvetoed-1, injexp)

  return injvetoed, injexp, prob

# =============================================================================
# Convert a data quality bit mask into segments
# =============================================================================

def _bits(i, n=8):
  """
    Convert integer bit mask into binary bits. Returns a list of 0s or 1s
    from 2^0 up to 2^n.

    Example:
    
    >>> _bits(295, n=8)
    [1, 1, 1, 0, 0, 1, 0, 0]
  """
  return [(0, 1)[int(i)>>j & 1] for j in xrange(n)]

def DQSegments(time, data, dq_key):

  """
    Returns a glue.segments.segmentlistdict of active segments for each bit
    in a dq_key.
  """

  segdict = segments.segmentlistdict()
  for key in dq_key:
    segdict[key] = segments.segmentlist()

  # convert DQ bits into segments
  for i in xrange(len(data)):
    binary = _bits(data[i], len(dq_key))
    seg = segments.segment(time[i], time[i]-1)
    for j, key in enumerate(dq_key):
      if binary[j] == 1:
        segdict[key].append(seg)

  segdict = segdict.coalesce()

  return segdict

