# Copyright (C) 2009  Kipp Cannon, Chad Hanna
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


import sys


from glue import iterutils
from glue import segments
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw.utils import search_summary as ligolw_search_summary
from glue.ligolw.utils import segments as ligolw_segments
from pylal import SnglInspiralUtils


#
# =============================================================================
#
#                                  Functions
#
# =============================================================================
#


def get_thinca_rings_by_available_instruments(connection, program_name = "thinca"):
  """
  Return the thinca rings from the database at the given connection.  The
  rings are returned as a glue.segments.segmentlistdict indexed by the
  set of instruments that were analyzed in that ring.

  Example:

  >>> seglists = get_thinca_rings_by_available_instruments(connection)
  >>> print seglists.keys()
  [frozenset(['H1', 'L1'])]
  """
  # extract raw rings indexed by available instrument set

  xmldoc = dbtables.get_xml(connection)
  seglists = segments.segmentlistdict()
  for row in map(dbtables.table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName).row_from_cols, connection.cursor().execute("""
SELECT
  search_summary.*
FROM
  search_summary
  JOIN process ON (
    process.process_id == search_summary.process_id
  )
WHERE
  process.program == ?
  """, (program_name,))):
    available_instruments = frozenset(row.get_ifos())
    try:
      seglists[available_instruments].append(row.get_out())
    except KeyError:
      seglists[available_instruments] = [row.get_out()]
  xmldoc.unlink()

  # remove rings that are exact duplicates on the assumption that there are
  # zero-lag and time-slide thinca jobs represented in the same document

  return segments.segmentlistdict((key, segments.segmentlist(sorted(set(value)))) for key, value in seglists.items())


def get_thinca_zero_lag_segments(connection, program_name = "thinca"):
  """
  Return the thinca rings from the database at the given connection.  The
  rings are returned as a coalesced glue.segments.segmentlistdict indexed
  by instrument.

  Example:

  >>> seglists = get_thinca_zero_lag_segments(connection)
  >>> print seglists.keys()
  ['H1', 'L1']

  This function is most useful if only zero-lag segments are desired
  because it allows for more convenient manipulation of the segment lists
  using the methods in glue.segments.  If information about background
  segments or the original ring boundaries is desired the data returned by
  get_thinca_rings_by_available_instruments() is required.
  """
  # extract the raw rings indexed by instrument

  xmldoc = dbtables.get_xml(connection)
  seglists = ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, program_name)
  xmldoc.unlink()

  # remove rings that are exact duplicates on the assumption that there are
  # zero-lag and time-slide thinca jobs represented in the same document

  seglists = segments.segmentlistdict((key, segments.segmentlist(set(value))) for key, value in seglists.items())

  # coalesce the remaining segments making sure we don't loose livetime in
  # the process

  durations_before = abs(seglists)
  seglists.coalesce()
  if abs(seglists) != durations_before:
    raise ValueError, "detected overlapping thinca rings"

  # done

  return seglists


def get_veto_segments(connection, name):
  """
  Return a coalesced glue.segments.segmentlistdict object containing the
  segments of the given name extracted from the database at the given
  connection.
  """
  xmldoc = dbtables.get_xml(connection)
  seglists = ligolw_segments.segmenttable_get_by_name(xmldoc, name).coalesce()
  xmldoc.unlink()
  return seglists


def get_background_offset_vectors(connection):
  """
  Return a list of the non-zero offset vectors extracted from the database
  at the given connection.  Each offset vector is returned as a dictionary
  mapping instrument name to offset.
  """
  xmldoc = dbtables.get_xml(connection)
  offset_vectors = [offsetvector for offsetvector in dbtables.table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).as_dict().values() if any(offsetvector.values())]
  xmldoc.unlink()
  return offset_vectors


def get_thinca_livetimes(ring_sets, veto_segments, offset_vectors, verbose = False):
  # FIXME:  somebody should document this
  livetimes = {}
  for available_instruments, rings in ring_sets.items():
    for on_instruments in (combo for m in range(2, len(available_instruments) + 1) for combo in iterutils.choices(sorted(available_instruments), m)):
      if verbose:
        print >>sys.stderr, "%s/%s" % (",".join(on_instruments), ",".join(sorted(available_instruments))),
      on_instruments = frozenset(on_instruments)
      if on_instruments not in livetimes:
        livetimes[on_instruments] = [0.0] * len(offset_vectors)
      for i, livetime in enumerate(SnglInspiralUtils.compute_thinca_livetime(on_instruments, available_instruments - on_instruments, rings, veto_segments, offset_vectors)):
        livetimes[on_instruments][i] += livetime
  return livetimes

