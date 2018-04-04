# Copyright (C) 2006  Duncan A. Brown
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

from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import search_summary as ligolw_search_summary
from pylal import git_version
from pylal import SnglInspiralUtils
from pylal import snglcluster
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

lsctables.LIGOTimeGPS = LIGOTimeGPS

__author__ = "Duncan Brown <dbrown@ligo.caltech.edu>"


#
# =============================================================================
#
#                                 Preparation
#
# =============================================================================
#

def get_tables(doc):
  snglinspiraltable = table.get_table(
    doc, lsctables.SnglInspiralTable.tableName)

  input_times = None
  output_times = None
  try:
    searchsummtable = table.get_table(
      doc, lsctables.SearchSummaryTable.tableName)
    input_times = searchsummtable.get_inlist().extent()
    output_times = searchsummtable.get_outlist().extent()
  except ValueError:
    pass
    
  return input_times, output_times, snglinspiraltable


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#

def append_process(doc, **kwargs):
  process = ligolw_process.append_process(
    doc, program = "ligolw_sicluster", version = git_version.verbose_msg,
    cvs_repository = "lscsoft", cvs_entry_time = git_version.date,
    comment = kwargs["comment"])

  ligolw_process.append_process_params(doc, process, 
    [("--cluster-window", "lstring", kwargs["cluster_window"])])
  if kwargs["snr_threshold"] > 0:
    ligolw_process.append_process_params(doc, process, 
      [("--snr-threshold", "lstring", kwargs["snr_threshold"])])
  if kwargs["sort_descending_snr"]:
    ligolw_process.append_process_params(doc, process, 
      [("--sort-descending-snr", "lstring", " ")])
  if kwargs["sort_ascending_snr"]:
    ligolw_process.append_process_params(doc, process, 
      [("--sort-ascending-snr", "lstring", " ")])

  return process


#
# =============================================================================
#
#                             Clustering Algorithm
#
# =============================================================================
#

def SnglInspiralCluster(a, b):
  """
  Replace a with a cluster constructed from a and b. 
  """
  if b.snr >= a.snr:
    return b
  else:
    return a

#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#

def ligolw_sicluster(doc, **kwargs):
  # Extract segments and tables
  inseg, outseg, snglinspiraltable = get_tables(doc)

  # Add process information
  try:
    process = append_process(doc, **kwargs)
  except ValueError:
    process = None

  # Delete all triggers below threshold
  if kwargs["snr_threshold"] > 0:
    thresh = float(kwargs["snr_threshold"])
    if kwargs["verbose"]:
      print >>sys.stderr, "discarding triggers with snr < %f ..." % \
        kwargs["snr_threshold"]
    for i in range(len(snglinspiraltable) - 1, -1, -1):
      if snglinspiraltable[i].snr <= thresh:
        del snglinspiraltable[i]

  # Cluster
  snglcluster.cluster_events(
    snglinspiraltable,
    testfunc = lambda a, b: SnglInspiralUtils.CompareSnglInspiral(a, b, twindow = kwargs["cluster_window"]),
    clusterfunc = SnglInspiralCluster,
    sortfunc = SnglInspiralUtils.CompareSnglInspiralByEndTime,
    bailoutfunc = lambda a, b: SnglInspiralUtils.CompareSnglInspiral(a, b, twindow = kwargs["cluster_window"]),
    verbose = kwargs["verbose"]
  )

  # Sort by signal-to-noise ratio
  if kwargs["sort_ascending_snr"] or kwargs["sort_descending_snr"]:
    if kwargs["verbose"]:
      print >>sys.stderr, "sorting by snr ..."
    snglinspiraltable.sort(SnglInspiralUtils.CompareSnglInspiralBySnr)
    if kwargs["sort_descending_snr"]:
      snglinspiraltable.reverse()

  # Add search summary information
  if process and inseg and outseg:
    ligolw_search_summary.append_search_summary(doc, process, inseg = inseg, outseg = outseg, 
      nevents = len(snglinspiraltable))
  if process:
    ligolw_process.set_process_end_time(process)

  return doc
