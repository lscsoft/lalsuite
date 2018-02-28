# Copyright (C) 2012  Duncan M. Macleod, Kipp Cannon
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


"""
Inspiral injection identification library. 

Contains code providing the capacity to search a list of multi_inspiral
candidates for events matching entries in a sim_inspiral list of
software injections, recording the matches as inspiral <--> injection
coincidences using the standard coincidence infrastructure. 
"""

import bisect
import sys

from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw.utils import process as ligolw_process

from pylal import git_version, ligolw_thinca, llwapp, \
                  SimInspiralUtils, MultiInspiralUtils
from lalburst import timeslides as ligolw_tisi
from pylal.xlal import tools
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

lsctables.CoincMapTable.RowType = lsctables.CoincMap = tools.CoincMap

__author__ = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__credits__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Speed Hacks
#
# =============================================================================
#


lsctables.LIGOTimeGPS = LIGOTimeGPS


def multi_inspiral___cmp__(self, other):
    # compare self's end time to the LIGOTimeGPS instance other
    return (cmp(self.end_time, other.seconds) or
            cmp(self.end_time_ns, other.nanoseconds))

lsctables.MultiInspiral.__cmp__ = multi_inspiral___cmp__


#
# =============================================================================
#
#                              Document Interface
#
# =============================================================================
#


MultiInspiralSICoincDef = lsctables.CoincDef(
                    search=u"inspiral", search_coinc_type=1,
                    description=u"sim_inspiral<-->multi_inspiral coincidences")

class DocContents(object):
    """A wrapper interface to the XML document.
    """
    def __init__(self, xmldoc, sbdef, process):
        #
        # store the process row
        #

        self.process = process

        #
        # locate the multi_inspiral and sim_inspiral tables
        #

        self.multiinspiraltable = table.get_table(
                                      xmldoc,
                                      lsctables.MultiInspiralTable.tableName)
        self.siminspiraltable = table.get_table(
                                    xmldoc,
                                    lsctables.SimInspiralTable.tableName)

        #
        # get out segment lists for programs that generated
        # triggers (currently only used for time_slide vector
        # construction)
        #

        search_summary = table.get_table(
                             xmldoc, lsctables.SearchSummaryTable.tableName)
        pids = set(self.multiinspiraltable.getColumnByName("process_id"))
        seglists = search_summary.get_out_segmentlistdict(pids)\
                                 .coalesce()

        #
        # construct the zero-lag time slide needed to cover the
        # instruments listed in all the triggers, then determine
        # its ID (or create it if needed)
        #
        # FIXME:  in the future, the sim_inspiral table should
        # indicate time slide at which the injection was done
        #

        self.tisi_id = ligolw_tisi.get_time_slide_id(
                           xmldoc,
                           dict.fromkeys(seglists, 0.0),
                           create_new=process)

        #
        # get coinc_definer row for sim_inspiral <--> multi_inspiral
        # coincs; this creates a coinc_definer table if the
        # document doesn't have one
        #

        self.sb_coinc_def_id = llwapp.get_coinc_def_id(
                                   xmldoc,
                                   sbdef.search,
                                   sbdef.search_coinc_type,
                                   create_new=True,
                                   description=sbdef.description)

        #
        # get coinc table, create one if needed
        #

        try:
            self.coinctable = table.get_table(xmldoc,
                                              lsctables.CoincTable.tableName)
        except ValueError:
            self.coinctable = lsctables.New(lsctables.CoincTable)
            xmldoc.childNodes[0].appendChild(self.coinctable)
        self.coinctable.sync_next_id()

        #
        # get coinc_map table, create one if needed
        #

        try:
            self.coincmaptable = table.get_table(
                                     xmldoc,
                                     lsctables.CoincMapTable.tableName)
        except ValueError:
            self.coincmaptable = lsctables.New(lsctables.CoincMapTable)
            xmldoc.childNodes[0].appendChild(self.coincmaptable)

        #
        # sort multi_inspiral table by end time
        #

        self.multiinspiraltable.sort(lambda a, b:
                                         cmp(a.end_time, b.end_time) or
                                         cmp(a.end_time_ns, b.end_time_ns))

    def inspirals_near_endtime(self, t, dt):
        """Return a list of the inspiral events whose peak times are
        within self.inspiral_end_time_window of t.
        """
        return self.multiinspiraltable[
                   bisect.bisect_left(self.multiinspiraltable, t - dt):
                   bisect.bisect_right(self.multiinspiraltable, t + dt)]

    def sort_triggers_by_id(self):
        """Sort the multi_inspiral table's rows by ID (tidy-up document
        for output).
        """
        self.multiinspiraltable.sort(lambda a, b: cmp(a.event_id, b.event_id))

    def new_coinc(self, coinc_def_id):
        """Construct a new coinc_event row attached to the given
        process, and belonging to the set of coincidences defined
        by the given coinc_def_id.
        """
        coinc = lsctables.Coinc()
        coinc.process_id = self.process.process_id
        coinc.coinc_def_id = coinc_def_id
        coinc.coinc_event_id = self.coinctable.get_next_id()
        coinc.time_slide_id = self.tisi_id
        coinc.set_instruments(None)
        coinc.nevents = 0
        coinc.likelihood = None
        self.coinctable.append(coinc)
        return coinc


#
# =============================================================================
#
#                           Add Process Information
#
# =============================================================================
#

process_program_name = "ligolw_inspinjfind"

def append_process(xmldoc, match_algorithm, time_window, loudest_by, comment):
    """Convenience wrapper for adding process metadata to the document.
    """
    process = ligolw_process.append_process(xmldoc,
                                    program=process_program_name,
                                    version=__version__,
                                    cvs_repository=u"lscsoft",
                                    cvs_entry_time=__date__, comment=comment)
    params = [(u"--match-algorithm", u"lstring", match_algorithm),
              (u"--time-window", u"lstring", time_window)]
    if loudest_by:
        params.append((u"--loudest-by", u"lstring", loudest_by))
    ligolw_process.append_process_params(xmldoc, process, params)
    return process


#
# =============================================================================
#
#                 Build sim_inspiral <--> multi_inspiral Coincidences
#
# =============================================================================
#


def find_multi_inspiral_matches(contents, sim, time_window, loudest_by=None):
    """Scan the inspiral table for triggers matching sim.
    """
    events = [inspiral for inspiral in
              contents.inspirals_near_endtime(sim.get_end(), time_window)]
    if loudest_by and len(events):
        if hasattr(lsctables.MultiInspiral, "get_%s" % loudest_by):
            rank = getattr(lsctables.MultiInspiral, "get_%s" % loudest_by)
        else:
            rank = lambda x: getattr(x, loudest_by)
        return sorted(events, key=lambda mi: rank(mi))[-1:]
    else:
        return events


def add_sim_inspiral_coinc(contents, sim, inspirals):
    """Create a coinc_event in the coinc table, and add arcs in the
    coinc_event_map table linking the sim_inspiral row and the list of
    multi_inspiral rows to the new coinc_event row.
    """
    coinc = contents.new_coinc(contents.sb_coinc_def_id)
    if inspirals:
        coinc.set_instruments(
            lsctables.instrument_set_from_ifos(inspirals[0].ifos))
    coinc.nevents = len(inspirals)
    coincmap = lsctables.CoincMap()
    coincmap.coinc_event_id = coinc.coinc_event_id
    coincmap.table_name = sim.simulation_id.table_name
    coincmap.event_id = sim.simulation_id
    contents.coincmaptable.append(coincmap)
    for event in inspirals:
        coincmap = lsctables.CoincMap()
        coincmap.coinc_event_id = coinc.coinc_event_id
        coincmap.table_name = event.event_id.table_name
        coincmap.event_id = event.event_id
        contents.coincmaptable.append(coincmap)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_miinjfind(xmldoc, process, search, time_window, loudest_by=None,
                     verbose=False):
    """Parse an XML document and find coincidences between entries
    in sim_inspiral and multi_inspiral tables.
    """
    if verbose:
        sys.stderr.write("indexing ...\n")

    sbdef = {"inspiral": MultiInspiralSICoincDef}[search]

    contents = DocContents(xmldoc=xmldoc, sbdef=sbdef, process=process)
    N = len(contents.siminspiraltable)

    #
    # Find sim_inspiral <--> multi_inspiral coincidences.
    #

    if verbose:
        sys.stderr.write("constructing %s:\n" % sbdef.description)
    for n, sim in enumerate(contents.siminspiraltable):
        if verbose:
            sys.stderr.write("\t%.1f%%\r" % (100.0 * n / N))
            sys.stderr.flush()
        inspirals = find_multi_inspiral_matches(contents, sim, time_window,
                                                loudest_by=loudest_by)
        if inspirals:
            add_sim_inspiral_coinc(contents, sim, inspirals)
    if verbose:
        sys.stderr.write("\t100.0%\n")

    #
    # Restore the original event order.
    #

    if verbose:
        sys.stderr.write("finishing ...\n")
    contents.sort_triggers_by_id()

    #
    # Done.
    #

    return xmldoc
