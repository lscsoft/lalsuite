# Copyright (C) 2012  Matthew West
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
#                                    Preamble
#
# =============================================================================
#

"""

"""

import itertools
import math
from optparse import OptionParser
import re
import sys
import os
try:
    any
    all
except NameError:
    # Python < 2.5
    from glue.iterutils import any, all

from glue import iterutils
from glue import lal
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue import git_version

__author__ = "Matt West <matthew.west@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date



#
# =============================================================================
#
#                         Depopulate ligolw_xml tables 
#
# =============================================================================
#

def depopulate_sngl_inspiral(xmldoc, verbose = False):
    """
    This function takes the lists of event_ids from the sngl_inspiral and coinc_event_map 
    tables and determine the difference, such that the newlist contains only  non-coinc 
    single-ifo triggers. Then it remove these non-coinc triggers from the sngl_inspiral table.
    """
    sngls_tbl = lsctables.table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
    sngls_tbl_eid = sngls_tbl.getColumnByName("event_id")

    coinc_map_tbl = lsctables.table.get_table(xmldoc, lsctables.CoincMapTable.tableName)

    if len(coinc_map_tbl) == 0:
        if verbose:
            print >> sys.stderr, "This file lacks any coincident events. All %i single-ifo "\
                % len(sngls_tbl_eid) + "inspiral triggers have been removed."
        del sngls_tbl[:]
    else:
        coinc_map_tbl_eid = set(coinc_map_tbl.getColumnByName("event_id"))
        non_coincs = set(sngls_tbl_eid) - coinc_map_tbl_eid
        if verbose:
            print >> sys.stderr, "%i single-ifo inspiral triggers not associated with a "\
                % len(non_coincs) + "coincident event are being removed."

        coinc_sngls_tbl = xmldoc.childNodes[0].insertBefore( lsctables.New(lsctables.SnglInspiralTable), sngls_tbl)
        if verbose:
            print >> sys.stderr, "Creating mapping dictionary."
        # Create a mapping dictionary
        sngls_tbl_dict = {}
        for idx, eid in enumerate(sngls_tbl_eid):
            sngls_tbl_dict[eid] = idx
        if verbose:
            print >> sys.stderr, "Creating new table."
        # Create the new coincident sngls table
        for idx, event_id in enumerate(coinc_map_tbl_eid):
            coinc_sngls_tbl.insert( idx, sngls_tbl[sngls_tbl_dict[event_id]] )
        xmldoc.childNodes[0].removeChild(sngls_tbl)
    if verbose:
        print >> sys.stderr, "Done removing non-coincident events."


def remove_process_rows(xmldoc, process_ids, verbose = False):
    proc_tbl = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
    proc_param_tbl = table.get_table(xmldoc, lsctables.ProcessParamsTable.tableName)

    # Remove rows in process table whose process_id are found in sd_pids
    len_p_tbl = len(proc_tbl)
    for i, row in enumerate(proc_tbl[::-1]):
        if row.process_id in process_ids:
            del proc_tbl[len_p_tbl-i-1]
    if verbose:
        print >> sys.stderr, "\t%i rows in the process table have been removed" \
        %(len_p_tbl - len(proc_tbl))

    # Remove rows in process_params table whose process_id are found in sd_pids
    len_pp_tbl = len(proc_param_tbl)
    for i, row in enumerate(proc_param_tbl[::-1]):
        if row.process_id in process_ids:
            del proc_param_tbl[len_pp_tbl-i-1]
    if verbose:
        print >> sys.stderr, "\t%i rows in the process param table have been removed" \
        %(len_pp_tbl - len(proc_param_tbl))


def drop_segment_tables(xmldoc, verbose = False):
    """
    Drop the segment, segment_definer & segment_summary tables from the
    xmldoc. In addition, remove the rows in the process & process_params
    tables that have process_ids found in the segment_definer table.
    """
    seg_tbl = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
    seg_sum_tbl = table.get_table(xmldoc, lsctables.SegmentSumTable.tableName)
    seg_def_tbl = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
    # determine the unique process_ids for the segment tables
    sd_pids = set(seg_def_tbl.getColumnByName("process_id"))

    if verbose:
        print >> sys.stderr, "Depopulate process tables of segment process_ids"
    remove_process_rows(xmldoc, sd_pids, verbose=verbose)

    # remove segment, segment_definer & segment_summary tables from xmldoc
    xmldoc.childNodes[0].removeChild(seg_tbl)
    seg_tbl.unlink()
    xmldoc.childNodes[0].removeChild(seg_sum_tbl)
    seg_sum_tbl.unlink()
    xmldoc.childNodes[0].removeChild(seg_def_tbl)
    seg_def_tbl.unlink()
    if verbose:
        print >> sys.stderr, "segment, segment-definer & segment-summary tables dropped from xmldoc"


def drop_vetodef_table(xmldoc, verbose = False):
    """
    Drop the veto_definer table from the xmldoc and remove the rows in the
    process & process_params tables that have process_ids found in the
    veto_definer table.
    """
    veto_def_tbl = table.get_table(xmldoc, lsctables.VetoDefTable.tableName)
    # determine the unique process_ids for the segment tables
    vd_pids = set(veto_def_tbl.getColumnByName("process_id"))

    if verbose:
        print >> sys.stderr, "Depopulate process tables of veto_definer process_ids"
    remove_process_rows(xmldoc, vd_pids, verbose=verbose)

    # remove veto_definer table from xmldoc
    xmldoc.childNodes[0].removeChild(veto_def_tbl)
    veto_def_tbl.unlink()
    if verbose:
        print >> sys.stderr, "Veto-Definer table dropped from xmldoc"


def depopulate_experiment_tables(xmldoc, verbose = False):
    """
    Removes entries from the experiment tables that do not have events
    or durations in them. In other words, if none of the rows in the
    experiment_summary table that are assoicated with a single experiment_id
    have neither an event in them nor any duration then all of the rows in the
    experiment_summary table associated with that experiment_id are deleted,
    as well as the corresponding row in the experiment table. If, however, just
    one of the rows in the experiment_summary table associated with an experiment_id
    have at least a 1 in their nevents column or at least 1 second in the
    durations column then nothing associated with that experiment_id are deleted.
    (In other words, if just one time slide in an experiment has just 1 event or
    is just 1 second long, then none of the time slides in that experiment are deleted,
    even if all of the other time slides in that experiment have nothing in them.)
    """

    if verbose:
        print >>sys.stderr, "Depopulating the experiment tables...",

    #
    # find the experiment and experiment summary table
    #

    try:
        experiment_table = table.get_table(xmldoc, lsctables.ExperimentTable.tableName)
    except ValueError:
        # no table --> no-op
        if verbose:
            print >>sys.stderr, "Cannot find the experiment table"
        return

    try:
        experiment_summ_table = table.get_table(xmldoc, lsctables.ExperimentSummaryTable.tableName)
    except ValueError:
        # no table --> no-op
        if verbose:
            print >>sys.stderr, "Cannot find the experiment_summary table"
        return

    del_eid_indices = []
    del_esid_indices = []

    for mm, erow in enumerate(experiment_table):
        this_eid = erow.experiment_id
        es_index_list = []
        for nn, esrow in enumerate(experiment_summ_table):
            if esrow.experiment_id == this_eid and (esrow.duration or esrow.nevents):
                # something in this experiment, go on to next experiment_id
                break
            if esrow.experiment_id == this_eid:
                es_index_list.append(nn)
            if nn == len(experiment_summ_table) - 1:
                # if get to here, nothing in that experiment, mark id and indices
                # for removal
                del_eid_indices.append(mm)
                del_esid_indices += es_index_list

    # delte all experiments who's eids fall in del_eid_indices
    del_eid_indices.sort(reverse = True)
    for eid_index in del_eid_indices:
        del experiment_table[eid_index]
    # delete all experiment_summaries whose esids fall in del_esid_indices
    del_esid_indices.sort(reverse = True)
    for esid_index in del_esid_indices:
        del experiment_summ_table[esid_index]

    if verbose:
        print >> sys.stderr, "removed %i empty experiment(s) from the experiment table and %i " + \
            "associated time slides from the experiment_summary table." \
            %( len(del_eid_indices), len(del_esid_indices) )


#
# =============================================================================
#
#                   experiment and experiment_summ tables 
#
# =============================================================================
#

def get_experiment_times(xmldoc):
    """
    Use the start & end-times stored in the segment_summary table to define 
    the experiment times.  This presumes that the vetoes file has been added
    to the file being analyzed.
    """
    segment_summary_tbl = lsctables.table.get_table(xmldoc, lsctables.SegmentSumTable.tableName)
    expr_start_time = min(segment_summary_tbl.getColumnByName("start_time"))
    expr_end_time = max(segment_summary_tbl.getColumnByName("end_time"))

    return expr_start_time, expr_end_time

def populate_experiment_table(
    xmldoc,
    search_group,
    trigger_program,
    lars_id,
    instruments,
    comments = None,
    add_inst_subsets = False,
    verbose = False
):
    """
    Populate the experiment table using the given entries. If
    add_inst_subsets is set to True, will write additional
    entries for every possible sub-combination of the given instrument
    set. Returns a dictionary of experiment_ids keyed by the instrument
    set.

    @xmldoc: xmldoc to get/write table to
    @lars_id: lars_id of the experiment
    @search_group: lsc group that performed the experiment (e.g., cbc)
    @trigger_program: name of the program that performed the analysis
        (e.g., inspiral, ringdown, etc.)
    @comments: any desired comments
    @add_inst_subsets: will write an entry for every possible subset
        of @instruments
    @verbose: be verbose
    """

    if verbose:
        print >> sys.stderr, "\tPopulating the Experiment table..."

    # find the experiment table or create one if needed
    try:
        expr_table = table.get_table(xmldoc, lsctables.ExperimentTable.tableName)
    except ValueError:
        expr_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentTable))

    # determine experiment start and end times
    expr_start_time, expr_end_time = get_experiment_times(xmldoc)

    # write entry to the experiment table for the given instruments if it doesn't already exist
    experiment_ids = {}

    for nn in range(len(instruments), 1, -1):
        if not add_inst_subsets and (nn != len(instruments)):
            break
        # add every possible sub-combination of the instrument set if
        # they're not already in the table
        for sub_combo in iterutils.choices( list(instruments), nn ):
            if frozenset(sub_combo) not in experiment_ids:
                experiment_ids[frozenset(sub_combo)] = expr_table.write_new_expr_id(
                    search_group,
                    trigger_program,
                    lars_id,
                    sub_combo,
                    expr_start_time,
                    expr_end_time,
                    comments = comments
                )

    return experiment_ids

def get_experiment_type(xmldoc, time_slide_dict):
    """
    Determine the which experiment type(s) the coincident triggers in this
    file belong to.  It uses information from the inspiral files stored in
    the process params table to decide if the triggers come from playground
    time or are from an injection run.  If the time_slide_dict has more than
    one entry, then the triggers are from a slide run.
    """
    # get the param column of the process_params table
    process_params_tbl = lsctables.table.get_table(xmldoc, lsctables.ProcessParamsTable.tableName)
    pp_value = set(process_params_tbl.getColumnByName("value"))
    pp_param = set(process_params_tbl.getColumnByName("param"))

    zero_lag_dict = dict([dict_entry for dict_entry in time_slide_dict.items() if not any( dict_entry[1].values() )])

    # determine experiment type(s)
    datatypes = ['all_data', 'exclude_play', 'playground']

    usertags = set(['FULL_DATA','PLAYGROUND'])
    if len(zero_lag_dict):
        if ('--injection-file' in pp_param) and not (pp_value & usertags):
            datatypes = ['simulation']
        elif 'PLAYGROUND' in pp_value:
            datatypes = ['playground']

    if len(time_slide_dict) > len(zero_lag_dict):
        datatypes += ['slide']

    return datatypes


def populate_experiment_summ_table(
    xmldoc,
    experiment_id,
    time_slide_dict,
    veto_def_name,
    verbose = False
):
    """
    Populate the experiment_summ_table using an experiment_id, a
    veto_def_name, and a list of time_slide ids.

    @xmldoc: xmldoc to get/write table to
    @experiment_id: experiment_id to be added to the table.
    @veto_def_name: veto_def_name to be added to the table.
    @time_slide_dict: time_slide table as dictionary; used to set time_slide_id
        column and figure out whether or not is zero-lag. Can either be the result
        of lsctables.time_slide_table.as_dict or any dictionary having same format.
    """

    if verbose:
        print >> sys.stderr, "\tPopulating the Experiment Summary table..."

    # find the experiment_summary table or create one if needed
    try:
        expr_summ_table = table.get_table(xmldoc, lsctables.ExperimentSummaryTable.tableName)
    except ValueError:
        expr_summ_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentSummaryTable))

    # populate the experiment_summary table
    datatypes = get_experiment_type(xmldoc, time_slide_dict)

    for type in datatypes:
        for slide_id in time_slide_dict:
            if type == "slide" and not any( time_slide_dict[slide_id].values() ):
                continue
            if type != "slide" and any( time_slide_dict[slide_id].values() ):
                continue
            expr_summ_table.write_experiment_summ(
                experiment_id,
                slide_id,
                veto_def_name,
                type,
                sim_proc_id = None
            )
            if type != "slide":
                break

def get_on_instruments(xmldoc, trigger_program):
    process_tbl = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
    instruments = set([])
    for row in process_tbl:
        if row.program == trigger_program:
            instruments.add(row.ifos)
    return instruments


def generate_experiment_tables(xmldoc, **cmdline_opts):
    """
    Create or adds entries to the experiment table and experiment_summ
    table using instruments pulled from the search summary table and
    offsets pulled from the time_slide table.
    """

    if cmdline_opts["verbose"]:
        print >> sys.stderr, "Populating the experiment and experiment_summary tables using " + \
            "search_summary and time_slide tables..."

    # Get the instruments that were on
    instruments = get_on_instruments(xmldoc, cmdline_opts["trigger_program"])

    # find the experiment & experiment_summary table or create one if needed
    try:
        table.get_table(xmldoc, lsctables.ExperimentSummaryTable.tableName)
    except ValueError:
        xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentSummaryTable))

    try:
        table.get_table(xmldoc, lsctables.ExperimentTable.tableName)
    except ValueError:
        xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentTable))

    # Populate the experiment table
    experiment_ids = populate_experiment_table(
        xmldoc,
        cmdline_opts["search_group"],
        cmdline_opts["trigger_program"],
        cmdline_opts["lars_id"],
        instruments,
        comments = cmdline_opts["comment"],
        add_inst_subsets = True,
        verbose = cmdline_opts["verbose"]
    )

    # Get the time_slide table as dict
    time_slide_dict = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).as_dict()

    # Populate the experiment_summary table
    for instruments in experiment_ids:
        populate_experiment_summ_table(
            xmldoc,
            experiment_ids[instruments],
            time_slide_dict,
            cmdline_opts["vetoes_name"],
            verbose = cmdline_opts["verbose"]
        )


def populate_experiment_map(xmldoc, veto_def_name, verbose = False):
    from glue.pipeline import s2play as is_in_playground

    #
    # find the experiment_map table or create one if needed
    #
    if verbose:
        print >> sys.stderr, "\tMapping coinc events to experiment_summary table..."

    try:
        expr_map_table = table.get_table(xmldoc, lsctables.ExperimentMapTable.tableName)
    except ValueError:
        expr_map_table = xmldoc.childNodes[0].appendChild(lsctables.New(lsctables.ExperimentMapTable))

    #
    # find the coinc_event table
    #

    coinc_event_table = table.get_table(xmldoc, lsctables.CoincTable.tableName)

    #
    # Index the coinc_inspiral table as a dictionary
    #

    coinc_index = dict((row.coinc_event_id, row) for row in table.get_table(xmldoc, lsctables.CoincInspiralTable.tableName))

    #
    # Get the time_slide_table as dict
    #

    time_slide_dict = table.get_table(xmldoc, lsctables.TimeSlideTable.tableName).as_dict()

    #
    # find the experiment & experiment summary tables
    #

    expr_table = table.get_table(xmldoc, lsctables.ExperimentTable.tableName)
    expr_summ_table = table.get_table(xmldoc, lsctables.ExperimentSummaryTable.tableName)

    #
    # determine what experiment datatypes are in this file
    #

    datatypes = get_experiment_type(xmldoc, time_slide_dict)

    #
    # cycle through the coincs in the coinc_inspiral table
    #
    for coinc in coinc_event_table:

        #
        # get the experiment and experiment_summ_id for this coinc
        #

        for expr in expr_table:
            if expr.instruments == coinc.instruments:
                expr_id =  expr.experiment_id

        in_playground = is_in_playground( coinc_index[coinc.coinc_event_id].end_time ) 
        for type in datatypes:
            # make sure that all_data triggers also get mapped to the right
            # all_data sub-type: exclude_play or playground
            if in_playground & (type == "exclude_play"):
                continue
            elif (not in_playground) & (type == "playground"):
                continue
            # if doing zerolag and slides together, make sure the triggers are mapped correctly
            elif type == "slide" and not any( time_slide_dict[coinc.time_slide_id].values() ):
                continue
            elif type != "slide" and any( time_slide_dict[coinc.time_slide_id].values() ):
                continue
            else:
                # map the coinc to an experiment
                expr_map = lsctables.ExperimentMap()
                expr_map.coinc_event_id = coinc.coinc_event_id

                expr_map.experiment_summ_id = expr_summ_table.get_expr_summ_id(
                    expr_id,
                    coinc.time_slide_id,
                    veto_def_name,
                    type,
                    sim_proc_id = None
                )
                if not expr_map.experiment_summ_id:
                    raise ValueError, "%s experiment_summ_id could not be found with %s" \
                    %( type, ','.join([ str(expr_id), str(coinc.time_slide_id), veto_def_name ]))

                # map the experiment
                expr_map_table.append(expr_map)
                # Increment number of events in nevents column by 1
                expr_summ_table.add_nevents( expr_map.experiment_summ_id, 1 )


