#!/usr/bin/env python
#
# Copyright (C) 2013  Leo Singer
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
"""
LIGO-LW convenience functions.
"""
__author__ = "Leo Singer <leo.singer@ligo.org>"


# Python standard library imports.
import itertools
import operator

# LIGO-LW XML imports.
from pylal import ligolw_inspinjfind
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw import table as ligolw_table
from pylal import ligolw_thinca
from glue.ligolw import lsctables


def get_temlate_bank_f_low(xmldoc):
    """Determine the low frequency cutoff from a template bank file,
    whether the template bank was produced by lalapps_tmpltbank or
    lalapps_cbc_sbank. bayestar_sim_to_tmpltbank does not have a command
    line for the low-frequency cutoff; instead, it is recorded in a row
    of the search_summvars table."""
    try:
        template_bank_f_low, = ligolw_process.get_process_params(xmldoc,
            'tmpltbank', '--low-frequency-cutoff')
    except ValueError:
        try:
            template_bank_f_low, = ligolw_process.get_process_params(xmldoc,
                'lalapps_cbc_sbank', '--flow')
        except ValueError:
            try:
                search_summvars_table = ligolw_table.get_table(xmldoc,
                    lsctables.SearchSummVarsTable.tableName)
                template_bank_f_low, = (search_summvars.value
                    for search_summvars in search_summvars_table
                    if search_summvars.name == 'low-frequency cutoff')
            except ValueError:
                raise ValueError("Could not determine low-frequency cutoff")
    return template_bank_f_low


def sim_and_sngl_inspirals_for_xmldoc(xmldoc):
    """Retrieve (as a generator) all of the
    (sim_inspiral, (sngl_inspiral, sngl_inspiral, ... sngl_inspiral) tuples from
    found coincidences in a LIGO-LW XML document."""

    # Look up necessary tables.
    coinc_table = ligolw_table.get_table(xmldoc, lsctables.CoincTable.tableName)
    coinc_def_table = ligolw_table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
    coinc_map_table = ligolw_table.get_table(xmldoc, lsctables.CoincMapTable.tableName)

    # Look up coinc_def ids.
    sim_coinc_def_id = coinc_def_table.get_coinc_def_id(
        ligolw_inspinjfind.InspiralSCExactCoincDef.search,
        ligolw_inspinjfind.InspiralSCExactCoincDef.search_coinc_type,
        create_new=False)

    def events_for_coinc_event_id(coinc_event_id):
        for coinc_map in coinc_map_table:
            if coinc_map.coinc_event_id == coinc_event_id:
                for row in ligolw_table.get_table(xmldoc, coinc_map.table_name):
                    column_name = coinc_map.event_id.column_name
                    if getattr(row, column_name) == coinc_map.event_id:
                        yield coinc_map.event_id, row

    # Loop over all coinc_event <-> sim_inspiral coincs.
    for sim_coinc in coinc_table:

        # If this is not a coinc_event <-> sim_inspiral coinc, skip it.
        if sim_coinc.coinc_def_id != sim_coinc_def_id:
            continue

        # Locate the sim_inspiral and coinc events.
        sim_inspiral = None
        coinc = None
        for event_id, event in events_for_coinc_event_id(sim_coinc.coinc_event_id):
            if event_id.table_name == ligolw_table.StripTableName(lsctables.SimInspiralTable.tableName):
                if sim_inspiral is not None:
                    raise RuntimeError("Found more than one matching sim_inspiral entry")
                sim_inspiral = event
            elif event_id.table_name == ligolw_table.StripTableName(lsctables.CoincTable.tableName):
                if coinc is not None:
                    raise RuntimeError("Found more than one matching coinc entry")
                coinc = event
            else:
                raise RuntimeError("Did not expect coincidence to contain an event of type '%s'" % event_id.table_name)

        sngl_inspirals = tuple(event
            for event_id, event in events_for_coinc_event_id(coinc.coinc_event_id))

        yield sim_inspiral, sngl_inspirals


def coinc_and_sngl_inspirals_for_xmldoc(xmldoc):
    """Retrieve (as a generator) all of the
    (sngl_inspiral, sngl_inspiral, ... sngl_inspiral) tuples from coincidences
    in a LIGO-LW XML document."""

    # Look up necessary tables.
    coinc_table = ligolw_table.get_table(xmldoc, lsctables.CoincTable.tableName)
    coinc_def_table = ligolw_table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
    coinc_map_table = ligolw_table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
    sngl_inspiral_table = ligolw_table.get_table(xmldoc, lsctables.SnglInspiralTable.tableName)

    # Look up coinc_def id.
    sngl_sngl_coinc_def_ids = set(row.coinc_def_id for row in coinc_def_table
        if (row.search, row.search_coinc_type) ==
        (ligolw_thinca.InspiralCoincDef.search,
        ligolw_thinca.InspiralCoincDef.search_coinc_type))

    # Indices to speed up lookups by ID.
    key = operator.attrgetter('coinc_event_id')
    coinc_maps_by_coinc_event_id = dict((coinc_event_id, tuple(coinc_maps))
        for coinc_event_id, coinc_maps
        in itertools.groupby(sorted(coinc_map_table, key=key), key=key))
    sngl_inspirals_by_event_id = dict((sngl_inspiral.event_id, sngl_inspiral)
        for sngl_inspiral in sngl_inspiral_table)

    # Loop over all sngl_inspiral <-> sngl_inspiral coincs.
    for coinc in coinc_table:
        if coinc.coinc_def_id in sngl_sngl_coinc_def_ids:
            coinc_maps = coinc_maps_by_coinc_event_id[coinc.coinc_event_id]
            yield coinc, tuple(sngl_inspirals_by_event_id[coinc_map.event_id]
                for coinc_map in coinc_maps)


def psd_filenames_by_process_id_for_xmldoc(xmldoc):
    """Retrieve a dictionary mapping process_ids to reference PSD filenames."""
    return dict((process_param.process_id, process_param.value)
        for process_param
        in ligolw_table.get_table(xmldoc, lsctables.ProcessParamsTable.tableName)
        if process_param.param == '--reference-psd')
