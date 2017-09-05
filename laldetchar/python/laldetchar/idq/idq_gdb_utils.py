# Copyright (C) 2014 Reed Essick
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## \defgroup laldetchar_py_idq_idq_gdb_utils methods used for GraceDB interface
## \ingroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import idq_gdb_utils
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>)

"""Module with utility functions used for generating iDQ input to GraceDB.
"""

from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

## \addtogroup laldetchar_py_idq_idq_gdb_utils
# @{

# Utility functions used for generating iDQ input to GraceDB.

import glob

import numpy as np
import re as re

from laldetchar.idq import event
from laldetchar.idq import ovl

from laldetchar.idq import idq

from laldetchar.idq import idq_tables
from glue.ligolw import ligolw
from laldetchar.idq import idq_tables_dbutils
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw import table

#===================================================================================================

def get_glitch_times(glitch_xmlfiles):
    """
    Returns list of (gps,gps_ns) tuples for all events contained in the glitch_xmlfiles.
    """
    # load files in database
    connection, cursor = idq_tables_dbutils.load_xml_files_into_database(glitch_xmlfiles)

    # get table names from the database
    tablenames = dbtables.get_table_names(connection)
    if not tablenames:
        print "No tables were found in the database."
        return []

    # check if glitch table is present
    if not idq_tables.IDQGlitchTable.tableName in tablenames:
        print "No glitch table is found in database."
        print "Can not perform requested query."
        return []

    data = cursor.execute("""SELECT gps, gps_ns FROM """ + \
        idq_tables.IDQGlitchTable.tableName).fetchall()
    # close database
    connection.close()
    return data

def get_glitch_ovl_snglburst_summary_info(glitch_xmlfiles, glitch_columns, ovl_columns, snglburst_columns):
    """
    Generates summary info table for glitch events stored in glitch_xmlfiles.
    Returns list of (ifo, gps, gps_ns, rank, fap, ovl_channel, trig_type, trig_snr) tuples.
    Each tuple in the list corresponds to a glitch event.
    """
    # load files in database
    connection, cursor = idq_tables_dbutils.load_xml_files_into_database(glitch_xmlfiles)

    # get glitch gps times and ovl channels
    data = idq_tables_dbutils.get_get_glitch_ovl_sngburst_data(\
        connection, cursor, glitch_columns, ovl_columns, snglburst_columns)

    # close database
    connection.close()
    return data

def get_glitch_ovl_channels(glitch_xmlfiles):
    """
    Gets ovl channels for glitch events from glitch_xmlfiles.
    Returns list of (gps_seconds, gps_nanonsecons, ovl_channel) tuples.
    Each tuple in the list corresponds to a glitch event.
    """
    # load files in database
    connection, cursor = idq_tables_dbutils.load_xml_files_into_database(glitch_xmlfiles)

    # get glitch gps times and ovl channels
    data = idq_tables_dbutils.get_glitch_ovl_data(connection, cursor, \
        ['gps', 'gps_ns'], ['aux_channel'])
    # close database
    connection.close()
    return data

def extract_ovl_vconfigs( rank_frames, channame, traindir, start, end, metric='eff/dt' ):
    """
    returns a dictionary mapping active vconfigs to segments
    does NOT include "none" channel
    """
    vconfigs = []
    for rnkfr in rank_frames:
        trained, calib = idq.extract_timeseries_ranges( rnkfr )
        classifier = idq.extract_fap_name( rnkfr ) 

        vetolist = glob.glob( "%s/%d_%d/ovl/ovl/*vetolist.eval"%(traindir, trained[0], trained[1]) )        
        if len(vetolist) != 1:
            raise ValueError( "trouble finding a single vetolist file for : %s"%rnkfr )
        vetolist=vetolist[0]
        v = event.loadstringtable( vetolist )

        rankmap = { 0:[(None, None, None, None, 0, 0)] }

        for line in v:
            metric_exp = float(line[ovl.vD['metric_exp']])
            if metric == 'eff/dt':
                rnk = ovl.effbydt_to_rank( metric_exp )
            elif metric == 'vsig':
                rnk = ovl.vsig_to_rank( metric_exp )
            elif metric == 'useP': 
                rnk = ovl.useP_to_rank( metric_exp )
            else:
                raise ValueError("metric=%s not understood"%metric)
            if rankmap.has_key(rnk):
                rankmap[rnk].append( (line[ovl.vD['vchan']], float(line[ovl.vD['vthr']]), float(line[ovl.vD['vwin']]), metric, metric_exp, rnk ))
            else:
                rankmap[rnk] = [(line[ovl.vD['vchan']], float(line[ovl.vD['vthr']]), float(line[ovl.vD['vwin']]), metric, metric_exp, rnk )]

        for key, value in rankmap.items():
            rankmap[key] = tuple(value)

        t, ts = idq.combine_gwf( [rnkfr], [channame])
        t = t[0]
        truth = (start <= t)*(t <= end)
        t = t[truth]
        ts = ts[0][truth]
        if not len(ts):
            continue

        configs = rankmap[ts[0]]
        segStart = t[0]
        for T, TS in zip(t, ts):
            if rankmap[TS] != configs:
                vconfigs.append( (configs, [segStart, T] ) )
                segStart = T
                configs = rankmap[TS]
            else:
                pass 
        vconfigs.append( (configs, [segStart, T+t[1]-t[0]] ) )

    configs = {}
    for vconfig, seg in vconfigs:
        if configs.has_key( vconfig ):
            configs[vconfig].append( seg )
        else:
            configs[vconfig] = [ seg ]
    for key, value in configs.items():
        value = event.andsegments( [event.fixsegments( value ), [[start,end]] ] )
        if event.livetime( value ):
            configs[key] = event.fixsegments( value )
        else:
            raise ValueError("somehow picked up a config with zero livetime...")

    return vconfigs, configs, {"vchan":0, "vthr":1, "vwin":2, "metric":3, "metric_exp":4, "rank":5}

##@}
