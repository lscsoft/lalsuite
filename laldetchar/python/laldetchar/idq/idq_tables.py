# Copyright (C) 2014 Reed Essick, Ruslan Vaulin
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

## \addtogroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import idq_tables
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>)

"""Module to keep definitions of various tables for iDQ pipeline."""

from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import types as ligolwtypes
from glue.ligolw import ilwd
from glue.ligolw import lsctables
from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>, Ruslan Vaulin <ruslan.vaulin@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date


## \addtogroup laldetchar_py_idq_idq_tables
# @{

# contains definitions for iDQ xml table classes and objects

# define types of coincidences within iDQ pipeline

IDQCoincDef = {
    'idq_glitch<-->sngl_burst': ('idq',0),
    'idq_glitch<-->ovl_data': ('idq', 1)
} 


def coinc_to_ovl_data(xmldoc):
    """Function returns list of (idq_glitch_object, ovl_data_object) tuples 
       where objects in the tuple are mapped to each other via coinc tables.
	"""

    # get necessary tables from xmldoc
    coinc_def_table = table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
    coinc_table = table.get_table(xmldoc, lsctables.CoincTable.tableName)
    coinc_map_table = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
    idq_glitch_table = table.get_table(xmldoc, lsctables.IDQGlitchTable.tableName)
    ovl_data_table = table.get_table(xmldoc, lsctables.OVLDataTable.tableName)
    
    # get coinc_def_ids
    #coinc_def_id = coinc_def_table.get_coinc_def_id(
    #            search = IDQCoincDef['idq_glitch<-->ovl_data'][0],
    #            search_coinc_type = IDQCoincDef['idq_glitch<-->ovl_data'][1],
    #            create_new = False,
    #            description = 'idq_glitch<-->ovl_data')
    coinc_def_ids = [row.coinc_def_id for row in coinc_def_table if row.description == 'idq_glitch<-->ovl_data']
        
        
    
    # use this id to get all coinc_event ids
    ovl_coinc_ids = [coinc.coinc_event_id for coinc in coinc_table if coinc.coinc_def_id in coinc_def_ids]

    # convert idq_glitch and ovl tables into dictionaries for a quick lookup
    glitches = dict([(glitch.event_id, glitch) for glitch in idq_glitch_table])
    ovl_data = dict([(row.event_id, row) for row in ovl_data_table])

    
    # create dictionary of connected events in coinc_event_map. 
    # We can not assume any specific order of rows in the table.
    connected_events_dict = {}
    for row in coinc_map_table:
        try: connected_events_dict[row.coinc_event_id].append(row)
        except: connected_events_dict[row.coinc_event_id] = [row]
    
    glitch_table_name =  table.StripTableName(lsctables.IDQGlitchTable.tableName)
    ovl_data_table_name = table.StripTableName(lsctables.OVLDataTable.tableName)
     
    glitch_ovl_pairs = []
    for coinc_id in ovl_coinc_ids:
        # get connectected events for this id
        connected_events = connected_events_dict[coinc_id]
        if len(connected_events) == 2:
            for event in connected_events:
                if event.table_name == glitch_table_name:
                    glitch_event = glitches[event.event_id]
                elif  event.table_name == ovl_data_table_name:
                    ovl_event = ovl_data[event.event_id]
                else:
                    print event.table_name
                    raise ValueError("Event is not the row of either " + \
                                     glitch_table_name + \
                                     " or " + ovl_data_table_name
                                     )
        else:
            raise Exception("Glitch-OVL coincidence must contain exactly 2 events. "\
                            + str(len(connected_events))+ " events are found instead."
                            )
        glitch_ovl_pairs.append((glitch_event, ovl_event))
    return glitch_ovl_pairs
            
        
    
def coinc_to_triggers(xmldoc, trigger_types):
    """ Function returns list of glitch-trigger(s) coincident events. 
        Coincident event in the list is represented by tuple (glitch_object, [trigger1, trigger2, ...]).
        trigger_types is the list of trigger type names corresponding to "search" column of the sngl_burst table
    """


    # get necessary tables from xmldoc
    coinc_def_table = table.get_table(xmldoc, lsctables.CoincDefTable.tableName)
    coinc_table = table.get_table(xmldoc, lsctables.CoincTable.tableName)
    coinc_map_table = table.get_table(xmldoc, lsctables.CoincMapTable.tableName)
    idq_glitch_table = table.get_table(xmldoc, lsctables.IDQGlitchTable.tableName)
    sngl_burst_table = table.get_table(xmldoc, lsctables.SnglBurstTable.tableName)
    
    
    # get coinc_def_id
    #coinc_def_id = coinc_def_table.get_coinc_def_id(
    #            search = IDQCoincDef['idq_glitch<-->sngl_burst'][0],
    #            search_coinc_type = IDQCoincDef['idq_glitch<-->sngl_burst'][1],
    #            create_new = False,
    #            description = 'idq_glitch<-->sngl_burst')
    coinc_def_ids = [row.coinc_def_id for row in coinc_def_table if row.description == 'idq_glitch<-->sngl_burst']
    
    # use this id to get all coinc_event ids
    trig_coinc_ids = [coinc.coinc_event_id for coinc in coinc_table if coinc.coinc_def_id in coinc_def_ids]
    

    # convert idq_glitch and sngl_burst tables into dictionaries for a quick lookup
    glitches = dict([ (glitch.event_id, glitch) for glitch in idq_glitch_table])
    triggers = dict([ (row.event_id, row) for row in sngl_burst_table if row.search in trigger_types])

    # create dictionary of connected events using coinc_event_map. 
    # We can not assume any specific order of rows in the table.
    connected_events_dict = {}
    for row in coinc_map_table:
        try: connected_events_dict[row.coinc_event_id].append(row)
        except: connected_events_dict[row.coinc_event_id] = [row]
    
    glitch_table_name = table.StripTableName(lsctables.IDQGlitchTable.tableName)
    sngl_burst_table_name = table.StripTableName(lsctables.SnglBurstTable.tableName)
    
    glitch_trig_tuples = []
    for coinc_id in trig_coinc_ids:
        # get connectected events for this id
        connected_events = connected_events_dict[coinc_id]
        connected_trigs = []
        if len(connected_events) >= 2:
            for event in connected_events:
                if event.table_name == glitch_table_name:
                    glitch_event = glitches[event.event_id]
                elif  event.table_name == sngl_burst_table_name:
                    try: connected_trigs.append(triggers[event.event_id])
                    except: # no trigger with that id, it is probably of different type.
                        pass
                else:
                    raise ValueError("Event is not the row of either " + \
                                     glitch_table_name + " or " + sngl_burst_table_name
                                     )
        else:
            raise Exception("Glitch-Triggers coincidences must contain at least 2 events. " \
                            + str(len(connected_events))+ " events are found instead."
                            )
        glitch_trig_tuples.append((glitch_event, connected_trigs))
    return glitch_trig_tuples

    
    
    
    
    
    
	
	
	
	
	
	
	


#
# =============================================================================
#
#                              glitch:table
#
# =============================================================================
#






GlitchID = ilwd.get_ilwdchar_class(u"glitch", u"event_id")

class GlitchTable(table.Table):

    tableName = 'glitch:table'
    validcolumns = {
        'event_id': 'ilwd:char',
        'ifo': 'lstring',
        'gps': 'int_4s',
        'gps_ns': 'int_4s',
        'rank': 'real_4',
        'fap': 'real_4',
        'likelihood': 'real_4',
        }
    next_id = GlitchID(0)

class Glitch(object):

    __slots__ = GlitchTable.validcolumns.keys()


GlitchTable.RowType = Glitch


#
# =============================================================================
#
#                              ovl_data:table
#
# =============================================================================
#

OVLDataID = ilwd.get_ilwdchar_class(u"ovl_data", u"event_id")

class OVLDataTable(table.Table):

    tableName = 'ovl_data:table'
    validcolumns = {
        'event_id': 'ilwd:char',
        'ifo': 'lstring',
        'aux_channel': 'lstring',
        'veto_thr': 'real_4',
        'veto_win': 'real_4',
        }
    next_id = OVLDataID(0)

class OVLData(object):

    __slots__ = OVLDataTable.validcolumns.keys()


OVLDataTable.RowType = OVLData

##@}

