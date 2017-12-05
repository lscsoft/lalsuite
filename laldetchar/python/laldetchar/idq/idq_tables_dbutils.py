
# Copyright (C) 2014 Ruslan Vaulin
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

## \defgroup laldetchar_py_idq_idq_tables_dbutils sqlite utils for iDQ Tables
## \ingroup laldetchar_py_idq
## Synopsis
# ~~~
# from laldetchar.idq import idq_tables_utils
# ~~~
# \author Ruslan Vaulin (<ruslan.vaulin@ligo.org>)

"""
Module with utility functions for manipulating idq tables using sqlite interface. 
"""

from laldetchar.idq import idq_tables
from glue.ligolw import ligolw
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw.utils import ligolw_sqlite
from glue.ligolw import dbtables
from glue import segments, segmentsUtils
from glue import lal
import sqlite3

## \addtogroup laldetchar_py_idq_idq_tables_dbutils
# @{

# idq content handler 
# atribute connection must be set to sqlite.connection object 
# before using this class e.g. for loading tables into database from xml files
# atribute connection should be set to None, after completion of operation(s)
# on the database.
class idq_content_handler(ligolw.LIGOLWContentHandler):
    connection = None

#    def set_connection(self, connection):   
#        self.connection = connection
#    def unset_connection(self):
#        self.connection = None
#    return idq_content_handler
#    def __init__(self, connection, *args):
#        super(idq_content_handler, self).__init__(*args)
#        self.connection = connection
        



def load_xml_files_into_database(xmlfiles, verbose = False):
    """
    Initializes in-memory sqlite data base and loads xmlfiles into it
    Returns connection and curusor objects 
    """
    # initialize sqlite connection and cursor object
    connection = sqlite3.connect(":memory:")
    cursor = connection.cursor()
    # set connection in content handler
    idq_content_handler.connection = connection
    # load files into database
    if verbose:
        print "Loading files into data base ..."
    for filename in xmlfiles:
        ligolw_sqlite.insert_from_url(filename, \
            contenthandler = idq_tables.use_in_db(idq_content_handler), \
            verbose = verbose)
    # unset connection in content handler
    idq_content_handler.connection = None

    return connection, cursor

def remove_redundant_entries(connection, cursor, verbose = False):
    """
    Clean up idq tables by removing duplicate or redundant rows. 
    """
    # get table names from the database
    tablenames = dbtables.get_table_names(connection)
    if not tablenames:
        print "No tables were found in the database. Nothing to update."
        return 

    # clean up process  and process params tables
    if (lsctables.ProcessTable.tableName in tablenames) \
        and (lsctables.ProcessParamsTable.tableName in tablenames):
        if verbose:
            print "Removing redundant entries from process and process_params tables ..."
        # get the lowest process_id from process table
        min_process_id = cursor.execute('''SELECT MIN(process_id) FROM process''').fetchone()[0]
        # delete redundant rows from process and process_params table that corresponds to different process ids
        cursor.execute('''DELETE FROM ''' + lsctables.ProcessParamsTable.tableName +
            ''' WHERE process_id != ?''',(min_process_id,))
        cursor.execute('''DELETE FROM ''' + lsctables.ProcessTable.tableName +
            ''' WHERE process_id != ?''',(min_process_id,))

        # get all tables that use process_id 
        prc_id_tables = cursor.execute("SELECT name FROM sqlite_master WHERE type == 'table' AND sql LIKE '%process_id%'").fetchall()

        if verbose:
            print "Reassigning process_ids in all tables in the database ..."
        # loop over these tables and set all process_ids to min_process_id
        for (table_name,) in prc_id_tables:
            cursor.execute('''UPDATE ''' + table_name + ''' SET process_id = ?''', (min_process_id,))
        connection.commit()
            
    # clean up coinc_definer and coinc event tables  
    if ( lsctables.CoincDefTable.tableName in tablenames ) and \
        (lsctables.CoincTable.tableName in tablenames ):
        if verbose:
            print "Removing redundant entries from coinc_definer and coinc_event tables ..."
        # get lowest (coinc_def_ids, search_coinc_type, description ) from coinc_def table
        # corresponding to unique search_coinc_type and description
        min_coinc_def_ids  = cursor.execute(''' SELECT MIN(coinc_def_id), search_coinc_type, description 
            FROM coinc_definer GROUP BY search_coinc_type, description ''').fetchall()
        
        if verbose:
            print "Updating coinc_def_id(s) in coinc_event table ..."
        # loop over them and update coinc_event table replacing redundant coinc_def ids
        for (min_coinc_def_id, min_search_coinc_type, min_description) in min_coinc_def_ids:  
                cursor.execute('''UPDATE ''' + lsctables.CoincTable.tableName +
                    ''' SET coinc_def_id = ? 
                    WHERE coinc_def_id IN  
                    (  
                    SELECT coinc_def_id FROM ''' + lsctables.CoincDefTable.tableName +
                    ''' WHERE ( search_coinc_type = ? ) AND ( description = ? ) 
                    )''', (min_coinc_def_id, min_search_coinc_type, min_description, ))
        if verbose:
            print "Removing unused entries from coinc_definer table ..."
        # remove redundant and already unused records from coinc_def table
        cursor.execute('''DELETE FROM ''' + lsctables.CoincDefTable.tableName +
            ''' WHERE coinc_def_id NOT IN 
            ( 
            SELECT MIN(coinc_def_id) FROM '''  + lsctables.CoincDefTable.tableName +
            ''' GROUP BY search_coinc_type, description 
            ) ''')
        connection.commit()
    if verbose:
        print "Finished updating tables"
        
def gpstime_inside_segment(gps_seconds, gps_nanoseconds, segment_string):
    """
    Convenience function to be used within sqlite query.
    It checks if a given gps in is inside the segment or not.
    Returns boolen.
    """
    # convert segment string into segment
    segment = segmentsUtils.from_range_strings([segment_string], boundtype=lal.LIGOTimeGPS)[0]
    return lal.LIGOTimeGPS(gps_seconds, gps_nanoseconds) in segment

def delete_glitch_events_in_segment(connection, cursor, segment):
    """
    Removes glitch events whose gps times fall within segment. 
    All rows from the other tables in database connected to these glitch
    events via coinc_map tables are also deteled.
    @param connection UNDOCUMENTED
    @param cursor UNDOCUMENTED
    @param segment is glue segment
    """
    # get table names from the database
    tablenames = dbtables.get_table_names(connection)
    if not tablenames:
        print "No tables were found in the database. Nothing to do."
        return 
    # check if glitch table is present
    if not idq_tables.IDQGlitchTable.tableName in tablenames:
        print "No glitch table is found in database. Nothing to do"
        return 
    # create sqlite function for testing if glitch gos time is inside the segment("md5", 1, md5sum)
    connection.create_function('gpstime_inside_segment', 3, gpstime_inside_segment)
    
    # convert segment inot segment string
    segment_string = segmentsUtils.to_range_strings([segment])[0]
        
    # check if coinc_def table is present
    if lsctables.CoincDefTable.tableName in tablenames:
        # get all distinct coinc types from coinc_def table 
        coinc_descriptions = cursor.execute('''SELECT DISTINCT description FROM ''' + 
            lsctables.CoincDefTable.tableName).fetchall()
        # loop over coinc descriptions
        for (description,)  in coinc_descriptions:
            tables_in_coinc = idq_tables.CoincDefToTableNames[description]
            # work only with coincs that involve glitch table
            if idq_tables.IDQGlitchTable.tableName in tables_in_coinc:
                # remove glitch table, we will work with it later
                connected_tables = [t for t in tables_in_coinc if t != idq_tables.IDQGlitchTable.tableName]
                # loop over other tables in the coinc
                for t in connected_tables:
                    # delete entries connected to glitch events that should be removed
                    cursor.execute('''DELETE FROM ''' + t + \
                    ''' WHERE event_id IN (SELECT t.event_id FROM ''' + t + \
                    ''' AS t JOIN coinc_event_map AS c1 ON t.event_id == c1.event_id''' + \
                    ''' JOIN coinc_event_map AS c2 ON c1.coinc_event_id == c2.coinc_event_id''' + \
                    ''' JOIN ''' + idq_tables.IDQGlitchTable.tableName + \
                    ''' AS g ON c2.event_id == g.event_id WHERE gpstime_inside_segment(g.gps, g.gps_ns, ?))''',\
                     (segment_string,))
                                        
    # delete entries from coinc_event table corresponding to deleted events
    cursor.execute('''DELETE FROM  coinc_event WHERE coinc_event_id IN''' + \
    ''' (SELECT c1.coinc_event_id FROM coinc_event_map AS c1 JOIN ''' + \
    idq_tables.IDQGlitchTable.tableName + \
    ''' AS g ON c1.event_id == g.event_id WHERE gpstime_inside_segment(g.gps, g.gps_ns, ?))''',\
    (segment_string,))
    
    # delete entries from coinc_map table corresponding to deleted events
    cursor.execute('''DELETE FROM  coinc_event_map WHERE event_id IN''' + \
    ''' (SELECT c1.event_id FROM coinc_event_map AS c1''' + \
    ''' JOIN coinc_event_map AS c2 ON c1.coinc_event_id == c2.coinc_event_id ''' + \
    ''' JOIN ''' + idq_tables.IDQGlitchTable.tableName +\
    ''' AS g ON c2.event_id == g.event_id WHERE gpstime_inside_segment(g.gps, g.gps_ns, ?))''',\
    (segment_string,))
    
                        
    # delete entries in glitch table
    cursor.execute('''DELETE FROM ''' + idq_tables.IDQGlitchTable.tableName +
        ''' WHERE gpstime_inside_segment(gps, gps_ns, ?)''', (segment_string,))
    
    # commit everythings
    connection.commit()

    
    
        
def delete_glitch_events_in_segmentlist(connection, cursor, segmentlist):
    """
    Remove glitch events that fall inside of a segment in the segmentlist. 
    @param connection UNDOCUMENTED
    @param cursor UNDOCUMENTED
    @param segmentlist is an instance of glue.segments.segmentlist class
    """
    # loop over segments
    for segment in segmentlist:
        delete_glitch_events_in_segment(connection, cursor, segment)


def get_glitch_ovl_data(connection, cursor, glitch_columns, ovl_columns):
    """
    Connects glitch and ovl tables, and returns requested data, 
    as a list of tuples with values for requested columns in the order of
    [(glitch_columns, ovl_columns), ...].
    @param connection UNDOCUMENTED
    @param cursor UNDOCUMENTED
    @param glitch_columns is the list of requested  glitch table column names,
    @param ovl_columns  is the list of requested ovl table column names.
    """
    # get table names from the database
    tablenames = dbtables.get_table_names(connection)
    if not tablenames:
        print "No tables were found in the database."
        return []
    # check if glitch table and ovl tables are present
    if not ( (idq_tables.IDQGlitchTable.tableName) in tablenames\
        and (idq_tables.OVLDataTable.tableName in tablenames) ):
        print "No glitch or ovl table is found in database."
        print "Need both to perform requested query."
        return []
    
    # construct feature string with requested columns for SELECT statement 
    glitch_features = ', '.join(['g.' + name for name in glitch_columns])
    ovl_features = ', '.join(['o.' + name for name in ovl_columns])
    features_string = glitch_features + ', ' + ovl_features
    # do the query 
    
    data = cursor.execute('''SELECT ''' + features_string +  ''' FROM ''' + \
        idq_tables.IDQGlitchTable.tableName + \
        ''' AS g JOIN coinc_event_map AS c1 ON g.event_id == c1.event_id''' + \
        ''' JOIN coinc_event_map AS c2 ON c1.coinc_event_id == c2.coinc_event_id''' + \
        ''' JOIN '''  + idq_tables.OVLDataTable.tableName + \
        ''' AS o ON o.event_id == c2.event_id''').fetchall()
    
    
    return data
    
    
def get_glitch_snglburst_data(connection, cursor, glitch_columns, snglburst_columns):
    """
    Connects glitch and snglburst tables, and returns requested data, 
    as a list of tuples with values for requested columns in the order of
    [(glitch_columns, snglburst_columns), ...].
    @param connection UNDOCUMENTED
    @param cursor UNDOCUMENTED
    @param glitch_columns is the list of requested  glitch table column names,
    @param snglburst_columns  is the list of requested snglburst table column names.
    """
    # get table names from the database
    tablenames = dbtables.get_table_names(connection)
    if not tablenames:
        print "No tables were found in the database."
        return []
    # check if glitch table and ovl tables are present
    if not ( (idq_tables.IDQGlitchTable.tableName in tablenames)\
        and (lsctables.SnglBurstTable.tableName in tablenames) ):
        print "No glitch or snglburst table is found in database."
        print "Need both to perform requested query."
        return []
    
    # construct feature string with requested columns for SELECT statement 
    glitch_features = ', '.join(['g.' + name for name in glitch_columns])
    snglburst_features = ', '.join(['b.' + name for name in snglburst_columns])
    features_string = glitch_features + ', ' + snglburst_features
    # do the query 
    
    data = cursor.execute('''SELECT ''' + features_string +  ''' FROM ''' + \
        idq_tables.IDQGlitchTable.tableName + \
        ''' AS g JOIN coinc_event_map AS c1 ON g.event_id == c1.event_id''' + \
        ''' JOIN coinc_event_map AS c2 ON c1.coinc_event_id == c2.coinc_event_id''' + \
        ''' JOIN '''  + lsctables.SnglBurstTable.tableName + \
        ''' AS b ON b.event_id == c2.event_id''').fetchall()
    
    return data
    
def get_get_glitch_ovl_sngburst_data(connection, cursor, glitch_columns, ovl_columns, snglburst_columns):
    """
    Connects glitch, ovl and snglburst tables, and returns requested data, 
    as a list of tuples with values for requested columns in the order of
    [(glitch_columns, ovl_columns, snglburst_columns), ...].
    @param connection UNDOCUMENTED
    @param cursor UNDOCUMENTED
    @param glitch_columns is the list of requested  glitch table column names,
    @param ovl_columns  is the list of requested ovl table column names,
    @param snglburst_columns  is the list of requested snglburst table column names.
    """
    # get table names from the database
    tablenames = dbtables.get_table_names(connection)
    if not tablenames:
        print "No tables were found in the database."
        return []
    # check if glitch table, ovl and snglburst tables are present
    if not ( (idq_tables.IDQGlitchTable.tableName in tablenames)\
        and (idq_tables.OVLDataTable.tableName in tablenames)\
        and (lsctables.SnglBurstTable.tableName in tablenames) ):
        print "Not all of glitch, ovl and snglburst tables are found in database."
        print "Need all three to perform requested query."
        return []
    
    # construct feature string with requested columns for SELECT statement 
    glitch_features = ', '.join(['g.' + name for name in glitch_columns])
    ovl_features = ', '.join(['o.' + name for name in ovl_columns])
    snglburst_features = ', '.join(['b.' + name for name in snglburst_columns])
    features_string = glitch_features + ', ' + ovl_features + ', ' + snglburst_features
    # do the query 
     
    data = cursor.execute('''SELECT ''' + features_string +  ''' FROM ''' + \
        idq_tables.IDQGlitchTable.tableName + \
        ''' AS g JOIN coinc_event_map AS c1 ON g.event_id == c1.event_id''' + \
        ''' JOIN coinc_event_map AS c2 ON c1.coinc_event_id == c2.coinc_event_id''' + \
        ''' JOIN '''  + idq_tables.OVLDataTable.tableName + \
        ''' AS o ON o.event_id == c2.event_id''' + \
        ''' JOIN coinc_event_map AS c3 ON g.event_id == c3.event_id''' + \
        ''' JOIN coinc_event_map AS c4 ON c4.coinc_event_id == c3.coinc_event_id''' + \
        ''' JOIN '''  + lsctables.SnglBurstTable.tableName + \
        ''' AS b ON b.event_id == c4.event_id''').fetchall()
    
    return data
   
##@}
