#!/usr/bin/env python
# Id$
#
# Copyright (C) 2009  Larne Pekowsky
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

import sys
import os
import re

from glue.segments import segment, segmentlist
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.segmentdb import query_engine
from glue.ligolw import types as ligolwtypes



#
# =============================================================================
#
#                     Routines to set up backends
#
# =============================================================================
#


def get_all_files_in_range(dirname, starttime, endtime, pad=64):
    """Returns all files in dirname and all its subdirectories whose
    names indicate that they contain segments in the range starttime
    to endtime"""
    
    ret = []

    # Maybe the user just wants one file...
    if os.path.isfile(dirname):
        if re.match('.*-[0-9]*-[0-9]*\.xml$', dirname):
            return [dirname]
        else:
            return ret

    first_four_start = starttime / 100000
    first_four_end   = endtime   / 100000

    for filename in os.listdir(dirname):
        if re.match('.*-[0-9]{5}$', filename):
            dirtime = int(filename[-5:])
            if dirtime >= first_four_start and dirtime <= first_four_end:
                ret += get_all_files_in_range(os.path.join(dirname,filename), starttime, endtime, pad=pad)
        elif re.match('.*-[0-9]{4}$', filename):
            dirtime = int(filename[-4:])
            if dirtime >= first_four_start and dirtime <= first_four_end:
                ret += get_all_files_in_range(os.path.join(dirname,filename), starttime, endtime, pad=pad)
        elif re.match('.*-[0-9]*-[0-9]*\.xml$', filename):
            file_time = int(filename.split('-')[-2])
            if file_time >= (starttime-pad) and file_time <= (endtime+pad):
                ret.append(os.path.join(dirname,filename))
        else:
            # Keep recursing, we may be looking at directories of
            # ifos, each of which has directories with times
            ret += get_all_files_in_range(os.path.join(dirname,filename), starttime, endtime, pad=pad)

    return ret



def setup_database(database_location):
    """ 1. Determine protocol"""
    try:
	# When no protocol is given:
    	if database_location.find('://') == -1:
           msg = "Error: Please specify protocol in your --segment-url argument in the format PROTOCOL://HOST"
	   msg +="\nFor example: --segment-url https://segdb.ligo.caltech.edu"
    	   msg += "\nSupported protocols include: http, https, ldbd, ldbdi"
    	   msg += "\nRun with --help for usage"
	   raise ValueError(msg)

	# When wrong protocol is given:
	protocol = database_location[:database_location.find('://')].lower()
	if protocol not in ("http","https","ldbd","ldbdi"): 
	   msg = "Error: protocol %s not supported" % protocol
           msg += "\nPlease specify correct protocol in your --segment-url argument in the format PROTOCOL://HOST"
    	   msg += "\nSupported protocols include: http, https, ldbd, ldbdi"
	   msg += "\nRun with --help for usage"
	   raise ValueError(msg)
    except ValueError, e:
      	sys.stderr.write('%s\n' % str(e))
      	sys.exit(1)

    """ 2. Determine host and port"""
    host_and_port = database_location[(len(protocol)+3):]
    if host_and_port.find(':') < 0:
       # if no port number given, set default port respectively:
       host = host_and_port
       if protocol == 'http':
         port = 80
       elif protocol == 'https':
         port = 443
       elif protocol == 'ldbd':
         port = 30015
       elif protocol == 'ldbdi':
         port = 30016
    else:
       host, portString = host_and_port.split(':')
       port = int(portString)

    if port == 30020:
       sys.stderr.write("Error: PORT 30020 no longer provide segment database service\n")
       sys.exit(1)

    """ 3. Set up connection to LDBD(W)Server """
    client = None
    identity = "/DC=org/DC=doegrids/OU=Services/CN=ldbd/"

    if protocol.startswith('http'):
      from glue import LDBDWClient
      if protocol == "https":
         identity += host
      else:
         identity = None
      try:
        client = LDBDWClient.LDBDClient(host,port,protocol,identity)
      except Exception, e:
        sys.stderr.write("Unable to connect to LDBD Server at %s://%s:%d \n" % (protocol,host, port) + str(e))
        sys.exit(1)

    elif protocol.startswith('ldbd'):
      from glue import LDBDClient
      if protocol == "ldbd":
         identity += host
         from glue import gsiserverutils
      else:
         identity = None
      try:
        client = LDBDClient.LDBDClient(host,port,identity)
      except Exception, e:
        sys.stderr.write("Unable to connect to LDBD Server at %s://%s:%d\n" % (protocol,host, port) + str(e))
        try:
          if gsiserverutils.checkCredentials():
             sys.stderr.write("Got the following error : \n" + str(e))
             sys.stderr.write("Run wiht --help for usage\n")
        except UnboundLocalError:
          pass
        sys.exit(1)

    else:
        raise ValueError( "invalid url for segment database" )
    

    return client


#
# =============================================================================
#
#        Routines to find segment information in databases/XML docs
#
# =============================================================================
#



def query_segments(engine, table, segdefs):
    # each segdef is a list containing:
    #     ifo, name, version, start_time, end_time, start_pad, end_pad


    # The trivial case: if there's nothing to do, return no time
    if len(segdefs) == 0:
        return [ segmentlist([]) ]

    #
    # For the sake of efficiency we query the database for all the segdefs at once
    # This constructs a clause that looks for one
    #
    def make_clause(table, segdef):
        ifo, name, version, start_time, end_time, start_pad, end_pad = segdef

        sql = " (segment_definer.ifos = '%s' " % ifo
        sql += "AND segment_definer.name = '%s' " % name
        sql += "AND segment_definer.version = %s " % version
        sql += "AND NOT (%d > %s.end_time OR %s.start_time > %d)) " % (start_time, table, table, end_time)

        return sql

    clauses = [make_clause(table, segdef) for segdef in segdefs]

    sql  = 'SELECT segment_definer.ifos, segment_definer.name, segment_definer.version, '
    sql += ' %s.start_time, %s.end_time ' % (table, table)
    sql += ' FROM segment_definer, %s '   % table
    sql += ' WHERE %s.segment_def_id = segment_definer.segment_def_id AND ' % table

    if engine.__class__ == query_engine.LdbdQueryEngine:
        sql += " %s.segment_def_cdb = segment_definer.creator_db AND " % table
    
    sql += '( ' + ' OR '.join(clauses) + ' )'

    rows = engine.query(sql)

    #
    # The result of a query will be rows of the form
    #    ifo, name, version, start_time, end_time
    #
    # We want to associate each returned row with the segdef it belongs to so that
    # we can apply the correct padding.
    #
    # If segdefs were uniquely spcified by (ifo, name, version) this would
    # be easy, but it may happen that we're looking for the same segment definer
    # at multiple disjoint times.  In particular this can happen if the user
    # didn't specify a version number; in that case we might have version 2
    # of some flag defined over multiple disjoint segment_definers.
    #
    results = []

    for segdef in segdefs:
        ifo, name, version, start_time, end_time, start_pad, end_pad = segdef

        search_span      = segment(start_time, end_time)
        search_span_list = segmentlist([search_span])

        # See whether the row belongs to the current segdef.  Name, ifo and version must match
        # and the padded segment must overlap with the range of the segdef.
        def matches(row):
            return ( row[0].strip() == ifo and row[1] == name and int(row[2]) == int(version)
                     and search_span.intersects(segment(row[3] + start_pad, row[4] + start_pad)) )

        # Add the padding.  Segments may extend beyond the time of interest, chop off the excess.
        def pad_and_truncate(row_start, row_end):
            tmp = segmentlist([segment(row_start + start_pad, row_end + end_pad)])
            # No coalesce needed as a list with a single segment is already coalesced
            tmp &= search_span_list

            # The intersection is guaranteed to be non-empty if the row passed match()
            # PR 2969: The above comment is incorrect.  Negative padding may cause
            # an empty intersection.
            if len(tmp) == 0:
                return segment(0,0)
            else:
                return tmp[0]

        # Build a segment list from the returned segments, padded and trunctated.  The segments will
        # not necessarily be disjoint, if the padding crosses gaps.  They are also not gauranteed to
        # be in order, since there's no ORDER BY in the query.  So the list needs to be coalesced
        # before arithmatic can be done with it.
        result  = segmentlist( [pad_and_truncate(row[3], row[4]) for row in rows if matches(row)] ).coalesce()

        # This is not needed: since each of the segments are constrained to be within the search
        # span the whole list must be as well.
        # result &= search_span_list
        
        results.append(result)

    return results


def expand_version_number(engine, segdef):
    ifo, name, version, start_time, end_time, start_pad, end_pad = segdef

    if version != '*':
        return [segdef]

    # Start looking at the full interval
    intervals = segmentlist([segment(start_time, end_time)])

    # Find the maximum version number
    sql  = "SELECT max(version) FROM segment_definer "
    sql += "WHERE  segment_definer.ifos = '%s' " % ifo
    sql += "AND   segment_definer.name = '%s' " % name

    rows    = engine.query(sql)
    try:
        version = len(rows[0]) and rows[0][0] or 1
    except:
        version = None

    results = []

    while version > 0:
        for interval in intervals:
            segs = query_segments(engine, 'segment_summary', [(ifo, name, version, interval[0], interval[1], 0, 0)])

            for seg in segs[0]:
                results.append( (ifo, name, version, seg[0], seg[1], 0, 0) )

        intervals.coalesce()
        intervals -= segs[0]

        version -= 1

    return results




def find_segments(doc, key, use_segment_table = True):
    key_pieces = key.split(':')
    while len(key_pieces) < 3:
        key_pieces.append('*')

    filter_func = lambda x: str(x.ifos) == key_pieces[0] and (str(x.name) == key_pieces[1] or key_pieces[1] == '*') and (str(x.version) == key_pieces[2] or key_pieces[2] == '*') 

    # Find all segment definers matching the critieria
    seg_def_table = table.get_table(doc, lsctables.SegmentDefTable.tableName)
    seg_defs      = filter(filter_func, seg_def_table)
    seg_def_ids   = map(lambda x: str(x.segment_def_id), seg_defs)

    # Find all segments belonging to those definers
    if use_segment_table:
        seg_table     = table.get_table(doc, lsctables.SegmentTable.tableName)
        seg_entries   = filter(lambda x: str(x.segment_def_id) in seg_def_ids, seg_table)
    else:
        seg_sum_table = table.get_table(doc, lsctables.SegmentSumTable.tableName)
        seg_entries   = filter(lambda x: str(x.segment_def_id) in seg_def_ids, seg_sum_table)

    # Combine into a segmentlist
    ret = segmentlist(map(lambda x: segment(x.start_time, x.end_time), seg_entries))

    ret.coalesce()

    return ret

#
# =============================================================================
#
#                      General utilities
#
# =============================================================================
#
def ensure_segment_table(connection):
    """Ensures that the DB represented by connection posses a segment table.
    If not, creates one and prints a warning to stderr"""

    count = connection.cursor().execute("SELECT count(*) FROM sqlite_master WHERE name='segment'").fetchone()[0]

    if count == 0:
        sys.stderr.write("WARNING: None of the loaded files contain a segment table\n")
        theClass  = lsctables.TableByName['segment']
        statement = "CREATE TABLE IF NOT EXISTS segment (" + ", ".join(map(lambda key: "%s %s" % (key, ligolwtypes.ToSQLiteType[theClass.validcolumns[key]]), theClass.validcolumns)) + ")"

        connection.cursor().execute(statement)



# =============================================================================
#
#                    Routines to write data to XML documents
#
# =============================================================================
#

def add_to_segment_definer(xmldoc, proc_id, ifo, name, version, comment=''):
    try:
        seg_def_table = table.get_table(xmldoc, lsctables.SegmentDefTable.tableName)
    except:
        seg_def_table = lsctables.New(lsctables.SegmentDefTable, columns = ["process_id", "segment_def_id", "ifos", "name", "version", "comment"])
        xmldoc.childNodes[0].appendChild(seg_def_table)

    seg_def_id                     = seg_def_table.get_next_id()
    segment_definer                = lsctables.SegmentDef()
    segment_definer.process_id     = proc_id
    segment_definer.segment_def_id = seg_def_id
    segment_definer.ifos           = ifo
    segment_definer.name           = name
    segment_definer.version        = version
    segment_definer.comment        = comment

    seg_def_table.append(segment_definer)

    return seg_def_id



def add_to_segment(xmldoc, proc_id, seg_def_id, sgmtlist):
    try:
        segtable = table.get_table(xmldoc, lsctables.SegmentTable.tableName)
    except:
        segtable = lsctables.New(lsctables.SegmentTable, columns = ["process_id", "segment_def_id", "segment_id", "start_time", "start_time_ns", "end_time", "end_time_ns"])
        xmldoc.childNodes[0].appendChild(segtable)

    for seg in sgmtlist:
        segment                = lsctables.Segment()
        segment.process_id     = proc_id
        segment.segment_def_id = seg_def_id
        segment.segment_id     = segtable.get_next_id()
        segment.start_time     = seg[0]
        segment.start_time_ns  = 0
        segment.end_time       = seg[1]
        segment.end_time_ns    = 0

        segtable.append(segment)


def add_to_segment_summary(xmldoc, proc_id, seg_def_id, sgmtlist, comment=''):
    try:
        seg_sum_table = table.get_table(xmldoc, lsctables.SegmentSumTable.tableName)
    except:
        seg_sum_table = lsctables.New(lsctables.SegmentSumTable, columns = ["process_id", "segment_def_id", "segment_sum_id", "start_time", "start_time_ns", "end_time", "end_time_ns", "comment"])
        xmldoc.childNodes[0].appendChild(seg_sum_table)

    for seg in sgmtlist:
        segment_sum                = lsctables.SegmentSum()
        segment_sum.process_id     = proc_id
        segment_sum.segment_def_id = seg_def_id
        segment_sum.segment_sum_id = seg_sum_table.get_next_id()
        segment_sum.start_time     = seg[0]
        segment_sum.start_time_ns  = 0
        segment_sum.end_time       = seg[1]
        segment_sum.end_time_ns    = 0
        segment_sum.comment        = comment

        seg_sum_table.append(segment_sum)


def add_segment_info(doc, proc_id, segdefs, segments, segment_summaries):

    for i in range(len(segdefs)):
        ifo, name, version, start_time, end_time, start_pad, end_pad = segdefs[i]

        seg_def_id = add_to_segment_definer(doc, proc_id, ifo, name, version)

        add_to_segment_summary(doc, proc_id, seg_def_id, segment_summaries[i])

        if segments:
            add_to_segment(doc, proc_id, seg_def_id, segments[i])

#
# =============================================================================
#
#                      Routines that should be obsolete
#
# =============================================================================
#

def build_segment_list(engine, gps_start_time, gps_end_time, ifo, segment_name, version = None, start_pad = 0, end_pad = 0):
    """Optains a list of segments for the given ifo, name and version between the
    specified times.  If a version is given the request is straightforward and is
    passed on to build_segment_list_one.  Otherwise more complex processing is
    performed (not yet implemented)"""
    if version is not None:
        return build_segment_list_one(engine, gps_start_time, gps_end_time, ifo, segment_name, version, start_pad, end_pad)

    # This needs more sophisticated logic, for the moment just return the latest
    # available version
    sql  = "SELECT max(version) FROM segment_definer "
    sql += "WHERE  segment_definer.ifos = '%s' " % ifo
    sql += "AND   segment_definer.name = '%s' " % segment_name

    rows = engine.query(sql)
    version = len(rows[0]) and rows[0][0] or 1

    return build_segment_list_one(engine, gps_start_time, gps_end_time, ifo, segment_name, version, start_pad, end_pad)


def build_segment_list_one(engine, gps_start_time, gps_end_time, ifo, segment_name, version = None, start_pad = 0, end_pad = 0):
    """Builds a list of segments satisfying the given criteria """
    seg_result = segmentlist([])
    sum_result = segmentlist([])

    # Is there any way to get segment and segement summary in one query?
    # Maybe some sort of outer join where we keep track of which segment
    # summaries we've already seen.
    sql = "SELECT segment_summary.start_time, segment_summary.end_time "
    sql += "FROM segment_definer, segment_summary "
    sql += "WHERE segment_summary.segment_def_id = segment_definer.segment_def_id "
    sql += "AND   segment_definer.ifos = '%s' " % ifo
    if engine.__class__ == query_engine.LdbdQueryEngine:
       sql += "AND segment_summary.segment_def_cdb = segment_definer.creator_db "
    sql += "AND   segment_definer.name = '%s' " % segment_name
    sql += "AND   segment_definer.version = %s " % version
    sql += "AND NOT (%s > segment_summary.end_time OR segment_summary.start_time > %s)" % (gps_start_time, gps_end_time)

    rows = engine.query(sql)

    for sum_start_time, sum_end_time in rows:
        sum_start_time = (sum_start_time < gps_start_time) and gps_start_time or sum_start_time
        sum_end_time = (sum_end_time > gps_end_time) and gps_end_time or sum_end_time

        sum_result |= segmentlist([segment(sum_start_time, sum_end_time)])

    # We can't use queries paramaterized with ? since the ldbd protocol doesn't support it...
    sql = "SELECT segment.start_time + %d, segment.end_time + %d " % (start_pad, end_pad)
    sql += "FROM segment, segment_definer "
    sql += "WHERE segment.segment_def_id = segment_definer.segment_def_id "

    if engine.__class__ == query_engine.LdbdQueryEngine:
       sql += "AND segment.segment_def_cdb = segment_definer.creator_db "
    sql += "AND   segment_definer.ifos = '%s' " % ifo
    sql += "AND   segment_definer.name = '%s' " % segment_name
    sql += "AND   segment_definer.version = %s " % version
    sql += "AND NOT (%s > segment.end_time OR segment.start_time > %s)" % (gps_start_time, gps_end_time)

    rows = engine.query(sql)
    
    for seg_start_time, seg_end_time in rows:
        seg_start_time = (seg_start_time < gps_start_time) and gps_start_time or seg_start_time
        seg_end_time = (seg_end_time > gps_end_time) and gps_end_time or seg_end_time

        seg_result |= segmentlist([segment(seg_start_time, seg_end_time)])

    engine.close()

    return sum_result, seg_result



def run_query_segments(doc, proc_id, engine, gps_start_time, gps_end_time, included_segments_string, excluded_segments_string = None, write_segments = True, start_pad = 0, end_pad = 0):
    """Runs a segment query.  This was originally part of ligolw_query_segments, but now is also
    used by ligolw_segments_from_cats.

    The write_segments option is provided so callers can coalesce segments obtained over
    sever invocations (as segments_from_cats does).
    """
    
    if write_segments:
        all_ifos = {}

        for ifo, segment_name, version in split_segment_ids(included_segments_string.split(',')):
            all_ifos[ifo] = True


        new_seg_def_id = add_to_segment_definer(doc, proc_id, ''.join(all_ifos.keys()), 'result', 0)
        add_to_segment_summary(doc, proc_id, new_seg_def_id, [[gps_start_time, gps_end_time]])

    result = segmentlist([])

    for ifo, segment_name, version in split_segment_ids(included_segments_string.split(',')):
        sum_segments, seg_segments = build_segment_list(engine, gps_start_time, gps_end_time, ifo, segment_name, version, start_pad, end_pad)
        seg_def_id                 = add_to_segment_definer(doc, proc_id, ifo, segment_name, version)

        add_to_segment_summary(doc, proc_id, seg_def_id, sum_segments)

        # and accumulate segments 
        result |= seg_segments

    # Excluded segments are not required
    if excluded_segments_string:
        excluded_segments = segmentlist([])

        for ifo, segment_name, version in split_segment_ids(excluded_segments_string.split(',')):
            sum_segments, seg_segments = build_segment_list(engine, gps_start_time, gps_end_time, ifo, segment_name, version)
            excluded_segments |= seg_segments

        result = result - excluded_segments

    result.coalesce()
    
    # Add the segments
    if write_segments:
        add_to_segment(doc, proc_id, new_seg_def_id, result)

    return result



def split_segment_ids(segment_ids):
    """Given an array of strings of the form ifo:name and
    ifo:name:version, returns an array of tuples of the form (ifo,
    name, version) where version may be None"""
    
    def split_segment_id(segment_id):
        temp = segment_id.split(':')
        if len(temp) == 2:
            temp.append(None)
        elif temp[2] == '*':
            temp[2] = None
        else:
            temp[2] = int(temp[2])
            
        return temp

    return map(split_segment_id, segment_ids)


