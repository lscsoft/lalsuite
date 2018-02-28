#!/usr/bin/env python

import sys
import os

from glue.segmentdb import query_engine
from glue.segmentdb import segmentdb_utils
from glue.segments import segment, segmentlist


def get_manually(gps_start_time, gps_end_time):
    db_location = os.environ['S6_SEGMENT_SERVER']
    segment_connection   = segmentdb_utils.setup_database(db_location)
    engine = query_engine.LdbdQueryEngine(segment_connection)

    # 1. Get v1 science segments

    sql = "SELECT segment.start_time, segment.end_time "
    sql += "FROM segment_definer, segment "
    sql += "WHERE segment.segment_def_id = segment_definer.segment_def_id "
    sql += "AND   segment_definer.ifos = 'H1' " 
    sql += "AND   segment.segment_def_cdb = segment_definer.creator_db "
    sql += "AND   segment_definer.name = 'DMT-SCIENCE' "
    sql += "AND   segment_definer.version = 1 " 
    sql += "AND NOT (%s > segment.end_time OR segment.start_time > %s)" % (gps_start_time, gps_end_time)

    v1_science_segments = segmentlist([segment(row[0], row[1]) for row in engine.query(sql)]).coalesce()

    # 2. Get v2 science summaries

    sql = "SELECT segment_summary.start_time, segment_summary.end_time "
    sql += "FROM segment_definer, segment_summary "
    sql += "WHERE segment_summary.segment_def_id = segment_definer.segment_def_id "
    sql += "AND   segment_definer.ifos = 'H1' " 
    sql += "AND   segment_summary.segment_def_cdb = segment_definer.creator_db "
    sql += "AND   segment_definer.name = 'DMT-SCIENCE' "
    sql += "AND   segment_definer.version = 2 " 
    sql += "AND NOT (%s > segment_summary.end_time OR segment_summary.start_time > %s)" % (gps_start_time, gps_end_time)

    v2_science_summaries = segmentlist([segment(row[0], row[1]) for row in engine.query(sql)]).coalesce()

    # 1. Get v2 science segments

    sql = "SELECT segment.start_time, segment.end_time "
    sql += "FROM segment_definer, segment "
    sql += "WHERE segment.segment_def_id = segment_definer.segment_def_id "
    sql += "AND   segment_definer.ifos = 'H1' " 
    sql += "AND   segment.segment_def_cdb = segment_definer.creator_db "
    sql += "AND   segment_definer.name = 'DMT-SCIENCE' "
    sql += "AND   segment_definer.version = 2 " 
    sql += "AND NOT (%s > segment.end_time OR segment.start_time > %s)" % (gps_start_time, gps_end_time)

    v2_science_segments = segmentlist([segment(row[0], row[1]) for row in engine.query(sql)]).coalesce()

    result = (v1_science_segments - v2_science_summaries) + v2_science_segments

    result.coalesce()

    result &= segmentlist([segment(gps_start_time, gps_end_time)])

    return result

def get_from_segment_query(gps_start_time, gps_end_time):
    pipe = os.popen('ligolw_segment_query --database --query-segments --include-segments H1:DMT-SCIENCE -s %d -e %d | ligolw_print -t segment -c start_time -c end_time' % (gps_start_time, gps_end_time))

    lines  = [l.strip().split(',') for l in pipe]
    result = segmentlist([segment(int(l[0]), int(l[1])) for l in lines])

    result.coalesce()

    return result


if __name__ == '__main__':
    # Brennan's illustration of the bug
    gps_start_time = 931035615 
    gps_end_time   = 934148055

    failed = False

    if get_manually(gps_start_time, gps_end_time) != get_from_segment_query(gps_start_time, gps_end_time):
        failed = True
        print >>sys.stderr, "Test 1 failed"

    # Tomoki's illustration of the bug
    gps_start_time = 932774415 
    gps_end_time   = 933379215

    if get_manually(gps_start_time, gps_end_time) != get_from_segment_query(gps_start_time, gps_end_time):
        failed = True
        print >>sys.stderr, "Test 2 failed"

    if not failed:
        print "All tests succeeded"

