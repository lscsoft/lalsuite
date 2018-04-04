#
# Copyright (C) 2006  Kipp C. Cannon
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

import sys
import os
import time
import socket
import pwd

try:
  import DB2
except:
  pass

try:
  from glue import gpstime
  from glue import segments
except ImportError:
  raise ImportError("Error, unable to import modules from glue. Check that glue is correctly installed and in your PYTHONPATH.")

#================================================================================
__author__ = "Ping Wei <piwei@physics.syr.edu>"
from glue import git_version
__date__ = git_version.date
__version__ = git_version.id
__src__ = "$Source$"
#================================================================================

def coalesce_seg(database, start_time, end_time):
  ret = 0            #assume execution successufl

  try:
    st = int(start_time)
    et = int(end_time)
    db = str(database.strip())
    #-------------------------------------------------------------------
    # Set up environment and get needed values
    #-------------------------------------------------------------------
    # Set up connection to the database
    dbconn = DB2.connect(dsn=db, uid='', pwd='', autoCommit=True)
    curs = dbconn.cursor()

    # create a new process_id
    sql = "select hex(GENERATE_UNIQUE()) from sysibm.sysdummy1"
    curs.execute(sql)
    hex_procid = curs.fetchone()[0]
    process_id = 'x' + '\'' + hex_procid + '\''

    # determine the local creator_db
    sql = "SELECT DEFAULT FROM SYSCAT.COLUMNS WHERE "
    sql += "TABNAME = 'PROCESS' AND COLNAME = 'CREATOR_DB'"
    curs.execute(sql)
    creator_db = int(curs.fetchone()[0])

  
    # prepare values for the new row to be inserted into the process table
    program = os.path.abspath(sys.argv[0])
    node = socket.gethostname()
    username = pwd.getpwuid(os.getuid()).pw_name
    unix_procid = os.getpid()
    proc_start_time = gpstime.GpsSecondsFromPyUTC(time.time())
    end_time = None
    jobid = 0
    domain = 'coalesce_local'

    # insert new row into process table
    sql = "INSERT INTO process "
    sql += "(program, is_online, node, username, unix_procid, start_time, jobid, domain, process_id, creator_db) "
    sql += "VALUES ('%s', 0, '%s', '%s', %d, %d, %d, '%s',%s, %d)" % (program, node, username, unix_procid, proc_start_time, jobid, domain, process_id, creator_db)
    curs.execute(sql)

    # get the BLOB process_id for later reference
    sql = "SELECT BLOB(process_id) from process where hex(process_id)='%s' " % hex_procid
    curs.execute(sql)
    blob_procid = curs.fetchone()[0]


    #========================================================================
    #
    #                                Main
    #
    #========================================================================
    # Algorithm:
    # 1. Find distinct version 1 segment type from segment_summary table witnin start_time, end_time range 
    # 2. Find segments and intervals to coalesce
    # 3. Coalesce segments and intervals
    # 4. Insert coaleseced segments back in to the database
    # 5. Delete uncoalesced segments and intervals from the database


    # 1. Find distinct segment types matching our criteria from segment_summary within the specified time range
    sql  = "SELECT distinct(hex(segment_summary.segment_def_id)) FROM segment_summary, segment_definer, process "
    sql += "WHERE segment_summary.segment_def_id=segment_definer.segment_def_id "
    sql += "AND segment_summary.segment_def_cdb=segment_definer.creator_db "
    sql += "AND segment_summary.process_id=process.process_id "
    sql += "AND segment_summary.creator_db=process.creator_db "
    # Removed next line so that all segments are coalesced: this will be slower up front but faster for queries and the long run
    #sql += "AND ((segment_definer.name like 'DMT-%' and segment_definer.version=1) or (process.ifos='V1' and process.program='SegOnline')) "
    sql += "AND segment_summary.start_time <=%d " % et
    sql += "AND segment_summary.end_time >= %d " % st
    curs.execute(sql)
    def_ids = curs.fetchall()
    if not def_ids:
      data_existence = 0
    else:
      data_existence = 1

    # loop in the segment types to fetch, coalesce, insert and delete
    for d in def_ids:
      # get the BLOB segment_def_id for later use 
      sql = "SELECT BLOB(segment_def_id), ifos, name, version, creator_db " 
      sql += "FROM segment_definer " 
      sql += "WHERE hex(segment_def_id) = '%s' " % d[0]

      curs.execute(sql)
      result = curs.fetchone()
      blob_defid = result[0]
      ifos = result[1].strip() 
      name = result[2]
      ver = result[3]
      def_cdb = result[4]

      # 2. Find segments and intervals to coalesce
      # get the segment start_time, end_time to coalesce, and according primary key to delete
      try:
        curs.execute("drop view seg_view")
      except:
        pass
      sql = "CREATE view seg_view (st,et,seg_id) AS "
      sql += "SELECT start_time,end_time, segment_id from segment "
      sql += "WHERE hex(segment_def_id) = '%s' " % d[0]
      sql += "AND segment.start_time <=%d " % et
      sql += "AND segment.end_time >= %d " % st
      sys.stdout.write("Selecting segments to coalesce for %s version:%d %s ... \n" % (ifos,ver, name))
      curs.execute(sql)

      curs.execute("SELECT st,et from seg_view")
      seg_bf_cos = curs.fetchall()   # get the segments to coalesce

      # get the summary start_time, end_time to coalesce, and according primary key to delete
      try:
        curs.execute("drop view sum_view")
      except:
        pass
      sql = "CREATE view sum_view (st,et,sum_id) AS "
      sql += "SELECT start_time,end_time, segment_sum_id from segment_summary "
      sql += "WHERE hex(segment_def_id) = '%s' " % d[0]
      sql += "AND segment_summary.start_time <=%d " % et
      sql += "AND segment_summary.end_time >= %d " % st
      curs.execute(sql)

      curs.execute("SELECT st,et from sum_view")
      sum_bf_cos = curs.fetchall()   # get the summarys to coalesce

      # 3. Coalesce segments and intervals
      sys.stdout.write("Coalescing segments ... \n")
      segs = segments.segmentlist([]) 
      sums = segments.segmentlist([]) 
      for bf in seg_bf_cos:
        seg = segments.segment(int(bf[0]), int(bf[1]))
        segs.append(seg) 
      for bf in sum_bf_cos:
        sum = segments.segment(int(bf[0]), int(bf[1]))
        sums.append(sum) 

      segs.coalesce()
      sums.coalesce()


      # 4. Insert coaleseced segments back in to the database
      # insert coalesced segs into segment table
      insert_list = []
      for s in segs:
        # generate unique id for insertion
        curs.execute("VALUES BLOB(GENERATE_UNIQUE())")
        prim_id = curs.fetchone()[0]
        # generate a list of values to insert using executemany()
        insert_list.append((prim_id, creator_db, s[0], s[1], blob_defid, def_cdb, blob_procid))

      sql = "INSERT INTO segment "
      sql += "(segment_id, creator_db, start_time, end_time, segment_def_id, segment_def_cdb, process_id) "
      sql += "VALUES (?,?,?,?,?,?,?) "
      sys.stdout.write("Inserting coalesced segments back in ... \n")
      curs.executemany(sql, insert_list)

      # insert coalesced sums into segment_summary table
      insert_list = []
      for s in sums:
        # generate unique id for insertion
        curs.execute("VALUES BLOB(GENERATE_UNIQUE())")
        prim_id = curs.fetchone()[0]
        # generate a list of values to insert using executemany()
        insert_list.append((prim_id, creator_db, s[0], s[1], blob_defid, def_cdb, blob_procid))
      sql = "INSERT INTO segment_summary "
      sql += "(segment_sum_id, creator_db, start_time, end_time, segment_def_id, segment_def_cdb, process_id) "
      sql += "VALUES (?,?,?,?,?,?,?) "
      curs.executemany(sql, insert_list)

      # 5. Delete uncoalesced segments and intervals from the database
      sys.stdout.write("Deleting un-coaleseced segments ... \n\n")
      sql = "DELETE FROM segment "
      sql += "WHERE segment_id in (select seg_id from seg_view) "
      sql += "AND process_id != %s " % process_id
      curs.execute(sql)

      sql = "DELETE FROM segment_summary "
      sql += "WHERE segment_sum_id in (select sum_id from sum_view) "
      sql += "AND process_id != %s " % process_id
      curs.execute(sql)

    # update end_time in process table
    sql = "update process set end_time=%d where hex(process_id)='%s' " % (gpstime.GpsSecondsFromPyUTC(time.time()),hex_procid)
    curs.execute(sql)
  
    try:  
      curs.execute("drop view seg_view")
      curs.execute("drop view sum_view")
    except:
      pass
    curs.close()

  except Exception as e:
    ret = str(e)
    sys.stdout.write("%s\n" % ret)

  return ret,data_existence
