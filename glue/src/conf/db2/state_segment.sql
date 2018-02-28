-- This file is part of the Grid LSC User Environment (GLUE)
-- 
-- GLUE is free software: you can redistribute it and/or modify it under the
-- terms of the GNU General Public License as published by the Free Software
-- Foundation, either version 3 of the License, or (at your option) any later
-- version.
-- 
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
-- FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
-- details.
-- 
-- You should have received a copy of the GNU General Public License along
-- with this program.  If not, see <http://www.gnu.org/licenses/>.

CREATE VIEW state_segment (
-- Create a view from the segment, segment_def_map, segment_lfn_map and 
-- segment_definer tables that will give us a virtual table that the 
-- segment publishing script can insert data into.

-- The creator database of the entry in the segment table
      creator_db,

-- Unique process ID of the process which defined this segment
      process_id,

-- Unique segment ID of this segment
      segment_id,

-- Segment start and end times, in GPS seconds.
      start_time,
      start_time_ns,
      end_time,
      end_time_ns,

-- Science segment number
      segnum,
      
-- Segment definer ID which explains the meaning of this segment 
      segment_def_id,
      segment_def_cdb,
      
-- ID of the LFN of the frame from which segment was derived
      lfn_id,
      lfn_cdb
) AS
SELECT 
-- Elements of the view from the individual tables
      segment.creator_db,
      segment.process_id,
      segment.segment_id,
      segment.start_time,
      segment.start_time_ns,
      segment.end_time,
      segment.end_time_ns,
      segment.segnum,
      segment_definer.segment_def_id,
      segment_definer.creator_db,
      segment_lfn_map.lfn_id,
      segment_lfn_map.creator_db
FROM
      segment,
      segment_def_map,
      segment_definer,
      segment_lfn_map
WHERE

-- releate the segment to the segment_definer
      segment.segment_id = segment_def_map.segment_id 
      AND segment.creator_db = segment_def_map.segment_cdb
      AND segment_def_map.segment_def_id = segment_definer.segment_def_id 
      AND segment_def_map.segment_def_cdb = segment_definer.creator_db
      
-- relate the segment to the lfn id
      AND segment.segment_id = segment_lfn_map.segment_id
      AND segment.creator_db = segment_lfn_map.segment_cdb

-- only select state segments from the segment_definer table
      AND segment_definer.state_vec_major IS NOT NULL 
      AND segment_definer.state_vec_minor IS NOT NULL
@


-- Create a trigger which will allow inserts into the state_segment view
-- by creating the appropriate entries in the segment, segment_def_map,
-- and segment_lfn_map tables
CREATE TRIGGER state_segment_i instead OF INSERT ON state_segment
      referencing new AS n 
      FOR each ROW MODE db2sql
      BEGIN atomic

-- Check that a segment with the same start time, end time and segment_definer
-- as the segment we are trying to insert does not already exist in the
-- segment table. Raise and error if it does.
      VALUES ( 
        CASE WHEN ( (SELECT count(s.segment_id) FROM 
            segment AS s, segment_def_map AS m WHERE
            m.segment_def_id = n.segment_def_id 
            AND m.segment_def_cdb = n.segment_def_cdb
            AND s.start_time = n.start_time 
            AND s.start_time_ns = n.start_time_ns 
            AND s.end_time = n.end_time 
            AND s.end_time_ns = n.end_time_ns 
            AND s.segment_id = m.segment_id
            AND s.creator_db = m.segment_cdb) > 0 ) 
        THEN raise_error ( '70001', 'state_segment rows must be unique' ) 
        ELSE 0 END );
        
-- Create an entry in the segment table with the active flag set to 1
      INSERT INTO segment
            (process_id,segment_id,active,
            start_time,start_time_ns,end_time,end_time_ns,segnum)
      VALUES 
            (n.process_id,n.segment_id,1,
            n.start_time,n.start_time_ns,n.end_time,n.end_time_ns,n.segnum);

-- Create an entry in the segment_def_map table which maps to the
-- segment_def_id that the user has provided for this segment
      INSERT INTO segment_def_map
            (seg_def_map_id,
              process_id,segment_id,segment_def_id,segment_def_cdb)
      VALUES 
            (GENERATE_UNIQUE(),
              n.process_id,n.segment_id,n.segment_def_id,
              CASE WHEN n.segment_def_cdb IS NOT NULL
                THEN n.segment_def_cdb
                ELSE 1 -- segment_def_cdb
              END
              );

-- Create an entry in the lfm_map table which maps the lfn_id which the
-- user has provided to this segment
      INSERT INTO segment_lfn_map
            (seg_lfn_map_id,process_id,segment_id,lfn_id,lfn_cdb)
      VALUES 
            (GENERATE_UNIQUE(),n.process_id,n.segment_id,n.lfn_id,
              CASE WHEN n.lfn_cdb IS NOT NULL
                THEN n.lfn_cdb
                ELSE 1 -- lfn_cdb
              END
              );

    END
@
