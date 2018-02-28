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

CREATE TABLE segment_lfn_map
(
-- Create a map between segments and logical file names (LFNs)

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Process ID of the program which created this map
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Unique ID for this map
      seg_lfn_map_id     CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the segment
      segment_cdb        INTEGER NOT NULL WITH DEFAULT 1,
      segment_id         CHAR(13) FOR BIT DATA NOT NULL,

-- ID of the logiacal file name
      lfn_cdb            INTEGER NOT NULL WITH DEFAULT 1,
      lfn_id             CHAR(13) FOR BIT DATA NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT seglfnmap_pk
      PRIMARY KEY (seg_lfn_map_id, creator_db),

-- This map must be created by some process
      CONSTRAINT seglfnmap_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE,

-- It must point to a segment
      CONSTRAINT seglfnmap_fk_sid
      FOREIGN KEY (segment_id, segment_cdb)
          REFERENCES segment(segment_id, creator_db),

-- It must point to a definition
      CONSTRAINT seglfnmap_fk_sdid
      FOREIGN KEY (lfn_id, lfn_cdb)
          REFERENCES lfn(lfn_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on segment ID
CREATE INDEX seglfnmap_sid on segment_lfn_map(segment_id, segment_cdb)
;
-- Create an index based on lfn ID
CREATE INDEX seglfnmap_lfnid on segment_lfn_map(lfn_id, lfn_cdb)
;
