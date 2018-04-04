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

CREATE TABLE segment_definer
(
-- List of processes which define segments.  Note that multiple processes
-- can define segments in the same segment_group.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROGRAM WHICH IS DEFINING SEGMENTS
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- INFORMATION ABOUT THE SEGMENT
-- Unique identification string for this segment definition
      segment_def_id     CHAR(13) FOR BIT DATA NOT NULL,

-- The name of the run which this segment corresponds to (S2, A4, M7, etc.)
      run                CHAR(4) NOT NULL,
-- Interferometer(s) for which these segments are meaningful
      ifos               CHAR(12) NOT NULL,
-- Descriptive name for this group of segments (e.g. 'Science', 'Dust')
      name               VARCHAR(128) NOT NULL,
-- Version number for segment group (to allow re-evaluation)
      version            INTEGER NOT NULL,
-- Optional user comment about this segment_group
      comment            VARCHAR(255),

-- OPTIONAL ADDITIONAL STATE VECTOR DATA
-- These should be null if the information does not come from the state vector
      state_vec_major    INTEGER,
      state_vec_minor    INTEGER,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT segdef_pk
      PRIMARY KEY (segment_def_id, creator_db),

      CONSTRAINT segdef_sk
      UNIQUE (run, ifos, name, version),

      CONSTRAINT segdef_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index for quicker scanning of a given segment type
CREATE INDEX segdef_cind ON segment_definer(name, version) CLUSTER
;
-- Create an index based on the run name
CREATE INDEX segdef_run ON segment_definer(run)
;
-- Create an index based on the keys that state_segment needs
CREATE INDEX segdef_ind_value ON segment_definer(ifos,state_vec_major,state_vec_minor)
;
