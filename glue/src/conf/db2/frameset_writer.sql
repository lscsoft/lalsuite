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

CREATE TABLE frameset_writer
(
-- List of processes which create framesets.  Note that multiple processes
-- can write framesets in the same frameset_group.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROGRAM WHICH IS WRITING FRAMESETS
-- Program name
      program            VARCHAR(64) NOT NULL,
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- INFORMATION ABOUT THE FRAMESETS BEING WRITTEN
-- Base name for this group of framesets (e.g. 'H2.F')
      frameset_group     CHAR(48) NOT NULL,
-- Source of this data.  Use 'DAQS' for original raw data, otherwise the
-- name of the frameset_group from which this new frameset_group is derived.
-- If the new frameset_group is derived from multiple frameset_groups, list
-- them all, separated by spaces.
      data_source        VARCHAR(240) NOT NULL,
-- Interferometer(s) with information in the frames
      ifos               CHAR(12),

-- Optional user remark about this frameset_group
      comment             VARCHAR(240),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT fswriter_pk
      PRIMARY KEY (frameset_group, process_id, creator_db),

      CONSTRAINT fswriter_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Also create a clustering index for quicker scanning
CREATE INDEX fswriter_ind_fsgrp ON frameset_writer(frameset_group) CLUSTER;
