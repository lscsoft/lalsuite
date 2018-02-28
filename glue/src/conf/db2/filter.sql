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

CREATE TABLE filter
(
-- Table of filter instances used by GDS and astrophysics-search programs.
-- Note that this table should contain an entry for each invocation of the
-- program (similar to the process table).  It is also possible for a single
-- process to use multiple filters with the same name but different
-- parameters.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROGRAM WHICH INITIALIZES THIS FILTER
-- Program name
      program            VARCHAR(64) NOT NULL,
-- Program start time (GPS seconds)
      start_time         INTEGER NOT NULL,
-- Process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- INFORMATION ABOUT THIS FILTER
-- Filter name (e.g. "FFT")
      filter_name        CHAR(64) NOT NULL,
-- Unique identifier for this invocation of the filter
      filter_id          CHAR(13) FOR BIT DATA NOT NULL,
-- Parameter set identifier.  Permits an association between multiple
-- invocations of a filter which use the same set of input parameters.
-- Probably not filled initially, but updated later.
      param_set          INTEGER,

-- User comment which describes the filter
      comment            VARCHAR(240),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT filter_pk
      PRIMARY KEY (filter_id, creator_db),

-- Foreign key relationship to process table.  The 'ON DELETE CASCADE'
-- modifier means that if a row in the process table is deleted, then
-- all its associated filters are deleted too.
      CONSTRAINT filter_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on process ID
CREATE INDEX filter_ind_pid ON filter(process_id)
;
-- Create an index based on program name
CREATE INDEX filter_ind_prog ON filter(program)
;
-- Create an index based on filter name
CREATE INDEX filter_ind_name ON filter(filter_name)
;
-- Create an index based on filter comment
CREATE INDEX filter_ind_comm ON filter(comment)
;
