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

CREATE TABLE summ_value
(
-- Table to record a value about a particular time interval

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique ID
      summ_value_id      CHAR(13) FOR BIT DATA NOT NULL,

-- INFORMATION ABOUT THE PROCESS WHICH RECORDED THE VALUE
-- Program name
      program            VARCHAR(64) NOT NULL,
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- TIME INTERVAL FROM WHICH THIS VALUE WAS CALCULATED
-- Group name for frameset which determined this time interval, if any
      frameset_group     VARCHAR(48),
-- segment type which determined this time interval, if any
      segment_def_cdb    INTEGER NOT NULL WITH DEFAULT 1,
      segment_def_id     CHAR(13) FOR BIT DATA,
-- Start and end times (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,

-- THE SUMMARY VALUE
-- Site or interferometer to which this applies (H0, H1, H2, L0, L1)
      ifo                CHAR(2) NOT NULL,
-- Descriptive name
      name               VARCHAR(128) NOT NULL,
-- The value itself (must be a real number)
      value              REAL,
-- Optional uncertainty on the value
      error              REAL,
-- An optional 4-byte integer value or bitmask
      intvalue           INTEGER,
-- Optional comment
      comment            VARCHAR(80),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
  
      CONSTRAINT summ_value_pk
      PRIMARY KEY (summ_value_id, creator_db),

      CONSTRAINT summval_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX summval_ind_time ON summ_value(start_time, end_time)
;
-- Create an index based on program name
CREATE INDEX summval_ind_prog ON summ_value(program, start_time, name)
;
-- Create an index based on process_id
CREATE INDEX summval_ind_pid ON summ_value(process_id)
;
-- Create an index based on frameset_group
CREATE INDEX summval_ind_fsg ON summ_value(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summval_ind_sgrp ON summ_value(segment_def_id)
;
