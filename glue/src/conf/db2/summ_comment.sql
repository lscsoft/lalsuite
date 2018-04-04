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

CREATE TABLE summ_comment
(
-- Table to attach a comment to a particular time interval

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH ADDED THE COMMENT (may be null)
-- Program name
      program            VARCHAR(64),
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA,

-- ORIGIN OF THE COMMENT
-- Name of person who made comment
      submitter          VARCHAR(48) NOT NULL,

-- TIME INTERVAL TO WHICH THIS COMMENT APPLIES
-- Group name for frameset which determined this time interval, if any
      frameset_group     VARCHAR(48),
-- Group and version of segment which determined this time interval, if any
      segment_def_cdb    INTEGER NOT NULL WITH DEFAULT 1, 
      segment_def_id     CHAR(13) FOR BIT DATA,
-- Start and end times (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,

-- COMMENT AND ASSOCIATED INFO
-- Interferometer or site to which the comment applies (if appropriate)
      ifo                CHAR(2),
-- The comment itself
      text               VARCHAR(1000) NOT NULL,
-- Unique identifier for this comment (needed for primary key)
      summ_comment_id    CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT summcomm_pk
      PRIMARY KEY (summ_comment_id, creator_db),

-- Note that process_id is allowed to be null, in which case no check is made.
      CONSTRAINT summcomm_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time the comment refers to
CREATE INDEX summcomm_ind_time ON summ_comment(start_time, end_time)
;
-- Create an index based on frameset_group
CREATE INDEX summcomm_ind_fsg ON summ_comment(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summcomm_ind_sgrp ON summ_comment(segment_def_id)
;
