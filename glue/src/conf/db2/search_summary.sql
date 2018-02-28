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

CREATE TABLE search_summary
(
-- This table contains information about the search proved by the dso
-- First version distributed by Duncan, 14 Feb 2002
-- Modified by Peter, 21 Feb 2002

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH RAN THIS SEARCH
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- BASIC INFORMATION ABOUT THE SHARED OBJECT
-- Shared object name
      shared_object      VARCHAR(64) NOT NULL,
-- LALwrapper CVS tag
      lalwrapper_cvs_tag VARCHAR(128) NOT NULL,
-- LAL CVS tag
      lal_cvs_tag        VARCHAR(128) NOT NULL,
-- User comment
      comment            VARCHAR(240),

-- Interferometers used
      ifos               VARCHAR(12),

-- INFORMATION ABOUT THE DATA RECIVED BY THE SHARED OBJECT
-- The start and stop time of the data passed to the search by the datacond
-- Input data start and end times, in GPS seconds and nanoseconds.
      in_start_time      INTEGER NOT NULL,
      in_start_time_ns   INTEGER NOT NULL,
      in_end_time        INTEGER NOT NULL,
      in_end_time_ns     INTEGER NOT NULL,

-- INFORMATION ABOUT THE DATA FILTERED BY THE SHARED OBJECT
-- The start and stop times of the data that are used to create events
-- Output data start and end times, in GPS seconds and nanoseconds.
      out_start_time     INTEGER NOT NULL,
      out_start_time_ns  INTEGER NOT NULL,
      out_end_time       INTEGER NOT NULL,
      out_end_time_ns    INTEGER NOT NULL,

-- SUMMARY VARIABLES WHICH ARE MORE-OR-LESS COMMON TO ALL SEARCHES
-- Number of events generated
      nevents            INTEGER,
-- Number of Beowulf nodes used
      nnodes             INTEGER,
-- What else??

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- The "primary key" must be unique, so the definition below implies that there
-- can be at most one search_summary entry per process_id entry
      CONSTRAINT s_summary_pk
      PRIMARY KEY (process_id, creator_db),

-- Require this to correspond to an entry in the process table
      CONSTRAINT s_summary_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index based on shared object name
CREATE INDEX s_summary_cind ON search_summary(shared_object) CLUSTER
;
-- Create an index based on input time
CREATE INDEX s_summary_ind_tin ON search_summary(in_start_time, in_end_time)
;
-- Create an index based on output time
CREATE INDEX s_summary_ind_tout ON search_summary(out_start_time, out_end_time)
;
-- Create an index based on interferometers
CREATE INDEX s_summar_ind_ifos ON search_summary(ifos)
;
