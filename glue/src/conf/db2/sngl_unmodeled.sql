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

CREATE TABLE sngl_unmodeled
(
-- Event table for searches for "unmodeled" sources, i.e. sources for which
-- the waveform is unknown.  There will probably be several algorithms in
-- use, so this table includes the algorithm name.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS EVENT
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,
-- Brief keyword to identify the search technique, e.g. "power"
      search             VARCHAR(24) NOT NULL,
-- Filter identifier (indicates type of filter, plus parameters).  May be null
      filter_id          CHAR(13) FOR BIT DATA,
-- Interferometer
      ifo                CHAR(2) NOT NULL,
-- Channel that was analyzed
      channel            VARCHAR(64),

-- TIME OF THE EVENT
-- The start time of this event (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
-- The time duration of this event (seconds)
      duration           REAL NOT NULL,

-- PROPERTIES OF THE EVENT
-- Absolute signal amplitude (fractional strain)
      amplitude          REAL NOT NULL,
-- Signal to noise ratio.
      snr                REAL,
-- Confidence variable
      confidence         REAL,
-- Note: additional properties may be recorded in the sngl_unmodeled_v table.

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT s_unmod_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT s_unmod_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE,

-- Note that filter_id is allowed to be null, in which case no check is made.
      CONSTRAINT s_unmod_fk_filt
      FOREIGN KEY (filter_id, creator_db)
          REFERENCES filter(filter_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index based on search algorithm
CREATE INDEX s_unmod_cind ON sngl_unmodeled(search) CLUSTER
;
-- Create an index based on search algorithm plus time and snr
CREATE INDEX s_unmod_ind_sea ON sngl_unmodeled(search, start_time, snr)
;
-- Create an index based on time and snr
CREATE INDEX s_unmod_ind_tim ON sngl_unmodeled(start_time, snr)
;
-- Create an index based on process_id
CREATE INDEX s_unmod_ind_pid ON sngl_unmodeled(process_id, start_time)
;
-- Create an index based on filter_id
CREATE INDEX s_unmod_ind_fid ON sngl_unmodeled(filter_id, start_time)
;
-- Create an SQL trigger so that if a sngl_unmodeled entry is deleted, any
-- associated sngl_datasource, sngl_transdata, and/or sngl_mime entries
-- are deleted too.
-- Must be done this way because there is no foreign-key relationship.
-- Run script sngl_unmodeled_tr_del.sql to create delete trigger for 
-- sngl_datasource, sngl_transdata, and sngl_mime records.
