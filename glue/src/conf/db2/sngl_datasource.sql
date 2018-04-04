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

CREATE TABLE sngl_datasource
(
-- Pointer to specific data which prompted a single_interferometer trigger.
-- That is, this indicates what data the program was looking at when it
-- decided to generate the trigger.  In general, this datasource time
-- interval will be longer than the duration of the transient it contains.
-- Note that there can be only one sngl_datasource entry per trigger, but
-- it can list multiple channels.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Table with event to which this applies (gds_trigger, sngl_inspiral, etc.)
      event_table        CHAR(18) NOT NULL,
-- Trigger/event identifier
      event_id           CHAR(13) FOR BIT DATA NOT NULL,
-- Site or interferometer to which this applies (H0, H1, H2, L0, L1)
      ifo                CHAR(2) NOT NULL,

-- Source of this data.  Use 'DAQS' for original raw data, otherwise the
-- name of the frameset_group read.  If multiple frameset_groups were read,
-- list them all, separated by spaces.
      data_source        VARCHAR(240),
-- Channel(s) the data come from.  If more than one, list them all,
-- separated by spaces.
      channels           VARCHAR(512) NOT NULL,

-- The beginning of the time interval for the data (GPS seconds/nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
-- The end of the time interval for the data (GPS seconds/nanoseconds)
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT s_dsource_pk
      PRIMARY KEY (event_id, creator_db),

-- Foreign key relationship to process table.
      CONSTRAINT s_dsource_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE

-- We cannot set up a foreign key based on event_id since the parent
-- table varies.  Instead, we set up a trigger to do the equivalent check.
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index based on event_table
CREATE INDEX s_dsource_cind ON sngl_datasource(event_table) CLUSTER
;
-- Create an index based on process_id
CREATE INDEX s_dsource_ind_pid ON sngl_datasource(process_id)
;
