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

CREATE TABLE summ_spectrum
(
-- Table to contain a summary spectrum derived from the data in a specific
-- time interval.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,
      
-- Unique ID
      summ_spectrum_id   CHAR(13) FOR BIT DATA NOT NULL,

-- INFORMATION ABOUT THE PROCESS WHICH PRODUCED THIS SPECTRUM
-- Program name
      program            VARCHAR(64) NOT NULL,
-- Unique process ID
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- TIME INTERVAL FROM WHICH THIS SPECTRUM WAS DERIVED
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
-- Number of frames actually used to create spectrum
      frames_used        INTEGER,
-- Start and end frequencies (stop frequency is derived from these)
      start_frequency    DOUBLE NOT NULL,
      delta_frequency    DOUBLE NOT NULL,
      mimetype        	 VARCHAR(64) NOT NULL,
	  
-- CHANNEL (OR PSEUDO-CHANNEL) NAME
-- The channel/pseudo-channel name should indicate the interferometer or site
      channel            VARCHAR(240) NOT NULL,

-- SPECTRUM DATA AND ASSOCIATED INFO
-- Spectrum type (descriptive name)
      spectrum_type      VARCHAR(128) NOT NULL,
-- The spectrum itself is stored in a Binary Large OBject (BLOB).
-- We specify COMPACT since we do not expect this ever to be updated.
      spectrum           BLOB(1M) COMPACT NOT NULL,
-- Length of the spectrum (in bytes)
      spectrum_length    INTEGER NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT summspect_pk
      PRIMARY KEY (summ_spectrum_id, creator_db),

      CONSTRAINT summspect_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index based on spectrum_type
CREATE INDEX summspect_cind ON summ_spectrum(spectrum_type) CLUSTER
;
-- Create an index based on channel name
CREATE INDEX summspect_chan ON summ_spectrum(channel)
;
-- Create an index based on time
CREATE INDEX summspect_ind_time ON summ_spectrum(start_time, end_time)
;
-- Create an index based on process_id
CREATE INDEX summspect_ind_pid ON summ_spectrum(process_id)
;
-- Create an index based on frameset_group
CREATE INDEX summspect_ind_fsg ON summ_spectrum(frameset_group)
;
-- Create an index based on segment_group
CREATE INDEX summspect_ind_sgrp ON summ_spectrum(segment_def_id)
;
