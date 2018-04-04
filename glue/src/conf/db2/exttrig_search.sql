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

CREATE TABLE exttrig_search
(
-- Lists individual segment-results associated with external triggers
-- for use in statistical analysis

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- INTERFEROMETERS USED FOR THIS EVENT
-- Interferometers used (Each IFO has 2-chars ID separated by ' ')
      ifos               VARCHAR(18),

-- INFORMATION ABOUT THE EXTERNAL TRIGGER
-- Type of external trigger (grb, supernova, cosmic ray, etc.)
      exttrig_type       VARCHAR(32) NOT NULL,
-- ID of this trigger (e.g. the "event_id" from the external_trigger table)
      exttrig_id         VARCHAR(40),
-- Start time of external trigger (GPS seconds, nanoseconds)
      start_time         INTEGER,
      start_time_ns      INTEGER,
-- Right ascension of source
      source_ra          REAL,
-- Error in right ascension of source
      source_ra_err      REAL,
-- Declination of source
      source_dec         REAL,
-- Error in declination of source
      source_dec_err     REAL,

-- INFORMATION ABOUT THIS SEGMENT
-- 1 for on-source, 0 for off-source
      segment_type       INTEGER NOT NULL,
-- Duration of the segment in seconds
      duration           REAL,
-- Delay between start of external trigger and start of segment @IFO1
-- (as ordered in IFOS)
      ref_delay          REAL,
-- Delay between start of segment @IFO1 and start of segment @IFO2, etc.
-- (as ordered in IFOS)
      lag_1_2            REAL,
      lag_1_3            REAL,
      lag_1_4            REAL,

-- RESULTS OF THE CALCULATION
-- Amplitude of, e.g., cross-correlation statistic
      amplitude          REAL,
-- Signal-to-noise ratio (if appropriate)
      snr                REAL,
-- Percent of data which is considered "bad data"
      percent_bad        REAL,
-- Variance in data segment at each IFO
      var_1              REAL,
      var_2              REAL,
      var_3              REAL,
      var_4              REAL,

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT exttrig_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT exttrig_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX exttrig_ind_time ON exttrig_search(start_time)
;
-- Create an index based on process_id
CREATE INDEX exttrig_ind_pid ON exttrig_search(process_id)
;
