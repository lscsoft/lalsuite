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

CREATE TABLE sngl_burst
(
-- Event table for single-interferometer burst-event search.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS EVENT
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,
-- Filter identifier (indicates type of filter, plus parameters).  May be null
      filter_id          CHAR(13) FOR BIT DATA,
-- Interferometer
      ifo                CHAR(2) NOT NULL,
-- Brief keyword to identify the search technique, e.g. "template" or "FCT"
      search             VARCHAR(24) NOT NULL,
-- Channel that was analyzed
      channel            VARCHAR(64),

-- TIME OF THE EVENT - THESE QUANTITIES ARE METHOD-DEFINED
-- The start time of this burst event (in GPS seconds and nanoseconds)
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
-- The stop  time of this burst event (in GPS seconds and nanoseconds)
      stop_time          INTEGER,
      stop_time_ns       INTEGER,
-- The time duration of this burst event (seconds), by definition this is the stop time minus the start time
      duration           REAL NOT NULL,

-- PROPERTIES OF THE EVENT - THESE QUANTITIES ARE METHOD-DEFINED
-- Low frequency (start frequency) of event (Hz)
      flow               REAL,
-- High frequency (end frequency) of event (Hz)
      fhigh              REAL,
-- Center of frequency band in which observation is made (Hz)
      central_freq       REAL,
--  Range of frequency observed (Hz), by definition fhigh minus flow
      bandwidth          REAL,
-- Absolute signal amplitude (fractional strain)
      amplitude          REAL NOT NULL,
-- Signal to noise ratio
      snr                REAL,
-- Confidence variable
      confidence         REAL,
-- Chi-squared statistic, and number of degrees of freedom
      chisq              REAL,
      chisq_dof          REAL,
-- time-frequency volume of the event, i.e, the number of pixels
      tfvolume           REAL,
-- strength of the event in calibrated strain root-sum-square (in-band)
      hrss              REAL,
-- time shift of the time series
      time_lag           REAL,

-- Properties of the event filled in by a GENERAL BURST PARAMETER  
-- estimation code - these quantities are NOT METHOD-DEPENDED
-- PRECISE DEFINITION IS PENDING
-- The peak time of this burst event (in GPS seconds and nanoseconds)
      peak_time          INTEGER,
      peak_time_ns       INTEGER,
-- Peak frequency of the event (Hz)
      peak_frequency     REAL,
-- Peak calibrated strain of the event
      peak_strain        REAL,
-- Error in peak time (in seconds)
      peak_time_error    REAL,
-- Error in peak frequency of the event (Hz)
      peak_frequency_error REAL,
-- Error in peak calibrated strain of the event
      peak_strain_error  REAL,

-- Properties of the (M)ost (S)ignificant pixel of the event.
-- These quantities are method-depended and they are useful for
-- trigger generators that want to capture more information about
-- the distribution of power in the event cluster.
-- Quantities are the same as in the TIME/PROPERTIES OF THE EVENT
-- sections above
      ms_start_time      INTEGER,
      ms_start_time_ns   INTEGER,
      ms_stop_time       INTEGER,
      ms_stop_time_ns    INTEGER,
      ms_duration        REAL,
      ms_flow            REAL,
      ms_fhigh           REAL,
      ms_bandwidth       REAL,
      ms_hrss 	         REAL,
      ms_snr             REAL,
      ms_confidence      REAL,

-- This is a set of three name value pairs that get populated by particular
-- trigger generators when they don't fit into the general structure
-- outlined above
      param_one_name	 VARCHAR(64),
      param_one_value	 DOUBLE,
      param_two_name	 VARCHAR(64),
      param_two_value	 DOUBLE,
      param_three_name	 VARCHAR(64),
      param_three_value	 DOUBLE,

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
      
      CONSTRAINT s_burst_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT s_burst_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE,

-- Note that filter_id is allowed to be null, in which case no check is made.
      CONSTRAINT s_burst_fk_filt
      FOREIGN KEY (filter_id, creator_db)
          REFERENCES filter(filter_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX s_burst_ind_tim ON sngl_burst(start_time, snr)
;
-- Create an index based on process_id
CREATE INDEX s_burst_ind_pid ON sngl_burst(process_id, start_time)
;
-- Create an index based on filter_id
CREATE INDEX s_burst_ind_fid ON sngl_burst(filter_id, start_time)
;
-- Create an SQL trigger so that if a sngl_burst entry is deleted, any
-- associated sngl_datasource, sngl_transdata, and/or sngl_mime entries
-- are deleted too.
-- Must be done this way because there is no foreign-key relationship.
-- Run script sngl_burst_tr_del.sql to create delete trigger for 
-- sngl_datasource, sngl_transdata, and sngl_mime records.
