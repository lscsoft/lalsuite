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

CREATE TABLE multi_inspiral
(
-- Event table for multi-interferometer binary-inspiral search.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS EVENT
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,
-- Filter identifier (indicates type of filter, plus parameters).  May be null
      filter_id          CHAR(13) FOR BIT DATA,
-- Interferometers used for this search
      ifos               CHAR(12) NOT NULL,
-- Brief keyword to identify the search technique, e.g. "Template" or "FCT"
      search             VARCHAR(16) NOT NULL,

-- TIME OF THE EVENT
-- The coalescence time of this inspiral event (GPS seconds and nanoseconds)
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,
-- Time of impulse that would generate a peak filter output at this time
-- (end_time-impulse_time is a fixed quantity for a given template)
      impulse_time       INTEGER,
      impulse_time_ns    INTEGER,

-- PROPERTIES OF THE EVENT
-- Absolute strain amplitude of fitted template at 100 Hz
      amplitude          REAL NOT NULL,
-- Effective distance, assuming optimally oriented binary (in megaparsecs)
      eff_distance       REAL,
-- Coalescence phase angle (radians)
      coa_phase          REAL,

-- Mass of the larger object (in solar mass units)
      mass1              REAL,
-- Mass of the smaller (or equal) object
      mass2              REAL,
-- Chirp mass [ (m1+m2)*eta^(3/5) ] and eta [ m1*m2/(m1+m2)^2 ]
      mchirp             REAL,
      eta                REAL,
-- Chirp-time parameters
      tau0               REAL,
      tau2               REAL,
      tau3               REAL,
      tau4               REAL,
      tau5               REAL,
-- ttotal is defined as tau0+tau2-tau3+tau4-tau5
      ttotal             REAL,

-- Amplitude signal-to-noise ratio
      snr                REAL,
-- Chi-squared statistic, and number of degrees of freedom
      chisq              REAL,
      chisq_dof          INTEGER,
-- Variance of filter output (nominally normalized snr) near time of event
      sigmasq            REAL,

-- Direction of Hanford-to-Livingston ray at time of event
-- (i.e. the central axis of the cone on which the source lies)
      ligo_axis_ra       REAL,
      ligo_axis_dec      REAL,
-- Wave arrival angle with respect to Hanford-to-Livingston ray, and error
      ligo_angle         REAL,
      ligo_angle_sig     REAL,

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT m_inspiral_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT m_inspiral_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE,

-- Note that filter_id is allowed to be null, in which case no check is made.
      CONSTRAINT m_inspiral_fk_filt
      FOREIGN KEY (filter_id, creator_db)
          REFERENCES filter(filter_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create a clustering index based on search algorithm
CREATE INDEX m_inspiral_cind ON multi_inspiral(search) CLUSTER
;
-- Create an index based on search algorithm plus time and snr
CREATE INDEX m_inspiral_ind_sea ON multi_inspiral(search, end_time, snr)
;
-- Create an index based on time and snr
CREATE INDEX m_inspiral_ind_tim ON multi_inspiral(end_time, snr)
;
-- Create an index based on process_id
CREATE INDEX m_inspiral_ind_pid ON multi_inspiral(process_id, end_time)
;
-- Create an index based on filter_id
CREATE INDEX m_inspiral_ind_fid ON multi_inspiral(filter_id, end_time)
;
