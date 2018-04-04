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

CREATE TABLE external_trigger
(
-- This table contains information about External triggers 
-- valuable for LIGO such as:
-- Gamma-Ray Bursts, Supernovae, Bar candidates, etc.
-- For further information, please 
-- contact Szabolcs Marka (smarka@ligo.caltech.edu)
-- or visit www.ligo.caltech.edu/~smarka/ExternalTriggers/index.html

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- Event ID (From the alert/parser itself, not assigned by LDAS)
      event_id           VARCHAR(40) NOT NULL,

-- Trigger information
      notice_type        VARCHAR(40) NOT NULL,
      trigger_type       VARCHAR(40) NOT NULL,
      trigger_number     INTEGER,
      trigger_source     VARCHAR(20),
      distributor        VARCHAR(20),

-- Time of the event
--  in GPS seconds/nanoseconds
--  duration (in seconds) 
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      duration           REAL NOT NULL,
      date_of_max        VARCHAR(40),
      time_of_discovery  VARCHAR(40),

-- Size of the event
--  intensity is in mCrab units
      amplitude          REAL NOT NULL,
      intensity          REAL,
      magnitude          REAL,
      trigger_mag        INTEGER NOT NULL,

-- Event/detector sensitivity properties
--  Frequency/Bandwith in Hz
--  Counting rate in counts/sec
      frequency1         REAL NOT NULL,
      frequency2         REAL,
      bandwidth1         REAL NOT NULL,
      bandwidth2         REAL,
      signal_to_noise    REAL NOT NULL,
      local_StoN         REAL,
      confidence         REAL NOT NULL,
      rate               INTEGER,
      max_size           REAL,

-- Event location and error, in decimal degrees
      trigger_error      REAL NOT NULL,
      right_ascension    REAL NOT NULL,
      ra_center          REAL,
      ra_z               REAL,
      ra_error           REAL,
      declination        REAL NOT NULL,
      dec_center         REAL,
      dec_z              REAL,
      dec_error          REAL,

-- Angular Distance to Sun and Moon, in degrees
      sun_position       REAL,
      sun_distance       REAL,
      moon_position      REAL,
      moon_distance      REAL,

-- Host galaxy properties
--  Location in decimal degrees
--  Distance in Mpc
--  Velocity in Km/seconds
--  Photometry data in magnitudes and/or W/m2
--  Apparent size in arcminutes x arcminutes
--  Object offset from galactic center in decimal arcseconds
      galaxy             VARCHAR(40),
      g_type             VARCHAR(20),
      g_ra               REAL,
      g_dec              REAL,
      g_magnitude        REAL,
      g_luminosity       REAL,
      g_distance         REAL,
      g_pos_angle        REAL,
      g_velocity         REAL,
      g_redshift         REAL,
      g_apparent_size    VARCHAR(20),
      g_offset_ew        REAL,
      g_offset_ns        REAL,

-- Comments
      comment            VARCHAR(256),
      comment1           VARCHAR(128),
      comment2           VARCHAR(128),
      comment3           VARCHAR(128),
      comment4           VARCHAR(128),
      comment5           VARCHAR(128),
      comment6           VARCHAR(80),
      comment7           VARCHAR(80),
      comment8           VARCHAR(80),
      comment9           VARCHAR(80),
      comment10          VARCHAR(80),
      comment11          VARCHAR(80),
      comment12          VARCHAR(80),
      comment13          VARCHAR(80),
      comment14          VARCHAR(80),
      comment15          VARCHAR(80),

-- Pointers to extra information/resources (HTML/WEB etc.)
      web                VARCHAR(80),
      web1               VARCHAR(80),
      web2               VARCHAR(80),
      web3               VARCHAR(80),
      web4               VARCHAR(80),
      web5               VARCHAR(80),

-- Extra fields
      int1               INTEGER,
      int2               INTEGER,
      int3               INTEGER,
      real1              REAL,
      real2              REAL,
      real3              REAL
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX ext_trig_time ON external_trigger(start_time)
;
