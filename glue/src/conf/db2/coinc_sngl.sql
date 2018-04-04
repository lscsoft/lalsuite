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

CREATE TABLE coinc_sngl
(
-- List of approximate coincidences between single-interferometer events
-- from different interferometers and/or of different event types.
-- Currently, only includes LIGO interferometers plus information about
-- an associated gamma-ray burst, if any.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS ENTRY
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Unique identifier for this coincidence
      coinc_id           CHAR(13) FOR BIT DATA NOT NULL,

-- Single-interferometer events which make up this coincidence (NULL if no
-- match).  Note that we must specify the "creator_db" for each individual
-- event, since it may not be the same as the database on which this
-- coincidence record is being created.
      h1_inspiral_id     CHAR(13) FOR BIT DATA,
      h1_inspiral_cdb    INTEGER,
      h1_burst_id        CHAR(13) FOR BIT DATA,
      h1_burst_cdb       INTEGER,
      h1_ringdown_id     CHAR(13) FOR BIT DATA,
      h1_ringdown_cdb    INTEGER,
      h1_unmodeled_id    CHAR(13) FOR BIT DATA,
      h1_unmodeled_cdb   INTEGER,
      h1_gdstrig_id      CHAR(13) FOR BIT DATA,
      h1_gdstrig_cdb     INTEGER,

      h2_inspiral_id     CHAR(13) FOR BIT DATA,
      h2_inspiral_cdb    INTEGER,
      h2_burst_id        CHAR(13) FOR BIT DATA,
      h2_burst_cdb       INTEGER,
      h2_ringdown_id     CHAR(13) FOR BIT DATA,
      h2_ringdown_cdb    INTEGER,
      h2_unmodeled_id    CHAR(13) FOR BIT DATA,
      h2_unmodeled_cdb   INTEGER,
      h2_gdstrig_id      CHAR(13) FOR BIT DATA,
      h2_gdstrig_cdb     INTEGER,

      l1_inspiral_id     CHAR(13) FOR BIT DATA,
      l1_inspiral_cdb    INTEGER,
      l1_burst_id        CHAR(13) FOR BIT DATA,
      l1_burst_cdb       INTEGER,
      l1_ringdown_id     CHAR(13) FOR BIT DATA,
      l1_ringdown_cdb    INTEGER,
      l1_unmodeled_id    CHAR(13) FOR BIT DATA,
      l1_unmodeled_cdb   INTEGER,
      l1_gdstrig_id      CHAR(13) FOR BIT DATA,
      l1_gdstrig_cdb     INTEGER,

-- Interferometer coincidence time in GPS seconds/nanoseconds
      coinc_time         INTEGER NOT NULL,
      coinc_time_ns      INTEGER NOT NULL,

-- Variable describing the quality of the coincidence
      coinc_quality      REAL NOT NULL,

-- Direction of Hanford-to-Livingston ray at time of event
-- (i.e. the central axis of the cone on which the source lies)
      ligo_axis_ra       REAL,
      ligo_axis_dec      REAL,
-- Wave arrival angle with respect to Hanford-to-Livingston ray, and error
      ligo_angle         REAL,
      ligo_angle_sig     REAL,

-- Gamma-ray burst event, if any
      grb_id             VARCHAR(64),
      grb_time           INTEGER,
      grb_time_ns        INTEGER,
-- Location of gamma-ray burst in the sky, if applicable
      grb_sky_ra         REAL,
      grb_sky_dec        REAL,

-- Place to indicate any other non-LIGO event in coincidence, if any
      other_external     VARCHAR(80),

-- PHYSICAL PROPERTIES FOR EVENT
-- Masses of inspiral objects
      inspiral_mass1     REAL,
      inspiral_mass2     REAL,
-- Q value for ringdown
      ringdown_q         REAL,
-- Fundamental ringdown frequency
      ringdown_freq      REAL,
-- Black hole mass from ringdown
      ringdown_mass      REAL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT coincs_pk
      PRIMARY KEY (coinc_id, creator_db),

      CONSTRAINT coincs_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db),

-- Foreign-key constraints are checked only for non-NULL values
      CONSTRAINT coincs_fk_h1in
      FOREIGN KEY (h1_inspiral_id, h1_inspiral_cdb)
          REFERENCES sngl_inspiral(event_id, creator_db),
      CONSTRAINT coincs_fk_h1bu
      FOREIGN KEY (h1_burst_id, h1_burst_cdb)
          REFERENCES sngl_burst(event_id, creator_db),
      CONSTRAINT coincs_fk_h1ri
      FOREIGN KEY (h1_ringdown_id, h1_ringdown_cdb)
          REFERENCES sngl_ringdown(event_id, creator_db),
      CONSTRAINT coincs_fk_h1un
      FOREIGN KEY (h1_unmodeled_id, h1_unmodeled_cdb)
          REFERENCES sngl_unmodeled(event_id, creator_db),
      CONSTRAINT coincs_fk_h1gds
      FOREIGN KEY (h1_gdstrig_id, h1_gdstrig_cdb)
          REFERENCES gds_trigger(event_id, creator_db),

      CONSTRAINT coincs_fk_h2in
      FOREIGN KEY (h2_inspiral_id, h2_inspiral_cdb)
          REFERENCES sngl_inspiral(event_id, creator_db),
      CONSTRAINT coincs_fk_h2bu
      FOREIGN KEY (h2_burst_id, h2_burst_cdb)
          REFERENCES sngl_burst(event_id, creator_db),
      CONSTRAINT coincs_fk_h2ri
      FOREIGN KEY (h2_ringdown_id, h2_ringdown_cdb)
          REFERENCES sngl_ringdown(event_id, creator_db),
      CONSTRAINT coincs_fk_h2un
      FOREIGN KEY (h2_unmodeled_id, h2_unmodeled_cdb)
          REFERENCES sngl_unmodeled(event_id, creator_db),
      CONSTRAINT coincs_fk_h2gds
      FOREIGN KEY (h2_gdstrig_id, h2_gdstrig_cdb)
          REFERENCES gds_trigger(event_id, creator_db),

      CONSTRAINT coincs_fk_l1in
      FOREIGN KEY (l1_inspiral_id, l1_inspiral_cdb)
          REFERENCES sngl_inspiral(event_id, creator_db),
      CONSTRAINT coincs_fk_l1bu
      FOREIGN KEY (l1_burst_id, l1_burst_cdb)
          REFERENCES sngl_burst(event_id, creator_db),
      CONSTRAINT coincs_fk_l1ri
      FOREIGN KEY (l1_ringdown_id, l1_ringdown_cdb)
          REFERENCES sngl_ringdown(event_id, creator_db),
      CONSTRAINT coincs_fk_l1un
      FOREIGN KEY (l1_unmodeled_id, l1_unmodeled_cdb)
          REFERENCES sngl_unmodeled(event_id, creator_db),
      CONSTRAINT coincs_fk_l1gds
      FOREIGN KEY (l1_gdstrig_id, l1_gdstrig_cdb)
          REFERENCES gds_trigger(event_id, creator_db)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX coincs_ind_time ON coinc_sngl(coinc_time)
;
-- Create an index based on process_id
CREATE INDEX coincs_ind_pid ON coinc_sngl(process_id, coinc_time)
;
