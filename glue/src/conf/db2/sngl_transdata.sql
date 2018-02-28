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

CREATE TABLE sngl_transdata
(
-- Record of transformed data upon which trigger decision was based.  There
-- may be multiple entries for a particular trigger instance, if desired.

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
-- One word or a few words to indicate what is stored in this table
      transdata_name     VARCHAR(32) NOT NULL,

-- Dimensionality of the data array (1 or 2)
      dimensions         INTEGER NOT NULL,

-- X-axis parameters
-- number of bins
      x_bins             INTEGER NOT NULL,
-- starting and ending values of x-axis
      x_start            DOUBLE NOT NULL,
      x_end              DOUBLE NOT NULL,
-- Units of x axis (e.g. 'GPS seconds', 'Hz', etc.)
      x_units            CHAR(16) NOT NULL,

-- Y-axis parameters (if the data array is 2-dimensional)
-- number of bins
      y_bins             INTEGER,
-- starting and ending values of y-axis
      y_start            DOUBLE,
      y_end              DOUBLE,
-- Units of y axis (e.g. 'GPS seconds', 'Hz', etc.)
      y_units            CHAR(16),

-- Parameters for the data in the array
      data_type          CHAR(16) NOT NULL,
      data_units         CHAR(16) NOT NULL,

-- Transformed data itself is stored in a Binary Large OBject (BLOB).
-- We specify COMPACT since we do not expect this ever to be updated.
      transdata          BLOB(1M) COMPACT NOT NULL,
-- Length of data, in bytes
      transdata_length   INTEGER NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT s_trans_pk
      PRIMARY KEY (event_id, creator_db, transdata_name),

-- Foreign key relationship to process table.
      CONSTRAINT s_trans_fk_pid
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
CREATE INDEX s_trans_cind ON sngl_transdata(event_table) CLUSTER
;
-- Create an index based on process_id
CREATE INDEX s_trans_ind_pid ON sngl_transdata(process_id)
;
