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

CREATE TABLE search_summvars
(
-- This table contains search-specific summary variables in the form of
-- name/value pairs.  Any given search can create an arbitrary number of
-- entries in this table.
-- Created by Peter, 21 Feb 2002

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH RAN THIS SEARCH
-- Process which generated this event
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Unique id for this row
      search_summvar_id  CHAR(13) FOR BIT DATA NOT NULL,

-- Name/value pairs.  The value can be either a string or a number (expressed
-- as a double-precision real number, even if the value is an integer).
-- To do this, two columns are provided; just fill the appropriate one and
-- leave the other one blank.
      name               VARCHAR(64) NOT NULL,
      string             VARCHAR(256),
      value              DOUBLE,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT s_summvar_pk
      PRIMARY KEY (search_summvar_id, creator_db),

-- Require this to correspond to an entry in the search_summary table
      CONSTRAINT s_summvar_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES search_summary(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on the process which created this entry
CREATE INDEX s_summvar_ind_pid ON search_summvars(process_id, creator_db)
;
-- Create an index based on name
CREATE INDEX s_summvar_ind_name ON search_summvars(name)
;
