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

CREATE TABLE runlist
(
-- This table contains the list of data runs.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Site (H for Hanford, L for Livingston, etc.)
-- If more than one DAQS system operates at a site, use different codes.
      site               CHAR(4) NOT NULL,

-- Run number
      run                INTEGER NOT NULL,

-- Time range (GPS seconds)
      start_time         INTEGER NOT NULL,
      end_time           INTEGER,

-- Optional additional information about this run
      run_type           VARCHAR(64),
      comment            VARCHAR(1024),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT runlist_pk
      PRIMARY KEY (site, run)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on site and time
CREATE INDEX runlist_ind_time ON runlist(site, start_time)
;
