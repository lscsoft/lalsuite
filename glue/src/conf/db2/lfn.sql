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

CREATE TABLE lfn
(
-- An entry in the process table may have one or more logial file names (LFNs)
-- associated with it. These LFNs can be used in a query to a replica location
-- service (RLS) to obtain the physical file names (PFNs) associated with the
-- entry in the process table. A process to LFN mapping is typically used to 
-- allow storage of information that is too large to be kept in the database
-- (e.g. raw trigger files)

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique process ID of the process which wrote this lfn entry
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Unique identifier for this LFN
      lfn_id             CHAR(13) FOR BIT DATA NOT NULL,

-- logical file name
      name               VARCHAR(255) NOT NULL,

-- optional information, if this file conforms to the frame file syntax
      comment            VARCHAR(240),
      start_time         INTEGER,
      end_time           INTEGER,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- LFN IDs must be unique
      CONSTRAINT lfn_pk
      PRIMARY KEY (lfn_id,creator_db),

      CONSTRAINT lfn_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
CREATE INDEX lfn_ind_name ON lfn(name)
;
