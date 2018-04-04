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

CREATE TABLE ligolw_mon
(
-- Table used to communicate with the DMT LIGOLwMon program.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- INFORMATION ABOUT THE PROCESS WHICH GENERATED THIS DATA
-- Process which generated this trigger
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- TIME OF THE TRIGGER
-- The time at which to display the trigger (GPS seconds and nanoseconds)
      time               INTEGER NOT NULL,
      time_ns            INTEGER NOT NULL,

-- PROPERTIES OF THE TRIGGER
-- amplitude is mandatory
      amplitude          REAL NOT NULL,
-- optional extra information about the event
      confidence         REAL,
      frequency          REAL,

-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL,

-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT ligolw_mon_pk
      PRIMARY KEY (event_id, creator_db),

      CONSTRAINT ligolw_mon_fk_pid
      FOREIGN KEY (process_id, creator_db)
          REFERENCES process(process_id, creator_db)
          ON DELETE CASCADE,
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on time
CREATE INDEX ligolw_mon_ind_t ON ligolw_mon(time)
;
-- Create an index based on time and time_ns
CREATE INDEX ligolw_mon_ind_tns ON ligolw_mon(time, time_ns)
;
