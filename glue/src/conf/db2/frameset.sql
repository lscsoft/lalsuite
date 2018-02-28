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

CREATE TABLE frameset
(
-- A "frameset" is a set of data frames contained in the same file.  It is
-- the smallest unit of raw data which can be read and analyzed, since
-- dictionary information stored at the beginning of the file is needed to 
-- interpret frames appearing later in the file.

-- Database which created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,

-- Unique process ID of the process which wrote this frameset
      process_id         CHAR(13) FOR BIT DATA NOT NULL,

-- Base name for the group to which this frameset belongs (e.g. 'H2.F')
      frameset_group     CHAR(48) NOT NULL,

-- INFORMATION ABOUT THE CONTENTS OF THIS FRAMESET
-- Database which created chanlist entry (which may be different from the
-- database creating this frameset entry)
      chanlist_cdb       INTEGER NOT NULL,
-- Channel set identifier
      chanlist_id        CHAR(13) FOR BIT DATA,
-- Frameset start and end times, in GPS seconds and nanoseconds.
-- Note that end_time is the time at the END of the last frame
-- included in this frameset.  Thus, these two times always differ
-- by at least the length of one frame.
      start_time         INTEGER NOT NULL,
      start_time_ns      INTEGER NOT NULL,
      end_time           INTEGER NOT NULL,
      end_time_ns        INTEGER NOT NULL,
-- Number of frames in frameset
      n_frames           INTEGER NOT NULL,
-- Number of missing frames within frameset
      missing_frames     INTEGER NOT NULL,
-- Size of the frameset, in bytes
      n_bytes            INTEGER NOT NULL,

-- FILENAME FOR THIS FRAMESET
-- This uniquely identifies the frameset.  No two different framesets
-- can have the same name (but if there are multiple copies of a
-- frameset in different places, they would have the same name).
-- Normally the name will consist of an interferometer code, a GPS time,
-- and a suffix which indicates the frame type, e.g. 'H2-628318531.F'.
      name               VARCHAR(80) NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

-- Note that (frameset_group,start_time) should be sufficient to uniquely
-- identify a frameset, but we include end_time in the primary key to
-- facilitate faster queries.
      CONSTRAINT frameset_pk
      PRIMARY KEY (frameset_group, start_time, end_time),

-- Also create a unique index for the frameset name
      CONSTRAINT frameset_uni_name
      UNIQUE (name),

      CONSTRAINT frameset_fk_fswrit
      FOREIGN KEY (frameset_group, process_id, creator_db)
          REFERENCES frameset_writer(frameset_group, process_id, creator_db),

      CONSTRAINT frameset_fk_chanli
      FOREIGN KEY (chanlist_cdb, chanlist_id)
          REFERENCES frameset_chanlist(creator_db, chanlist_id)
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
