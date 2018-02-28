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

CREATE TABLE sngl_mime
(
-- Table to store arbitrary binary data, identified by a MIME type,
-- associated with a given trigger/event.  There may be multiple entries
-- for a particular trigger instance, if desired.

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

-- ORIGIN OF THE DATA
-- Generator of the data (program name, 'video camera', etc.)
      origin             VARCHAR(64),
-- Original filename (optional)
      filename           VARCHAR(64),
-- Name of person who entered data into database
      submitter          VARCHAR(48),

-- BINARY DATA AND ASSOCIATED INFO
-- The data itself is stored in a Binary Large OBject (BLOB).
-- We specify COMPACT since we do not expect this ever to be updated.
      mimedata           BLOB(1M) COMPACT NOT NULL,
-- Length of the binary data (in bytes)
      mimedata_length    INTEGER NOT NULL,
-- MIME content-type (e.g. 'image/jpeg')
      mimetype           VARCHAR(64) NOT NULL,
-- Brief description of the contents, e.g. "power spectrum"
      descrip            VARCHAR(64),
-- Optional comment about this data
      comment            VARCHAR(240),
-- Unique identifier for this data (needed for primary key)
      sngl_mime_id       CHAR(13) FOR BIT DATA NOT NULL,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT s_mime_pk
      PRIMARY KEY (event_id, creator_db, sngl_mime_id),

-- Foreign key relationship to process table.
      CONSTRAINT s_mime_fk_pid
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
CREATE INDEX s_mime_cind ON sngl_mime(event_table) CLUSTER
;
-- Create an index based on descrip
CREATE INDEX s_mime_ind_desc ON sngl_mime(descrip)
;
-- Create an index based on process_id
CREATE INDEX s_mime_ind_pid ON sngl_mime(process_id)
;
