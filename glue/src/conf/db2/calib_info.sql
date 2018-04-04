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

CREATE TABLE calib_info
(
-- Table to store calibration information in binary form

-- Time interval for which this record is valid, in GPS seconds
      valid_start        INTEGER NOT NULL,
      valid_end          INTEGER NOT NULL,

-- Information about the origin of this calibration information
-- (e.g. name of the program which produced it)
      origin             VARCHAR(128),
-- Generation time in GPS seconds
      origin_time        INTEGER NOT NULL,

-- Information about the insertion of this entry into the database
-- Note: in an SQL query, you can use a 'YYYYMMDDhhmmss' date/time string (UTC)
-- as follows:     " ... where insert_time < timestamp('20020604181522') "
      insertion_time        TIMESTAMP NOT NULL WITH DEFAULT CURRENT TIMESTAMP,
-- Name of the person (or automated program) who inserted the entry into DB2
      submitter          VARCHAR(48) NOT NULL,

-- Channel name
      channel            VARCHAR(64) NOT NULL,
-- Physical units that channel counts are converted into,
-- e.g. "strain" or "m/s"
      units              VARCHAR(64) NOT NULL,
-- Calibration type bitmask (1=amplitude,2=offset,4=time delay,
--    8=transfer function, 16=pole/zero)
      caltype            INTEGER NOT NULL,
-- Byte string containing the binary data structure
      caldata            BLOB(1M) COMPACT NOT NULL,

-- Optional comment
      comment            VARCHAR(256),

      CONSTRAINT calinfo_pk
      PRIMARY KEY (channel, units, valid_start, valid_end, origin_time )
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
