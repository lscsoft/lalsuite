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

CREATE TABLE dbmdctest
(
  -- This table is used by Peter's database MDC scripts.  It should not be
  -- created automatically, but should be created manually when needed.

  seqid      INTEGER NOT NULL,
  tshort     SMALLINT,
  tint       INTEGER,
  tlong      BIGINT,
  treal      REAL,
  tdouble    DOUBLE,
  tchar      CHAR(20),
  tvarchar   VARCHAR(20),
  tcharbin   CHAR(13) FOR BIT DATA,
  tblob      BLOB(32),
  tbigblob   BLOB(1048576),

      CONSTRAINT dbmdctest_pk
      PRIMARY KEY (seqid)
);
