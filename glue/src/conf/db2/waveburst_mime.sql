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

CREATE TABLE waveburst_mime  
(
-- Table with smallest rectangle of pixels containing the cluster

-- Process which generated this entry
      process_id         CHAR(13) FOR BIT DATA NOT NULL, 
-- LIGO/LSC site that created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,
-- Unique identifier for this event
      event_id           CHAR(13) FOR BIT DATA NOT NULL, 
-- Description of BLOB format if different from t_size by f_size array
-- of REAL4 entries
      mimetype            VARCHAR(64), 
-- Dimensions of the rectangle
-- Number of time steps
      t_size             INTEGER, 
-- Number of frequency steps
      f_size             INTEGER, 
-- Rectangle contains original amplitudes of each pixel
      cluster_o          BLOB(102400) NOT LOGGED COMPACT, 
-- Rectangle contains percentile amplitudes of each pixel
      cluster_p          BLOB(102400) NOT LOGGED COMPACT, 
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
      
      CONSTRAINT wbm_fk_wb
      FOREIGN KEY (event_id, creator_db)
          REFERENCES waveburst(event_id,creator_db)
          ON DELETE CASCADE

)   
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;

