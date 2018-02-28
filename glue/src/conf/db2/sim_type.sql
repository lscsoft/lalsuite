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

CREATE TABLE sim_type
(
-- Classification of simulations. This table is not changed by DSO but by hand.
-- This table should migrate from one database to another together with data.

-- Simulations types are enumerated here automatically by DB2
      sim_type          INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY ( START WITH 0 , INCREMENT BY 1 , NO CACHE ), 
-- Simulation name
      name               VARCHAR(24) NOT NULL, 
-- Simulation description.
      comment            VARCHAR(64),
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,
      
      CONSTRAINT s_type_pk
      PRIMARY KEY (sim_type)
)   

-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;
-- Create an index based on name
CREATE INDEX s_type_ind_name ON sim_type(name)
;

