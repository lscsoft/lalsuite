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

CREATE TABLE sim_inst_params
(
-- Parameters of the simulation instance
	
-- Simulation instance id
      simulation_id               CHAR(13) FOR BIT DATA NOT NULL,
-- LIGO/LSC site that created this entry
      creator_db         INTEGER NOT NULL WITH DEFAULT 1,
-- Parameter name
      name               VARCHAR(24) NOT NULL, 
-- Parameter descrition
      comment            VARCHAR(64),
-- Parameter value
      value              DOUBLE,
      
-- Insertion time (automatically assigned by the database)
      insertion_time     TIMESTAMP WITH DEFAULT CURRENT TIMESTAMP,

      CONSTRAINT sim_in_par_fk_si
      FOREIGN KEY (simulation_id, creator_db)
          REFERENCES sim_inst(simulation_id, creator_db)
          ON DELETE CASCADE
)   
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;


