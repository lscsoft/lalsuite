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

CREATE TABLE sim_type_params
(
-- Parameters for the simulation types from sim_type table

-- Foreign key refering to sim_type.type
      sim_type          INTEGER NOT NULL, 
-- Parameter name
      name               VARCHAR(24) NOT NULL, 
-- Parameter description
      comment            VARCHAR(64), 
-- Parameter value
      value              DOUBLE,

      CONSTRAINT sim_type_par_pk
      PRIMARY KEY (sim_type, name),

      CONSTRAINT sim_type_par_fk_st 
      FOREIGN KEY (sim_type)
          REFERENCES sim_type(sim_type)
          ON DELETE CASCADE
)
-- The following line is needed for this table to be replicated to other sites
DATA CAPTURE CHANGES
;

