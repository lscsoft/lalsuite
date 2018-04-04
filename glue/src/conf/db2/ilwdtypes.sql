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

CREATE TABLE ilwdtypes
(
-- This is a test table to test insertion of ilwd types.

-- char_s  int_4s  lstring char_u  int_4u  real_4
-- int_2s  int_8s  real_8 int_2u  int_8u 

	type_char_s		CHAR(20),
    type_char_v     VARCHAR(64),
    type_lstring    VARCHAR(255),
	
--	integers  
	type_int_2s		SMALLINT WITH DEFAULT 0,
	type_int_2u		SMALLINT WITH DEFAULT 0,
	type_int_4s		INTEGER WITH DEFAULT 0,
	type_int_4u		INTEGER WITH DEFAULT 0,
	type_int_8s		BIGINT WITH DEFAULT 0,
	type_int_8u		BIGINT WITH DEFAULT 0,
	
--	unique Id or blob
	type_char_u    	CHAR(13) FOR BIT DATA,
	type_blob		BLOB(512K) COMPACT,
	type_clob		CLOB(512K) COMPACT,
	
--	double
	type_real_4		REAL WITH DEFAULT 0,
	type_real_8		DOUBLE WITH DEFAULT 0
);

grant select on table ilwdtypes to user ldasdbro; 	


	
	

     
