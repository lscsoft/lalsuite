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

echo table process;
select * from process where username in ('U_953671340','U_953671341');
echo table process_params;
select * from process_params where param in ('Param_953671340', 'Param_953671341' );
echo table filter;
select * from filter where param_set in (953671340, 953671341);
echo table filter_params;
select * from filter_params where param in ('param_953671340', 'param_953671341');
echo table sngl_inspiral;
select * from sngl_inspiral where end_time in (637706554, 637706555, 637706556);
echo table sngl_datasource;
select * from sngl_datasource where start_time in (637706553, 637706554, 637706555);
