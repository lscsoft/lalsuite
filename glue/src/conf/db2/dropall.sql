-- 
-- delete all tables in proper order
--
-- to run, first connect to db2 via command line connect 
-- e.g. db2 connect to cit_1 user db2sol7s using <password>
--
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

echo dropping tables;

-- misc tables
echo table calib_info;
drop table calib_info;
echo table runlist;
drop table runlist;
echo table search_summvars;
drop table search_summvars;
echo table search_summary;
drop table search_summary;

-- multi events
echo table exttrig_search;
drop table exttrig_search;
echo table coinc_sngl;
drop table coinc_sngl;
echo table multi_burst;
drop table multi_burst;
echo table multi_inspiral;
drop table multi_inspiral;

-- single events 
echo table waveburst_mime;
drop table waveburst_mime;
echo table waveburst;
drop table waveburst;
echo table gds_trigger;
drop table gds_trigger;
echo table sngl_burst;
drop table sngl_burst;
echo table sngl_block;
drop table sngl_block;
echo drop table sngl_dperiodic;
drop table sngl_dperiodic;
echo drop table sngl_inspiral;
drop table sngl_inspiral;
echo drop table sngl_ringdown;
drop table sngl_ringdown;
echo drop table sngl_unmodeled;
drop table sngl_unmodeled;
echo drop table sngl_unmodeled_v;
drop table sngl_unmodeled_v;

echo drop table summ_mime;
drop table summ_mime;
echo drop table summ_comment;
drop table summ_comment;
echo drop table summ_spectrum;
drop table summ_spectrum;
echo drop table summ_csd;
drop table summ_csd;
echo drop table summ_statistics;
drop table summ_statistics;
echo drop table summ_value;
drop table summ_value;

echo table sngl_mime;
drop table sngl_mime;
echo table sngl_transdata;
drop table sngl_transdata;
echo table sngl_datasource;
drop table sngl_datasource;

echo table frameset_loc;
drop table frameset_loc;
echo table frameset;
drop table frameset;
echo table frameset_chanlist;
drop table frameset_chanlist;
echo table frameset_writer;
drop table frameset_writer;

echo table segment;
drop table segment;
echo table segment_definer;
drop table segment_definer;

-- simulation tables
echo table sim_inst_params;
drop table sim_inst_params;
echo table sim_inst;
drop table sim_inst;
echo table sim_type_params;
drop table sim_type_params;
echo table sim_type;
drop table sim_type;

-- filter and params
echo table filter;
drop table filter;
echo table filter_params;
drop table filter_params;
echo table process_params;
drop table process_params;

-- drop process table last
echo table process;
drop table process;



