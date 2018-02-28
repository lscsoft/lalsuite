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
-- 
echo deleting tables;

-- misc tables
echo table calib_info;
delete from calib_info;
echo table runlist;
delete from runlist;
echo table search_summvars;
delete from search_summvars;
echo table search_summary;
delete from search_summary;

-- multi events
echo table exttrig_search;
delete from exttrig_search;
echo table coinc_sngl;
delete from coinc_sngl;
echo table multi_burst;
delete from multi_burst;
echo table multi_inspiral;
delete from multi_inspiral;

-- single events 
echo table waveburst_mime;
delete from waveburst_mime;
echo table waveburst;
delete from waveburst;
echo table gds_trigger;
delete from gds_trigger;
echo table sngl_burst;
delete from sngl_burst;
echo delete from sngl_block;
delete from sngl_block;
echo delete from sngl_dperiodic;
delete from sngl_dperiodic;
echo delete from sngl_inspiral;
delete from sngl_inspiral;
echo delete from sngl_ringdown;
delete from sngl_ringdown;
echo delete from sngl_unmodeled;
delete from sngl_unmodeled;
echo delete from sngl_unmodeled_v;
delete from sngl_unmodeled_v;

echo delete from summ_mime;
delete from summ_mime;
echo delete from summ_comment;
delete from summ_comment;
echo delete from summ_spectrum;
delete from summ_spectrum;
echo delete from summ_csd;
delete from summ_csd;
echo delete from summ_statistics;
delete from summ_statistics;
echo delete from summ_value;
delete from summ_value;

echo deleting sngl_mime;
delete from sngl_mime;
echo table sngl_transdata;
delete from sngl_transdata;
echo deleting sngl_datasource;
delete from sngl_datasource;

echo table frameset_loc;
delete from frameset_loc;
echo table frameset;
delete from frameset;
echo table frameset_chanlist;
delete from frameset_chanlist;
echo table frameset_writer;
delete from frameset_writer;

echo table segment;
delete from segment;
echo table segment_definer;
delete from segment_definer;

-- simulation tables
echo table sim_inst_params;
delete from sim_inst_params;
echo table sim_inst;
delete from sim_inst;
echo table sim_type_params;
delete from sim_type_params;
echo table sim_type;
delete from sim_type;

-- filter and params
echo table filter;
delete from filter;
echo table filter_params;
delete from filter_params;
echo table process_params;
delete from process_params;

-- delete process table last
echo table process;
delete from process;
