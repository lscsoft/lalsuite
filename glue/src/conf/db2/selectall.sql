-- 
-- select count from all tables 
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

-- misc tables
echo select count(*) from calib_info;
select count(*) from calib_info for read only ;
echo select count(*) from runlist;
select count(*) from runlist for read only ;
echo select count(*) from search_summvars;
select count(*) from search_summvars for read only ;
echo select count(*) from search_summary;
select count(*) from search_summary for read only ;

-- multi events
echo select count(*) from exttrig_search;
select count(*) from exttrig_search for read only ;
echo select count(*) from coinc_sngl;
select count(*) from coinc_sngl for read only ;
echo select count(*) from multi_burst;
select count(*) from multi_burst for read only ;
echo select count(*) from multi_inspiral;
select count(*) from multi_inspiral for read only ;

-- single events 
echo select count(*) from waveburst_mime;
select count(*) from waveburst_mime for read only ;
echo select count(*) from waveburst;
select count(*) from waveburst for read only ;
echo select count(*) from gds_trigger;
select count(*) from gds_trigger for read only ;
echo select count(*) from sngl_burst;
select count(*) from sngl_burst for read only ;
echo select count(*) from sngl_block;
select count(*) from sngl_block for read only ;
echo select count(*) from sngl_dperiodic;
select count(*) from sngl_dperiodic for read only ;
echo select count(*) from sngl_inspiral;
select count(*) from sngl_inspiral for read only ;
echo select count(*) from sngl_ringdown;
select count(*) from sngl_ringdown for read only ;
echo select count(*) from sngl_unmodeled;
select count(*) from sngl_unmodeled for read only ;
echo select count(*) from sngl_unmodeled_v;
select count(*) from sngl_unmodeled_v for read only ;

echo select count(*) from summ_mime;
select count(*) from summ_mime for read only ;
echo select count(*) from summ_comment;
select count(*) from summ_comment for read only ;
echo select count(*) from summ_spectrum;
select count(*) from summ_spectrum for read only ;
echo select count(*) from summ_csd;
select count(*) from summ_csd for read only ;
echo select count(*) from summ_statistics;
select count(*) from summ_statistics for read only ;
echo select count(*) from summ_value;
select count(*) from summ_value for read only ;

echo select count(*) from sngl_mime;
select count(*) from sngl_mime for read only ;
echo select count(*) from sngl_transdata;
select count(*) from sngl_transdata for read only ;
echo select count(*) from sngl_datasource;
select count(*) from sngl_datasource for read only ;

echo select count(*) from frameset_loc;
select count(*) from frameset_loc for read only ;
echo select count(*) from frameset;
select count(*) from frameset for read only ;
echo select count(*) from frameset_chanlist;
select count(*) from frameset_chanlist for read only ;
echo select count(*) from frameset_writer;
select count(*) from frameset_writer for read only ;

echo select count(*) from segment;
select count(*) from segment for read only ;
echo select count(*) from segment_definer;
select count(*) from segment_definer for read only ;

-- simulation tables
echo select count(*) from sim_inst_params;
select count(*) from sim_inst_params for read only ;
echo select count(*) from sim_inst;
select count(*) from sim_inst for read only ;
echo select count(*) from sim_type_params;
select count(*) from sim_type_params for read only ;
echo select count(*) from sim_type;
select count(*) from sim_type for read only ;

-- filter and params
echo select count(*) from filter;
select count(*) from filter for read only ;
echo select count(*) from filter_params;
select count(*) from filter_params for read only ;
echo select count(*) from process_params;
select count(*) from process_params for read only ;

-- process table last
echo select count(*) from process;
select count(*) from process for read only ;



