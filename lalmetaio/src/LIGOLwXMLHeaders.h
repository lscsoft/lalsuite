/*
*  Copyright (C) 2007 Andres C. Rodriguez, Sukanta Bose, Alexander Dietz, Duncan Brown, Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Robert Adam Mercer, Saikat Ray-Majumder, Stephen Fairhurst, Xavier Siemens, Sean Seader, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOLwXMLHeaders.h
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <lal/LALStdlib.h>

/**
\author Brown, D. A.
\file
\ingroup lalmetaio

\brief This header provides provides <tt>\#define</tt>s for the common elements of LIGO light weight XML files.

\heading{Synopsis}
\code
#include <lal/LIGOLwXMLHeaders.h>
\endcode

It provides the XML header and footer, as well as table definitions for the various metadata tables.
It will need to be kept up to date with changes in the LIGO database table definitions.  The quantities which are
defined in this file are

<ul>
<li> LIGOLW_XML_HEADER</li>
<li> LIGOLW_XML_FOOTER</li>
<li> LIGOLW_XML_TABLE_FOOTER</li>
<li> LIGOLW_XML_PROCESS</li>
<li> PROCESS_ROW</li>
<li> LIGOLW_XML_PROCESS_PARAMS</li>
<li> PROCESS_PARAMS_ROW</li>
<li> LIGOLW_XML_SEARCH_SUMMARY</li>
<li> SEARCH_SUMMARY_ROW</li>
<li> LIGOLW_XML_SEARCH_SUMMVARS</li>
<li> SEARCH_SUMMVARS_ROW</li>
<li> LIGOLW_XML_SIM_RINGDOWN</li>
<li> SIM_RINGDOWN_ROW</li>
<li> LIGOLW_XML_SUMM_VALUE</li>
<li> SUMM_VALUE_ROW</li>
<li> LIGOLW_XML_SIM_INST_PARAMS</li>
<li> SIM_INST_PARAMS_ROW</li>
<li> LIGOLW_XML_STOCHASTIC</li>
<li> STOCHASTIC_ROW</li>
<li> LIGOLW_XML_STOCH_SUMM</li>
<li> STOCH_SUMM_ROW</li>
<li> LIGOLW_XML_EXT_TRIGGERS</li>
<li> EXT_TRIGGERS_ROW</li>
</ul>

*/

#ifndef _LIGOLWXMLHEADERS_H
#define _LIGOLWXMLHEADERS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#define PRINT_LIGOLW_XML_HEADER(fp) ( \
fputs( "<?xml version='1.0' encoding='utf-8' ?>\n", fp) == EOF || \
fputs( "<!DOCTYPE LIGO_LW SYSTEM \"http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt\">", fp ) == EOF || \
fputs( "<LIGO_LW>\n", fp ) == EOF )


#define PRINT_LIGOLW_XML_FOOTER(fp) \
fputs( "</LIGO_LW>", fp )

#define PRINT_LIGOLW_XML_TABLE_FOOTER(fp) ( \
fputs( "\n", fp ) == EOF || \
fputs( "      </Stream>\n", fp ) == EOF || \
fputs( "   </Table>\n", fp ) )

#define PRINT_LIGOLW_XML_PROCESS(fp) ( \
fputs( "   <Table Name=\"processgroup:process:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:program\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:version\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:cvs_repository\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:cvs_entry_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:is_online\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:node\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:username\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:unix_procid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:jobid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:domain\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:ifos\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"processgroup:process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"processgroup:process:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_ROW \
"         \"%s\",\"%s\",\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",\"process:process_id:0\""

#define PRINT_LIGOLW_XML_PROCESS_PARAMS(fp) ( \
fputs( "   <Table Name=\"process_paramsgroup:process_params:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:program\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:param\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:type\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_paramsgroup:process_params:value\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"process_paramsgroup:process_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_PARAMS_ROW \
"         \"%s\",\"process:process_id:0\",\"%s\",\"%s\",\"%s\""

#define PRINT_LIGOLW_XML_SEARCH_SUMMARY(fp) ( \
fputs( "   <Table Name=\"search_summarygroup:search_summary:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:shared_object\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:lalwrapper_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:lal_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:ifos\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:in_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:out_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:nevents\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summarygroup:search_summary:nnodes\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"search_summarygroup:search_summary:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMARY_ROW \
"         \"process:process_id:0\",\"standalone\",\"\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d"

#define PRINT_LIGOLW_XML_SEARCH_SUMMVARS(fp) (  \
fputs( "   <Table Name=\"search_summvarsgroup:search_summvars:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:name\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:string\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:value\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvarsgroup:search_summvars:search_summvar_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"search_summvarsgroup:search_summvars:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMVARS_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%22.16e,\"search_summvars:search_summvar_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SNGL_RINGDOWN(fp) ( \
fputs( "   <Table Name=\"sngl_ringdowngroup:sngl_ringdown:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:start_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:frequency\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:quality\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:mass\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:spin\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:epsilon\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:num_clust_trigs\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:ds2_H1H2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:ds2_H1L1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:ds2_H1V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:ds2_H2L1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:ds2_H2V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:ds2_L1V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:eff_dist\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:sigma_sq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"sngl_ringdowngroup:sngl_ringdown:event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_ringdowngroup:sngl_ringdown:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_RINGDOWN_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%22.16e,\"sngl_ringdown:event_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_RINGDOWN(fp) ( \
fputs( "   <Table Name=\"sim_ringdowngroup:sim_ringdown:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF ||  \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:coordinates\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:geocent_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:geocent_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:h_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:h_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:l_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:l_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:v_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:v_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:start_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:inclination\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:frequency\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:quality\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:mass\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:spin\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:epsilon\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:hrss\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:hrss_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:hrss_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:hrss_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_ringdowngroup:sim_ringdown:simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_ringdowngroup:sim_ringdown:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_RINGDOWN_ROW \
  "         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"sim_ringdown:simulation_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SUMM_VALUE(fp) ( \
fputs( "   <Table Name=\"summ_valuegroup:summ_value:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:program\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:ifo\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:value\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_valuegroup:summ_value:summ_value_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"summ_valuegroup:summ_value:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SUMM_VALUE_ROW \
"         \"%s\",\"process:process_id:0\",%d,%d,%d,%d,\"%s\",\"%s\",%e,\"%s\",\"summ_value:summ_value_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_INST_PARAMS(fp) ( \
fputs( "   <Table Name=\"sim_inst_paramsgroup:sim_inst_params:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sim_inst_paramsgroup:sim_inst_params:value\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_inst_paramsgroup:sim_inst_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INST_PARAMS_ROW \
"         \"sim_inst:simulation_id:0\",\"%s\",\"%s\",%22.16e"

#define PRINT_LIGOLW_XML_STOCHASTIC(fp) ( \
fputs( "   <Table Name=\"stochasticgroup:stochastic:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:duration\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:duration_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:cc_stat\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochasticgroup:stochastic:cc_sigma\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"stochasticgroup:stochastic:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCHASTIC_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_STOCH_SUMM(fp) ( \
fputs( "   <Table Name=\"stochsummgroup:stochsumm:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:y_opt\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"stochsummgroup:stochsumm:error\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"stochsummgroup:stochsumm:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCH_SUMM_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_EXT_TRIGGERS(fp) ( \
fputs( " <Table Name=\"external_trigger:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_alts\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_band\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_fluence\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_fluence_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_name\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_peak\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_peak_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:det_snr\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:email_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_dec_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_epoch\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_err_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_ra_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:start_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:start_time_ns\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_z\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_z_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_comments\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_id\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_sequence\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:notice_url\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_dec_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_fov_ra_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_loc_ele\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_loc_lat\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:obs_loc_long\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:ligo_fave_lho\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:ligo_fave_llo\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:ligo_delay\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_number_gcn\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_number_grb\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"external_trigger:event_status\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"external_trigger:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )


#define EXT_TRIGGERS_ROW \
"         \"process:process_id:0\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\",\"%s\",%d, %e, %e, \"%s\",\"%s\", %e , %e, %d, %d, \"%s\",%e, %e, \"%s\",\"%s\",\"%s\",%d, \"%s\",\"%s\",%e, %e,%e,%e, %e, %e, %e, %e, %e, %e, %d, \"%s\" , %d"

#define PRINT_LIGOLW_XML_FILTER(fp) ( \
fputs( "   <Table Name=\"filtergroup:filter:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"filtergroup:filter:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"filtergroup:filter:program\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"filtergroup:filter:start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"filtergroup:filter:filter_name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"filtergroup:filter:comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"filtergroup:filter:filter_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"filtergroup:filter:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define FILTER_ROW \
"         \"process:process_id:0\",\"%s\",%d,\"%s\",\"%s\",\"filter:filter_id:%" LAL_UINT8_FORMAT "\""

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOLWXMLHEADERS_H */
