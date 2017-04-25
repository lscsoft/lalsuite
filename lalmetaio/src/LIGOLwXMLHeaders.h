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
 * \author Brown, D. A.
 * \file
 * \ingroup lalmetaio_general
 *
 * \brief This header provides provides <tt>\#define</tt>s for the common elements of LIGO light weight XML files.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LIGOLwXMLHeaders.h>
 * \endcode
 *
 * It provides the XML header and footer, as well as table definitions for the various metadata tables.
 * It will need to be kept up to date with changes in the LIGO database table definitions.  The quantities which are
 * defined in this file are
 *
 * <ul>
 * <li> LIGOLW_XML_HEADER</li>
 * <li> LIGOLW_XML_FOOTER</li>
 * <li> LIGOLW_XML_TABLE_FOOTER</li>
 * <li> LIGOLW_XML_PROCESS</li>
 * <li> PROCESS_ROW</li>
 * <li> LIGOLW_XML_PROCESS_PARAMS</li>
 * <li> PROCESS_PARAMS_ROW</li>
 * <li> LIGOLW_XML_SEARCH_SUMMARY</li>
 * <li> SEARCH_SUMMARY_ROW</li>
 * <li> LIGOLW_XML_SEARCH_SUMMVARS</li>
 * <li> SEARCH_SUMMVARS_ROW</li>
 * <li> LIGOLW_XML_SIM_RINGDOWN</li>
 * <li> SIM_RINGDOWN_ROW</li>
 * <li> LIGOLW_XML_SUMM_VALUE</li>
 * <li> SUMM_VALUE_ROW</li>
 * <li> LIGOLW_XML_SIM_INST_PARAMS</li>
 * <li> SIM_INST_PARAMS_ROW</li>
 * <li> LIGOLW_XML_STOCHASTIC</li>
 * <li> STOCHASTIC_ROW</li>
 * <li> LIGOLW_XML_STOCH_SUMM</li>
 * <li> STOCH_SUMM_ROW</li>
 * <li> LIGOLW_XML_EXT_TRIGGERS</li>
 * <li> EXT_TRIGGERS_ROW</li>
 * </ul>
 *
 */

#ifndef _LIGOLWXMLHEADERS_H
#define _LIGOLWXMLHEADERS_H

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#define PRINT_LIGOLW_XML_HEADER(fp) (( \
fputs( "<?xml version='1.0' encoding='utf-8' ?>\n", fp) == EOF || \
fputs( "<!DOCTYPE LIGO_LW SYSTEM \"http://ldas-sw.ligo.caltech.edu/doc/ligolwAPI/html/ligolw_dtd.txt\">", fp ) == EOF || \
fputs( "<LIGO_LW>\n", fp ) == EOF ) ? EOF : 0)


#define PRINT_LIGOLW_XML_FOOTER(fp) \
fputs( "</LIGO_LW>", fp )

#define PRINT_LIGOLW_XML_TABLE_FOOTER(fp) ( \
fputs( "\n", fp ) == EOF || \
fputs( "      </Stream>\n", fp ) == EOF || \
fputs( "   </Table>\n", fp ) )

#define PRINT_LIGOLW_XML_PROCESS(fp) ( \
fputs( "   <Table Name=\"process:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"program\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"version\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cvs_repository\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"cvs_entry_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"is_online\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"node\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"username\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"unix_procid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"jobid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"domain\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ifos\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"process:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_ROW \
"         \"%s\",\"%s\",\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",\"process:process_id:0\""

#define PRINT_LIGOLW_XML_PROCESS_PARAMS(fp) ( \
fputs( "   <Table Name=\"process_params:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"program\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"param\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"type\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"value\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"process_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_PARAMS_ROW \
"         \"%s\",\"process:process_id:0\",\"%s\",\"%s\",\"%s\""

#define PRINT_LIGOLW_XML_SEARCH_SUMMARY(fp) ( \
fputs( "   <Table Name=\"search_summary:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"shared_object\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"lalwrapper_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"lal_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ifos\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"in_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"in_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"in_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"in_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"out_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"out_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"out_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"out_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"nevents\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"nnodes\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"search_summary:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMARY_ROW \
"         \"process:process_id:0\",\"standalone\",\"\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d"

#define PRINT_LIGOLW_XML_SEARCH_SUMMVARS(fp) (  \
fputs( "   <Table Name=\"search_summvars:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"name\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"string\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"value\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"search_summvar_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"search_summvars:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMVARS_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%22.16e,\"search_summvars:search_summvar_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SNGL_RINGDOWN(fp) ( \
fputs( "   <Table Name=\"sngl_ringdown:table\">\n", fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"start_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"frequency\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"quality\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mass\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"epsilon\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"num_clust_trigs\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ds2_H1H2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ds2_H1L1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ds2_H1V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ds2_H2L1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ds2_H2V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ds2_L1V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"sigma_sq\" Type=\"real_8\"/>\n", fp ) == EOF || \
fputs( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
fputs( "      <Stream Name=\"sngl_ringdown:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_RINGDOWN_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%22.16e,\"sngl_ringdown:event_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_RINGDOWN(fp) ( \
fputs( "   <Table Name=\"sim_ringdown:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF ||  \
fputs( "      <Column Name=\"waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"coordinates\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"geocent_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"geocent_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"h_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"h_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"l_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"l_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"v_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"v_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"inclination\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"frequency\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"quality\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"mass\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"spin\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"epsilon\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"hrss\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"hrss_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"hrss_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"hrss_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_ringdown:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_RINGDOWN_ROW \
  "         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"sim_ringdown:simulation_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SUMM_VALUE(fp) ( \
fputs( "   <Table Name=\"summ_value:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"program\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"value\" Type=\"real_4\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"summ_value_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"summ_value:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SUMM_VALUE_ROW \
"         \"%s\",\"process:process_id:0\",%d,%d,%d,%d,\"%s\",\"%s\",%e,\"%s\",\"summ_value:summ_value_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_INST_PARAMS(fp) ( \
fputs( "   <Table Name=\"sim_inst_params:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"value\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"sim_inst_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INST_PARAMS_ROW \
"         \"sim_inst:simulation_id:0\",\"%s\",\"%s\",%22.16e"

#define PRINT_LIGOLW_XML_STOCHASTIC(fp) ( \
fputs( "   <Table Name=\"stochastic:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"duration\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"duration_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cc_stat\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"cc_sigma\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"stochastic:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCHASTIC_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_STOCH_SUMM(fp) ( \
fputs( "   <Table Name=\"stochsumm:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"y_opt\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"error\" Type=\"real_8\"/>\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"stochsumm:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCH_SUMM_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_EXT_TRIGGERS(fp) ( \
fputs( " <Table Name=\"external_trigger:table\">\n" , fp ) == EOF || \
fputs( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_alts\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_band\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_fluence\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_fluence_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_name\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_peak\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_peak_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"det_snr\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"email_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_dec_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_epoch\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_err_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_ra_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"start_time_ns\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_z\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_z_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"notice_comments\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"notice_id\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"notice_sequence\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"notice_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"notice_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"notice_url\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"obs_fov_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"obs_fov_dec_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"obs_fov_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"obs_fov_ra_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"obs_loc_ele\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"obs_loc_lat\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"obs_loc_long\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ligo_fave_lho\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ligo_fave_llo\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"ligo_delay\" Type=\"real_4\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_number_gcn\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_number_grb\" Type=\"lstring\" />\n" , fp ) == EOF || \
fputs( "      <Column Name=\"event_status\" Type=\"int_4s\" />\n" , fp ) == EOF || \
fputs( "      <Stream Name=\"external_trigger:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )


#define EXT_TRIGGERS_ROW \
"         \"process:process_id:0\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\",\"%s\",%d, %e, %e, \"%s\",\"%s\", %e , %e, %d, %d, \"%s\",%e, %e, \"%s\",\"%s\",\"%s\",%d, \"%s\",\"%s\",%e, %e,%e,%e, %e, %e, %e, %e, %e, %e, %d, \"%s\" , %d"

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOLWXMLHEADERS_H */
