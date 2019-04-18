/*
*  Copyright (C) 2007 Andres C. Rodriguez, Sukanta Bose, Alexander Dietz, Duncan Brown, Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Robert Adam Mercer, Saikat Ray-Majumder, Anand Sengupta, Stephen Fairhurst, Xavier Siemens, Sean Seader, Thomas Cokelaer
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
 * File Name: LIGOLwXML.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A.
 * \file
 * \ingroup lalmetaio_general
 *
 * \brief Routines to write LIGO metadata database structures to LIGO lightweight XML files.
 *
 * ### Description ###
 *
 * The routine \c LALOpenLIGOLwXMLFile
 *
 * The routine \c LALCloseLIGOLwXMLFile prints the standard LIGO lightweight
 * XML footer, \c LIGOLW_XML_FOOTER given in LIGOLwXMLHeaders.h, and closes
 * the file stream pointed to by <tt>xml->fp</tt>.
 *
 * The routine \c LALBeginLIGOLwXMLTable prints the table header.  The type of
 * table to begin is specified by the \c table argument.  The appropriate
 * headers are again contained in LIGOLwXMLHeaders.h and contain the table name as
 * well as the names and data types of each of the columns in the table.  In
 * addition, it sets <tt>xml->first</tt> to 1 and <tt>xml->table</tt> to the requested
 * table.
 *
 * The routine \c LALEndLIGOLwXMLTable prints the table footer.  This is the
 * same for all tables, and given by \c LIGOLW_XML_TABLE_FOOTER in
 * LIGOLwXMLHeaders.h.  Additionally, <tt>xml->table</tt> is set to \c no_table.
 *
 * The routine \c LALWriteLIGOLwXMLTable writes the content of the xml table.
 * The type of table to be written is specified by \c table.  The contents of
 * the table should be stored as a linked list in <tt>tablePtr->table</tt>.  The data
 * is written using the row format for the specified table given in
 * LIGOLwXMLHeaders.h.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * <tt>fopen()</tt>
 * <tt>XLALFilePrintf()</tt>
 * <tt>fclose()</tt>
 *
 * ### Notes ###
 *
 * In order to change a table definition in LAL, changes must be made in
 * several places.  It is necessary to update the structure which is used to store
 * the information in memory as well as the reading and writing codes.  Below is a
 * list of all the files which must be updated.
 * <ul>
 * <li>  Update the LAL table definition in \ref LIGOMetadataTables.h</li>
 *
 * <li>  Update the LIGOLwXML writing code:</li>
 * <ol>
 * <li>  Change the table header written at to the LIGOLwXML file.  This is
 * \#defined in \ref LIGOLwXMLHeaders.h.  For example, to change the
 * \c sngl_inspiral table, you must edit \c LIGOLW_XML_SNGL_INSPIRAL.</li>
 *
 * <li> Change the row format of the LIGOLwXML file.  This is \#defined in
 * \ref LIGOLwXMLHeaders.h.  For example, to change the <tt> sngl_inspiral</tt>
 * table, you must edit \c SNGL_INSPIRAL_ROW.</li>
 *
 * <li> Change the XLALFilePrintf command which writes the table rows.  This is contained
 * in \ref LIGOLwXML.c.
 * </ol>
 *
 * <li> Update the LIGOLwXML reading code:</li>
 * <ol>
 * <li> Add/remove columns from the table directory of the table in question.
 * This is contained in \ref LIGOLwXMLRead.c, either in
 * \c LALCreateMetaTableDir or in the specific reading function.</li>
 *
 * <li> Check that all columns read in from the XML table are stored in memory.
 * This requires editing the table specific reading codes in
 * \ref LIGOLwXMLRead.c.</li>
 * </ol>
 * </ul>
 *
 */


#include <LIGOLwXMLlegacy.h>
#include <LIGOLwXMLHeaders.h>
#include <lal/LALVCSInfo.h>
#include <lal/LIGOMetadataTables.h>


/**\name Error Codes */ /*@{*/
#define LIGOLWXMLH_ENULL 1
#define LIGOLWXMLH_ENNUL 2
#define LIGOLWXMLH_EALOC 3
#define LIGOLWXMLH_EUTAB 4
#define LIGOLWXMLH_EOPEN 5
#define LIGOLWXMLH_ECLOS 6
#define LIGOLWXMLH_EBGNT 7
#define LIGOLWXMLH_ENTAB 8
#define LIGOLWXMLH_EENDT 8
#define LIGOLWXMLH_ETMSM 9
#define LIGOLWXMLH_ETNOP 10
#define LIGOLWXMLH_MSGENULL "Null pointer"
#define LIGOLWXMLH_MSGENNUL "Non-null pointer"
#define LIGOLWXMLH_MSGEALOC "Memory allocation error"
#define LIGOLWXMLH_MSGEUTAB "Unknown metadata table type"
#define LIGOLWXMLH_MSGEOPEN "Error opening XML file"
#define LIGOLWXMLH_MSGECLOS "Closing an XML file with an open table"
#define LIGOLWXMLH_MSGEBGNT "Begining a table without ending previous table"
#define LIGOLWXMLH_MSGENTAB "No table type specified"
#define LIGOLWXMLH_MSGEENDT "Ending a table without an beginning a table"
#define LIGOLWXMLH_MSGETMSM "Table type mismatch"
#define LIGOLWXMLH_MSGETNOP "Table not begun for writing"
/*@}*/

#define PRINT_LIGOLW_XML_PROCESS(fp) ( \
XLALFilePuts( "   <Table Name=\"process:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"program\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"version\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cvs_repository\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cvs_entry_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"is_online\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"node\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"username\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"unix_procid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"jobid\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"domain\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifos\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"process:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_ROW \
"         \"%s\",\"%s\",\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",\"process:process_id:0\""

#define PRINT_LIGOLW_XML_PROCESS_PARAMS(fp) ( \
XLALFilePuts( "   <Table Name=\"process_params:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"program\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"param\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"type\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"value\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"process_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define PROCESS_PARAMS_ROW \
"         \"%s\",\"process:process_id:0\",\"%s\",\"%s\",\"%s\""

#define PRINT_LIGOLW_XML_SEARCH_SUMMARY(fp) ( \
XLALFilePuts( "   <Table Name=\"search_summary:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"shared_object\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"lalwrapper_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"lal_cvs_tag\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"comment\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifos\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"in_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"in_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"in_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"in_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"out_start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"out_start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"out_end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"out_end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"nevents\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"nnodes\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"search_summary:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMARY_ROW \
"         \"process:process_id:0\",\"standalone\",\"\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d"

#define PRINT_LIGOLW_XML_SEARCH_SUMMVARS(fp) (  \
XLALFilePuts( "   <Table Name=\"search_summvars:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"name\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"string\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"value\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search_summvar_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"search_summvars:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SEARCH_SUMMVARS_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%22.16e,\"search_summvars:search_summvar_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SNGL_RINGDOWN(fp) ( \
XLALFilePuts( "   <Table Name=\"sngl_ringdown:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"frequency\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"quality\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"epsilon\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"num_clust_trigs\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ds2_H1H2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ds2_H1L1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ds2_H1V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ds2_H2L1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ds2_H2V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ds2_L1V1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigma_sq\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sngl_ringdown:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_RINGDOWN_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%22.16e,\"sngl_ringdown:event_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_RINGDOWN(fp) ( \
XLALFilePuts( "   <Table Name=\"sim_ringdown:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF ||  \
XLALFilePuts( "      <Column Name=\"waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coordinates\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"geocent_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"geocent_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v_start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v_start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"inclination\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"frequency\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"quality\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"epsilon\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"hrss\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"hrss_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"hrss_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"hrss_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sim_ringdown:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_RINGDOWN_ROW \
  "         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"sim_ringdown:simulation_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SUMM_VALUE(fp) ( \
XLALFilePuts( "   <Table Name=\"summ_value:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"program\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"value\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"summ_value_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"summ_value:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SUMM_VALUE_ROW \
"         \"%s\",\"process:process_id:0\",%d,%d,%d,%d,\"%s\",\"%s\",%e,\"%s\",\"summ_value:summ_value_id:%" LAL_UINT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_INST_PARAMS(fp) ( \
XLALFilePuts( "   <Table Name=\"sim_inst_params:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"name\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"comment\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"value\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sim_inst_params:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INST_PARAMS_ROW \
"         \"sim_inst:simulation_id:0\",\"%s\",\"%s\",%22.16e"

#define PRINT_LIGOLW_XML_STOCHASTIC(fp) ( \
XLALFilePuts( "   <Table Name=\"stochastic:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"duration\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"duration_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cc_stat\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cc_sigma\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"stochastic:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCHASTIC_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_STOCH_SUMM(fp) ( \
XLALFilePuts( "   <Table Name=\"stochsumm:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel_one\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel_two\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_min\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_max\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"y_opt\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"error\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"stochsumm:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define STOCH_SUMM_ROW \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.2f,%.2f,%22.16e,%22.16e"

#define PRINT_LIGOLW_XML_EXT_TRIGGERS(fp) ( \
XLALFilePuts( " <Table Name=\"external_trigger:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_alts\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_band\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_fluence\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_fluence_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_name\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_peak\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_peak_int\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"det_snr\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"email_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_dec_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_epoch\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_err_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_ra_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"start_time_ns\" Type=\"int_4s\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_z\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_z_err\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"notice_comments\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"notice_id\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"notice_sequence\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"notice_time\" Type=\"int_4s\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"notice_type\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"notice_url\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"obs_fov_dec\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"obs_fov_dec_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"obs_fov_ra\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"obs_fov_ra_width\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"obs_loc_ele\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"obs_loc_lat\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"obs_loc_long\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ligo_fave_lho\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ligo_fave_llo\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ligo_delay\" Type=\"real_4\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_number_gcn\" Type=\"int_4s\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_number_grb\" Type=\"lstring\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_status\" Type=\"int_4s\" />\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"external_trigger:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define EXT_TRIGGERS_ROW \
"         \"process:process_id:0\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\", \"%s\",\"%s\",%d, %e, %e, \"%s\",\"%s\", %e , %e, %d, %d, \"%s\",%e, %e, \"%s\",\"%s\",\"%s\",%d, \"%s\",\"%s\",%e, %e,%e,%e, %e, %e, %e, %e, %e, %e, %d, \"%s\" , %d"

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL(fp) ( \
XLALFilePuts( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"kappa\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chi\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau4\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau5\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha4\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha5\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha6\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"beta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma4\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma5\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma6\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma7\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma8\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"Gamma9\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1x\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1y\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1z\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2x\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2y\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2z\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW \
"         \"process:process_id:%ld\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,\"sngl_inspiral:event_id:%ld\""

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BNS(fp) ( \
XLALFilePuts( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mtotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BNS \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_SNGL_INSPIRAL_BCV(fp) ( \
XLALFilePuts( "   <Table Name=\"sngl_inspiral:table\">\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifo\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"channel\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"template_duration\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_distance\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq\" Type=\"real_8\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SNGL_INSPIRAL_ROW_BCV \
"         \"process:process_id:0\",\"%s\",\"%s\",\"%s\",%d,%d,%22.16e,%22.16e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%22.16e,%e"

#define PRINT_LIGOLW_XML_MULTI_INSPIRAL(fp) ( \
XLALFilePuts( "   <Table Name=\"multi_inspiral:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ifos\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"search\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"impulse_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amplitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h1h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chi\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"kappa\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"tau5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ttotal\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"trace_snr\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"snr_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_7\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_8\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_9\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_term_10\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_h1\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_h2\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_l\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_g\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_t\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sigmasq_v\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sngl_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bank_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sngl_bank_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_h1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_h2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cont_chisq_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"sngl_cont_chisq_dof\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ra\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"dec\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ligo_angle\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ligo_angle_sig\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"inclination\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"polarization\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"null_statistic\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"null_stat_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"null_stat_degen\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"event_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h2quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h2quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v1quad_re\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v1quad_im\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coh_snr_h1h2\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"cohSnrSqLocal\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"autoCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"crossCorrCohSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"autoCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"crossCorrNullSq\"  Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ampMetricEigenVal1\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"ampMetricEigenVal2\"  Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"time_slide:time_slide_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"multi_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define MULTI_INSPIRAL_ROW \
  "         \"process:process_id:0\",\"%s\",\"%s\",%d,%d,%22.16e,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%d,%e,%d,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"multi_inspiral:event_id:%" LAL_INT8_FORMAT "\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,\"time_slide:time_slide_id:%" LAL_INT8_FORMAT "\""

#define PRINT_LIGOLW_XML_SIM_INSPIRAL(fp) ( \
XLALFilePuts( "   <Table Name=\"sim_inspiral:table\">\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", fp ) == EOF ||  \
XLALFilePuts( "      <Column Name=\"waveform\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"geocent_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"geocent_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"h_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"l_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"g_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"t_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v_end_time\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"v_end_time_ns\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"source\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mass2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"mchirp\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"distance\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"longitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"latitude\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"inclination\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"coa_phase\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"polarization\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"psi3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha1\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha2\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha3\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha4\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha5\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"alpha6\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"beta\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin1z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2x\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2y\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"spin2z\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"theta0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"phi0\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_lower\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"f_final\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_h\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_g\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_t\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"numrel_mode_min\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"numrel_mode_max\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"numrel_data\" Type=\"lstring\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"amp_order\" Type=\"int_4s\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"taper\" Type=\"lstring\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"bandpass\" Type=\"int_4s\"/>\n", fp ) == EOF || \
XLALFilePuts( "      <Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , fp ) == EOF || \
XLALFilePuts( "      <Stream Name=\"sim_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", fp ) == EOF )

#define SIM_INSPIRAL_ROW \
"         \"process:process_id:%ld\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%22.16e,\"%s\",%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%d,\"%s\",%d,\"%s\",%d,\"sim_inspiral:simulation_id:%ld\""



void
LALOpenLIGOLwXMLFile (
    LALStatus          *status,
    LIGOLwXMLStream    *xml,
    const CHAR         *path
    )

{
  LIGOLwXMLStream *new;
  XLAL_PRINT_DEPRECATION_WARNING("XLALOpenLIGOLwXMLFile");
  INITSTATUS(status);
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( ! xml->fp, status, LIGOLWXMLH_ENNUL, LIGOLWXMLH_MSGENNUL );

  new = XLALOpenLIGOLwXMLFile( path );
  if ( ! new )
  {
    ABORT( status, LIGOLWXMLH_EOPEN, LIGOLWXMLH_MSGEOPEN );
  }

  *xml = *new;
  XLALFree(new);

  RETURN( status );
}


void
LALCloseLIGOLwXMLFile (
    LALStatus          *status,
    LIGOLwXMLStream    *xml
    )

{
  LIGOLwXMLStream *copy;
  XLAL_PRINT_DEPRECATION_WARNING("XLALCloseLIGOLwXMLFile");
  /* print the xml footer and close the file handle */
  INITSTATUS(status);
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  /* make an XLALFree()'able copy */
  copy = XLALMalloc(sizeof(*copy));
  if ( ! copy )
  {
    ABORT( status, LIGOLWXMLH_ECLOS, LIGOLWXMLH_MSGECLOS );
  }
  *copy = *xml;
  if ( XLALCloseLIGOLwXMLFile( copy ) )
  {
    XLALFree( copy );
    ABORT( status, LIGOLWXMLH_ECLOS, LIGOLWXMLH_MSGECLOS );
  }
  xml->fp = NULL;
  RETURN( status );
}


void
LALBeginLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml,
    MetadataTableType    table
    )

{
  /* print the header for the xml table */
  INITSTATUS(status);
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  if ( xml->table != no_table )
  {
    ABORT( status, LIGOLWXMLH_EBGNT, LIGOLWXMLH_MSGEBGNT );
  }

  switch( table )
  {
    case no_table:
      ABORT( status, LIGOLWXMLH_ENTAB, LIGOLWXMLH_MSGENTAB );
      break;
    case process_table:
      (void)PRINT_LIGOLW_XML_PROCESS( xml->fp );
      break;
    case process_params_table:
      (void)PRINT_LIGOLW_XML_PROCESS_PARAMS( xml->fp );
      break;
    case search_summary_table:
      (void)PRINT_LIGOLW_XML_SEARCH_SUMMARY( xml->fp );
      break;
    case search_summvars_table:
      (void)PRINT_LIGOLW_XML_SEARCH_SUMMVARS( xml->fp );
      break;
    case sngl_inspiral_table:
      (void)PRINT_LIGOLW_XML_SNGL_INSPIRAL( xml->fp );
      break;
    case sngl_inspiral_table_bns:
      (void)PRINT_LIGOLW_XML_SNGL_INSPIRAL_BNS( xml->fp );
      break;
    case sngl_inspiral_table_bcv:
      (void)PRINT_LIGOLW_XML_SNGL_INSPIRAL_BCV( xml->fp );
      break;
    case sngl_ringdown_table:
      (void)PRINT_LIGOLW_XML_SNGL_RINGDOWN( xml->fp );
      break;
    case multi_inspiral_table:
      (void)PRINT_LIGOLW_XML_MULTI_INSPIRAL( xml->fp );
      break;
    case sim_inspiral_table:
      (void)PRINT_LIGOLW_XML_SIM_INSPIRAL( xml->fp );
      break;
    case sim_ringdown_table:
      (void)PRINT_LIGOLW_XML_SIM_RINGDOWN( xml->fp );
      break;
    case summ_value_table:
      (void)PRINT_LIGOLW_XML_SUMM_VALUE( xml->fp );
      break;
    case sim_inst_params_table:
      (void)PRINT_LIGOLW_XML_SIM_INST_PARAMS( xml->fp );
      break;
    case stochastic_table:
      (void)PRINT_LIGOLW_XML_STOCHASTIC( xml->fp );
      break;
    case ext_triggers_table:
      (void)PRINT_LIGOLW_XML_EXT_TRIGGERS( xml->fp );
      break;
    default:
      ABORT( status, LIGOLWXMLH_EUTAB, LIGOLWXMLH_MSGEUTAB );
  }
  xml->first = 1;
  xml->rowCount = 0;
  xml->table = table;
  RETURN( status );
}


void
LALEndLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml
    )

{
  /* print the header for the xml table */
  INITSTATUS(status);
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  if ( xml->table == no_table )
  {
    ABORT( status, LIGOLWXMLH_EENDT, LIGOLWXMLH_MSGEENDT );
  }
  (void)PRINT_LIGOLW_XML_TABLE_FOOTER( xml->fp );
  xml->table = no_table;
  RETURN( status );
}

/* macro to print a comma on subsequent table rows */
#define FIRST_TABLE_ROW \
  if ( xml->first ) \
{ \
  xml->first = 0; \
} else \
{ \
  XLALFilePrintf( xml->fp, ",\n" ); \
}


void
LALWriteLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml,
    MetadataTable        tablePtr,
    MetadataTableType    table
    )

{
  /* print contents of the database struct into the xml table */
  INITSTATUS(status);
  ASSERT( xml, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  ASSERT( xml->fp, status, LIGOLWXMLH_ENULL, LIGOLWXMLH_MSGENULL );
  if ( xml->table == no_table )
  {
    ABORT( status, LIGOLWXMLH_ETNOP, LIGOLWXMLH_MSGETNOP );
  }
  if ( xml->table != table )
  {
    ABORT( status, LIGOLWXMLH_ETMSM, LIGOLWXMLH_MSGETMSM );
  }
  switch( table )
  {
    case no_table:
      ABORT( status, LIGOLWXMLH_ENTAB, LIGOLWXMLH_MSGENTAB );
      break;
    case process_table:
      while( tablePtr.processTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, PROCESS_ROW,
              tablePtr.processTable->program,
              tablePtr.processTable->version,
              tablePtr.processTable->cvs_repository,
              tablePtr.processTable->cvs_entry_time.gpsSeconds,
              tablePtr.processTable->comment,
              tablePtr.processTable->is_online,
              tablePtr.processTable->node,
              tablePtr.processTable->username,
              tablePtr.processTable->unix_procid,
              tablePtr.processTable->start_time.gpsSeconds,
              tablePtr.processTable->end_time.gpsSeconds,
              tablePtr.processTable->jobid,
              tablePtr.processTable->domain,
              tablePtr.processTable->ifos
              );
        tablePtr.processTable = tablePtr.processTable->next;
        ++(xml->rowCount);
      }
      break;
    case process_params_table:
      while( tablePtr.processParamsTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, PROCESS_PARAMS_ROW,
              tablePtr.processParamsTable->program,
              tablePtr.processParamsTable->param,
              tablePtr.processParamsTable->type,
              tablePtr.processParamsTable->value
              );
        tablePtr.processParamsTable = tablePtr.processParamsTable->next;
        ++(xml->rowCount);
      }
      break;
    case search_summary_table:
      while( tablePtr.searchSummaryTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SEARCH_SUMMARY_ROW,
              lalVCSInfo.vcsTag,
              tablePtr.searchSummaryTable->comment,
              tablePtr.searchSummaryTable->ifos,
              tablePtr.searchSummaryTable->in_start_time.gpsSeconds,
              tablePtr.searchSummaryTable->in_start_time.gpsNanoSeconds,
              tablePtr.searchSummaryTable->in_end_time.gpsSeconds,
              tablePtr.searchSummaryTable->in_end_time.gpsNanoSeconds,
              tablePtr.searchSummaryTable->out_start_time.gpsSeconds,
              tablePtr.searchSummaryTable->out_start_time.gpsNanoSeconds,
              tablePtr.searchSummaryTable->out_end_time.gpsSeconds,
              tablePtr.searchSummaryTable->out_end_time.gpsNanoSeconds,
              tablePtr.searchSummaryTable->nevents,
              tablePtr.searchSummaryTable->nnodes
              );
        tablePtr.searchSummaryTable = tablePtr.searchSummaryTable->next;
        ++(xml->rowCount);
      }
      break;
    case search_summvars_table:
      while( tablePtr.searchSummvarsTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SEARCH_SUMMVARS_ROW,
              tablePtr.searchSummvarsTable->name,
              tablePtr.searchSummvarsTable->string,
              tablePtr.searchSummvarsTable->value,
              xml->rowCount
              );
        tablePtr.searchSummvarsTable = tablePtr.searchSummvarsTable->next;
        ++(xml->rowCount);
      }
      break;
    case sngl_inspiral_table:
      while( tablePtr.snglInspiralTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SNGL_INSPIRAL_ROW,
              tablePtr.snglInspiralTable->process_id,
              tablePtr.snglInspiralTable->ifo,
              tablePtr.snglInspiralTable->search,
              tablePtr.snglInspiralTable->channel,
              tablePtr.snglInspiralTable->end.gpsSeconds,
              tablePtr.snglInspiralTable->end.gpsNanoSeconds,
              tablePtr.snglInspiralTable->end_time_gmst,
              tablePtr.snglInspiralTable->impulse_time.gpsSeconds,
              tablePtr.snglInspiralTable->impulse_time.gpsNanoSeconds,
              tablePtr.snglInspiralTable->template_duration,
              tablePtr.snglInspiralTable->event_duration,
              tablePtr.snglInspiralTable->amplitude,
              tablePtr.snglInspiralTable->eff_distance,
              tablePtr.snglInspiralTable->coa_phase,
              tablePtr.snglInspiralTable->mass1,
              tablePtr.snglInspiralTable->mass2,
              tablePtr.snglInspiralTable->mchirp,
              tablePtr.snglInspiralTable->mtotal,
              tablePtr.snglInspiralTable->eta,
              tablePtr.snglInspiralTable->kappa,
              tablePtr.snglInspiralTable->chi,
              tablePtr.snglInspiralTable->tau0,
              tablePtr.snglInspiralTable->tau2,
              tablePtr.snglInspiralTable->tau3,
              tablePtr.snglInspiralTable->tau4,
              tablePtr.snglInspiralTable->tau5,
              tablePtr.snglInspiralTable->ttotal,
              tablePtr.snglInspiralTable->psi0,
              tablePtr.snglInspiralTable->psi3,
              tablePtr.snglInspiralTable->alpha,
              tablePtr.snglInspiralTable->alpha1,
              tablePtr.snglInspiralTable->alpha2,
              tablePtr.snglInspiralTable->alpha3,
              tablePtr.snglInspiralTable->alpha4,
              tablePtr.snglInspiralTable->alpha5,
              tablePtr.snglInspiralTable->alpha6,
              tablePtr.snglInspiralTable->beta,
              tablePtr.snglInspiralTable->f_final,
              tablePtr.snglInspiralTable->snr,
              tablePtr.snglInspiralTable->chisq,
              tablePtr.snglInspiralTable->chisq_dof,
              tablePtr.snglInspiralTable->bank_chisq,
              tablePtr.snglInspiralTable->bank_chisq_dof,
              tablePtr.snglInspiralTable->cont_chisq,
              tablePtr.snglInspiralTable->cont_chisq_dof,
              tablePtr.snglInspiralTable->sigmasq,
	      tablePtr.snglInspiralTable->rsqveto_duration,
	      tablePtr.snglInspiralTable->Gamma[0],
	      tablePtr.snglInspiralTable->Gamma[1],
	      tablePtr.snglInspiralTable->Gamma[2],
	      tablePtr.snglInspiralTable->Gamma[3],
	      tablePtr.snglInspiralTable->Gamma[4],
	      tablePtr.snglInspiralTable->Gamma[5],
	      tablePtr.snglInspiralTable->Gamma[6],
	      tablePtr.snglInspiralTable->Gamma[7],
	      tablePtr.snglInspiralTable->Gamma[8],
	      tablePtr.snglInspiralTable->Gamma[9],
              tablePtr.snglInspiralTable->spin1x,
              tablePtr.snglInspiralTable->spin1y,
              tablePtr.snglInspiralTable->spin1z,
              tablePtr.snglInspiralTable->spin2x,
              tablePtr.snglInspiralTable->spin2y,
              tablePtr.snglInspiralTable->spin2z,
              tablePtr.snglInspiralTable->event_id );
        tablePtr.snglInspiralTable = tablePtr.snglInspiralTable->next;
        ++(xml->rowCount);
      }
      break;
    case sngl_inspiral_table_bns:
      while( tablePtr.snglInspiralTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SNGL_INSPIRAL_ROW_BNS,
              tablePtr.snglInspiralTable->ifo,
              tablePtr.snglInspiralTable->search,
              tablePtr.snglInspiralTable->channel,
              tablePtr.snglInspiralTable->end.gpsSeconds,
              tablePtr.snglInspiralTable->end.gpsNanoSeconds,
              tablePtr.snglInspiralTable->end_time_gmst,
              tablePtr.snglInspiralTable->template_duration,
              tablePtr.snglInspiralTable->eff_distance,
              tablePtr.snglInspiralTable->coa_phase,
              tablePtr.snglInspiralTable->mass1,
              tablePtr.snglInspiralTable->mass2,
              tablePtr.snglInspiralTable->mchirp,
              tablePtr.snglInspiralTable->mtotal,
              tablePtr.snglInspiralTable->eta,
              tablePtr.snglInspiralTable->tau0,
              tablePtr.snglInspiralTable->tau3,
              tablePtr.snglInspiralTable->ttotal,
              tablePtr.snglInspiralTable->f_final,
              tablePtr.snglInspiralTable->snr,
              tablePtr.snglInspiralTable->chisq,
              tablePtr.snglInspiralTable->chisq_dof,
              tablePtr.snglInspiralTable->bank_chisq,
              tablePtr.snglInspiralTable->bank_chisq_dof,
              tablePtr.snglInspiralTable->cont_chisq,
              tablePtr.snglInspiralTable->cont_chisq_dof,
              tablePtr.snglInspiralTable->sigmasq,
	      tablePtr.snglInspiralTable->rsqveto_duration );
        tablePtr.snglInspiralTable = tablePtr.snglInspiralTable->next;
        ++(xml->rowCount);
      }
      break;
    case sngl_inspiral_table_bcv:
      while( tablePtr.snglInspiralTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SNGL_INSPIRAL_ROW_BCV,
              tablePtr.snglInspiralTable->ifo,
              tablePtr.snglInspiralTable->search,
              tablePtr.snglInspiralTable->channel,
              tablePtr.snglInspiralTable->end.gpsSeconds,
              tablePtr.snglInspiralTable->end.gpsNanoSeconds,
              tablePtr.snglInspiralTable->end_time_gmst,
              tablePtr.snglInspiralTable->template_duration,
              tablePtr.snglInspiralTable->eff_distance,
              tablePtr.snglInspiralTable->coa_phase,
              tablePtr.snglInspiralTable->mchirp,
              tablePtr.snglInspiralTable->eta,
              tablePtr.snglInspiralTable->psi0,
              tablePtr.snglInspiralTable->psi3,
              tablePtr.snglInspiralTable->alpha,
              tablePtr.snglInspiralTable->f_final,
              tablePtr.snglInspiralTable->snr,
              tablePtr.snglInspiralTable->chisq,
              tablePtr.snglInspiralTable->chisq_dof,
              tablePtr.snglInspiralTable->bank_chisq,
              tablePtr.snglInspiralTable->bank_chisq_dof,
              tablePtr.snglInspiralTable->cont_chisq,
              tablePtr.snglInspiralTable->cont_chisq_dof,
              tablePtr.snglInspiralTable->sigmasq,
	      tablePtr.snglInspiralTable->rsqveto_duration );
        tablePtr.snglInspiralTable = tablePtr.snglInspiralTable->next;
        ++(xml->rowCount);
      }
      break;
    case sngl_ringdown_table:
      while( tablePtr.snglRingdownTable )
      {
        UINT8 id = xml->rowCount;
        if ( tablePtr.snglRingdownTable->event_id )
        {
          id = tablePtr.snglRingdownTable->event_id->id;
        }
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SNGL_RINGDOWN_ROW,
              tablePtr.snglRingdownTable->ifo,
              tablePtr.snglRingdownTable->channel,
              tablePtr.snglRingdownTable->start_time.gpsSeconds,
              tablePtr.snglRingdownTable->start_time.gpsNanoSeconds,
              tablePtr.snglRingdownTable->start_time_gmst,
              tablePtr.snglRingdownTable->frequency,
              tablePtr.snglRingdownTable->quality,
              tablePtr.snglRingdownTable->phase,
              tablePtr.snglRingdownTable->mass,
              tablePtr.snglRingdownTable->spin,
              tablePtr.snglRingdownTable->epsilon,
              tablePtr.snglRingdownTable->num_clust_trigs,
              tablePtr.snglRingdownTable->ds2_H1H2,
              tablePtr.snglRingdownTable->ds2_H1L1,
              tablePtr.snglRingdownTable->ds2_H1V1,
              tablePtr.snglRingdownTable->ds2_H2L1,
              tablePtr.snglRingdownTable->ds2_H2V1,
              tablePtr.snglRingdownTable->ds2_L1V1,
              tablePtr.snglRingdownTable->amplitude,
              tablePtr.snglRingdownTable->snr,
              tablePtr.snglRingdownTable->eff_dist,
              tablePtr.snglRingdownTable->sigma_sq,
              id
              );
        tablePtr.snglRingdownTable = tablePtr.snglRingdownTable->next;
        ++(xml->rowCount);
      }
      break;
    case multi_inspiral_table:
      while( tablePtr.multiInspiralTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, MULTI_INSPIRAL_ROW,
              tablePtr.multiInspiralTable->ifos,
              tablePtr.multiInspiralTable->search,
              tablePtr.multiInspiralTable->end_time.gpsSeconds,
              tablePtr.multiInspiralTable->end_time.gpsNanoSeconds,
              tablePtr.multiInspiralTable->end_time_gmst,
              tablePtr.multiInspiralTable->impulse_time.gpsSeconds,
              tablePtr.multiInspiralTable->impulse_time.gpsNanoSeconds,
              tablePtr.multiInspiralTable->amplitude,
              tablePtr.multiInspiralTable->distance,
              tablePtr.multiInspiralTable->eff_dist_h1,
              tablePtr.multiInspiralTable->eff_dist_h2,
              tablePtr.multiInspiralTable->eff_dist_l,
              tablePtr.multiInspiralTable->eff_dist_g,
              tablePtr.multiInspiralTable->eff_dist_t,
              tablePtr.multiInspiralTable->eff_dist_v,
              tablePtr.multiInspiralTable->eff_dist_h1h2,
              tablePtr.multiInspiralTable->coa_phase,
              tablePtr.multiInspiralTable->mass1,
              tablePtr.multiInspiralTable->mass2,
              tablePtr.multiInspiralTable->mchirp,
              tablePtr.multiInspiralTable->eta,
              tablePtr.multiInspiralTable->chi,
              tablePtr.multiInspiralTable->kappa,
              tablePtr.multiInspiralTable->tau0,
              tablePtr.multiInspiralTable->tau2,
              tablePtr.multiInspiralTable->tau3,
              tablePtr.multiInspiralTable->tau4,
              tablePtr.multiInspiralTable->tau5,
              tablePtr.multiInspiralTable->ttotal,
              tablePtr.multiInspiralTable->snr,
              tablePtr.multiInspiralTable->snr_dof,
              tablePtr.multiInspiralTable->chisq,
              tablePtr.multiInspiralTable->chisq_dof,
              tablePtr.multiInspiralTable->bank_chisq,
              tablePtr.multiInspiralTable->bank_chisq_dof,
              tablePtr.multiInspiralTable->cont_chisq,
              tablePtr.multiInspiralTable->cont_chisq_dof,
              tablePtr.multiInspiralTable->trace_snr,
              tablePtr.multiInspiralTable->snr_h1,
              tablePtr.multiInspiralTable->snr_h2,
              tablePtr.multiInspiralTable->snr_l,
              tablePtr.multiInspiralTable->snr_g,
              tablePtr.multiInspiralTable->snr_t,
              tablePtr.multiInspiralTable->snr_v,
              tablePtr.multiInspiralTable->amp_term_1,
              tablePtr.multiInspiralTable->amp_term_2,
              tablePtr.multiInspiralTable->amp_term_3,
              tablePtr.multiInspiralTable->amp_term_4,
              tablePtr.multiInspiralTable->amp_term_5,
              tablePtr.multiInspiralTable->amp_term_6,
              tablePtr.multiInspiralTable->amp_term_7,
              tablePtr.multiInspiralTable->amp_term_8,
              tablePtr.multiInspiralTable->amp_term_9,
              tablePtr.multiInspiralTable->amp_term_10,
              tablePtr.multiInspiralTable->sigmasq_h1,
              tablePtr.multiInspiralTable->sigmasq_h2,
              tablePtr.multiInspiralTable->sigmasq_l,
              tablePtr.multiInspiralTable->sigmasq_g,
              tablePtr.multiInspiralTable->sigmasq_t,
              tablePtr.multiInspiralTable->sigmasq_v,
              tablePtr.multiInspiralTable->chisq_h1,
              tablePtr.multiInspiralTable->chisq_h2,
              tablePtr.multiInspiralTable->chisq_l,
              tablePtr.multiInspiralTable->chisq_g,
              tablePtr.multiInspiralTable->chisq_t,
              tablePtr.multiInspiralTable->chisq_v,
              tablePtr.multiInspiralTable->sngl_chisq_dof,
              tablePtr.multiInspiralTable->bank_chisq_h1,
              tablePtr.multiInspiralTable->bank_chisq_h2,
              tablePtr.multiInspiralTable->bank_chisq_l,
              tablePtr.multiInspiralTable->bank_chisq_g,
              tablePtr.multiInspiralTable->bank_chisq_t,
              tablePtr.multiInspiralTable->bank_chisq_v,
              tablePtr.multiInspiralTable->sngl_bank_chisq_dof,
              tablePtr.multiInspiralTable->cont_chisq_h1,
              tablePtr.multiInspiralTable->cont_chisq_h2,
              tablePtr.multiInspiralTable->cont_chisq_l,
              tablePtr.multiInspiralTable->cont_chisq_g,
              tablePtr.multiInspiralTable->cont_chisq_t,
              tablePtr.multiInspiralTable->cont_chisq_v,
              tablePtr.multiInspiralTable->sngl_cont_chisq_dof,
              tablePtr.multiInspiralTable->ra,
              tablePtr.multiInspiralTable->dec,
              tablePtr.multiInspiralTable->ligo_angle,
	      tablePtr.multiInspiralTable->ligo_angle_sig,
	      tablePtr.multiInspiralTable->inclination,
              tablePtr.multiInspiralTable->polarization,
	      tablePtr.multiInspiralTable->null_statistic,
              tablePtr.multiInspiralTable->null_stat_h1h2,
              tablePtr.multiInspiralTable->null_stat_degen,
              tablePtr.multiInspiralTable->event_id->id,
              crealf(tablePtr.multiInspiralTable->h1quad),
              cimagf(tablePtr.multiInspiralTable->h1quad),
              crealf(tablePtr.multiInspiralTable->h2quad),
              cimagf(tablePtr.multiInspiralTable->h2quad),
              crealf(tablePtr.multiInspiralTable->l1quad),
              cimagf(tablePtr.multiInspiralTable->l1quad),
              crealf(tablePtr.multiInspiralTable->g1quad),
              cimagf(tablePtr.multiInspiralTable->g1quad),
              crealf(tablePtr.multiInspiralTable->t1quad),
              cimagf(tablePtr.multiInspiralTable->t1quad),
              crealf(tablePtr.multiInspiralTable->v1quad),
              cimagf(tablePtr.multiInspiralTable->v1quad),
	      tablePtr.multiInspiralTable->coh_snr_h1h2,
	      tablePtr.multiInspiralTable->cohSnrSqLocal,
	      tablePtr.multiInspiralTable->autoCorrCohSq,
	      tablePtr.multiInspiralTable->crossCorrCohSq,
	      tablePtr.multiInspiralTable->autoCorrNullSq,
	      tablePtr.multiInspiralTable->crossCorrNullSq,
	      tablePtr.multiInspiralTable->ampMetricEigenVal1,
	      tablePtr.multiInspiralTable->ampMetricEigenVal2,
              tablePtr.multiInspiralTable->time_slide_id->id
              );
        tablePtr.multiInspiralTable = tablePtr.multiInspiralTable->next;
        ++(xml->rowCount);
      }
      break;
    case sim_inspiral_table:
      {
      while( tablePtr.simInspiralTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SIM_INSPIRAL_ROW,
              tablePtr.simInspiralTable->process_id,
              tablePtr.simInspiralTable->waveform,
              tablePtr.simInspiralTable->geocent_end_time.gpsSeconds,
              tablePtr.simInspiralTable->geocent_end_time.gpsNanoSeconds,
              tablePtr.simInspiralTable->h_end_time.gpsSeconds,
              tablePtr.simInspiralTable->h_end_time.gpsNanoSeconds,
              tablePtr.simInspiralTable->l_end_time.gpsSeconds,
              tablePtr.simInspiralTable->l_end_time.gpsNanoSeconds,
              tablePtr.simInspiralTable->g_end_time.gpsSeconds,
              tablePtr.simInspiralTable->g_end_time.gpsNanoSeconds,
              tablePtr.simInspiralTable->t_end_time.gpsSeconds,
              tablePtr.simInspiralTable->t_end_time.gpsNanoSeconds,
              tablePtr.simInspiralTable->v_end_time.gpsSeconds,
              tablePtr.simInspiralTable->v_end_time.gpsNanoSeconds,
              tablePtr.simInspiralTable->end_time_gmst,
              tablePtr.simInspiralTable->source,
              tablePtr.simInspiralTable->mass1,
              tablePtr.simInspiralTable->mass2,
              tablePtr.simInspiralTable->mchirp,
              tablePtr.simInspiralTable->eta,
              tablePtr.simInspiralTable->distance,
              tablePtr.simInspiralTable->longitude,
              tablePtr.simInspiralTable->latitude,
              tablePtr.simInspiralTable->inclination,
              tablePtr.simInspiralTable->coa_phase,
              tablePtr.simInspiralTable->polarization,
              tablePtr.simInspiralTable->psi0,
              tablePtr.simInspiralTable->psi3,
              tablePtr.simInspiralTable->alpha,
              tablePtr.simInspiralTable->alpha1,
              tablePtr.simInspiralTable->alpha2,
              tablePtr.simInspiralTable->alpha3,
              tablePtr.simInspiralTable->alpha4,
              tablePtr.simInspiralTable->alpha5,
              tablePtr.simInspiralTable->alpha6,
              tablePtr.simInspiralTable->beta,
              tablePtr.simInspiralTable->spin1x,
              tablePtr.simInspiralTable->spin1y,
              tablePtr.simInspiralTable->spin1z,
              tablePtr.simInspiralTable->spin2x,
              tablePtr.simInspiralTable->spin2y,
              tablePtr.simInspiralTable->spin2z,
              tablePtr.simInspiralTable->theta0,
              tablePtr.simInspiralTable->phi0,
              tablePtr.simInspiralTable->f_lower,
              tablePtr.simInspiralTable->f_final,
              tablePtr.simInspiralTable->eff_dist_h,
              tablePtr.simInspiralTable->eff_dist_l,
              tablePtr.simInspiralTable->eff_dist_g,
              tablePtr.simInspiralTable->eff_dist_t,
              tablePtr.simInspiralTable->eff_dist_v,
	      tablePtr.simInspiralTable->numrel_mode_min,
	      tablePtr.simInspiralTable->numrel_mode_max,
	      tablePtr.simInspiralTable->numrel_data,
	      tablePtr.simInspiralTable->amp_order,
	      tablePtr.simInspiralTable->taper,
	      tablePtr.simInspiralTable->bandpass,
	      tablePtr.simInspiralTable->simulation_id
              );
        tablePtr.simInspiralTable = tablePtr.simInspiralTable->next;
        ++(xml->rowCount);
        }
      }
      break;
    case sim_ringdown_table:
      {
        while( tablePtr.simRingdownTable )
        {
          FIRST_TABLE_ROW
            XLALFilePrintf( xml->fp, SIM_RINGDOWN_ROW,
                tablePtr.simRingdownTable->waveform,
                tablePtr.simRingdownTable->coordinates,
                tablePtr.simRingdownTable->geocent_start_time.gpsSeconds,
                tablePtr.simRingdownTable->geocent_start_time.gpsNanoSeconds,
                tablePtr.simRingdownTable->h_start_time.gpsSeconds,
                tablePtr.simRingdownTable->h_start_time.gpsNanoSeconds,
                tablePtr.simRingdownTable->l_start_time.gpsSeconds,
                tablePtr.simRingdownTable->l_start_time.gpsNanoSeconds,
                tablePtr.simRingdownTable->v_start_time.gpsSeconds,
                tablePtr.simRingdownTable->v_start_time.gpsNanoSeconds,
                tablePtr.simRingdownTable->start_time_gmst,
                tablePtr.simRingdownTable->longitude,
                tablePtr.simRingdownTable->latitude,
                tablePtr.simRingdownTable->distance,
                tablePtr.simRingdownTable->inclination,
                tablePtr.simRingdownTable->polarization,
                tablePtr.simRingdownTable->frequency,
                tablePtr.simRingdownTable->quality,
                tablePtr.simRingdownTable->phase,
                tablePtr.simRingdownTable->mass,
                tablePtr.simRingdownTable->spin,
                tablePtr.simRingdownTable->epsilon,
                tablePtr.simRingdownTable->amplitude,
                tablePtr.simRingdownTable->eff_dist_h,
                tablePtr.simRingdownTable->eff_dist_l,
                tablePtr.simRingdownTable->eff_dist_v,
                tablePtr.simRingdownTable->hrss,
                tablePtr.simRingdownTable->hrss_h,
                tablePtr.simRingdownTable->hrss_l,
                tablePtr.simRingdownTable->hrss_v,
                xml->rowCount
                  );
          tablePtr.simRingdownTable = tablePtr.simRingdownTable->next;
          ++(xml->rowCount);
        }
      }
      break;
    case summ_value_table:
      while( tablePtr.summValueTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SUMM_VALUE_ROW,
              tablePtr.summValueTable->program,
              tablePtr.summValueTable->start_time.gpsSeconds,
              tablePtr.summValueTable->start_time.gpsNanoSeconds,
              tablePtr.summValueTable->end_time.gpsSeconds,
              tablePtr.summValueTable->end_time.gpsNanoSeconds,
              tablePtr.summValueTable->ifo,
              tablePtr.summValueTable->name,
              tablePtr.summValueTable->value,
              tablePtr.summValueTable->comment,
              xml->rowCount
              );
        tablePtr.snglInspiralTable = tablePtr.snglInspiralTable->next;
        ++(xml->rowCount);
      }
      break;
    case sim_inst_params_table:
      while( tablePtr.simInstParamsTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, SIM_INST_PARAMS_ROW,
              tablePtr.simInstParamsTable->name,
              tablePtr.simInstParamsTable->comment,
              tablePtr.simInstParamsTable->value
              );
        tablePtr.simInstParamsTable = tablePtr.simInstParamsTable->next;
        ++(xml->rowCount);
      }
      break;
    case stochastic_table:
      while( tablePtr.stochasticTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, STOCHASTIC_ROW,
              tablePtr.stochasticTable->ifo_one,
              tablePtr.stochasticTable->ifo_two,
              tablePtr.stochasticTable->channel_one,
              tablePtr.stochasticTable->channel_two,
              tablePtr.stochasticTable->start_time.gpsSeconds,
              tablePtr.stochasticTable->start_time.gpsNanoSeconds,
              tablePtr.stochasticTable->duration.gpsSeconds,
              tablePtr.stochasticTable->duration.gpsNanoSeconds,
              tablePtr.stochasticTable->f_min,
              tablePtr.stochasticTable->f_max,
              tablePtr.stochasticTable->cc_stat,
              tablePtr.stochasticTable->cc_sigma
              );
        tablePtr.stochasticTable = tablePtr.stochasticTable->next;
        ++(xml->rowCount);
      }
      break;
    case stoch_summ_table:
      while( tablePtr.stochSummTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, STOCH_SUMM_ROW,
              tablePtr.stochSummTable->ifo_one,
              tablePtr.stochSummTable->ifo_two,
              tablePtr.stochSummTable->channel_one,
              tablePtr.stochSummTable->channel_two,
              tablePtr.stochSummTable->start_time.gpsSeconds,
              tablePtr.stochSummTable->start_time.gpsNanoSeconds,
              tablePtr.stochSummTable->end_time.gpsSeconds,
              tablePtr.stochSummTable->end_time.gpsNanoSeconds,
              tablePtr.stochSummTable->f_min,
              tablePtr.stochSummTable->f_max,
              tablePtr.stochSummTable->y_opt,
              tablePtr.stochSummTable->error
              );
        tablePtr.stochSummTable = tablePtr.stochSummTable->next;
        ++(xml->rowCount);
      }
      break;
    case ext_triggers_table:
      while( tablePtr.extTriggerTable )
      {
        FIRST_TABLE_ROW
          XLALFilePrintf( xml->fp, EXT_TRIGGERS_ROW,
              tablePtr.extTriggerTable->det_alts,
              tablePtr.extTriggerTable->det_band,
              tablePtr.extTriggerTable->det_fluence,
              tablePtr.extTriggerTable->det_fluence_int,
              tablePtr.extTriggerTable->det_name,
              tablePtr.extTriggerTable->det_peak,
              tablePtr.extTriggerTable->det_peak_int,
              tablePtr.extTriggerTable->det_snr,
              tablePtr.extTriggerTable->email_time,
              tablePtr.extTriggerTable->event_dec,
              tablePtr.extTriggerTable->event_dec_err,
              tablePtr.extTriggerTable->event_epoch,
              tablePtr.extTriggerTable->event_err_type,
              tablePtr.extTriggerTable->event_ra,
              tablePtr.extTriggerTable->event_ra_err,
              tablePtr.extTriggerTable->start_time,
              tablePtr.extTriggerTable->start_time_ns,
              tablePtr.extTriggerTable->event_type,
              tablePtr.extTriggerTable->event_z,
              tablePtr.extTriggerTable->event_z_err,
              tablePtr.extTriggerTable->notice_comments,
              tablePtr.extTriggerTable->notice_id,
              tablePtr.extTriggerTable->notice_sequence,
              tablePtr.extTriggerTable->notice_time,
              tablePtr.extTriggerTable->notice_type,
              tablePtr.extTriggerTable->notice_url,
              tablePtr.extTriggerTable->obs_fov_dec,
              tablePtr.extTriggerTable->obs_fov_dec_width,
              tablePtr.extTriggerTable->obs_fov_ra,
              tablePtr.extTriggerTable->obs_fov_ra_width,
	      tablePtr.extTriggerTable->obs_loc_ele,
	      tablePtr.extTriggerTable->obs_loc_lat,
	      tablePtr.extTriggerTable->obs_loc_long,
	      tablePtr.extTriggerTable->ligo_fave_lho,
	      tablePtr.extTriggerTable->ligo_fave_llo,
	      tablePtr.extTriggerTable->ligo_delay,
	      tablePtr.extTriggerTable->event_number_gcn,
	      tablePtr.extTriggerTable->event_number_grb,
	      tablePtr.extTriggerTable->event_status
	    );
        tablePtr.extTriggerTable = tablePtr.extTriggerTable->next;
        ++(xml->rowCount);
      }
      break;
    default:
      ABORT( status, LIGOLWXMLH_EUTAB, LIGOLWXMLH_MSGEUTAB );
  }
  RETURN( status );
}
