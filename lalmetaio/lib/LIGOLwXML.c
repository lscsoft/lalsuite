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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
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
 * The routine \c XLALCreateLIGOLwXMLFileName creates a name for a  LIGO lightweight XML file that is in accordance with the specifications of document T050017.
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
 */


#include <lal/FileIO.h>
#include <lal/LALMalloc.h>
#include <lal/LALVCSInfo.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/XLALError.h>
#include <LIGOLwXMLHeaders.h>


/**
 * Open an XML file for writing.  The return value is a pointer to a new
 * LIGOLwXMLStream file handle or NULL on failure.
 */
LIGOLwXMLStream *
XLALOpenLIGOLwXMLFile (
    const char *path
)
{
  LIGOLwXMLStream *new;

  /* malloc a new XML file handle */

  new = XLALMalloc( sizeof( *new ) );
  if ( ! new )
    XLAL_ERROR_NULL( XLAL_EFUNC );

  /* fopen() the underlying C file */

  new->fp = XLALFileOpen( path, "w" );
  if ( ! new->fp )
  {
    XLALFree(new);
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }

  /* initialize the table flag */

  new->table = no_table;

  /* write the XML header */

  if ( PRINT_LIGOLW_XML_HEADER( new->fp ) < 0 )
  {
    XLALFileClose( new->fp );
    XLALFree( new );
    XLAL_ERROR_NULL( XLAL_EIO );
  }

  /* done */

  return new;
}


/**
 * Close an XML stream.  On failure the stream is left in an undefined
 * state, and has not been free()'ed.  Sorry.
 */

int
XLALCloseLIGOLwXMLFile (
  LIGOLwXMLStream *xml
)
{
  if ( xml )
  {
    if ( xml->table != no_table)
      /* trying to close the file in the middle of a table */
      XLAL_ERROR(XLAL_EFAILED);
    if ( PRINT_LIGOLW_XML_FOOTER( xml->fp ) < 0 )
      /* can't write XML footer */
      XLAL_ERROR( XLAL_EIO );
    if ( XLALFileClose( xml->fp ) )
      /* fclose() on the underlying C file failed */
      XLAL_ERROR( XLAL_EFUNC );
  }

  XLALFree( xml );

  return 0;
}


/**
 * Write a process table to an XML file.
 */


int XLALWriteLIGOLwXMLProcessTable(
	LIGOLwXMLStream *xml,
	const ProcessTable *process
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"process:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:program\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:version\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:cvs_repository\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:cvs_entry_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:comment\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:is_online\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:node\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:username\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:unix_procid\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:jobid\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:domain\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:ifos\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"process:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; process; process = process->next) {
		if(XLALFilePrintf(xml->fp, "%s\"%s\",\"%s\",\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",%ld",
			row_head,
			process->program,
			process->version,
			process->cvs_repository,
			process->cvs_entry_time.gpsSeconds,
			process->comment,
			process->is_online,
			process->node,
			process->username,
			process->unix_procid,
			process->start_time.gpsSeconds,
			process->end_time.gpsSeconds,
			process->jobid,
			process->domain,
			process->ifos,
			process->process_id
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}


/**
 * Write a process_params table to an XML file.
 */


int XLALWriteLIGOLwXMLProcessParamsTable(
	LIGOLwXMLStream *xml,
	const ProcessParamsTable *process_params
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"process_params:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process_params:program\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process_params:param\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process_params:type\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process_params:value\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"process_params:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; process_params; process_params = process_params->next) {
		if(XLALFilePrintf(xml->fp, "%s\"%s\",%ld,\"%s\",\"%s\",\"%s\"",
			row_head,
			process_params->program,
			process_params->process_id,
			process_params->param,
			process_params->type,
			process_params->value
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}


/**
 * Write a search_summary table to an XML file.
 */


int XLALWriteLIGOLwXMLSearchSummaryTable(
	LIGOLwXMLStream *xml,
	const SearchSummaryTable *search_summary
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"search_summary:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:shared_object\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:lalwrapper_cvs_tag\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:lal_cvs_tag\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:comment\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:ifos\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:in_start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:in_start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:in_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:in_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:out_start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:out_start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:out_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:out_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:nevents\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search_summary:nnodes\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"search_summary:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; search_summary; search_summary = search_summary->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,\"standalone\",\"\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
			row_head,
			search_summary->process_id,
			lalVCSInfo.vcsTag,
			search_summary->comment,
			search_summary->ifos,
			search_summary->in_start_time.gpsSeconds,
			search_summary->in_start_time.gpsNanoSeconds,
			search_summary->in_end_time.gpsSeconds,
			search_summary->in_end_time.gpsNanoSeconds,
			search_summary->out_start_time.gpsSeconds,
			search_summary->out_start_time.gpsNanoSeconds,
			search_summary->out_end_time.gpsSeconds,
			search_summary->out_end_time.gpsNanoSeconds,
			search_summary->nevents,
			search_summary->nnodes
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}


/**
 * Write a sngl_burst table to an XML file.
 */


int XLALWriteLIGOLwXMLSnglBurstTable(
	LIGOLwXMLStream *xml,
	const SnglBurst *sngl_burst
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sngl_burst:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:ifo\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:search\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:channel\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:peak_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:peak_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:duration\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:central_freq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:bandwidth\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:amplitude\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:snr\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:confidence\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:chisq\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:chisq_dof\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_burst:event_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"sngl_burst:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; sngl_burst; sngl_burst = sngl_burst->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.16g,%.16g,%ld",
			row_head,
			sngl_burst->process_id,
			sngl_burst->ifo,
			sngl_burst->search,
			sngl_burst->channel,
			sngl_burst->start_time.gpsSeconds,
			sngl_burst->start_time.gpsNanoSeconds,
			sngl_burst->peak_time.gpsSeconds,
			sngl_burst->peak_time.gpsNanoSeconds,
			sngl_burst->duration,
			sngl_burst->central_freq,
			sngl_burst->bandwidth,
			sngl_burst->amplitude,
			sngl_burst->snr,
			sngl_burst->confidence,
			sngl_burst->chisq,
			sngl_burst->chisq_dof,
			sngl_burst->event_id
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}

int XLALWriteLIGOLwXMLSnglInspiralTable(
	LIGOLwXMLStream *xml,
	const SnglInspiralTable *sngl_inspiral
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sngl_inspiral:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:ifo\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:search\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:channel\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:end_time_gmst\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:impulse_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:impulse_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:template_duration\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:event_duration\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:amplitude\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:eff_distance\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:coa_phase\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:mass1\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:mass2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:mchirp\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:mtotal\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:eta\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:kappa\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:chi\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:tau0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:tau2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:tau3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:tau4\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:tau5\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:ttotal\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:psi0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:psi3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:alpha\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:alpha1\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:alpha2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:alpha3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:alpha4\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:alpha5\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:alpha6\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:beta\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:f_final\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:snr\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:chisq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:chisq_dof\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:bank_chisq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:bank_chisq_dof\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:cont_chisq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:cont_chisq_dof\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:sigmasq\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:rsqveto_duration\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma1\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma4\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma5\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma6\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma7\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma8\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:Gamma9\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:spin1x\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:spin1y\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:spin1z\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:spin2x\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:spin2y\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:spin2z\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sngl_inspiral:event_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"sngl_inspiral:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; sngl_inspiral; sngl_inspiral = sngl_inspiral->next) {
		if( XLALFilePrintf(xml->fp,"%s%ld,\"%s\",\"%s\",\"%s\",%d,%d,%.16g,%d,%d,%.16g,%.16g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%u,%.8g,%u,%.8g,%u,%.16g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%ld",
			   row_head,
			   sngl_inspiral->process_id,
			   sngl_inspiral->ifo,
			   sngl_inspiral->search,
			   sngl_inspiral->channel,
			   sngl_inspiral->end.gpsSeconds,
			   sngl_inspiral->end.gpsNanoSeconds,
			   sngl_inspiral->end_time_gmst,
			   sngl_inspiral->impulse_time.gpsSeconds,
			   sngl_inspiral->impulse_time.gpsNanoSeconds,
			   sngl_inspiral->template_duration,
			   sngl_inspiral->event_duration,
			   sngl_inspiral->amplitude,
			   sngl_inspiral->eff_distance,
			   sngl_inspiral->coa_phase,
			   sngl_inspiral->mass1,
			   sngl_inspiral->mass2,
			   sngl_inspiral->mchirp,
			   sngl_inspiral->mtotal,
			   sngl_inspiral->eta,
			   sngl_inspiral->kappa,
			   sngl_inspiral->chi,
			   sngl_inspiral->tau0,
			   sngl_inspiral->tau2,
			   sngl_inspiral->tau3,
			   sngl_inspiral->tau4,
			   sngl_inspiral->tau5,
			   sngl_inspiral->ttotal,
			   sngl_inspiral->psi0,
			   sngl_inspiral->psi3,
			   sngl_inspiral->alpha,
			   sngl_inspiral->alpha1,
			   sngl_inspiral->alpha2,
			   sngl_inspiral->alpha3,
			   sngl_inspiral->alpha4,
			   sngl_inspiral->alpha5,
			   sngl_inspiral->alpha6,
			   sngl_inspiral->beta,
			   sngl_inspiral->f_final,
			   sngl_inspiral->snr,
			   sngl_inspiral->chisq,
			   sngl_inspiral->chisq_dof,
			   sngl_inspiral->bank_chisq,
			   sngl_inspiral->bank_chisq_dof,
			   sngl_inspiral->cont_chisq,
			   sngl_inspiral->cont_chisq_dof,
			   sngl_inspiral->sigmasq,
			   sngl_inspiral->rsqveto_duration,
			   sngl_inspiral->Gamma[0],
			   sngl_inspiral->Gamma[1],
			   sngl_inspiral->Gamma[2],
			   sngl_inspiral->Gamma[3],
			   sngl_inspiral->Gamma[4],
			   sngl_inspiral->Gamma[5],
			   sngl_inspiral->Gamma[6],
			   sngl_inspiral->Gamma[7],
			   sngl_inspiral->Gamma[8],
			   sngl_inspiral->Gamma[9],
			   sngl_inspiral->spin1x,  
			   sngl_inspiral->spin1y,
			   sngl_inspiral->spin1z,
			   sngl_inspiral->spin2x,
			   sngl_inspiral->spin2y,
			   sngl_inspiral->spin2z,
			   sngl_inspiral->event_id ) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */
	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */
	return 0;
}




/**
 * Write a sim_burst table to an XML file.
 */


int XLALWriteLIGOLwXMLSimBurstTable(
	LIGOLwXMLStream *xml,
	const SimBurst *sim_burst
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sim_burst:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:waveform\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:ra\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:dec\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:psi\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:time_geocent_gps\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:time_geocent_gps_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:time_geocent_gmst\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:duration\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:frequency\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:bandwidth\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:q\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:pol_ellipse_angle\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:pol_ellipse_e\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:amplitude\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:hrss\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:egw_over_rsquared\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:waveform_number\" Type=\"int_8u\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:time_slide_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sim_burst:simulation_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"sim_burst:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; sim_burst; sim_burst = sim_burst->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,\"%s\",%.16g,%.16g,%.16g,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%lu,%ld,%ld",
			row_head,
			sim_burst->process_id,
			sim_burst->waveform,
			sim_burst->ra,
			sim_burst->dec,
			sim_burst->psi,
			sim_burst->time_geocent_gps.gpsSeconds,
			sim_burst->time_geocent_gps.gpsNanoSeconds,
			sim_burst->time_geocent_gmst,
			sim_burst->duration,
			sim_burst->frequency,
			sim_burst->bandwidth,
			sim_burst->q,
			sim_burst->pol_ellipse_angle,
			sim_burst->pol_ellipse_e,
			sim_burst->amplitude,
			sim_burst->hrss,
			sim_burst->egw_over_rsquared,
			sim_burst->waveform_number,
			sim_burst->time_slide_id,
			sim_burst->simulation_id
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}


/**
 * Write a time_slide table to an XML file.
 */


int XLALWriteLIGOLwXMLTimeSlideTable(
	LIGOLwXMLStream *xml,
	const TimeSlide *time_slide
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"time_slide:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:time_slide_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:instrument\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:offset\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"time_slide:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; time_slide; time_slide = time_slide->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,%ld,\"%s\",%.16g",
			row_head,
			time_slide->process_id,
			time_slide->time_slide_id,
			time_slide->instrument,
			time_slide->offset
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}


/**
 * Write a segment table to an XML file.
 */

int XLALWriteLIGOLwXMLSegmentTable(
	LIGOLwXMLStream *xml,
	const SegmentTable *segment_table
)
{
	const char *row_head = "\n\t\t\t";

	if(xml->table != no_table) {
		XLALPrintError("a table is still open");
		XLAL_ERROR(XLAL_EFAILED);
	}

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"segment:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment:creator_db\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment:segment_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment:start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment:start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment:end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment:end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment_definer:segment_def_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment:segment_def_cdb\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"segment:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; segment_table; segment_table = segment_table->next) {
		if(XLALFilePrintf(xml->fp, "%s%d,%ld,%ld,%d,%d,%d,%d,%ld,%d",
			row_head,
			segment_table->creator_db,
			segment_table->process_id,
			segment_table->segment_id,
			segment_table->start_time.gpsSeconds,
			segment_table->start_time.gpsNanoSeconds,
			segment_table->end_time.gpsSeconds,
			segment_table->end_time.gpsNanoSeconds,
			segment_table->segment_def_id,
			segment_table->segment_def_cdb
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}

/**
 * Creates a XML filename accordingly to document T050017
 */

int XLALCreateLIGODataFileName(
        char* filename,
        size_t size,
        const char* dataSource,
        const char* dataDescription,
        const LIGOTimeGPS* gpsStartTime,
        const LIGOTimeGPS* gpsEndTime,
        const char* extension
)
{
     INT4 gpsDuration;

     /* check input structures */
     if (!filename || !dataSource || !dataDescription ||
	 !gpsStartTime || !gpsEndTime || !extension)
          XLAL_ERROR(XLAL_EFAULT);

     /* check the correctnes of the input strings */
     if ( strchr(dataSource, '-') || strchr(dataDescription, '-'))
     {
          filename = NULL;
          XLALPrintError("the input character strings contain invalid"
			 " dashes ('-').");
          XLAL_ERROR(XLAL_EINVAL);
      }

      /* calculate the GPS duration */
      gpsDuration = gpsEndTime->gpsSeconds - gpsStartTime->gpsSeconds;
      if (gpsEndTime->gpsNanoSeconds > 0) ++gpsDuration;

      /* and here put it all together */
      snprintf( filename, size, "%s-%s-%d-%d.%s",
		   dataSource, dataDescription, gpsStartTime->gpsSeconds,
		   gpsDuration, extension );

      /* return success */
      return 0;
}
