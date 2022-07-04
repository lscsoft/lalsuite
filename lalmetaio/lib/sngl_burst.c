/*
 * sngl_burst.c
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <metaio.h>


#include <lal/Date.h>
#include <lal/LALMalloc.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/XLALError.h>


/**
 * Create a SnglBurst structure.
 */
SnglBurst *XLALCreateSnglBurst(void)
{
	SnglBurst *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	new->process_id = new->event_id = -1;
	memset(new->ifo, 0, sizeof(new->ifo));
	memset(new->search, 0, sizeof(new->search));
	memset(new->channel, 0, sizeof(new->channel));
	XLALGPSSet(&new->start_time, 0, 0);
	XLALGPSSet(&new->peak_time, 0, 0);
	new->duration = XLAL_REAL4_FAIL_NAN;
	new->central_freq = XLAL_REAL4_FAIL_NAN;
	new->bandwidth = XLAL_REAL4_FAIL_NAN;
	new->amplitude = XLAL_REAL4_FAIL_NAN;
	new->snr = XLAL_REAL4_FAIL_NAN;
	new->confidence = XLAL_REAL4_FAIL_NAN;
	new->chisq = XLAL_REAL8_FAIL_NAN;
	new->chisq_dof = XLAL_REAL8_FAIL_NAN;

	return new;
}


/**
 * Free a SnglBurst.
 */
void XLALDestroySnglBurst(SnglBurst *event)
{
	XLALFree(event);
}


/**
 * Free a SnglBurst linked list.
 */
void XLALDestroySnglBurstTable(SnglBurst *head)
{
	while(head) {
		SnglBurst *next = head->next;
		XLALDestroySnglBurst(head);
		head = next;
	}
}


/**
 * Read the sngl_burst table from a LIGO Light Weight XML file into a
 * linked list of SnglBurst structures.
 */
SnglBurst *XLALSnglBurstTableFromLIGOLw(
	const char *filename
)
{
	static const char table_name[] = "sngl_burst";
	int miostatus;
	SnglBurst *head = NULL;
	SnglBurst **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int ifo;
		int search;
		int channel;
		int start_time;
		int start_time_ns;
		int peak_time;
		int peak_time_ns;
		int duration;
		int central_freq;
		int bandwidth;
		int amplitude;
		int snr;
		int confidence;
		int chisq;
		int chisq_dof;
		int event_id;
	} column_pos;

	/* open the file and find table */

	if(MetaioOpenFile(&env, filename)) {
		XLALPrintError("%s(): error opening \"%s\": %s\n", __func__, filename, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}
	if(MetaioOpenTableOnly(&env, table_name)) {
		MetaioAbort(&env);
		XLALPrintError("%s(): cannot find %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* find columns */

	XLALClearErrno();
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_INT_8S, 1);
	column_pos.ifo = XLALLIGOLwFindColumn(&env, "ifo", METAIO_TYPE_LSTRING, 1);
	column_pos.search = XLALLIGOLwFindColumn(&env, "search", METAIO_TYPE_LSTRING, 1);
	column_pos.channel = XLALLIGOLwFindColumn(&env, "channel", METAIO_TYPE_LSTRING, 1);
	column_pos.start_time = XLALLIGOLwFindColumn(&env, "start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.start_time_ns = XLALLIGOLwFindColumn(&env, "start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.peak_time = XLALLIGOLwFindColumn(&env, "peak_time", METAIO_TYPE_INT_4S, 1);
	column_pos.peak_time_ns = XLALLIGOLwFindColumn(&env, "peak_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.duration = XLALLIGOLwFindColumn(&env, "duration", METAIO_TYPE_REAL_4, 1);
	column_pos.central_freq = XLALLIGOLwFindColumn(&env, "central_freq", METAIO_TYPE_REAL_4, 1);
	column_pos.bandwidth = XLALLIGOLwFindColumn(&env, "bandwidth", METAIO_TYPE_REAL_4, 1);
	column_pos.amplitude = XLALLIGOLwFindColumn(&env, "amplitude", METAIO_TYPE_REAL_4, 1);
	column_pos.snr = XLALLIGOLwFindColumn(&env, "snr", METAIO_TYPE_REAL_4, 1);
	column_pos.confidence = XLALLIGOLwFindColumn(&env, "confidence", METAIO_TYPE_REAL_4, 1);
	column_pos.chisq = XLALLIGOLwFindColumn(&env, "chisq", METAIO_TYPE_REAL_8, 1);
	column_pos.chisq_dof = XLALLIGOLwFindColumn(&env, "chisq_dof", METAIO_TYPE_REAL_8, 1);
	column_pos.event_id = XLALLIGOLwFindColumn(&env, "event_id", METAIO_TYPE_INT_8S, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		SnglBurst *row = XLALCreateSnglBurst();

		if(!row) {
			XLALDestroySnglBurstTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
		if(strlen(env.ligo_lw.table.elt[column_pos.ifo].data.lstring.data) >= sizeof(row->ifo) ||
		strlen(env.ligo_lw.table.elt[column_pos.search].data.lstring.data) >= sizeof(row->search) ||
		strlen(env.ligo_lw.table.elt[column_pos.channel].data.lstring.data) >= sizeof(row->channel)) {
			XLALDestroySnglBurstTable(head);
			MetaioAbort(&env);
			XLALPrintError("%s(): failure reading %s table: string too long\n", __func__, table_name);
			XLAL_ERROR_NULL(XLAL_EIO);
		}
		strncpy(row->ifo, env.ligo_lw.table.elt[column_pos.ifo].data.lstring.data, sizeof(row->ifo) - 1);
		strncpy(row->search, env.ligo_lw.table.elt[column_pos.search].data.lstring.data, sizeof(row->search) - 1);
		strncpy(row->channel, env.ligo_lw.table.elt[column_pos.channel].data.lstring.data, sizeof(row->channel) - 1);
		XLALGPSSet(&row->start_time, env.ligo_lw.table.elt[column_pos.start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.start_time_ns].data.int_4s);
		XLALGPSSet(&row->peak_time, env.ligo_lw.table.elt[column_pos.peak_time].data.int_4s, env.ligo_lw.table.elt[column_pos.peak_time_ns].data.int_4s);
		row->duration = env.ligo_lw.table.elt[column_pos.duration].data.real_4;
		row->central_freq = env.ligo_lw.table.elt[column_pos.central_freq].data.real_4;
		row->bandwidth = env.ligo_lw.table.elt[column_pos.bandwidth].data.real_4;
		row->amplitude = env.ligo_lw.table.elt[column_pos.amplitude].data.real_4;
		row->snr = env.ligo_lw.table.elt[column_pos.snr].data.real_4;
		row->confidence = env.ligo_lw.table.elt[column_pos.confidence].data.real_4;
		row->chisq = env.ligo_lw.table.elt[column_pos.chisq].data.real_8;
		row->chisq_dof = env.ligo_lw.table.elt[column_pos.chisq_dof].data.real_8;
		row->event_id = env.ligo_lw.table.elt[column_pos.event_id].data.int_8s;
	}
	if(miostatus < 0) {
		XLALDestroySnglBurstTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroySnglBurstTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
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

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sngl_burst:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"ifo\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"channel\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"peak_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"peak_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"duration\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"central_freq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"bandwidth\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"amplitude\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"snr\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"confidence\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"chisq\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"chisq_dof\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"event_id\" Type=\"int_8s\"/>\n", xml->fp);
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
