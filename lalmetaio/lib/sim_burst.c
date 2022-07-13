/*
 * sim_burst.c
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
 * Create a SimBurst structure.
 */
SimBurst *XLALCreateSimBurst(void)
{
	SimBurst *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	new->process_id = new->time_slide_id = new->simulation_id = -1;
	memset(new->waveform, 0, sizeof(new->waveform));
	new->ra = XLAL_REAL8_FAIL_NAN;
	new->dec = XLAL_REAL8_FAIL_NAN;
	new->psi = XLAL_REAL8_FAIL_NAN;
	XLALGPSSet(&new->time_geocent_gps, 0, 0);
	new->time_geocent_gmst = XLAL_REAL8_FAIL_NAN;
	new->duration = XLAL_REAL8_FAIL_NAN;
	new->frequency = XLAL_REAL8_FAIL_NAN;
	new->bandwidth = XLAL_REAL8_FAIL_NAN;
	new->q = XLAL_REAL8_FAIL_NAN;
	new->pol_ellipse_angle = XLAL_REAL8_FAIL_NAN;
	new->pol_ellipse_e= XLAL_REAL8_FAIL_NAN;
	new->amplitude = XLAL_REAL8_FAIL_NAN;
	new->hrss = XLAL_REAL8_FAIL_NAN;
	new->egw_over_rsquared = XLAL_REAL8_FAIL_NAN;
	new->waveform_number = 0;

	return new;
}


/**
 * Destroy a SimBurst structure.
 */
void XLALDestroySimBurst(SimBurst *sim_burst)
{
	XLALFree(sim_burst);
}


/**
 * Destroy a SimBurst linked list.
 */
void XLALDestroySimBurstTable(SimBurst *head)
{
	while(head) {
		SimBurst *next = head->next;
		XLALDestroySimBurst(head);
		head = next;
	}
}


/**
 * Read the sim_burst table from a LIGO Light Weight XML file into a linked
 * list of SimBurst structures.
 */
SimBurst *XLALSimBurstTableFromLIGOLw(
	const char *filename
)
{
	static const char table_name[] = "sim_burst";
	int miostatus;
	SimBurst *head = NULL;
	SimBurst **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int waveform;
		int ra;
		int dec;
		int psi;
		int time_geocent_gps;
		int time_geocent_gps_ns;
		int time_geocent_gmst;
		int duration;
		int frequency;
		int bandwidth;
		int q;
		int pol_ellipse_angle;
		int pol_ellipse_e;
		int amplitude;
		int hrss;
		int egw_over_rsquared;
		int waveform_number;
		int time_slide_id;
		int simulation_id;
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
	column_pos.waveform = XLALLIGOLwFindColumn(&env, "waveform", METAIO_TYPE_LSTRING, 1);
	column_pos.ra = XLALLIGOLwFindColumn(&env, "ra", METAIO_TYPE_REAL_8, 0);
	column_pos.dec = XLALLIGOLwFindColumn(&env, "dec", METAIO_TYPE_REAL_8, 0);
	column_pos.psi = XLALLIGOLwFindColumn(&env, "psi", METAIO_TYPE_REAL_8, 0);
	column_pos.time_geocent_gps = XLALLIGOLwFindColumn(&env, "time_geocent_gps", METAIO_TYPE_INT_4S, 1);
	column_pos.time_geocent_gps_ns = XLALLIGOLwFindColumn(&env, "time_geocent_gps_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.time_geocent_gmst = XLALLIGOLwFindColumn(&env, "time_geocent_gmst", METAIO_TYPE_REAL_8, 0);
	column_pos.duration = XLALLIGOLwFindColumn(&env, "duration", METAIO_TYPE_REAL_8, 0);
	column_pos.frequency = XLALLIGOLwFindColumn(&env, "frequency", METAIO_TYPE_REAL_8, 0);
	column_pos.bandwidth = XLALLIGOLwFindColumn(&env, "bandwidth", METAIO_TYPE_REAL_8, 0);
	column_pos.q = XLALLIGOLwFindColumn(&env, "q", METAIO_TYPE_REAL_8, 0);
	column_pos.pol_ellipse_angle = XLALLIGOLwFindColumn(&env, "pol_ellipse_angle", METAIO_TYPE_REAL_8, 0);
	column_pos.pol_ellipse_e = XLALLIGOLwFindColumn(&env, "pol_ellipse_e", METAIO_TYPE_REAL_8, 0);
	column_pos.amplitude = XLALLIGOLwFindColumn(&env, "amplitude", METAIO_TYPE_REAL_8, 0);
	column_pos.hrss = XLALLIGOLwFindColumn(&env, "hrss", METAIO_TYPE_REAL_8, 0);
	column_pos.egw_over_rsquared = XLALLIGOLwFindColumn(&env, "egw_over_rsquared", METAIO_TYPE_REAL_8, 0);
	column_pos.waveform_number = XLALLIGOLwFindColumn(&env, "waveform_number", METAIO_TYPE_INT_8U, 0);
	column_pos.time_slide_id = XLALLIGOLwFindColumn(&env, "time_slide_id", METAIO_TYPE_INT_8S, 1);
	column_pos.simulation_id = XLALLIGOLwFindColumn(&env, "simulation_id", METAIO_TYPE_INT_8S, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		SimBurst *row = XLALCreateSimBurst();

		if(!row) {
			XLALDestroySimBurstTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* populate the columns */

		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
		if(strlen(env.ligo_lw.table.elt[column_pos.waveform].data.lstring.data) >= sizeof(row->waveform)) {
			XLALDestroySimBurst(row);
			XLALDestroySimBurstTable(head);
			MetaioAbort(&env);
			XLALPrintError("%s(): failure reading %s table: string too long\n", __func__, table_name);
			XLAL_ERROR_NULL(XLAL_EIO);
		}
		strncpy(row->waveform, env.ligo_lw.table.elt[column_pos.waveform].data.lstring.data, sizeof(row->waveform) - 1);
		if(column_pos.ra >= 0)
			row->ra = env.ligo_lw.table.elt[column_pos.ra].data.real_8;
		if(column_pos.dec >= 0)
			row->dec = env.ligo_lw.table.elt[column_pos.dec].data.real_8;
		if(column_pos.psi >= 0)
			row->psi = env.ligo_lw.table.elt[column_pos.psi].data.real_8;
		XLALGPSSet(&row->time_geocent_gps, env.ligo_lw.table.elt[column_pos.time_geocent_gps].data.int_4s, env.ligo_lw.table.elt[column_pos.time_geocent_gps_ns].data.int_4s);
		if(column_pos.time_geocent_gmst >= 0)
			row->time_geocent_gmst = env.ligo_lw.table.elt[column_pos.time_geocent_gmst].data.real_8;
		row->time_slide_id = env.ligo_lw.table.elt[column_pos.time_slide_id].data.int_8s;
		row->simulation_id = env.ligo_lw.table.elt[column_pos.simulation_id].data.int_8s;

		if(!strcmp(row->waveform, "StringCusp")) {
			if(column_pos.duration < 0 || column_pos.frequency < 0 || column_pos.amplitude < 0) {
				XLALDestroySimBurst(row);
				XLALDestroySimBurstTable(head);
				MetaioAbort(&env);
				XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
				XLAL_ERROR_NULL(XLAL_EIO);
			}
			row->duration = env.ligo_lw.table.elt[column_pos.duration].data.real_8;
			row->frequency = env.ligo_lw.table.elt[column_pos.frequency].data.real_8;
			row->amplitude = env.ligo_lw.table.elt[column_pos.amplitude].data.real_8;
		} else if(!strcmp(row->waveform, "SineGaussian")||!strcmp(row->waveform, "SineGaussianF")) {
			if(column_pos.duration < 0 || column_pos.frequency < 0 || column_pos.bandwidth < 0 || column_pos.q < 0 || column_pos.pol_ellipse_angle < 0 || column_pos.pol_ellipse_e < 0 || column_pos.hrss < 0) {
				XLALDestroySimBurst(row);
				XLALDestroySimBurstTable(head);
				MetaioAbort(&env);
				XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
				XLAL_ERROR_NULL(XLAL_EIO);
			}
			row->duration = env.ligo_lw.table.elt[column_pos.duration].data.real_8;
			row->frequency = env.ligo_lw.table.elt[column_pos.frequency].data.real_8;
			row->bandwidth = env.ligo_lw.table.elt[column_pos.bandwidth].data.real_8;
			row->q = env.ligo_lw.table.elt[column_pos.q].data.real_8;
			row->pol_ellipse_angle = env.ligo_lw.table.elt[column_pos.pol_ellipse_angle].data.real_8;
			row->pol_ellipse_e = env.ligo_lw.table.elt[column_pos.pol_ellipse_e].data.real_8;
			row->hrss = env.ligo_lw.table.elt[column_pos.hrss].data.real_8;
		} else if(!strcmp(row->waveform, "Gaussian")) {
			if(column_pos.duration < 0 || column_pos.hrss < 0) {
				XLALDestroySimBurst(row);
				XLALDestroySimBurstTable(head);
				MetaioAbort(&env);
				XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
				XLAL_ERROR_NULL(XLAL_EIO);
			}
			row->duration = env.ligo_lw.table.elt[column_pos.duration].data.real_8;
			row->hrss = env.ligo_lw.table.elt[column_pos.hrss].data.real_8;
		} else if(!strcmp(row->waveform, "BTLWNB")) {
			if(column_pos.duration < 0 || column_pos.frequency < 0 || column_pos.bandwidth < 0 || column_pos.pol_ellipse_angle < 0 || column_pos.pol_ellipse_e < 0 || column_pos.egw_over_rsquared < 0 || column_pos.waveform_number < 0) {
				XLALDestroySimBurst(row);
				XLALDestroySimBurstTable(head);
				MetaioAbort(&env);
				XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
				XLAL_ERROR_NULL(XLAL_EIO);
			}
			row->duration = env.ligo_lw.table.elt[column_pos.duration].data.real_8;
			row->frequency = env.ligo_lw.table.elt[column_pos.frequency].data.real_8;
			row->bandwidth = env.ligo_lw.table.elt[column_pos.bandwidth].data.real_8;
			row->pol_ellipse_angle = env.ligo_lw.table.elt[column_pos.pol_ellipse_angle].data.real_8;
			row->pol_ellipse_e = env.ligo_lw.table.elt[column_pos.pol_ellipse_e].data.real_8;
			row->egw_over_rsquared = env.ligo_lw.table.elt[column_pos.egw_over_rsquared].data.real_8;
			row->waveform_number = env.ligo_lw.table.elt[column_pos.waveform_number].data.int_8u;
		} else if(!strcmp(row->waveform, "Impulse")) {
			if(column_pos.amplitude < 0) {
				XLALDestroySimBurst(row);
				XLALDestroySimBurstTable(head);
				MetaioAbort(&env);
				XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
				XLAL_ERROR_NULL(XLAL_EIO);
			}
			row->amplitude = env.ligo_lw.table.elt[column_pos.amplitude].data.real_8;
		} else {
			/* unrecognized waveform */
			XLALDestroySimBurst(row);
			XLALDestroySimBurstTable(head);
			MetaioAbort(&env);
			XLALPrintError("%s(): unrecognized waveform \"%s\" in %s table\n", __func__, row->waveform, table_name);
			XLAL_ERROR_NULL(XLAL_EIO);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;
	}
	if(miostatus < 0) {
		XLALDestroySimBurstTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroySimBurstTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
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

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sim_burst:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"waveform\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"ra\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"dec\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"psi\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_geocent_gps\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_geocent_gps_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_geocent_gmst\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"duration\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"frequency\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"bandwidth\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"q\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"pol_ellipse_angle\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"pol_ellipse_e\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"amplitude\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"hrss\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"egw_over_rsquared\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"waveform_number\" Type=\"int_8u\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:time_slide_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"simulation_id\" Type=\"int_8s\"/>\n", xml->fp);
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
 * Assign simulation_id values to the entries in a sim_burst linked list.
 * All sim_burst rows in the list will be blamed on the given process_id,
 * and assigned simulation_ids in order starting with the given
 * simulation_id.  The return value is the next simulation_id after the
 * last one assigned to a row in the list.
 */
long XLALSimBurstAssignIDs(
	SimBurst *sim_burst,
	long process_id,
	long time_slide_id,
	long simulation_id)
{
	for(; sim_burst; sim_burst = sim_burst->next) {
		sim_burst->process_id = process_id;
		sim_burst->time_slide_id = time_slide_id;
		sim_burst->simulation_id = simulation_id++;
	}
	return simulation_id;
}
