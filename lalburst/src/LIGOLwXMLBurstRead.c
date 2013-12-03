/*
 * LIGOLwXMLBurstRead.c
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
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <metaio.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataBurstUtils.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/Date.h>

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
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_ILWD_CHAR, 1);
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
	column_pos.event_id = XLALLIGOLwFindColumn(&env, "event_id", METAIO_TYPE_ILWD_CHAR, 1);

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

		if((row->process_id = XLALLIGOLwParseIlwdChar(&env, column_pos.process_id, "process", "process_id")) < 0) {
			XLALDestroySnglBurstTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
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
		if((row->event_id = XLALLIGOLwParseIlwdChar(&env, column_pos.event_id, "sngl_burst", "event_id")) < 0) {
			XLALDestroySnglBurstTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
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
 * Read the sim_burst table from a LIGO Light Weight XML file into a linked
 * list of SimBurst structures.  If start is not NULL, then only rows whose
 * geocentre peak times are \f$\ge\f$ the given GPS time will be loaded, similarly
 * if end is not NULL.
 */
SimBurst *XLALSimBurstTableFromLIGOLw(
	const char *filename,
	const LIGOTimeGPS *start,
	const LIGOTimeGPS *end
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
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_ILWD_CHAR, 1);
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
	column_pos.time_slide_id = XLALLIGOLwFindColumn(&env, "time_slide_id", METAIO_TYPE_ILWD_CHAR, 1);
	column_pos.simulation_id = XLALLIGOLwFindColumn(&env, "simulation_id", METAIO_TYPE_ILWD_CHAR, 1);

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

		if((row->process_id = XLALLIGOLwParseIlwdChar(&env, column_pos.process_id, "process", "process_id")) < 0) {
			XLALDestroySimBurst(row);
			XLALDestroySimBurstTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
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
		if((row->time_slide_id = XLALLIGOLwParseIlwdChar(&env, column_pos.time_slide_id, "time_slide", "time_slide_id")) < 0) {
			XLALDestroySimBurst(row);
			XLALDestroySimBurstTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		if((row->simulation_id = XLALLIGOLwParseIlwdChar(&env, column_pos.simulation_id, "sim_burst", "simulation_id")) < 0) {
			XLALDestroySimBurst(row);
			XLALDestroySimBurstTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

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
		} else if(!strcmp(row->waveform, "SineGaussian")) {
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
		} else if(!strcmp(row->waveform, "BTLWNB")) {
			if(column_pos.duration < 0 || column_pos.frequency < 0 || column_pos.bandwidth < 0 || column_pos.egw_over_rsquared < 0 || column_pos.waveform_number < 0) {
				XLALDestroySimBurst(row);
				XLALDestroySimBurstTable(head);
				MetaioAbort(&env);
				XLALPrintError("%s(): failure reading %s table: missing required column\n", __func__, table_name);
				XLAL_ERROR_NULL(XLAL_EIO);
			}
			row->duration = env.ligo_lw.table.elt[column_pos.duration].data.real_8;
			row->frequency = env.ligo_lw.table.elt[column_pos.frequency].data.real_8;
			row->bandwidth = env.ligo_lw.table.elt[column_pos.bandwidth].data.real_8;
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

		/* if outside accepted time window, discard */

		if((start && XLALGPSDiff(start, &row->time_geocent_gps) > 0) || (end && XLALGPSDiff(end, &row->time_geocent_gps) < 0)) {
			XLALDestroySimBurst(row);
			continue;
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
