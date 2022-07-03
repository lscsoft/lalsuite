/*
 * sngl_ringdown.c
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
 * Create a SnglRingdownTable structure.
 */
SnglRingdownTable *XLALCreateSnglRingdownTableRow(
	const ProcessTable *process
)
{
	SnglRingdownTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	memset(new->ifo, 0, sizeof(new->ifo));
	memset(new->channel, 0, sizeof(new->channel));
	if(process) {
		/*new->process_id = process->process_id;*/
	} else {
		/*new->process_id = -1; */
	}

	return new;
}


/**
 * Destroy a SnglRingdownTable structure.
 */
void XLALDestroySnglRingdownTableRow(
	SnglRingdownTable *row
)
{
	XLALFree(row);
}


/**
 * Destroy a SnglRingdownTable linked list.
 */
void XLALDestroySnglRingdownTable(
	SnglRingdownTable *head
)
{
	while(head) {
		SnglRingdownTable *next = head->next;
		XLALDestroySnglRingdownTableRow(head);
		head = next;
	}
}


/**
 * Read the sngl_ringdown table from a LIGO Light Weight XML file into a
 * linked list of SnglRingdownTable structures.
 */
SnglRingdownTable *XLALSnglRingdownTableFromLIGOLw (
	const char *filename
)
{
	static const char table_name[] = "sngl_ringdown";
	int miostatus;
	SnglRingdownTable *head = NULL;
	SnglRingdownTable **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int ifo;
		int channel;
		int start_time;
		int start_time_ns;
		int start_time_gmst;
		int frequency;
		int quality;
		int phase;
		int mass;
		int spin;
		int epsilon;
		int num_clust_trigs;
		int ds2_H1H2;
		int ds2_H1L1;
		int ds2_H1V1;
		int ds2_H2L1;
		int ds2_H2V1;
		int ds2_L1V1;
		int amplitude;
		int snr;
		int eff_dist;
		int sigma_sq;
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
	column_pos.channel = XLALLIGOLwFindColumn(&env, "channel", METAIO_TYPE_LSTRING, 1);
	column_pos.start_time = XLALLIGOLwFindColumn(&env, "start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.start_time_ns = XLALLIGOLwFindColumn(&env, "start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.start_time_gmst = XLALLIGOLwFindColumn(&env, "start_time_gmst", METAIO_TYPE_REAL_8, 1);
	column_pos.frequency = XLALLIGOLwFindColumn(&env, "frequency", METAIO_TYPE_REAL_4, 1);
	column_pos.quality = XLALLIGOLwFindColumn(&env, "quality", METAIO_TYPE_REAL_4, 1);
	column_pos.phase = XLALLIGOLwFindColumn(&env, "phase", METAIO_TYPE_REAL_4, 1);
	column_pos.mass = XLALLIGOLwFindColumn(&env, "mass", METAIO_TYPE_REAL_4, 1);
	column_pos.spin = XLALLIGOLwFindColumn(&env, "spin", METAIO_TYPE_REAL_4, 1);
	column_pos.epsilon = XLALLIGOLwFindColumn(&env, "epsilon", METAIO_TYPE_REAL_4, 1);
	column_pos.num_clust_trigs = XLALLIGOLwFindColumn(&env, "num_clust_trigs", METAIO_TYPE_INT_4S, 1);
	column_pos.ds2_H1H2 = XLALLIGOLwFindColumn(&env, "ds2_H1H2", METAIO_TYPE_REAL_4, 1);
	column_pos.ds2_H1L1 = XLALLIGOLwFindColumn(&env, "ds2_H1L1", METAIO_TYPE_REAL_4, 1);
	column_pos.ds2_H1V1 = XLALLIGOLwFindColumn(&env, "ds2_H1V1", METAIO_TYPE_REAL_4, 1);
	column_pos.ds2_H2L1 = XLALLIGOLwFindColumn(&env, "ds2_H2L1", METAIO_TYPE_REAL_4, 1);
	column_pos.ds2_H2V1 = XLALLIGOLwFindColumn(&env, "ds2_H2V1", METAIO_TYPE_REAL_4, 1);
	column_pos.ds2_L1V1 = XLALLIGOLwFindColumn(&env, "ds2_L1V1", METAIO_TYPE_REAL_4, 1);
	column_pos.amplitude = XLALLIGOLwFindColumn(&env, "amplitude", METAIO_TYPE_REAL_4, 1);
	column_pos.snr = XLALLIGOLwFindColumn(&env, "snr", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist = XLALLIGOLwFindColumn(&env, "eff_dist", METAIO_TYPE_REAL_4, 1);
	column_pos.sigma_sq = XLALLIGOLwFindColumn(&env, "sigma_sq", METAIO_TYPE_REAL_4, 1);
	column_pos.event_id = XLALLIGOLwFindColumn(&env, "event_id", METAIO_TYPE_INT_8S, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
        }

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		SnglRingdownTable *row = XLALCreateSnglRingdownTableRow(NULL);

		if(!row) {
			XLALDestroySnglRingdownTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		/*row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;*/
		strncpy(row->ifo, env.ligo_lw.table.elt[column_pos.ifo].data.lstring.data, sizeof(row->ifo) - 1);
		strncpy(row->channel, env.ligo_lw.table.elt[column_pos.channel].data.lstring.data, sizeof(row->channel) - 1);
		XLALGPSSet(&row->start_time, env.ligo_lw.table.elt[column_pos.start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.start_time_ns].data.int_4s);
		row->start_time_gmst = env.ligo_lw.table.elt[column_pos.start_time_gmst].data.real_8;
		row->frequency = env.ligo_lw.table.elt[column_pos.frequency].data.real_4;
		row->quality = env.ligo_lw.table.elt[column_pos.quality].data.real_4;
		row->phase = env.ligo_lw.table.elt[column_pos.phase].data.real_4;
		row->mass = env.ligo_lw.table.elt[column_pos.mass].data.real_4;
		row->spin = env.ligo_lw.table.elt[column_pos.spin].data.real_4;
		row->epsilon = env.ligo_lw.table.elt[column_pos.epsilon].data.real_4;
		row->num_clust_trigs = env.ligo_lw.table.elt[column_pos.num_clust_trigs].data.int_4s;
		row->ds2_H1H2 = env.ligo_lw.table.elt[column_pos.ds2_H1H2].data.real_4;
		row->ds2_H1L1 = env.ligo_lw.table.elt[column_pos.ds2_H1L1].data.real_4;
		row->ds2_H1V1 = env.ligo_lw.table.elt[column_pos.ds2_H1V1].data.real_4;
		row->ds2_H2L1 = env.ligo_lw.table.elt[column_pos.ds2_H2L1].data.real_4;
		row->ds2_H2V1 = env.ligo_lw.table.elt[column_pos.ds2_H2V1].data.real_4;
		row->ds2_L1V1 = env.ligo_lw.table.elt[column_pos.ds2_L1V1].data.real_4;
		row->amplitude = env.ligo_lw.table.elt[column_pos.amplitude].data.real_4;
		row->snr = env.ligo_lw.table.elt[column_pos.snr].data.real_4;
		row->eff_dist = env.ligo_lw.table.elt[column_pos.eff_dist].data.real_4;
		row->sigma_sq = env.ligo_lw.table.elt[column_pos.sigma_sq].data.real_4;
		row->event_id = env.ligo_lw.table.elt[column_pos.event_id].data.int_8s;
	}
	if(miostatus < 0) {
		XLALDestroySnglRingdownTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroySnglRingdownTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}
