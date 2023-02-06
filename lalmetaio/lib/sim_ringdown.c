/*
 * sim_ringdown.c
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
 * Create a SimRingdownTable structure.
 */
SimRingdownTable *XLALCreateSimRingdownTableRow(
	const ProcessTable *process
)
{
	SimRingdownTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	memset(new->waveform, 0, sizeof(new->waveform));
	memset(new->coordinates, 0, sizeof(new->coordinates));
	if(process) {
		/*new->process_id = process->process_id;*/
	} else {
		/*new->process_id = -1; */
	}

	return new;
}


/**
 * Destroy a SimRingdownTable structure.
 */
SimRingdownTable *XLALDestroySimRingdownTableRow(
	SimRingdownTable *row
)
{
	SimRingdownTable *next = row ? row->next : NULL;
	XLALFree(row);
	return next;
}


/**
 * Destroy a SimRingdownTable linked list.
 */
void XLALDestroySimRingdownTable(
	SimRingdownTable *head
)
{
	while(head)
		head = XLALDestroySimRingdownTableRow(head);
}


/**
 * Read the sim_ringdown table from a LIGO Light Weight XML file into a
 * linked list of SimRingdownTable structures.
 */
SimRingdownTable *XLALSimRingdownTableFromLIGOLw(
	const char *filename
)
{
	static const char table_name[] = "sim_ringdown";
	int miostatus;
	SimRingdownTable *head = NULL;
	SimRingdownTable **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int waveform;
		int coordinates;
		int geocent_start_time;
		int geocent_start_time_ns;
		int h_start_time;
		int h_start_time_ns;
		int l_start_time;
		int l_start_time_ns;
		int v_start_time;
		int v_start_time_ns;
		int start_time_gmst;
		int longitude;
		int latitude;
		int distance;
		int inclination;
		int polarization;
		int frequency;
		int quality;
		int phase;
		int mass;
		int spin;
		int epsilon;
		int amplitude;
		int eff_dist_h;
		int eff_dist_l;
		int eff_dist_v;
		int hrss;
		int hrss_h;
		int hrss_l;
		int hrss_v;
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
	column_pos.coordinates = XLALLIGOLwFindColumn(&env, "coordinates", METAIO_TYPE_LSTRING, 1);
	column_pos.geocent_start_time = XLALLIGOLwFindColumn(&env, "geocent_start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.geocent_start_time_ns = XLALLIGOLwFindColumn(&env, "geocent_start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.h_start_time = XLALLIGOLwFindColumn(&env, "h_start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.h_start_time_ns = XLALLIGOLwFindColumn(&env, "h_start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.l_start_time = XLALLIGOLwFindColumn(&env, "l_start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.l_start_time_ns = XLALLIGOLwFindColumn(&env, "l_start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.v_start_time = XLALLIGOLwFindColumn(&env, "v_start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.v_start_time_ns = XLALLIGOLwFindColumn(&env, "v_start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.start_time_gmst = XLALLIGOLwFindColumn(&env, "start_time_gmst", METAIO_TYPE_REAL_8, 1);
	column_pos.longitude = XLALLIGOLwFindColumn(&env, "longitude", METAIO_TYPE_REAL_4, 1);
	column_pos.latitude = XLALLIGOLwFindColumn(&env, "latitude", METAIO_TYPE_REAL_4, 1);
	column_pos.distance = XLALLIGOLwFindColumn(&env, "distance", METAIO_TYPE_REAL_4, 1);
	column_pos.inclination = XLALLIGOLwFindColumn(&env, "inclination", METAIO_TYPE_REAL_4, 1);
	column_pos.polarization = XLALLIGOLwFindColumn(&env, "polarization", METAIO_TYPE_REAL_4, 1);
	column_pos.frequency = XLALLIGOLwFindColumn(&env, "frequency", METAIO_TYPE_REAL_4, 1);
	column_pos.quality = XLALLIGOLwFindColumn(&env, "quality", METAIO_TYPE_REAL_4, 1);
	column_pos.phase = XLALLIGOLwFindColumn(&env, "phase", METAIO_TYPE_REAL_4, 1);
	column_pos.mass = XLALLIGOLwFindColumn(&env, "mass", METAIO_TYPE_REAL_4, 1);
	column_pos.spin = XLALLIGOLwFindColumn(&env, "spin", METAIO_TYPE_REAL_4, 1);
	column_pos.epsilon = XLALLIGOLwFindColumn(&env, "epsilon", METAIO_TYPE_REAL_4, 1);
	column_pos.amplitude = XLALLIGOLwFindColumn(&env, "amplitude", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_h = XLALLIGOLwFindColumn(&env, "eff_dist_h", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_l = XLALLIGOLwFindColumn(&env, "eff_dist_l", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_v = XLALLIGOLwFindColumn(&env, "eff_dist_v", METAIO_TYPE_REAL_4, 1);
	column_pos.hrss = XLALLIGOLwFindColumn(&env, "hrss", METAIO_TYPE_REAL_4, 1);
	column_pos.hrss_h = XLALLIGOLwFindColumn(&env, "hrss_h", METAIO_TYPE_REAL_4, 1);
	column_pos.hrss_l = XLALLIGOLwFindColumn(&env, "hrss_l", METAIO_TYPE_REAL_4, 1);
	column_pos.hrss_v = XLALLIGOLwFindColumn(&env, "hrss_v", METAIO_TYPE_REAL_4, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
        }

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		SimRingdownTable *row = XLALCreateSimRingdownTableRow(NULL);

		if(!row) {
			XLALDestroySimRingdownTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		/*row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;*/
		strncpy(row->waveform, env.ligo_lw.table.elt[column_pos.waveform].data.lstring.data, sizeof(row->waveform) - 1);
		strncpy(row->coordinates, env.ligo_lw.table.elt[column_pos.coordinates].data.lstring.data, sizeof(row->coordinates) - 1);
		XLALGPSSet(&row->geocent_start_time, env.ligo_lw.table.elt[column_pos.geocent_start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.geocent_start_time_ns].data.int_4s);
		XLALGPSSet(&row->h_start_time, env.ligo_lw.table.elt[column_pos.h_start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.h_start_time_ns].data.int_4s);
		XLALGPSSet(&row->l_start_time, env.ligo_lw.table.elt[column_pos.l_start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.l_start_time_ns].data.int_4s);
		XLALGPSSet(&row->v_start_time, env.ligo_lw.table.elt[column_pos.v_start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.v_start_time_ns].data.int_4s);
		row->start_time_gmst = env.ligo_lw.table.elt[column_pos.start_time_gmst].data.real_8;
		row->longitude = env.ligo_lw.table.elt[column_pos.longitude].data.real_4;
		row->latitude = env.ligo_lw.table.elt[column_pos.latitude].data.real_4;
		row->distance = env.ligo_lw.table.elt[column_pos.distance].data.real_4;
		row->inclination = env.ligo_lw.table.elt[column_pos.inclination].data.real_4;
		row->polarization = env.ligo_lw.table.elt[column_pos.polarization].data.real_4;
		row->frequency = env.ligo_lw.table.elt[column_pos.frequency].data.real_4;
		row->quality = env.ligo_lw.table.elt[column_pos.quality].data.real_4;
		row->phase = env.ligo_lw.table.elt[column_pos.phase].data.real_4;
		row->mass = env.ligo_lw.table.elt[column_pos.mass].data.real_4;
		row->spin = env.ligo_lw.table.elt[column_pos.spin].data.real_4;
		row->epsilon = env.ligo_lw.table.elt[column_pos.epsilon].data.real_4;
		row->amplitude = env.ligo_lw.table.elt[column_pos.amplitude].data.real_4;
		row->eff_dist_h = env.ligo_lw.table.elt[column_pos.eff_dist_h].data.real_4;
		row->eff_dist_l = env.ligo_lw.table.elt[column_pos.eff_dist_l].data.real_4;
		row->eff_dist_v = env.ligo_lw.table.elt[column_pos.eff_dist_v].data.real_4;
		row->hrss = env.ligo_lw.table.elt[column_pos.hrss].data.real_4;
		row->hrss_h = env.ligo_lw.table.elt[column_pos.hrss_h].data.real_4;
		row->hrss_l = env.ligo_lw.table.elt[column_pos.hrss_l].data.real_4;
		row->hrss_v = env.ligo_lw.table.elt[column_pos.hrss_v].data.real_4;
	}
	if(miostatus < 0) {
		XLALDestroySimRingdownTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroySimRingdownTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}


/**
 * Write a sim_ringdown table to an XML file.
 */

int XLALWriteLIGOLwXMLSimRingdownTable(
	LIGOLwXMLStream *xml,
	const SimRingdownTable *sim_ringdown
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sim_ringdown:table\">\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"ilwd:char\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"waveform\" Type=\"lstring\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"coordinates\" Type=\"lstring\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"geocent_start_time\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"geocent_start_time_ns\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"h_start_time\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"h_start_time_ns\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"l_start_time\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"l_start_time_ns\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"v_start_time\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"v_start_time_ns\" Type=\"int_4s\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"start_time_gmst\" Type=\"real_8\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"longitude\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"latitude\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"distance\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"inclination\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"polarization\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"frequency\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"quality\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"phase\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mass\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"epsilon\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"amplitude\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_h\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"hrss\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"hrss_h\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"hrss_l\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"hrss_v\" Type=\"real_4\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Column Name=\"simulation_id\" Type=\"ilwd:char\"/>\n" , xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"sim_ringdown:table\" Type=\"Local\" Delimiter=\",\">\n", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; sim_ringdown; sim_ringdown = sim_ringdown->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%ld",
			row_head,
			0,	/* process_id */
			sim_ringdown->waveform,
			sim_ringdown->coordinates,
			sim_ringdown->geocent_start_time.gpsSeconds,
			sim_ringdown->geocent_start_time.gpsNanoSeconds,
			sim_ringdown->h_start_time.gpsSeconds,
			sim_ringdown->h_start_time.gpsNanoSeconds,
			sim_ringdown->l_start_time.gpsSeconds,
			sim_ringdown->l_start_time.gpsNanoSeconds,
			sim_ringdown->v_start_time.gpsSeconds,
			sim_ringdown->v_start_time.gpsNanoSeconds,
			sim_ringdown->start_time_gmst,
			sim_ringdown->longitude,
			sim_ringdown->latitude,
			sim_ringdown->distance,
			sim_ringdown->inclination,
			sim_ringdown->polarization,
			sim_ringdown->frequency,
			sim_ringdown->quality,
			sim_ringdown->phase,
			sim_ringdown->mass,
			sim_ringdown->spin,
			sim_ringdown->epsilon,
			sim_ringdown->amplitude,
			sim_ringdown->eff_dist_h,
			sim_ringdown->eff_dist_l,
			sim_ringdown->eff_dist_v,
			sim_ringdown->hrss,
			sim_ringdown->hrss_h,
			sim_ringdown->hrss_l,
			sim_ringdown->hrss_v,
			sim_ringdown->simulation_id
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
