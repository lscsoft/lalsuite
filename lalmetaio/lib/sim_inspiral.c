/*
 * sim_inspiral.c
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
 * Create a SimInspiralTable structure.
 */
SimInspiralTable *XLALCreateSimInspiralTableRow(
	const ProcessTable *process
)
{
	SimInspiralTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	memset(new->waveform, 0, sizeof(new->waveform));
	memset(new->source, 0, sizeof(new->source));
	memset(new->numrel_data, 0, sizeof(new->numrel_data));
	memset(new->taper, 0, sizeof(new->taper));
	if(process)
		new->process_id = process->process_id;
	else
		new->process_id = -1;	/* impossible */
	new->simulation_id = -1;	/* impossible */

	return new;
}


/**
 * Destroy a SimInspiralTable structure.
 */
void XLALDestroySimInspiralTableRow(
	SimInspiralTable *row
)
{
	XLALFree(row);
}


/**
 * Destroy a SimInspiralTable linked list.
 */
void XLALDestroySimInspiralTable(
	SimInspiralTable *head
)
{
	while(head) {
		SimInspiralTable *next = head->next;
		XLALDestroySimInspiralTableRow(head);
		head = next;
	}
}


/**
 * Assign simulation_id values to the entries in a SimInspiral linked list.
 * All SimInspiral rows in the list will be blamed on the given process_id,
 * and assigned sequential simulation_ids starting with the given
 * simulation_id.  The return value is the next simulation_id after the
 * last one assigned to a row in the list.
 */
long
XLALSimInspiralAssignIDs(
	SimInspiralTable *head,
	long process_id,
	long simulation_id
)
{
	for(; head; head = head->next) {
		head->process_id = process_id;
		head->simulation_id = simulation_id++;
	}
	return simulation_id;
}


/**
 * Read the sim_inspiral table from a LIGO Light Weight XML file into a
 * linked list of SimInspiralTable structures.
 */
SimInspiralTable *XLALSimInspiralTableFromLIGOLw (
	const char *filename
)
{
	static const char table_name[] = "sim_inspiral";
	int miostatus;
	SimInspiralTable *head = NULL;
	SimInspiralTable **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int waveform;
		int geocent_end_time;
		int geocent_end_time_ns;
		int h_end_time;
		int h_end_time_ns;
		int l_end_time;
		int l_end_time_ns;
		int g_end_time;
		int g_end_time_ns;
		int t_end_time;
		int t_end_time_ns;
		int v_end_time;
		int v_end_time_ns;
		int end_time_gmst;
		int source;
		int mass1;
		int mass2;
		int eta;
		int distance;
		int longitude;
		int latitude;
		int inclination;
		int coa_phase;
		int polarization;
		int psi0;
		int psi3;
		int alpha;
		int alpha1;
		int alpha2;
		int alpha3;
		int alpha4;
		int alpha5;
		int alpha6;
		int beta;
		int spin1x;
		int spin1y;
		int spin1z;
		int spin2x;
		int spin2y;
		int spin2z;
		int theta0;
		int phi0;
		int f_lower;
		int f_final;
		int mchirp;
		int eff_dist_h;
		int eff_dist_l;
		int eff_dist_g;
		int eff_dist_t;
		int eff_dist_v;
		int numrel_mode_min;
		int numrel_mode_max;
		int numrel_data;
		int amp_order;
		int taper;
		int bandpass;
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
	column_pos.geocent_end_time = XLALLIGOLwFindColumn(&env, "geocent_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.geocent_end_time_ns = XLALLIGOLwFindColumn(&env, "geocent_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.h_end_time = XLALLIGOLwFindColumn(&env, "h_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.h_end_time_ns = XLALLIGOLwFindColumn(&env, "h_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.l_end_time = XLALLIGOLwFindColumn(&env, "l_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.l_end_time_ns = XLALLIGOLwFindColumn(&env, "l_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.g_end_time = XLALLIGOLwFindColumn(&env, "g_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.g_end_time_ns = XLALLIGOLwFindColumn(&env, "g_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.t_end_time = XLALLIGOLwFindColumn(&env, "t_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.t_end_time_ns = XLALLIGOLwFindColumn(&env, "t_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.v_end_time = XLALLIGOLwFindColumn(&env, "v_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.v_end_time_ns = XLALLIGOLwFindColumn(&env, "v_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.end_time_gmst = XLALLIGOLwFindColumn(&env, "end_time_gmst", METAIO_TYPE_REAL_8, 1);
	column_pos.source = XLALLIGOLwFindColumn(&env, "source", METAIO_TYPE_LSTRING, 1);
	column_pos.mass1 = XLALLIGOLwFindColumn(&env, "mass1", METAIO_TYPE_REAL_4, 1);
	column_pos.mass2 = XLALLIGOLwFindColumn(&env, "mass2", METAIO_TYPE_REAL_4, 1);
	column_pos.eta = XLALLIGOLwFindColumn(&env, "eta", METAIO_TYPE_REAL_4, 1);
	column_pos.distance = XLALLIGOLwFindColumn(&env, "distance", METAIO_TYPE_REAL_4, 1);
	column_pos.longitude = XLALLIGOLwFindColumn(&env, "longitude", METAIO_TYPE_REAL_4, 1);
	column_pos.latitude = XLALLIGOLwFindColumn(&env, "latitude", METAIO_TYPE_REAL_4, 1);
	column_pos.inclination = XLALLIGOLwFindColumn(&env, "inclination", METAIO_TYPE_REAL_4, 1);
	column_pos.coa_phase = XLALLIGOLwFindColumn(&env, "coa_phase", METAIO_TYPE_REAL_4, 1);
	column_pos.polarization = XLALLIGOLwFindColumn(&env, "polarization", METAIO_TYPE_REAL_4, 1);
	column_pos.psi0 = XLALLIGOLwFindColumn(&env, "psi0", METAIO_TYPE_REAL_4, 1);
	column_pos.psi3 = XLALLIGOLwFindColumn(&env, "psi3", METAIO_TYPE_REAL_4, 1);
	column_pos.alpha = XLALLIGOLwFindColumn(&env, "alpha", METAIO_TYPE_REAL_4, 1);
	column_pos.alpha1 = XLALLIGOLwFindColumn(&env, "alpha1", METAIO_TYPE_REAL_4, 1);
	column_pos.alpha2 = XLALLIGOLwFindColumn(&env, "alpha2", METAIO_TYPE_REAL_4, 1);
	column_pos.alpha3 = XLALLIGOLwFindColumn(&env, "alpha3", METAIO_TYPE_REAL_4, 1);
	column_pos.alpha4 = XLALLIGOLwFindColumn(&env, "alpha4", METAIO_TYPE_REAL_4, 1);
	column_pos.alpha5 = XLALLIGOLwFindColumn(&env, "alpha5", METAIO_TYPE_REAL_4, 1);
	column_pos.alpha6 = XLALLIGOLwFindColumn(&env, "alpha6", METAIO_TYPE_REAL_4, 1);
	column_pos.beta = XLALLIGOLwFindColumn(&env, "beta", METAIO_TYPE_REAL_4, 1);
	column_pos.spin1x = XLALLIGOLwFindColumn(&env, "spin1x", METAIO_TYPE_REAL_4, 1);
	column_pos.spin1y = XLALLIGOLwFindColumn(&env, "spin1y", METAIO_TYPE_REAL_4, 1);
	column_pos.spin1z = XLALLIGOLwFindColumn(&env, "spin1z", METAIO_TYPE_REAL_4, 1);
	column_pos.spin2x = XLALLIGOLwFindColumn(&env, "spin2x", METAIO_TYPE_REAL_4, 1);
	column_pos.spin2y = XLALLIGOLwFindColumn(&env, "spin2y", METAIO_TYPE_REAL_4, 1);
	column_pos.spin2z = XLALLIGOLwFindColumn(&env, "spin2z", METAIO_TYPE_REAL_4, 1);
	column_pos.theta0 = XLALLIGOLwFindColumn(&env, "theta0", METAIO_TYPE_REAL_4, 1);
	column_pos.phi0 = XLALLIGOLwFindColumn(&env, "phi0", METAIO_TYPE_REAL_4, 1);
	column_pos.f_lower = XLALLIGOLwFindColumn(&env, "f_lower", METAIO_TYPE_REAL_4, 1);
	column_pos.f_final = XLALLIGOLwFindColumn(&env, "f_final", METAIO_TYPE_REAL_4, 1);
	column_pos.mchirp = XLALLIGOLwFindColumn(&env, "mchirp", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_h = XLALLIGOLwFindColumn(&env, "eff_dist_h", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_l = XLALLIGOLwFindColumn(&env, "eff_dist_l", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_g = XLALLIGOLwFindColumn(&env, "eff_dist_g", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_t = XLALLIGOLwFindColumn(&env, "eff_dist_t", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_dist_v = XLALLIGOLwFindColumn(&env, "eff_dist_v", METAIO_TYPE_REAL_4, 1);
	column_pos.numrel_mode_min = XLALLIGOLwFindColumn(&env, "numrel_mode_min", METAIO_TYPE_INT_4S, 1);
	column_pos.numrel_mode_max = XLALLIGOLwFindColumn(&env, "numrel_mode_max", METAIO_TYPE_INT_4S, 1);
	column_pos.numrel_data = XLALLIGOLwFindColumn(&env, "numrel_data", METAIO_TYPE_LSTRING, 1);
	column_pos.amp_order = XLALLIGOLwFindColumn(&env, "amp_order", METAIO_TYPE_INT_4S, 1);
	column_pos.taper = XLALLIGOLwFindColumn(&env, "taper", METAIO_TYPE_LSTRING, 1);
	column_pos.bandpass = XLALLIGOLwFindColumn(&env, "bandpass", METAIO_TYPE_INT_4S, 1);
	column_pos.simulation_id = XLALLIGOLwFindColumn(&env, "simulation_id", METAIO_TYPE_INT_8S, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		SimInspiralTable *row = XLALCreateSimInspiralTableRow(NULL);

		if(!row) {
			XLALDestroySimInspiralTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
		strncpy(row->waveform, env.ligo_lw.table.elt[column_pos.waveform].data.lstring.data, sizeof(row->waveform) - 1);
		strncpy(row->source, env.ligo_lw.table.elt[column_pos.source].data.lstring.data, sizeof(row->source) - 1);
		XLALGPSSet(&row->geocent_end_time, env.ligo_lw.table.elt[column_pos.geocent_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.geocent_end_time_ns].data.int_4s);
		XLALGPSSet(&row->h_end_time, env.ligo_lw.table.elt[column_pos.h_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.h_end_time_ns].data.int_4s);
		XLALGPSSet(&row->l_end_time, env.ligo_lw.table.elt[column_pos.l_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.l_end_time_ns].data.int_4s);
		XLALGPSSet(&row->g_end_time, env.ligo_lw.table.elt[column_pos.g_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.g_end_time_ns].data.int_4s);
		XLALGPSSet(&row->t_end_time, env.ligo_lw.table.elt[column_pos.t_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.t_end_time_ns].data.int_4s);
		XLALGPSSet(&row->v_end_time, env.ligo_lw.table.elt[column_pos.v_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.v_end_time_ns].data.int_4s);
		row->end_time_gmst = env.ligo_lw.table.elt[column_pos.end_time_gmst].data.real_8;
		row->mass1 = env.ligo_lw.table.elt[column_pos.mass1].data.real_4;
		row->mass2 = env.ligo_lw.table.elt[column_pos.mass2].data.real_4;
		row->eta = env.ligo_lw.table.elt[column_pos.eta].data.real_4;
		row->distance = env.ligo_lw.table.elt[column_pos.distance].data.real_4;
		row->longitude = env.ligo_lw.table.elt[column_pos.longitude].data.real_4;
		row->latitude = env.ligo_lw.table.elt[column_pos.latitude].data.real_4;
		row->inclination = env.ligo_lw.table.elt[column_pos.inclination].data.real_4;
		row->coa_phase = env.ligo_lw.table.elt[column_pos.coa_phase].data.real_4;
		row->polarization = env.ligo_lw.table.elt[column_pos.polarization].data.real_4;
		row->psi0 = env.ligo_lw.table.elt[column_pos.psi0].data.real_4;
		row->psi3 = env.ligo_lw.table.elt[column_pos.psi3].data.real_4;
		row->alpha = env.ligo_lw.table.elt[column_pos.alpha].data.real_4;
		row->alpha1 = env.ligo_lw.table.elt[column_pos.alpha1].data.real_4;
		row->alpha2 = env.ligo_lw.table.elt[column_pos.alpha2].data.real_4;
		row->alpha3 = env.ligo_lw.table.elt[column_pos.alpha3].data.real_4;
		row->alpha4 = env.ligo_lw.table.elt[column_pos.alpha4].data.real_4;
		row->alpha5 = env.ligo_lw.table.elt[column_pos.alpha5].data.real_4;
		row->alpha6 = env.ligo_lw.table.elt[column_pos.alpha6].data.real_4;
		row->beta = env.ligo_lw.table.elt[column_pos.beta].data.real_4;
		row->spin1x = env.ligo_lw.table.elt[column_pos.spin1x].data.real_4;
		row->spin1y = env.ligo_lw.table.elt[column_pos.spin1y].data.real_4;
		row->spin1z = env.ligo_lw.table.elt[column_pos.spin1z].data.real_4;
		row->spin2x = env.ligo_lw.table.elt[column_pos.spin2x].data.real_4;
		row->spin2y = env.ligo_lw.table.elt[column_pos.spin2y].data.real_4;
		row->spin2z = env.ligo_lw.table.elt[column_pos.spin2z].data.real_4;
		row->theta0 = env.ligo_lw.table.elt[column_pos.theta0].data.real_4;
		row->phi0 = env.ligo_lw.table.elt[column_pos.phi0].data.real_4;
		row->f_lower = env.ligo_lw.table.elt[column_pos.f_lower].data.real_4;
		row->f_final = env.ligo_lw.table.elt[column_pos.f_final].data.real_4;
		row->mchirp = env.ligo_lw.table.elt[column_pos.mchirp].data.real_4;
		row->eff_dist_h = env.ligo_lw.table.elt[column_pos.eff_dist_h].data.real_4;
		row->eff_dist_l = env.ligo_lw.table.elt[column_pos.eff_dist_l].data.real_4;
		row->eff_dist_g = env.ligo_lw.table.elt[column_pos.eff_dist_g].data.real_4;
		row->eff_dist_t = env.ligo_lw.table.elt[column_pos.eff_dist_t].data.real_4;
		row->eff_dist_v = env.ligo_lw.table.elt[column_pos.eff_dist_v].data.real_4;
		row->numrel_mode_min = env.ligo_lw.table.elt[column_pos.numrel_mode_min].data.int_4s;
		row->numrel_mode_max = env.ligo_lw.table.elt[column_pos.numrel_mode_max].data.int_4s;
		strncpy(row->numrel_data, env.ligo_lw.table.elt[column_pos.numrel_data].data.lstring.data, sizeof(row->numrel_data) - 1);
		row->amp_order = env.ligo_lw.table.elt[column_pos.amp_order].data.int_4s;
		strncpy(row->taper, env.ligo_lw.table.elt[column_pos.taper].data.lstring.data, sizeof(row->taper) - 1);
		row->bandpass = env.ligo_lw.table.elt[column_pos.bandpass].data.int_4s;
		row->simulation_id = env.ligo_lw.table.elt[column_pos.simulation_id].data.int_8s;
	}
	if(miostatus < 0) {
		XLALDestroySimInspiralTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroySimInspiralTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}


/**
 * Write a sim_inspiral table to an XML file.
 */

int XLALWriteLIGOLwXMLSimInspiralTable(
	LIGOLwXMLStream *xml,
	const SimInspiralTable *sim_inspiral
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sim_inspiral:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"waveform\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"geocent_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"geocent_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"h_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"h_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"l_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"l_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"g_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"g_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"t_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"t_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"v_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"v_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"source\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mass1\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mass2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mchirp\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eta\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"distance\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"longitude\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"latitude\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"inclination\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"coa_phase\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"polarization\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"psi0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"psi3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"alpha\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"alpha1\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"alpha2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"alpha3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"alpha4\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"alpha5\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"alpha6\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"beta\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin1x\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin1y\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin1z\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin2x\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin2y\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin2z\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"theta0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"phi0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"f_lower\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"f_final\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_h\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_l\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_g\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_t\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_dist_v\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"numrel_mode_min\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"numrel_mode_max\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"numrel_data\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"amp_order\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"taper\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"bandpass\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"simulation_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"sim_inspiral:table\" Type=\"Local\" Delimiter=\",\">\n", xml->fp);

	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; sim_inspiral; sim_inspiral = sim_inspiral->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.16g,\"%s\",%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%16g,%d,%d,\"%s\",%d,\"%s\",%d,%ld",
					row_head,
					sim_inspiral->process_id,
					sim_inspiral->waveform,
					sim_inspiral->geocent_end_time.gpsSeconds,
					sim_inspiral->geocent_end_time.gpsNanoSeconds,
					sim_inspiral->h_end_time.gpsSeconds,
					sim_inspiral->h_end_time.gpsNanoSeconds,
					sim_inspiral->l_end_time.gpsSeconds,
					sim_inspiral->l_end_time.gpsNanoSeconds,
					sim_inspiral->g_end_time.gpsSeconds,
					sim_inspiral->g_end_time.gpsNanoSeconds,
					sim_inspiral->t_end_time.gpsSeconds,
					sim_inspiral->t_end_time.gpsNanoSeconds,
					sim_inspiral->v_end_time.gpsSeconds,
					sim_inspiral->v_end_time.gpsNanoSeconds,
					sim_inspiral->end_time_gmst,
					sim_inspiral->source,
					sim_inspiral->mass1,
					sim_inspiral->mass2,
					sim_inspiral->mchirp,
					sim_inspiral->eta,
					sim_inspiral->distance,
					sim_inspiral->longitude,
					sim_inspiral->latitude,
					sim_inspiral->inclination,
					sim_inspiral->coa_phase,
					sim_inspiral->polarization,
					sim_inspiral->psi0,
					sim_inspiral->psi3,
					sim_inspiral->alpha,
					sim_inspiral->alpha1,
					sim_inspiral->alpha2,
					sim_inspiral->alpha3,
					sim_inspiral->alpha4,
					sim_inspiral->alpha5,
					sim_inspiral->alpha6,
					sim_inspiral->beta,
					sim_inspiral->spin1x,
					sim_inspiral->spin1y,
					sim_inspiral->spin1z,
					sim_inspiral->spin2x,
					sim_inspiral->spin2y,
					sim_inspiral->spin2z,
					sim_inspiral->theta0,
					sim_inspiral->phi0,
					sim_inspiral->f_lower,
					sim_inspiral->f_final,
					sim_inspiral->eff_dist_h,
					sim_inspiral->eff_dist_l,
					sim_inspiral->eff_dist_g,
					sim_inspiral->eff_dist_t,
					sim_inspiral->eff_dist_v,
					sim_inspiral->numrel_mode_min,
					sim_inspiral->numrel_mode_max,
					sim_inspiral->numrel_data,
					sim_inspiral->amp_order,
					sim_inspiral->taper,
					sim_inspiral->bandpass,
					sim_inspiral->simulation_id
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
