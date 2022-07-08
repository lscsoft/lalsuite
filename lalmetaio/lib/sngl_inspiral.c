/*
 * sngl_inspiral.c
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
 * Create a SnglInspiralTable structure.
 */
SnglInspiralTable *XLALCreateSnglInspiralTableRow(const ProcessTable *process)
{
	SnglInspiralTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	memset(new->ifo, 0, sizeof(new->ifo));
	memset(new->search, 0, sizeof(new->search));
	memset(new->channel, 0, sizeof(new->channel));
	if(process)
		new->process_id = process->process_id;
	else
		new->process_id = -1;	/* impossible */
	new->event_id = -1;	/* impossible */

	return new;
}


/**
 * Destroy a SnglInspiralTable structure.
 */
SnglInspiralTable *XLALDestroySnglInspiralTableRow(SnglInspiralTable *row)
{
	SnglInspiralTable *next = row ? row->next : NULL;
	XLALFree(row);
	return next;
}


/**
 * Destroy a SnglInspiralTable linked list.
 */
void XLALDestroySnglInspiralTable(SnglInspiralTable *head)
{
	while(head)
		head = XLALDestroySnglInspiralTableRow(head);
}


/**
 * Read the sngl_inspiral table from a LIGO Light Weight XML file into a
 * linked list of SnglInspiralTable structures.
 */
SnglInspiralTable *XLALSnglInspiralTableFromLIGOLw (
	const char *filename
)
{
	static const char table_name[] = "sngl_inspiral";
	int miostatus;
	SnglInspiralTable *head = NULL;
	SnglInspiralTable **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int ifo;
		int search;
		int channel;
		int end_time;
		int end_time_ns;
		int end_time_gmst;
		int impulse_time;
		int impulse_time_ns;
		int template_duration;
		int event_duration;
		int amplitude;
		int eff_distance;
		int coa_phase;
		int mass1;
		int mass2;
		int mchirp;
		int mtotal;
		int eta;
		int tau0;
		int tau2;
		int tau3;
		int tau4;
		int tau5;
		int ttotal;
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
		int f_final;
		int snr;
		int chisq;
		int chisq_dof;
		int bank_chisq;
		int bank_chisq_dof;
		int cont_chisq;
		int cont_chisq_dof;
		int sigmasq;
		int rsqveto_duration;
		int Gamma0;
		int Gamma1;
		int Gamma2;
		int Gamma3;
		int Gamma4;
		int Gamma5;
		int Gamma6;
		int Gamma7;
		int Gamma8;
		int Gamma9;
		int kappa;
		int chi;
		int spin1x;
		int spin1y;
		int spin1z;
		int spin2x;
		int spin2y;
		int spin2z;
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
	column_pos.end_time = XLALLIGOLwFindColumn(&env, "end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.end_time_ns = XLALLIGOLwFindColumn(&env, "end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.end_time_gmst = XLALLIGOLwFindColumn(&env, "end_time_gmst", METAIO_TYPE_REAL_8, 1);
	column_pos.impulse_time = XLALLIGOLwFindColumn(&env, "impulse_time", METAIO_TYPE_INT_4S, 1);
	column_pos.impulse_time_ns = XLALLIGOLwFindColumn(&env, "impulse_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.template_duration = XLALLIGOLwFindColumn(&env, "template_duration", METAIO_TYPE_REAL_8, 1);
	column_pos.event_duration = XLALLIGOLwFindColumn(&env, "event_duration", METAIO_TYPE_REAL_8, 1);
	column_pos.amplitude = XLALLIGOLwFindColumn(&env, "amplitude", METAIO_TYPE_REAL_4, 1);
	column_pos.eff_distance = XLALLIGOLwFindColumn(&env, "eff_distance", METAIO_TYPE_REAL_4, 1);
	column_pos.coa_phase = XLALLIGOLwFindColumn(&env, "coa_phase", METAIO_TYPE_REAL_4, 1);
	column_pos.mass1 = XLALLIGOLwFindColumn(&env, "mass1", METAIO_TYPE_REAL_4, 1);
	column_pos.mass2 = XLALLIGOLwFindColumn(&env, "mass2", METAIO_TYPE_REAL_4, 1);
	column_pos.mchirp = XLALLIGOLwFindColumn(&env, "mchirp", METAIO_TYPE_REAL_4, 1);
	column_pos.mtotal = XLALLIGOLwFindColumn(&env, "mtotal", METAIO_TYPE_REAL_4, 1);
	column_pos.eta = XLALLIGOLwFindColumn(&env, "eta", METAIO_TYPE_REAL_4, 1);
	column_pos.tau0 = XLALLIGOLwFindColumn(&env, "tau0", METAIO_TYPE_REAL_4, 1);
	column_pos.tau2 = XLALLIGOLwFindColumn(&env, "tau2", METAIO_TYPE_REAL_4, 1);
	column_pos.tau3 = XLALLIGOLwFindColumn(&env, "tau3", METAIO_TYPE_REAL_4, 1);
	column_pos.tau4 = XLALLIGOLwFindColumn(&env, "tau4", METAIO_TYPE_REAL_4, 1);
	column_pos.tau5 = XLALLIGOLwFindColumn(&env, "tau5", METAIO_TYPE_REAL_4, 1);
	column_pos.ttotal = XLALLIGOLwFindColumn(&env, "ttotal", METAIO_TYPE_REAL_4, 1);
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
	column_pos.f_final = XLALLIGOLwFindColumn(&env, "f_final", METAIO_TYPE_REAL_4, 1);
	column_pos.snr = XLALLIGOLwFindColumn(&env, "snr", METAIO_TYPE_REAL_4, 1);
	column_pos.chisq = XLALLIGOLwFindColumn(&env, "chisq", METAIO_TYPE_REAL_4, 1);
	column_pos.chisq_dof = XLALLIGOLwFindColumn(&env, "chisq_dof", METAIO_TYPE_INT_4S, 1);
	column_pos.bank_chisq = XLALLIGOLwFindColumn(&env, "bank_chisq", METAIO_TYPE_REAL_4, 1);
	column_pos.bank_chisq_dof = XLALLIGOLwFindColumn(&env, "bank_chisq_dof", METAIO_TYPE_INT_4S, 1);
	column_pos.cont_chisq = XLALLIGOLwFindColumn(&env, "cont_chisq", METAIO_TYPE_REAL_4, 1);
	column_pos.cont_chisq_dof = XLALLIGOLwFindColumn(&env, "cont_chisq_dof", METAIO_TYPE_INT_4S, 1);
	column_pos.sigmasq = XLALLIGOLwFindColumn(&env, "sigmasq", METAIO_TYPE_REAL_8, 1);
	column_pos.rsqveto_duration = XLALLIGOLwFindColumn(&env, "rsqveto_duration", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma0 = XLALLIGOLwFindColumn(&env, "Gamma0", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma1 = XLALLIGOLwFindColumn(&env, "Gamma1", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma2 = XLALLIGOLwFindColumn(&env, "Gamma2", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma3 = XLALLIGOLwFindColumn(&env, "Gamma3", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma4 = XLALLIGOLwFindColumn(&env, "Gamma4", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma5 = XLALLIGOLwFindColumn(&env, "Gamma5", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma6 = XLALLIGOLwFindColumn(&env, "Gamma6", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma7 = XLALLIGOLwFindColumn(&env, "Gamma7", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma8 = XLALLIGOLwFindColumn(&env, "Gamma8", METAIO_TYPE_REAL_4, 1);
	column_pos.Gamma9 = XLALLIGOLwFindColumn(&env, "Gamma9", METAIO_TYPE_REAL_4, 1);
	column_pos.kappa = XLALLIGOLwFindColumn(&env, "kappa", METAIO_TYPE_REAL_4, 1);
	column_pos.chi = XLALLIGOLwFindColumn(&env, "chi", METAIO_TYPE_REAL_4, 1);
	column_pos.spin1x = XLALLIGOLwFindColumn(&env, "spin1x", METAIO_TYPE_REAL_4, 1);
	column_pos.spin1y = XLALLIGOLwFindColumn(&env, "spin1y", METAIO_TYPE_REAL_4, 1);
	column_pos.spin1z = XLALLIGOLwFindColumn(&env, "spin1z", METAIO_TYPE_REAL_4, 1);
	column_pos.spin2x = XLALLIGOLwFindColumn(&env, "spin2x", METAIO_TYPE_REAL_4, 1);
	column_pos.spin2y = XLALLIGOLwFindColumn(&env, "spin2y", METAIO_TYPE_REAL_4, 1);
	column_pos.spin2z =XLALLIGOLwFindColumn(&env, "spin2z", METAIO_TYPE_REAL_4, 1);
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

		SnglInspiralTable *row = XLALCreateSnglInspiralTableRow(NULL);

		if(!row) {
			XLALDestroySnglInspiralTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
		strncpy(row->ifo, env.ligo_lw.table.elt[column_pos.ifo].data.lstring.data, sizeof(row->ifo) - 1);
		strncpy(row->search, env.ligo_lw.table.elt[column_pos.search].data.lstring.data, sizeof(row->search) - 1);
		strncpy(row->channel, env.ligo_lw.table.elt[column_pos.channel].data.lstring.data, sizeof(row->channel) - 1);
		XLALGPSSet(&row->end, env.ligo_lw.table.elt[column_pos.end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.end_time_ns].data.int_4s);
		row->end_time_gmst = env.ligo_lw.table.elt[column_pos.end_time_gmst].data.real_8;
		XLALGPSSet(&row->impulse_time, env.ligo_lw.table.elt[column_pos.impulse_time].data.int_4s, env.ligo_lw.table.elt[column_pos.impulse_time_ns].data.int_4s);
		row->template_duration = env.ligo_lw.table.elt[column_pos.template_duration].data.real_8;
		row->event_duration = env.ligo_lw.table.elt[column_pos.event_duration].data.real_8;
		row->amplitude = env.ligo_lw.table.elt[column_pos.amplitude].data.real_4;
		row->eff_distance = env.ligo_lw.table.elt[column_pos.eff_distance].data.real_4;
		row->coa_phase = env.ligo_lw.table.elt[column_pos.coa_phase].data.real_4;
		row->mass1 = env.ligo_lw.table.elt[column_pos.mass1].data.real_4;
		row->mass2 = env.ligo_lw.table.elt[column_pos.mass2].data.real_4;
		row->mchirp = env.ligo_lw.table.elt[column_pos.mchirp].data.real_4;
		row->mtotal = env.ligo_lw.table.elt[column_pos.mtotal].data.real_4;
		row->eta = env.ligo_lw.table.elt[column_pos.eta].data.real_4;
		row->tau0 = env.ligo_lw.table.elt[column_pos.tau0].data.real_4;
		row->tau2 = env.ligo_lw.table.elt[column_pos.tau2].data.real_4;
		row->tau3 = env.ligo_lw.table.elt[column_pos.tau3].data.real_4;
		row->tau4 = env.ligo_lw.table.elt[column_pos.tau4].data.real_4;
		row->tau5 = env.ligo_lw.table.elt[column_pos.tau5].data.real_4;
		row->ttotal = env.ligo_lw.table.elt[column_pos.ttotal].data.real_4;
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
		row->f_final = env.ligo_lw.table.elt[column_pos.f_final].data.real_4;
		row->snr = env.ligo_lw.table.elt[column_pos.snr].data.real_4;
		row->chisq = env.ligo_lw.table.elt[column_pos.chisq].data.real_4;
		row->chisq_dof = env.ligo_lw.table.elt[column_pos.chisq_dof].data.int_4s;
		row->bank_chisq = env.ligo_lw.table.elt[column_pos.bank_chisq].data.real_4;
		row->bank_chisq_dof = env.ligo_lw.table.elt[column_pos.bank_chisq_dof].data.int_4s;
		row->cont_chisq = env.ligo_lw.table.elt[column_pos.cont_chisq].data.real_4;
		row->cont_chisq_dof = env.ligo_lw.table.elt[column_pos.cont_chisq_dof].data.int_4s;
		row->sigmasq = env.ligo_lw.table.elt[column_pos.sigmasq].data.real_8;
		row->rsqveto_duration = env.ligo_lw.table.elt[column_pos.rsqveto_duration].data.real_4;
		row->Gamma[0] = env.ligo_lw.table.elt[column_pos.Gamma0].data.real_4;
		row->Gamma[1] = env.ligo_lw.table.elt[column_pos.Gamma1].data.real_4;
		row->Gamma[2] = env.ligo_lw.table.elt[column_pos.Gamma2].data.real_4;
		row->Gamma[3] = env.ligo_lw.table.elt[column_pos.Gamma3].data.real_4;
		row->Gamma[4] = env.ligo_lw.table.elt[column_pos.Gamma4].data.real_4;
		row->Gamma[5] = env.ligo_lw.table.elt[column_pos.Gamma5].data.real_4;
		row->Gamma[6] = env.ligo_lw.table.elt[column_pos.Gamma6].data.real_4;
		row->Gamma[7] = env.ligo_lw.table.elt[column_pos.Gamma7].data.real_4;
		row->Gamma[8] = env.ligo_lw.table.elt[column_pos.Gamma8].data.real_4;
		row->Gamma[9] = env.ligo_lw.table.elt[column_pos.Gamma9].data.real_4;
		row->kappa = env.ligo_lw.table.elt[column_pos.kappa].data.real_4;
		row->chi = env.ligo_lw.table.elt[column_pos.chi].data.real_4;
		row->spin1x = env.ligo_lw.table.elt[column_pos.spin1x].data.real_4;
		row->spin1y = env.ligo_lw.table.elt[column_pos.spin1y].data.real_4;
		row->spin1z = env.ligo_lw.table.elt[column_pos.spin1z].data.real_4;
		row->spin2x = env.ligo_lw.table.elt[column_pos.spin2x].data.real_4;
		row->spin2y = env.ligo_lw.table.elt[column_pos.spin2y].data.real_4;
		row->spin2z = env.ligo_lw.table.elt[column_pos.spin2z].data.real_4;
		row->event_id = env.ligo_lw.table.elt[column_pos.event_id].data.int_8s;
	}
	if(miostatus < 0) {
		XLALDestroySnglInspiralTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroySnglInspiralTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}


/**
 * Write a linked list of SnglInspiralTable structures to a the
 * sngl_inspiral table in a LIGO Light Weight XML file.
 */
int XLALWriteLIGOLwXMLSnglInspiralTable(
	LIGOLwXMLStream *xml,
	const SnglInspiralTable *sngl_inspiral
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"sngl_inspiral:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"ifo\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"search\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"channel\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"end_time_gmst\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"impulse_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"impulse_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"template_duration\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"event_duration\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"amplitude\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eff_distance\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"coa_phase\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mass1\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mass2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mchirp\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"mtotal\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"eta\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"kappa\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"chi\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"tau0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"tau2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"tau3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"tau4\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"tau5\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"ttotal\" Type=\"real_4\"/>\n", xml->fp);
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
	XLALFilePuts("\t\t<Column Name=\"f_final\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"snr\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"chisq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"chisq_dof\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"bank_chisq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"bank_chisq_dof\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"cont_chisq\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"cont_chisq_dof\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"sigmasq\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"rsqveto_duration\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma0\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma1\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma2\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma3\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma4\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma5\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma6\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma7\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma8\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"Gamma9\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin1x\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin1y\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin1z\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin2x\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin2y\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"spin2z\" Type=\"real_4\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"event_id\" Type=\"int_8s\"/>\n", xml->fp);
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
