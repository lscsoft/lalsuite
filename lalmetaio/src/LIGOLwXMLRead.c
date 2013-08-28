/*
 * Copyright (C) 2007 Andres C. Rodriguez, Alexander Dietz, Duncan Brown,
 * Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Robert
 * Adam Mercer, Saikat Ray-Majumder, Anand Sengupta, Stephen Fairhurst,
 * Xavier Siemens, Craig Robinson , Sean Seader, Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOLwXMLRead.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A. and Fairhurst, S.
 * \file
 * \ingroup lalmetaio
 *
 * \brief Routines to write LIGO metadata database structures to LIGO lightweight XML files.
 *
 * ### Description ###
 *
 * The routine \c LALSnglInspiralTableFromLIGOLw reads in a
 * \c sngl_inspiral table from the LIGOLwXML file specified in \c fileName.
 * It returns the number of triggers read in and \c eventHead provides a
 * pointer to the head of a linked list of \c SnglInspiralTables containing the
 * events.  It will return all events between the \c startEvent and
 * \c stopEvent; if these are set to 0 and -1 respectively, all events are
 * returned.
 *
 * The routine \c InspiralTmpltBankFromLIGOLw reads in a \c sngl_inspiral
 * table from the LIGOLwXML file specified in \c fileName. It returns the
 * number of templates read in and \c bankHead provides a pointer to the head
 * of a linked list of \c InspiralTemplates containing the templates read in.
 * It will return all events between the \c startTmplt and \c stopTmplt; if
 * these are set to 0 and -1 respectively, all events are returned.  Although a
 * \c sngl_inspiral table is read in, only those entries relevant for an
 * InspiralTemplate are read in and stored.
 *
 * The routine \c SimInspiralTableFromLIGOLw reads in a \c sim_inspiral
 * table from the LIGOLwXML file specified in \c fileName.  It returns the
 * number of rows read in and \c SimHead provides a pointer to the head of a
 * linked list of \c SimInspiralTables containing the events.  Additionally, a
 * \c startTime and \c endTime are specified.  Only simulated events
 * occuring between these times are returned.  If the \c endTime is set to
 * zero, then all events are returned.
 *
 * The routine \c XLALSearchSummaryTableFromLIGOLw reads in a
 * \c search_summary table from the LIGOLwXML file specified in
 * \c fileName.  It returns a pointer to the head of a linked list of
 * \c SearchSummaryTables.
 *
 * The routine \c SummValueTableFromLIGOLw reads in a \c summ_value
 * table from the LIGOLwXML file specified in \c fileName.  It returns the
 * number of rows read in and \c sumHead provides a pointer to the head of a
 * linked list of \c SummValueTables.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * Functions in the Metaio library:
 * <ul>
 * <li> \c MetaioFindColumn
 * </li><li> \c MetaioGetRow
 * </li><li> \c MetaioOpenTable
 * </li><li> \c MetaioClose
 * </li></ul>
 *
 * ### Notes ###
 *
 * %% Any relevant notes.
 *
 */


#include <string.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>

/**
 * Test a LIGO Light Weight XML file for the presence of a specific table.
 * Returns > 0 if the document contains the table, 0 if the document does
 * not contain the table, and < 0 on error.
 *
 * BUGS:
 *
 * - This function can't tell the difference between a missing table and an
 * unparseable document.  This is a limitation in libmetaio.
 *
 * - This function parses the entire file to determine if the table is
 * present, which is slow.
 *
 * - This entire approach to XML I/O is the wrong way to go.  What's needed
 * is a "load document" function, and a "save document" function.  DO NOT
 * attempt to write such functions by using this function to test for
 * every possible table one-by-one and loading the ones that are found.
 * Put the time into writing a proper XML I/O layer!!
 */
int XLALLIGOLwHasTable(const char *filename, const char *table_name)
{
	struct MetaioParseEnvironment env;
	int has_table;

	/*
	 * open the file and find table
	 */

	if(MetaioOpenFile(&env, filename)) {
		XLALPrintError("%s(): error opening \"%s\": %s\n", __func__, filename, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR(XLAL_EIO);
	}

	/*
	 * find table.  note:  parse errors are interpreted as "table is
	 * missing".  metaio provides no other mechanism for testing for
	 * the presence of a table.
	 */

	has_table = !MetaioOpenTableOnly(&env, table_name);
	/* FIXME:  when we can rely on newer versions of libmetaio use this
	 * function instead of what follows */
	/*MetaioClearErrno(&env);*/
	env.mierrno = 0;

	/*
	 * close
	 */

	if(MetaioClose(&env)) {
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR(XLAL_EIO);
	}

	/*
	 * done
	 */

	return has_table;
}


/**
 * Convenience wrapper for MetaioFindColumn(), translating to XLAL-style
 * error reporting and printing useful error messages on failure.  Returns
 * the integer index of the column, or a negative integer if the column is
 * not found or has the wrong type.  If required is non-zero, then an XLAL
 * error is reported if the column is missing, but if required is zero then
 * no error is generated for missing columns.  When a column is found, it's
 * type is checked and an XLAL error is reported if it does not match the
 * requested type.  Passing METAIO_TYPE_UNKNOWN disables the column type
 * test.
 */
int XLALLIGOLwFindColumn(
	struct MetaioParseEnvironment *env,
	const char *name,
	enum METAIO_Type type,
	int required
)
{
	int pos = MetaioFindColumn(env, name);
	if(pos >= 0) {
		/* column was found, check type */
		if(type != METAIO_TYPE_UNKNOWN && env->ligo_lw.table.col[pos].data_type != type) {
			XLALPrintError("%s(): column \"%s\" has wrong type\n", __func__, name);
			XLAL_ERROR(XLAL_EDATA);
		}
	} else if(required) {
		/* required column is missing */
		XLALPrintError("%s(): missing required column \"%s\"\n", __func__, name);
		XLAL_ERROR(XLAL_EDATA);
	}
	return pos;
}


/**
 * Convenience function to extract the integer part of an ilwd:char ID
 * string with some error checking.  If either of ilwd_char_table_name or
 * ilwd_char_column_name is not NULL, then the corresponding portion of the
 * ilwd:char string must match it exactly.  The return value is the
 * recovered integer suffix or < 0 on failure.
 */


long long XLALLIGOLwParseIlwdChar(
	const struct MetaioParseEnvironment *env,
	int column_number,
	const char *ilwd_char_table_name,
	const char *ilwd_char_column_name
)
{
	char *fmt;
	const char *ilwd_char = env->ligo_lw.table.elt[column_number].data.lstring.data;
	long long id;

	/*
	 * 8 = 1 for the '\0', 2 for the ':' characters, and 5 for the
	 * "%%lld" string
	 */

	fmt = malloc(strlen(ilwd_char_table_name ? ilwd_char_table_name : "%*[^:]") + strlen(ilwd_char_column_name ? ilwd_char_column_name : "%*[^:]") + 8);
	if(!fmt)
		XLAL_ERROR(XLAL_ENOMEM);

	sprintf(fmt, "%s:%s:%%lld", ilwd_char_table_name ? ilwd_char_table_name : "%*[^:]", ilwd_char_column_name ? ilwd_char_column_name : "%*[^:]");

	if(sscanf(ilwd_char, fmt, &id) < 1) {
		free(fmt);
		XLALPrintError("%s(): invalid %s \"%s\" for %s\n", __func__, ilwd_char_column_name ? ilwd_char_column_name : "ID", ilwd_char, ilwd_char_table_name ? ilwd_char_table_name : "table");
		XLAL_ERROR(XLAL_EDATA);
	}

	free(fmt);

	return id;
}


/**
 * Read the process table from a LIGO Light Weight XML file into a linked
 * list of ProcessTable structures.
 */
ProcessTable *XLALProcessTableFromLIGOLw(
	const char *filename
)
{
	static const char table_name[] = "process";
	int miostatus;
	ProcessTable *head = NULL;
	ProcessTable **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int program;
		int version;
		int cvs_repository;
		int cvs_entry_time;
		int comment;
		int is_online;
		int node;
		int username;
		int unix_procid;
		int start_time;
		int end_time;
		int jobid;
		int domain;
		int ifos;
		int process_id;
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
	column_pos.program = XLALLIGOLwFindColumn(&env, "program", METAIO_TYPE_LSTRING, 1);
	column_pos.version = XLALLIGOLwFindColumn(&env, "version", METAIO_TYPE_LSTRING, 1);
	column_pos.cvs_repository = XLALLIGOLwFindColumn(&env, "cvs_repository", METAIO_TYPE_LSTRING, 1);
	column_pos.cvs_entry_time = XLALLIGOLwFindColumn(&env, "cvs_entry_time", METAIO_TYPE_INT_4S, 1);
	column_pos.comment = XLALLIGOLwFindColumn(&env, "comment", METAIO_TYPE_LSTRING, 1);
	column_pos.is_online = XLALLIGOLwFindColumn(&env, "is_online", METAIO_TYPE_INT_4S, 1);
	column_pos.node = XLALLIGOLwFindColumn(&env, "node", METAIO_TYPE_LSTRING, 1);
	column_pos.username = XLALLIGOLwFindColumn(&env, "username", METAIO_TYPE_LSTRING, 1);
	column_pos.unix_procid = XLALLIGOLwFindColumn(&env, "unix_procid", METAIO_TYPE_INT_4S, 1);
	column_pos.start_time = XLALLIGOLwFindColumn(&env, "start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.end_time = XLALLIGOLwFindColumn(&env, "end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.jobid = XLALLIGOLwFindColumn(&env, "jobid", METAIO_TYPE_INT_4S, 1);
	column_pos.domain = XLALLIGOLwFindColumn(&env, "domain", METAIO_TYPE_LSTRING, 1);
	column_pos.ifos = XLALLIGOLwFindColumn(&env, "ifos", METAIO_TYPE_LSTRING, 1);
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_ILWD_CHAR, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		ProcessTable *row = XLALCreateProcessTableRow();

		if(!row) {
			XLALDestroyProcessTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		strncpy(row->program, env.ligo_lw.table.elt[column_pos.program].data.lstring.data, sizeof(row->program) - 1);
		strncpy(row->version, env.ligo_lw.table.elt[column_pos.version].data.lstring.data, sizeof(row->version) - 1);
		strncpy(row->cvs_repository, env.ligo_lw.table.elt[column_pos.cvs_repository].data.lstring.data, sizeof(row->cvs_repository) - 1);
		XLALGPSSet(&row->cvs_entry_time, env.ligo_lw.table.elt[column_pos.cvs_entry_time].data.int_4s, 0);
		strncpy(row->comment, env.ligo_lw.table.elt[column_pos.comment].data.lstring.data, sizeof(row->comment) - 1);
		row->is_online = env.ligo_lw.table.elt[column_pos.is_online].data.int_4s;
		strncpy(row->node, env.ligo_lw.table.elt[column_pos.node].data.lstring.data, sizeof(row->node) - 1);
		strncpy(row->username, env.ligo_lw.table.elt[column_pos.username].data.lstring.data, sizeof(row->username) - 1);
		row->unix_procid = env.ligo_lw.table.elt[column_pos.unix_procid].data.int_4s;
		XLALGPSSet(&row->start_time, env.ligo_lw.table.elt[column_pos.start_time].data.int_4s, 0);
		XLALGPSSet(&row->end_time, env.ligo_lw.table.elt[column_pos.end_time].data.int_4s, 0);
		row->jobid = env.ligo_lw.table.elt[column_pos.jobid].data.int_4s;
		strncpy(row->domain, env.ligo_lw.table.elt[column_pos.domain].data.lstring.data, sizeof(row->domain) - 1);
		strncpy(row->ifos, env.ligo_lw.table.elt[column_pos.ifos].data.lstring.data, sizeof(row->ifos) - 1);
		if((row->process_id = XLALLIGOLwParseIlwdChar(&env, column_pos.process_id, "process", "process_id")) < 0) {
			XLALDestroyProcessTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
	}
	if(miostatus < 0) {
		XLALDestroyProcessTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroyProcessTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}


/**
 * Read the process_params table from a LIGO Light Weight XML file into a
 * linked list of ProcessParamsTable structures.
 */
ProcessParamsTable *XLALProcessParamsTableFromLIGOLw(
	const char *filename
)
{
	static const char table_name[] = "process_params";
	int miostatus;
	ProcessParamsTable *head = NULL;
	ProcessParamsTable **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int program;
		int process_id;
		int param;
		int type;
		int value;
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
	column_pos.program = XLALLIGOLwFindColumn(&env, "program", METAIO_TYPE_LSTRING, 1);
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_ILWD_CHAR, 1);
	column_pos.param = XLALLIGOLwFindColumn(&env, "param", METAIO_TYPE_LSTRING, 1);
	column_pos.type = XLALLIGOLwFindColumn(&env, "type", METAIO_TYPE_LSTRING, 1);
	column_pos.value = XLALLIGOLwFindColumn(&env, "value", METAIO_TYPE_LSTRING, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		ProcessParamsTable *row = XLALCreateProcessParamsTableRow(NULL);

		if(!row) {
			XLALDestroyProcessParamsTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		strncpy(row->program, env.ligo_lw.table.elt[column_pos.program].data.lstring.data, sizeof(row->program) - 1);
		if((row->process_id = XLALLIGOLwParseIlwdChar(&env, column_pos.process_id, "process", "process_id")) < 0) {
			XLALDestroyProcessParamsTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		strncpy(row->param, env.ligo_lw.table.elt[column_pos.param].data.lstring.data, sizeof(row->param) - 1);
		strncpy(row->type, env.ligo_lw.table.elt[column_pos.type].data.lstring.data, sizeof(row->type) - 1);
		strncpy(row->value, env.ligo_lw.table.elt[column_pos.value].data.lstring.data, sizeof(row->value) - 1);
	}
	if(miostatus < 0) {
		XLALDestroyProcessParamsTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroyProcessParamsTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}


/**
 * Read the time_slide table from a LIGO Light Weight XML file into a
 * linked list of TimeSlide structures.
 */

TimeSlide *
XLALTimeSlideTableFromLIGOLw (
    const char *filename
)
{
	static const char table_name[] = "time_slide";
	int miostatus;
	TimeSlide *head = NULL;
	TimeSlide **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int time_slide_id;
		int instrument;
		int offset;
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
	column_pos.time_slide_id = XLALLIGOLwFindColumn(&env, "time_slide_id", METAIO_TYPE_ILWD_CHAR, 1);
	column_pos.instrument = XLALLIGOLwFindColumn(&env, "instrument", METAIO_TYPE_LSTRING, 1);
	column_pos.offset = XLALLIGOLwFindColumn(&env, "offset", METAIO_TYPE_REAL_8, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		TimeSlide *row = XLALCreateTimeSlide();

		if(!row) {
			XLALDestroyTimeSlideTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		if((row->process_id = XLALLIGOLwParseIlwdChar(&env, column_pos.process_id, "process", "process_id")) < 0) {
			XLALDestroyTimeSlideTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		if((row->time_slide_id = XLALLIGOLwParseIlwdChar(&env, column_pos.time_slide_id, "time_slide", "time_slide_id")) < 0) {
			XLALDestroyTimeSlideTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		strncpy(row->instrument, env.ligo_lw.table.elt[column_pos.instrument].data.lstring.data, sizeof(row->instrument) - 1);
		row->offset = env.ligo_lw.table.elt[column_pos.offset].data.real_8;
	}
	if(miostatus < 0) {
		XLALDestroyTimeSlideTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroyTimeSlideTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}


/**
 * Read the search_summary table from a LIGO Light Weight XML file into a
 * linked list of SearchSummaryTable structures.
 */
SearchSummaryTable *XLALSearchSummaryTableFromLIGOLw(
	const char *filename
)
{
	static const char table_name[] = "search_summary";
	int miostatus;
	SearchSummaryTable *head = NULL;
	SearchSummaryTable **next = &head;
	struct MetaioParseEnvironment env;
	struct {
		int process_id;
		int shared_object;
		int lalwrapper_cvs_tag;
		int lal_cvs_tag;
		int comment;
		int ifos;
		int in_start_time;
		int in_start_time_ns;
		int in_end_time;
		int in_end_time_ns;
		int out_start_time;
		int out_start_time_ns;
		int out_end_time;
		int out_end_time_ns;
		int nevents;
		int nnodes;
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
	column_pos.shared_object = XLALLIGOLwFindColumn(&env, "shared_object", METAIO_TYPE_LSTRING, 1);
	column_pos.lalwrapper_cvs_tag = XLALLIGOLwFindColumn(&env, "lalwrapper_cvs_tag", METAIO_TYPE_LSTRING, 1);
	column_pos.lal_cvs_tag = XLALLIGOLwFindColumn(&env, "lal_cvs_tag", METAIO_TYPE_LSTRING, 1);
	column_pos.comment = XLALLIGOLwFindColumn(&env, "comment", METAIO_TYPE_LSTRING, 1);
	column_pos.ifos = XLALLIGOLwFindColumn(&env, "ifos", METAIO_TYPE_LSTRING, 1);
	column_pos.in_start_time = XLALLIGOLwFindColumn(&env, "in_start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.in_start_time_ns = XLALLIGOLwFindColumn(&env, "in_start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.in_end_time = XLALLIGOLwFindColumn(&env, "in_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.in_end_time_ns = XLALLIGOLwFindColumn(&env, "in_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.out_start_time = XLALLIGOLwFindColumn(&env, "out_start_time", METAIO_TYPE_INT_4S, 1);
	column_pos.out_start_time_ns = XLALLIGOLwFindColumn(&env, "out_start_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.out_end_time = XLALLIGOLwFindColumn(&env, "out_end_time", METAIO_TYPE_INT_4S, 1);
	column_pos.out_end_time_ns = XLALLIGOLwFindColumn(&env, "out_end_time_ns", METAIO_TYPE_INT_4S, 1);
	column_pos.nevents = XLALLIGOLwFindColumn(&env, "nevents", METAIO_TYPE_INT_4S, 1);
	column_pos.nnodes = XLALLIGOLwFindColumn(&env, "nnodes", METAIO_TYPE_INT_4S, 1);

	/* check for failure (== a required column is missing) */

	if(XLALGetBaseErrno()) {
		MetaioAbort(&env);
		XLALPrintError("%s(): failure reading %s table\n", __func__, table_name);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/* loop over the rows in the file */

	while((miostatus = MetaioGetRow(&env)) > 0) {
		/* create a new row */

		SearchSummaryTable *row = XLALCreateSearchSummaryTableRow(NULL);

		if(!row) {
			XLALDestroySearchSummaryTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

		/* append to linked list */

		*next = row;
		next = &(*next)->next;

		/* populate the columns */

		if((row->process_id = XLALLIGOLwParseIlwdChar(&env, column_pos.process_id, "process", "process_id")) < 0) {
			XLALDestroySearchSummaryTable(head);
			MetaioAbort(&env);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		/* FIXME:  structure definition does not include elements
		 * for these columns */
		/*strncpy(row->shared_object, env.ligo_lw.table.elt[column_pos.shared_object].data.lstring.data, sizeof(row->shared_object) - 1);*/
		/*strncpy(row->lalwrapper_cvs_tag, env.ligo_lw.table.elt[column_pos.lalwrapper_cvs_tag].data.lstring.data, sizeof(row->lalwrapper_cvs_tag) - 1);*/
		/*strncpy(row->lal_cvs_tag, env.ligo_lw.table.elt[column_pos.lal_cvs_tag].data.lstring.data, sizeof(row->lal_cvs_tag) - 1);*/
		strncpy(row->comment, env.ligo_lw.table.elt[column_pos.comment].data.lstring.data, sizeof(row->comment) - 1);
		strncpy(row->ifos, env.ligo_lw.table.elt[column_pos.ifos].data.lstring.data, sizeof(row->ifos) - 1);
		XLALGPSSet(&row->in_start_time, env.ligo_lw.table.elt[column_pos.in_start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.in_start_time_ns].data.int_4s);
		XLALGPSSet(&row->in_end_time, env.ligo_lw.table.elt[column_pos.in_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.in_end_time_ns].data.int_4s);
		XLALGPSSet(&row->out_start_time, env.ligo_lw.table.elt[column_pos.out_start_time].data.int_4s, env.ligo_lw.table.elt[column_pos.out_start_time_ns].data.int_4s);
		XLALGPSSet(&row->out_end_time, env.ligo_lw.table.elt[column_pos.out_end_time].data.int_4s, env.ligo_lw.table.elt[column_pos.out_end_time_ns].data.int_4s);
		row->nevents = env.ligo_lw.table.elt[column_pos.nevents].data.int_4s;
		row->nnodes = env.ligo_lw.table.elt[column_pos.nnodes].data.int_4s;
	}
	if(miostatus < 0) {
		XLALDestroySearchSummaryTable(head);
		MetaioAbort(&env);
		XLALPrintError("%s(): I/O error parsing %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* close file */

	if(MetaioClose(&env)) {
		XLALDestroySearchSummaryTable(head);
		XLALPrintError("%s(): error parsing document after %s table: %s\n", __func__, table_name, env.mierrmsg.data ? env.mierrmsg.data : "unknown reason");
		XLAL_ERROR_NULL(XLAL_EIO);
	}

	/* done */

	return head;
}
