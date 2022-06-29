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
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
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
 * \ingroup lalmetaio_general
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
 * The routine \c XLALSearchSummaryTableFromLIGOLw reads in a
 * \c search_summary table from the LIGOLwXML file specified in
 * \c fileName.  It returns a pointer to the head of a linked list of
 * \c SearchSummaryTables.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <metaio.h>

#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/XLALError.h>

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
	MetaioClearErrno(&env);

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
	unsigned int type,
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
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_INT_8S, 1);

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
		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
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
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_INT_8S, 1);
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
		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
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
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_INT_8S, 1);
	column_pos.time_slide_id = XLALLIGOLwFindColumn(&env, "time_slide_id", METAIO_TYPE_INT_8S, 1);
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

		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
		row->time_slide_id = env.ligo_lw.table.elt[column_pos.time_slide_id].data.int_8s;
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
	column_pos.process_id = XLALLIGOLwFindColumn(&env, "process_id", METAIO_TYPE_INT_8S, 1);
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

		row->process_id = env.ligo_lw.table.elt[column_pos.process_id].data.int_8s;
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


/**
 * Read the sim_ringdown table from a LIGO Light Weight XML file into a
 * linked list of SimRingdownTable structures.
 */
SimRingdownTable *XLALSimRingdownTableFromLIGOLw (
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
