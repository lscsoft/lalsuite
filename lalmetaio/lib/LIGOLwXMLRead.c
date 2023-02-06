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
 * \author Cannon, K. C. and Brown, D. A. and Fairhurst, S.
 * \file
 * \ingroup lalmetaio_general
 *
 * \brief Routines to read tabular data from LIGO lightweight XML files.
 *
 * ### Description ###
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
