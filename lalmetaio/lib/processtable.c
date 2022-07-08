/*
 * Copyright (C) 2007-2012 Duncan Brown, Jolien Creighton, Kipp Cannon,
 * Reinhard Prix, Bernd Machenschalk
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


#define _GNU_SOURCE


#include <config.h>

#include <ctype.h>
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>


#include <metaio.h>


#include <lal/Date.h>
#include <lal/LALMalloc.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/XLALError.h>


/**
 * Create a ProcessTable structure.
 */
ProcessTable *XLALCreateProcessTableRow(void)
{
	ProcessTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	memset(new->program, 0, sizeof(new->program));
	memset(new->version, 0, sizeof(new->version));
	memset(new->cvs_repository, 0, sizeof(new->cvs_repository));
	XLALGPSSet(&new->cvs_entry_time, 0, 0);
	memset(new->comment, 0, sizeof(new->comment));
	new->is_online = 0;
	memset(new->node, 0, sizeof(new->node));
	memset(new->username, 0, sizeof(new->username));
	XLALGPSSet(&new->start_time, 0, 0);
	XLALGPSSet(&new->end_time, 0, 0);
	new->jobid = 0;
	memset(new->domain, 0, sizeof(new->domain));
	new->unix_procid = 0;
	memset(new->ifos, 0, sizeof(new->ifos));
	new->process_id = -1;	/* impossible */

	return new;
}


/**
 * Destroy a ProcessTable structure.
 */
ProcessTable *XLALDestroyProcessTableRow(ProcessTable *row)
{
	ProcessTable *next = row ? row->next : NULL;
	XLALFree(row);
	return next;
}


/**
 * Destroy a ProcessTable linked list.
 */
void XLALDestroyProcessTable(ProcessTable *head)
{
	while(head)
		head = XLALDestroyProcessTableRow(head);
}


/**
 * Return the next available process ID.
 */
long XLALProcessTableGetNextID(ProcessTable *head)
{
	long highest = -1;
	for(; head; head = head->next)
		if(head->process_id > highest)
			highest = head->process_id;
	return highest + 1;
}


/**
 * Count and return the number of rows in the linked list.
 */
int XLALCountProcessTable(ProcessTable *head)
{
	int length;

	for(length = 0; head; head = head->next)
		length++;

	return length;
}


#ifndef _WIN32
/**
 * Parses a cvs keyword string in the form $keyword:value$, returning a
 * newly-malloc()'ed string containing just the value.  Leading and
 * trailing whitespace is removed from the value string.  Returns NULL on
 * parse or malloc() failure.
 */
static char *cvs_get_keyword_value(const char *cvs_string)
{
	const char *value_start;
	const char *value_end;
	char *value;

	/*
	 * check for input
	 */

	if(!cvs_string)
		return NULL;

	/*
	 * string must start with '$'
	 */

	value_start = strchr(cvs_string, '$');
	if(!value_start)
		return NULL;

	/*
	 * keyword ends at ':'
	 */

	value_end = strchr(value_start, ':');
	if(!value_end) {
		/*
		 * allow for unset keyword
		 */
		value_start = strchr(value_start, '$');
		if(!value_start)
			return NULL;
		else
			--value_start;
	} else {
		value_start = value_end;
	}

	/*
	 * skip leading white space
	 */

	while(isspace(*++value_start));
	if(!*value_start)
		return NULL;

	/*
	 * string must end with '$'
	 */

	value_end = strchr(value_start, '$');
	if(!value_end)
		return NULL;

	/*
	 * skip trailing white space
	 */

	while(isspace(*--value_end));
	value_end++;
	if(value_end - value_start < 0)
		value_end = value_start;

	/*
	 * extract value.  +1 for the '\0' to be added
	 */

	value = malloc((value_end - value_start + 1) * sizeof(*value));
	if(!value)
		return NULL;
	memcpy(value, value_start, (value_end - value_start) * sizeof(*value));
	value[value_end - value_start] = '\0';

	/*
	 * done
	 */

	return value;
}


/**
 * Populate a pre-allocated ProcessTable structure.
 */
int XLALPopulateProcessTable(
	ProcessTable *ptable,
	const char *program_name,
	const char *cvs_revision,
	const char *cvs_source,
	const char *cvs_date,
	long process_id
)
{
	char *cvs_keyword_value;
	uid_t uid;
	struct passwd *pw;
	struct tm utc;

	/*
	 * program name entry
	 */

	snprintf(ptable->program, LIGOMETA_PROGRAM_MAX, "%s", program_name);

	/*
	 * cvs version
	 */

	cvs_keyword_value = cvs_get_keyword_value(cvs_revision);
	if(!cvs_keyword_value) {
		XLALPrintError("%s(): cannot parse \"%s\"\n", __func__, cvs_revision);
		XLAL_ERROR(XLAL_EINVAL);
	}
	snprintf(ptable->version, LIGOMETA_VERSION_MAX, "%s", cvs_keyword_value);
	free(cvs_keyword_value);

	/*
	 * cvs repository
	 */

	cvs_keyword_value = cvs_get_keyword_value(cvs_source);
	if(!cvs_keyword_value) {
		XLALPrintError("%s(): cannot parse \"%s\"\n", __func__, cvs_source);
		XLAL_ERROR(XLAL_EINVAL);
	}
	snprintf(ptable->cvs_repository, LIGOMETA_CVS_REPOSITORY_MAX, "%s", cvs_keyword_value);
	free(cvs_keyword_value);

	/*
	 * cvs check-in time
	 */

	cvs_keyword_value = cvs_get_keyword_value(cvs_date);
	if(!cvs_keyword_value) {
		XLALPrintError("%s(): cannot parse \"%s\"\n", __func__, cvs_date);
		XLAL_ERROR(XLAL_EINVAL);
	}
	if(!strptime(cvs_keyword_value, "%Y/%m/%d %T", &utc)) {
		if(!strptime(cvs_keyword_value, "%Y-%m-%d %T", &utc)) {
			XLALPrintError("%s(): cannot parse \"%s\"\n", __func__, cvs_keyword_value);
			free(cvs_keyword_value);
			XLAL_ERROR(XLAL_EINVAL);
		}
	}
	free(cvs_keyword_value);
	XLALClearErrno();
	XLALGPSSet(&ptable->cvs_entry_time, XLALUTCToGPS(&utc), 0);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/*
	 * comment
	 */

	snprintf(ptable->comment, LIGOMETA_COMMENT_MAX, " ");

	/*
	 * online flag and domain
	 */

	ptable->is_online = 0;
	snprintf(ptable->domain, LIGOMETA_DOMAIN_MAX, "lalapps");

	/*
	 * unix process id, username, host, and process_id
	 */

	ptable->unix_procid = getpid();
	if(!ptable->unix_procid)
		ptable->unix_procid = getppid();
	if(gethostname(ptable->node, LIGOMETA_NODE_MAX) < 0) {
		perror("could not determine host name");
		XLAL_ERROR(XLAL_ESYS);
	}
	uid = geteuid();
	if(!(pw = getpwuid(uid)))
		snprintf(ptable->username, LIGOMETA_USERNAME_MAX, "%d", uid);
	else
		snprintf(ptable->username, LIGOMETA_USERNAME_MAX, "%s", pw->pw_name);
	ptable->process_id = process_id;

	/*
	 * done
	 */

	return 0;
}

#else				/* _WIN32 */

int XLALPopulateProcessTable(
	ProcessTable * ptable,
	const char *program_name,
	const char *cvs_revision,
	const char *cvs_source,
	const char *cvs_date,
	long process_id
)
{
	fprintf(stderr, "XLALPopulateProcessTable() not implemented for WIN32\n");
	return 1;
}
#endif				/* __WIN32 */


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
 * Write a process table to an XML file.
 */
int XLALWriteLIGOLwXMLProcessTable(
	LIGOLwXMLStream *xml,
	const ProcessTable *process
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"process:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"program\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"version\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"cvs_repository\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"cvs_entry_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"comment\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"is_online\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"node\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"username\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"unix_procid\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"jobid\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"domain\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"ifos\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"process:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; process; process = process->next) {
		if(XLALFilePrintf(xml->fp, "%s\"%s\",\"%s\",\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",%d,%d,%d,%d,\"%s\",\"%s\",%ld",
			row_head,
			process->program,
			process->version,
			process->cvs_repository,
			process->cvs_entry_time.gpsSeconds,
			process->comment,
			process->is_online,
			process->node,
			process->username,
			process->unix_procid,
			process->start_time.gpsSeconds,
			process->end_time.gpsSeconds,
			process->jobid,
			process->domain,
			process->ifos,
			process->process_id
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
