/*
 * process_params.c
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
 * Create a ProcessParamsTable structure.
 */
ProcessParamsTable *XLALCreateProcessParamsTableRow(
	const ProcessTable *process
)
{
	ProcessParamsTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	memset(new->program, 0, sizeof(new->program));
	if(process)
		new->process_id = process->process_id;
	else
		new->process_id = -1;
	memset(new->param, 0, sizeof(new->param));
	memset(new->type, 0, sizeof(new->type));
	memset(new->value, 0, sizeof(new->value));

	return new;
}


/**
 * Destroy a ProcessParamsTable structure.
 */
void XLALDestroyProcessParamsTableRow(
	ProcessParamsTable *row
)
{
	XLALFree(row);
}


/**
 * Destroy a ProcessParamsTable linked list.
 */
void XLALDestroyProcessParamsTable(
	ProcessParamsTable * head
)
{
	while(head) {
		ProcessParamsTable *next = head->next;
		XLALDestroyProcessParamsTableRow(head);
		head = next;
	}
}


/**
 * Count the number of rows in a ProcessParamsTable linked the list.
 */
int XLALCountProcessParamsTable(
	ProcessParamsTable *head
)
{
	int length;

	for(length = 0; head; head = head->next)
		length++;

	return (length);
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
 * Write a process_params table to an XML file.
 */
int XLALWriteLIGOLwXMLProcessParamsTable(
	LIGOLwXMLStream *xml,
	const ProcessParamsTable *process_params
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"process_params:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"program\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"param\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"type\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"value\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"process_params:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; process_params; process_params = process_params->next) {
		if(XLALFilePrintf(xml->fp, "%s\"%s\",%ld,\"%s\",\"%s\",\"%s\"", row_head, process_params->program, process_params->process_id, process_params->param, process_params->type, process_params->value) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}
