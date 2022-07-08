/*
 * search_summary.c
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
#include <lal/LALVCSInfo.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/XLALError.h>


/**
 * Create a SearchSummaryTable structure.
 */
SearchSummaryTable *XLALCreateSearchSummaryTableRow(
	const ProcessTable * process
)
{
	SearchSummaryTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	if(process)
		new->process_id = process->process_id;
	else
		new->process_id = -1;
	memset(new->comment, 0, sizeof(new->comment));
	XLALGPSSet(&new->in_start_time, 0, 0);
	XLALGPSSet(&new->in_end_time, 0, 0);
	XLALGPSSet(&new->out_start_time, 0, 0);
	XLALGPSSet(&new->out_end_time, 0, 0);
	new->nevents = -1;
	new->nnodes = -1;
	memset(new->ifos, 0, sizeof(new->ifos));

	return new;
}


/**
 * Destroy a SearchSummaryTable structure.
 */
SearchSummaryTable *XLALDestroySearchSummaryTableRow(
	SearchSummaryTable *row
)
{
	SearchSummaryTable *next = row ? row->next : NULL;
	XLALFree(row);
	return next;
}


/**
 * Destroy a SearchSummaryTable linked list.
 */
void XLALDestroySearchSummaryTable(
	SearchSummaryTable *head
)
{
	while(head)
		head = XLALDestroySearchSummaryTableRow(head);
}


int XLALCompareSearchSummaryByOutTime(
	const void *a,
	const void *b
)
{
	const SearchSummaryTable *aPtr = *((const SearchSummaryTable * const *) a);
	const SearchSummaryTable *bPtr = *((const SearchSummaryTable * const *) b);

	INT8 ta = 0;
	INT8 tb = 0;

	/* determine the out start times */
	ta = XLALGPSToINT8NS(&(aPtr->out_start_time));
	tb = XLALGPSToINT8NS(&(bPtr->out_start_time));

	if(ta > tb)
		return 1;
	else if(ta < tb)
		return -1;
	else {
		/* determine the out end times */
		ta = XLALGPSToINT8NS(&(aPtr->out_end_time));
		tb = XLALGPSToINT8NS(&(bPtr->out_end_time));

		if(ta > tb)
			return 1;
		else if(ta < tb)
			return -1;
		else
			return 0;
	}
}


int XLALTimeSortSearchSummary(
	SearchSummaryTable **summHead,
	int (*comparfunc)(const void *, const void *)
)
{
	INT4 i;
	INT4 numSumms = 0;
	SearchSummaryTable *thisSearchSumm = NULL;
	SearchSummaryTable **summHandle = NULL;

	if(!summHead)
		XLAL_ERROR(XLAL_EIO);

	/* count the number of summs in the linked list */
	for(thisSearchSumm = *summHead; thisSearchSumm; thisSearchSumm = thisSearchSumm->next)
		++numSumms;
	if(!numSumms)
		return 0;

	/* allocate memory for an array of ptrs to sort and populate array */
	summHandle = LALCalloc(numSumms, sizeof(SearchSummaryTable *));
	for(i = 0, thisSearchSumm = *summHead; i < numSumms; ++i, thisSearchSumm = thisSearchSumm->next)
		summHandle[i] = thisSearchSumm;

	/* qsort the array using the specified function */
	qsort(summHandle, numSumms, sizeof(summHandle[0]), comparfunc);

	/* re-link the linked list in the right order */
	thisSearchSumm = *summHead = summHandle[0];
	for(i = 1; i < numSumms; ++i, thisSearchSumm = thisSearchSumm->next)
		thisSearchSumm->next = summHandle[i];
	thisSearchSumm->next = NULL;

	/* free the internal memory */
	LALFree(summHandle);

	return 0;
}


/**
 * Scan through a linked list of search_summary row objects and copy rows
 * whose ifos match the given string into a new linked list, and return the
 * address of the first entry in the new list.  Returns NULL if no matching
 * rows are found, or if an error occured (out of memory).  Use XLAL
 * error functions to check for an error status code if it matters which
 * condition occured.
 */
SearchSummaryTable *XLALIfoScanSearchSummary(
	SearchSummaryTable *input,
	CHAR *ifos
)
{
	SearchSummaryTable *output = NULL;
	SearchSummaryTable **next = &output;

	for(; input; input = input->next) {
		if(strcmp(input->ifos, ifos))
			/* IFOs don't match:  move to next row */
			continue;

		/* duplicate this row */
		*next = LALMalloc(sizeof(**next));
		if(!*next) {
			XLALDestroySearchSummaryTable(output);
			XLAL_ERROR_NULL(XLAL_ENOMEM);
		}
		/* FIXME:  this works OK now, but this is not guaranteed to
		 * duplicate the object if changes are made.  watch out for
		 * updates to the definition in the future */
		memcpy(*next, input, sizeof(**next));
		next = &(*next)->next;
		/* unlink the copy from the old list */
		*next = NULL;
	}

	return output;
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
		/*strncpy(row->shared_object, env.ligo_lw.table.elt[column_pos.shared_object].data.lstring.data, sizeof(row->shared_object) - 1); */
		/*strncpy(row->lalwrapper_cvs_tag, env.ligo_lw.table.elt[column_pos.lalwrapper_cvs_tag].data.lstring.data, sizeof(row->lalwrapper_cvs_tag) - 1); */
		/*strncpy(row->lal_cvs_tag, env.ligo_lw.table.elt[column_pos.lal_cvs_tag].data.lstring.data, sizeof(row->lal_cvs_tag) - 1); */
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
 * Write a search_summary table to an XML file.
 */


int XLALWriteLIGOLwXMLSearchSummaryTable(
	LIGOLwXMLStream *xml,
	const SearchSummaryTable *search_summary
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"search_summary:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"shared_object\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"lalwrapper_cvs_tag\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"lal_cvs_tag\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"comment\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"ifos\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"in_start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"in_start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"in_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"in_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"out_start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"out_start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"out_end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"out_end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"nevents\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"nnodes\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"search_summary:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; search_summary; search_summary = search_summary->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,\"standalone\",\"\",\"%s\",\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d", row_head, search_summary->process_id, lalVCSInfo.vcsTag, search_summary->comment, search_summary->ifos, search_summary->in_start_time.gpsSeconds, search_summary->in_start_time.gpsNanoSeconds, search_summary->in_end_time.gpsSeconds, search_summary->in_end_time.gpsNanoSeconds, search_summary->out_start_time.gpsSeconds, search_summary->out_start_time.gpsNanoSeconds, search_summary->out_end_time.gpsSeconds, search_summary->out_end_time.gpsNanoSeconds, search_summary->nevents, search_summary->nnodes) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}
