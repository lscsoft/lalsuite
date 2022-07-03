/*
 * time_slide.c
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
 * Create a TimeSlide structure.
 */
TimeSlide *XLALCreateTimeSlide(
	void
)
{
	TimeSlide *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	new->process_id = -1;
	new->time_slide_id = -1;
	memset(new->instrument, 0, sizeof(new->instrument));
	new->offset = 0;

	return new;
}


/**
 * Destroy a TimeSlide structure.
 */
void XLALDestroyTimeSlide(
	TimeSlide * row
)
{
	XLALFree(row);
}


/**
 * Destroy a TimeSlide linked list.
 */
void XLALDestroyTimeSlideTable(
	TimeSlide * head
)
{
	while(head) {
		TimeSlide *next = head->next;
		XLALDestroyTimeSlide(head);
		head = next;
	}
}


/**
 * Find and return the address of the first element in the linked list of
 * TimeSlide objects whose time_slide_id and instrument name equal the
 * values given.  TimeSlide elements whose instrument pointer is NULL are
 * skipped.  Returns NULL if no matching row is found.  This version is for
 * a linked list of const pointers, and returns a const pointer.  See also
 * XLALTimeSlideGetByIDAndInstrument().
 */
const TimeSlide *XLALTimeSlideConstGetByIDAndInstrument(
	const TimeSlide *time_slide,
	long time_slide_id,
	const char *instrument
)
{
	for(; time_slide && (time_slide->time_slide_id != time_slide_id || strcmp(time_slide->instrument, instrument)); time_slide = time_slide->next);
	return time_slide;
}


/**
 * Find and return the address of the first element in the linked list of
 * TimeSlide objects whose time_slide_id and instrument name equal the
 * values given.  TimeSlide elements whose instrument pointer is NULL are
 * skipped.  Returns NULL if no matching row is found.  This version is for
 * a linked list of non-const pointers, and returns a non-const pointer.
 * See also XLALTimeSlideConstGetByIDAndInstrument().  NOTE:  neither
 * version modifies the TimeSlide rows;  the two versions are identical,
 * they are provided to allow the const'ness to be passed through the
 * function.
 */
TimeSlide *XLALTimeSlideGetByIDAndInstrument(
	TimeSlide *time_slide,
	long time_slide_id,
	const char *instrument
)
{
	for(; time_slide && (time_slide->time_slide_id != time_slide_id || strcmp(time_slide->instrument, instrument)); time_slide = time_slide->next);
	return time_slide;
}


/**
 * Read the time_slide table from a LIGO Light Weight XML file into a
 * linked list of TimeSlide structures.
 */
TimeSlide *XLALTimeSlideTableFromLIGOLw(
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
 * Write a time_slide table to an XML file.
 */
int XLALWriteLIGOLwXMLTimeSlideTable(
	LIGOLwXMLStream *xml,
	const TimeSlide *time_slide
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"time_slide:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:time_slide_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:instrument\" Type=\"lstring\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide:offset\" Type=\"real_8\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"time_slide:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; time_slide; time_slide = time_slide->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,%ld,\"%s\",%.16g", row_head, time_slide->process_id, time_slide->time_slide_id, time_slide->instrument, time_slide->offset) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}
