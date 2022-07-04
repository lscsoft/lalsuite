/*
 * segment.c
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
 * Create a SegmentTable structure.
 */
SegmentTable *XLALCreateSegmentTableRow(
	const ProcessTable *process
)
{
	SegmentTable *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	new->next = NULL;
	if(process)
		new->process_id = process->process_id;
	else
		new->process_id = -1;
	new->creator_db = -1;
	new->segment_id = -1;
	new->segment_def_id = -1;
	new->segment_def_cdb = -1;
	XLALGPSSet(&new->start_time, 0, 0);
	XLALGPSSet(&new->end_time, 0, 0);

	return new;
}


/*
 * Destroy a SegmentTable structure.
 */
void XLALDestroySegmentTableRow(
	SegmentTable * row
)
{
	XLALFree(row);
}


/*
 * Destroy a SegmentTable linked list.
 */
void XLALDestroySegmentTable(
	SegmentTable * head
)
{
	while(head) {
		SegmentTable *next = head->next;
		XLALDestroySegmentTableRow(head);
		head = next;
	}
}


/**
 * Write a segment table to an XML file.
 */
int XLALWriteLIGOLwXMLSegmentTable(
	LIGOLwXMLStream *xml,
	const SegmentTable *segment_table
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"segment:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"process:process_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"start_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"start_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"end_time\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"end_time_ns\" Type=\"int_4s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"segment_definer:segment_def_id\" Type=\"int_8s\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Stream Name=\"segment:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; segment_table; segment_table = segment_table->next) {
		if(XLALFilePrintf(xml->fp, "%s%ld,%ld,%d,%d,%d,%d,%ld", row_head, segment_table->process_id, segment_table->segment_id, segment_table->start_time.gpsSeconds, segment_table->start_time.gpsNanoSeconds, segment_table->end_time.gpsSeconds, segment_table->end_time.gpsNanoSeconds, segment_table->segment_def_id) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}
