/*
*  Copyright (C) 2007 Duncan Brown
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: LIGOLwXML.h
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A.
 * \file
 * \ingroup lalmetaio_general
 * \brief Provides functions for writing the LIGO metdata database table structres to LIGO light weight XML files.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/LIGOLwXML.h>
 * \endcode
 *
 */

#ifndef _LIGOLWXML_H
#define _LIGOLWXML_H

#include <stdlib.h>
#include <lal/FileIO.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LIGOMetadataTables.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * This structure contains the file stream and current table type for
 * writing to LIGO lightweight XML files. It should not be manipulated
 * directly, but passed to the \c LIGOLwXML functions for their use.
 * <dl>
 * <dt>fp</dt><dd> The file stream pointer of the XML file.</dd>
 * <dt>first</dt><dd> Is this the first entry in the table.</dd>
 * <dt>rowCount</dt><dd> Counter for the number of rows in the current table.</dd>
 * <dt>table</dt><dd> The database table currently open.</dd>
 * </dl>
 *
 */
typedef struct
tagLIGOLwXMLStream
{
  LALFILE              *fp;
  INT4                  first;
  UINT8                 rowCount;
  MetadataTableType     table;
}
LIGOLwXMLStream;


LIGOLwXMLStream *
XLALOpenLIGOLwXMLFile (
    const char *path
    );

int
XLALCloseLIGOLwXMLFile (
    LIGOLwXMLStream *xml
    );

int XLALWriteLIGOLwXMLProcessTable(
	LIGOLwXMLStream *,
	const ProcessTable *
);

int XLALWriteLIGOLwXMLProcessParamsTable(
	LIGOLwXMLStream *,
	const ProcessParamsTable *
);

int XLALWriteLIGOLwXMLSearchSummaryTable(
	LIGOLwXMLStream *,
	const SearchSummaryTable *
);

int XLALWriteLIGOLwXMLSnglBurstTable(
	LIGOLwXMLStream *,
	const SnglBurst *
);

int XLALWriteLIGOLwXMLSnglInspiralTable(
	LIGOLwXMLStream *xml,
	const SnglInspiralTable *sngl_inspiral
);

int XLALWriteLIGOLwXMLSimBurstTable(
	LIGOLwXMLStream *,
	const SimBurst *
);

int XLALWriteLIGOLwXMLTimeSlideTable(
	LIGOLwXMLStream *,
	const TimeSlide *
);

int XLALWriteLIGOLwXMLSegmentTable(
	LIGOLwXMLStream *xml,
	const SegmentTable *segment_table
);

int XLALCreateLIGODataFileName(
        char* filename,
        size_t size,
        const char* dataSource,
        const char* dataDescription,
        const LIGOTimeGPS* startTime,
        const LIGOTimeGPS* endTime,
        const char* extension
);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOLIGOLWXML_H */
