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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
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
 * \ingroup lalmetaio
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

#include <stdio.h>
#include <lal/FileIO.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**\name Error Codes */ /*@{*/
#define LIGOLWXMLH_ENULL 1
#define LIGOLWXMLH_ENNUL 2
#define LIGOLWXMLH_EALOC 3
#define LIGOLWXMLH_EUTAB 4
#define LIGOLWXMLH_EOPEN 5
#define LIGOLWXMLH_ECLOS 6
#define LIGOLWXMLH_EBGNT 7
#define LIGOLWXMLH_ENTAB 8
#define LIGOLWXMLH_EENDT 8
#define LIGOLWXMLH_ETMSM 9
#define LIGOLWXMLH_ETNOP 10
#define LIGOLWXMLH_MSGENULL "Null pointer"
#define LIGOLWXMLH_MSGENNUL "Non-null pointer"
#define LIGOLWXMLH_MSGEALOC "Memory allocation error"
#define LIGOLWXMLH_MSGEUTAB "Unknown metadata table type"
#define LIGOLWXMLH_MSGEOPEN "Error opening XML file"
#define LIGOLWXMLH_MSGECLOS "Closing an XML file with an open table"
#define LIGOLWXMLH_MSGEBGNT "Begining a table without ending previous table"
#define LIGOLWXMLH_MSGENTAB "No table type specified"
#define LIGOLWXMLH_MSGEENDT "Ending a table without an beginning a table"
#define LIGOLWXMLH_MSGETMSM "Table type mismatch"
#define LIGOLWXMLH_MSGETNOP "Table not begun for writing"
/*@}*/

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

void
LALOpenLIGOLwXMLFile (
    LALStatus           *status,
    LIGOLwXMLStream     *xml,
    const CHAR          *path
    );

void
LALCloseLIGOLwXMLFile (
    LALStatus           *status,
    LIGOLwXMLStream     *xml
    );

void
LALBeginLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml,
    MetadataTableType    table
    );

void
LALEndLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml
    );

void
LALWriteLIGOLwXMLTable (
    LALStatus           *status,
    LIGOLwXMLStream     *xml,
    MetadataTable        tablePtr,
    MetadataTableType    table
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

int XLALWriteLIGOLwXMLTimeSlideSegmentMapTable(
	LIGOLwXMLStream *xml,
	const TimeSlideSegmentMapTable *time_slide_seg_map
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
