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
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LIGOLwXMLHV">
Author: Brown, D. A.
$Id$
</lalVerbatim>
<lalLaTeX>
\section{Header \texttt{LIGOLwXML.h}}
\label{s:LIGOLwXML.h}

Provides functions for writing the LIGO metdata database table structres to
LIGO light weight XML files.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LIGOLwXML.h>
\end{verbatim}

</lalLaTeX>
#endif

#ifndef _LIGOLWXML_H
#define _LIGOLWXML_H

#include <stdio.h>
#include <lal/FileIO.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>

NRCSID( LIGOLWXMLH, "$Id$" );

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

#if 0
<lalLaTeX>
\subsection*{Error conditions}
</lalLaTeX>
#endif
/* <lalErrTable> */
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
/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Structures}

\subsubsection*{Type \texttt{tagLIGOLwXMLStream}}
\idx[Type]{tagLIGOLwXMLStream}

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagLIGOLwXMLStream
{
  LALFILE              *fp;
  INT4                  first;
  UINT8                 rowCount;
  MetadataTableType     table;
}
LIGOLwXMLStream;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
This structure contains the file stream and current table type for
writing to LIGO lightweight XML files. It should not be manipulated
directly, but passed to the \verb+LIGOLwXML+ functions for their use.
\begin{description}
\item[\texttt{fp}] The file stream pointer of the XML file.
\item[\texttt{first}] Is this the first entry in the table.
\item[\texttt{rowCount}] Counter for the number of rows in the current table.
\item[\texttt{table}] The database table currently open.
\end{description}
</lalLaTeX>
#endif

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

int XLALWriteLIGOLwXMLSimBurstTable(
	LIGOLwXMLStream *,
	const SimBurst *
);

int XLALCreateLIGOLwXMLFileName(
        char* filename,
        const char* dataSource,
        const char* dataDescription,
        const LIGOTimeGPS* startTime,
        const LIGOTimeGPS* endTime,
        INT4 compress
);


#if 0
<lalLaTeX>
\vfill{\footnotesize\input{LIGOLwXMLHV}}
\newpage\input{LIGOLwXMLC}
</lalLaTeX> */
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LIGOLIGOLWXML_H */
