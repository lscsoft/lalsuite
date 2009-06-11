/*
*  Copyright (C) 2007 Jolien Creighton, Peter Shawhan
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

/* <lalVerbatim file="SegmentsIOHV">
Author: Peter Shawhan
Revision: $Id$
</lalVerbatim> */

/* <lalLaTeX>
\section{Header \texttt{SegmentsIO.h}}
\label{s:SegmentsIO.h}

Provides segment list reading and writing functions as part of the \texttt{lalsupport} library.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Segments.h>
#include <lal/SegmentsIO.h>
\end{verbatim}

\subsection*{Notes}
The baseline format of a segment list file is described at
\newline
\texttt{http://www.lsc-group.phys.uwm.edu/daswg/docs/technical/seglist\_format.html} .

\subsection*{Error codes}
</lalLaTeX> */

/* <lalErrTable> */
#define SEGMENTSIOH_ENULL 1
#define SEGMENTSIOH_MSGENULL "Null pointer passed to function"
#define SEGMENTSIOH_EINVAL 2
#define SEGMENTSIOH_MSGEINVAL "LALSegList structure was not properly initialized"
#define SEGMENTSIOH_EBADOPT 3
#define SEGMENTSIOH_MSGEBADOPT "Invalid option letter in options string"
#define SEGMENTSIOH_ENOFMT 4
#define SEGMENTSIOH_MSGENOFMT "No output format specified in options string"
#define SEGMENTSIOH_EOPENR 5
#define SEGMENTSIOH_MSGEOPENR "Error opening segment list file for reading"
#define SEGMENTSIOH_EOPENW 6
#define SEGMENTSIOH_MSGEOPENW "Error opening segment list file for writing"
#define SEGMENTSIOH_EFMT 7
#define SEGMENTSIOH_MSGEFMT "Segment list file is not in a recognized format"
#define SEGMENTSIOH_EPARSE 8
#define SEGMENTSIOH_MSGEPARSE "Parsing error while reading from file"
#define SEGMENTSIOH_EDOM 9
#define SEGMENTSIOH_MSGEDOM "GPS times do not represent a valid segment"
#define SEGMENTSIOH_EINT 10
#define SEGMENTSIOH_MSGEINT "Internal error in SegmentsIO module"
/* </lalErrTable> */

/* <lalLaTeX>
\vfill{\footnotesize\input{SegmentsIOHV}}
\newpage\input{SegmentsIOC}
\newpage\input{SegmentsIOTestC}
</lalLaTeX> */

#ifndef _SEGMENTSIO_H
#define _SEGMENTSIO_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <lal/FileIO.h>
#include <lal/Segments.h>

NRCSID( SEGMENTSIOH, "$Id$" );

/* Function prototypes */

void
LALSegListRead( LALStatus *status, LALSegList *seglist, const CHAR *fileName, const CHAR *options );

void
LALSegListWrite( LALStatus *status, LALSegList *seglist, const CHAR *fileName, const CHAR *options );

#ifdef __cplusplus
}
#endif

#endif /* _SEGMENTSIO_H */
