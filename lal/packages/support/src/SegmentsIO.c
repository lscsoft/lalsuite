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

#if 0 /* autodoc block */
<lalVerbatim file="SegmentsIOCV">
Author: Peter Shawhan
Revision: $Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{SegmentsIO.c}}

Segment list input/output routines in the \texttt{lalsupport} library.

\subsection*{Prototypes}
\idx{LALSegListRead()}
\idx{LALSegListWrite()}
\input{SegmentsIOCP}

\subsection*{Description}

\subsubsection*{Reading a segment list file}

The function \verb+LALSegListRead()+ reads a segment list from a file and
appends the segments to the specified segment list.  The segment list must
previously have been initialized.  If it already has some segments in it,
then they are retained.

If the segment read from the file includes the optional 'id' at the
beginning of the line, that is recorded in the segment list.  Otherwise,
an 'id' of $0$ is recorded.  Any additional information about the segment
(appearing on the line after the stop time) is ignored.

The function syntax includes an 'options' argument,
but no reading options are currently implemented.  In the meantime, the
user may pass either a null pointer or a CHAR pointer to a null string.

\subsubsection*{Writing a segment list file}

The function \verb+LALSegListWrite()+ writes a segment list to a file with
the specified file name, following the standard format described at \newline
\texttt{http://www.lsc-group.phys.uwm.edu/daswg/docs/technical/seglist\_format.html}
.  If the file already exists, it is overwritten.

The function syntax includes an 'options' argument which determines how
the segments are written out.  This argument must be a non-null CHAR
pointer to a string.  The string must contain one or more of the following
lowercase letters:

\begin{tabular}{|c|l|}
\hline
Character  & Effect \\ \hline
  \texttt{a} & \parbox[t]{4.5in}{\sloppy
        Causes the segment list to be written out in ASCII format.
        {\bf Currently, this is the only format supported, so \texttt{a}
        \emph{must} appear in the options string.}  In the future,
        some other format ({\it e.g.}\ XML) could be supported.} \\ \hline
  \texttt{i} & \parbox[t]{4.5in}{\sloppy
        Write the 'id' of each segment to the file, appearing before the
        GPS start and end times on the line.  If this option is not requested,
        then the GPS start and end times will be the first things on the
        line.} \\ \hline
  \texttt{d} & \parbox[t]{4.5in}{\sloppy
    Include the duration (in seconds) of the segment on each line of the file,
        appearing after the GPS start and end times on the line.} \\ \hline
\end{tabular}

\vfill{\footnotesize\input{SegmentsIOCV}}

</lalLaTeX>
#endif /* autodoc block */

#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/Segments.h>
#include <lal/SegmentsIO.h>
#include <lal/StringInput.h>

NRCSID( SEGMENTSIOC, "$Id$" );

/* <lalVerbatim file="SegmentsIOCP"> */
void
LALSegListRead( LALStatus *status, LALSegList *seglist, const CHAR *fileName, const CHAR *options )
{ /* </lalVerbatim> */
  const CHAR *optc;
  FILE *fp;
  CHAR line[4096];
  CHAR *cptr, *newcptr;
  INT4 errorCode;
  CHAR errorString[256];
  LIGOTimeGPS start, end;
  LALSeg seg;
  INT4 segid;
  INT4 xstatus;

  INITSTATUS( status, "LALSegListRead", SEGMENTSIOC );
  ATTATCHSTATUSPTR( status );

  /*-- Check arguments for sanity.  Note that the options string pointer is
    allowed to be NULL. --*/
  if ( seglist == NULL ) {
    ABORT( status, SEGMENTSIOH_ENULL, SEGMENTSIOH_MSGENULL );
  }
  if ( fileName == NULL ) {
    ABORT( status, SEGMENTSIOH_ENULL, SEGMENTSIOH_MSGENULL );
  }
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    ABORT( status, SEGMENTSIOH_EINVAL, SEGMENTSIOH_MSGEINVAL );
  }

  /*-- Parse options string, if any --*/
  if ( options ) {
    for ( optc = options; *optc != '\0'; optc++ ) {
      switch (*optc) {
      default:
	ABORT( status, SEGMENTSIOH_EBADOPT, SEGMENTSIOH_MSGEBADOPT );
      }
    }
  }

  /*-- Open the file for reading --*/
  fp = LALOpenDataFile( fileName );
  if ( fp == NULL ) {
    ABORT( status, SEGMENTSIOH_EOPENR, SEGMENTSIOH_MSGEOPENR );
  }

  /* Initially, set the errorCode and errorString for "bad format".  Once we
     have successfully parsed a segment, we will change these to "parsing
     error" */
  errorCode = SEGMENTSIOH_EFMT;
  strcpy( errorString, SEGMENTSIOH_MSGEFMT );

  /*-- Loop over lines in the file --*/
  while ( fgets(line,sizeof(line),fp) ) {

    /* Strip away any comment on the line by null-terminating at the '#' */
    cptr = strchr( line, '#' );
    if ( cptr ) { *cptr = '\0'; }

    /* Find the first non-whitespace character on the line */
    cptr = line;
    while ( isspace(*cptr) ) { cptr++; }
    /* If we consumed the whole line, go on to the next line */
    if ( *cptr == '\0' ) { continue; }

    /* Parse the first item on the line, which could be start time or 'id' */
    /* Initially assume it is the start time */
    LALStringToGPS( status->statusPtr, &start, cptr, &newcptr );
    BEGINFAIL( status ) {
      LALFclose( fp );
    } ENDFAIL( status );

    /* If no characters were consumed, then the item was not parsable as a
       GPS time (or integer) */
    if ( newcptr == cptr ) {
      LALFclose( fp );
      ABORT( status, errorCode, errorString );
    }

    /*-- Update the pointer --*/
    cptr = newcptr;

    /* If the value just read is an integer with eight or fewer digits,
       then it is the 'id' and we need to read the start time */
    if ( start.gpsSeconds < 100000000 ) {

      segid = start.gpsSeconds;

      /* Find the next non-whitespace character on the line */
      while ( isspace(*cptr) ) { cptr++; }
      /* If we consumed the whole line, then there is an error */
      if ( *cptr == '\0' ) {
	LALFclose( fp );
	ABORT( status, errorCode, errorString );
      }

      /* Parse the start time */
      LALStringToGPS( status->statusPtr, &start, cptr, &newcptr );
      BEGINFAIL( status ) {
	LALFclose( fp );
      } ENDFAIL( status );

      /* If no characters were consumed, then the item was not parsable as a
	 GPS time */
      if ( newcptr == cptr ) {
	LALFclose( fp );
	ABORT( status, errorCode, errorString );
      }

      /* Update the pointer */
      cptr = newcptr;

    } else {
      segid = 0;
    }

    /* Finally, read the end time of the segment */

    /* Find the next non-whitespace character on the line */
    while ( isspace(*cptr) ) { cptr++; }
    /* If we consumed the whole line, then there is an error */
    if ( *cptr == '\0' ) {
      LALFclose( fp );
      ABORT( status, errorCode, errorString );
    }

    /* Parse the end time */
    LALStringToGPS( status->statusPtr, &end, cptr, &newcptr );
    BEGINFAIL( status ) {
      LALFclose( fp );
    } ENDFAIL( status );

    /* If no characters were consumed, then the item was not parsable as a
       GPS time */
    if ( newcptr == cptr ) {
      LALFclose( fp );
      ABORT( status, errorCode, errorString );
    }

    /*-- If we get here, then we successfully parsed the basic segment info */
    /*-- Ignore any annotations on the segment! --*/

    /*-- Make a segment from these times --*/
    xstatus = XLALSegSet( &seg, &start, &end, segid );
    if ( xstatus ) {
      LALFclose( fp );
      if ( xlalErrno == XLAL_EDOM ) {
	ABORT( status, SEGMENTSIOH_EDOM, SEGMENTSIOH_MSGEDOM );
      } else {
	ABORT( status, SEGMENTSIOH_EINT, SEGMENTSIOH_MSGEINT );
      }
    }

    /*-- Append this segment to the segment list --*/
    xstatus = XLALSegListAppend( seglist, &seg );
    if ( xstatus ) {
      LALFclose( fp );
      ABORT( status, SEGMENTSIOH_EINT, SEGMENTSIOH_MSGEINT );
    }

    /*-- Update the errorCode and errorMessage to reflect that we have read
      one good segment, so any future error is a "parsing error" --*/
    if ( errorCode == SEGMENTSIOH_EFMT ) {
      errorCode = SEGMENTSIOH_EPARSE;
      strcpy( errorString, SEGMENTSIOH_MSGEPARSE );
    }

  }

  /*-- Close the file --*/
  LALFclose( fp );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="SegmentsIOCP"> */
void
LALSegListWrite( LALStatus *status, LALSegList *seglist, const CHAR *fileName, const CHAR *options )
{ /* </lalVerbatim> */
  const CHAR *optc;
  /* Option flags */
  INT4 ascii=0, includeID=0, includeDuration=0;
  FILE *fp;
  LALSeg *segp;
  REAL8 duration;

  INITSTATUS( status, "LALSegListWrite", SEGMENTSIOC );

  /*-- Check arguments for sanity --*/
  if ( seglist == NULL ) {
    ABORT( status, SEGMENTSIOH_ENULL, SEGMENTSIOH_MSGENULL );
  }
  if ( fileName == NULL ) {
    ABORT( status, SEGMENTSIOH_ENULL, SEGMENTSIOH_MSGENULL );
  }
  if ( options == NULL ) {
    ABORT( status, SEGMENTSIOH_ENULL, SEGMENTSIOH_MSGENULL );
  }
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    ABORT( status, SEGMENTSIOH_EINVAL, SEGMENTSIOH_MSGEINVAL );
  }

  /*-- Parse options string --*/
  for ( optc = options; *optc != '\0'; optc++ ) {
    switch (*optc) {
    case 'a': ascii=1; break;
    case 'i': includeID=1; break;
    case 'd': includeDuration=1; break;
    default:
      ABORT( status, SEGMENTSIOH_EBADOPT, SEGMENTSIOH_MSGEBADOPT );
    }
  }

  /*-- Check options for self-consistency --*/
  /* Currently, ASCII is the only valid output format */
  if ( ! ascii ) {
    ABORT( status, SEGMENTSIOH_ENOFMT, SEGMENTSIOH_MSGENOFMT );
  }

  /*-- Open the file for writing --*/
  fp = LALFopen( fileName, "w" );
  if ( fp == NULL ) {
    ABORT( status, SEGMENTSIOH_EOPENW, SEGMENTSIOH_MSGEOPENW );
  }

  /*-- Loop over segments in the list, writing out each one --*/
  for ( segp=seglist->segs; segp<seglist->segs+seglist->length; segp++ ) {

    /* Write id if specified */
    if ( includeID ) {
      fprintf( fp, "%8d ", segp->id );
    }

    /* Write GPS start and end times, with a number of decimals appropriate
       to represent the times in the segment list */
    if ( seglist->dplaces == 0 ) {
      fprintf( fp, "%10d %10d",
	       segp->start.gpsSeconds, segp->end.gpsSeconds );
    } else if ( seglist->dplaces <= 3 ) {
      fprintf( fp, "%10d.%03d %10d.%03d",
	       segp->start.gpsSeconds, segp->start.gpsNanoSeconds/1000000,
	       segp->end.gpsSeconds, segp->end.gpsNanoSeconds/1000000 );
    } else if ( seglist->dplaces <= 6 ) {
      fprintf( fp, "%10d.%06d %10d.%06d",
	       segp->start.gpsSeconds, segp->start.gpsNanoSeconds/1000,
	       segp->end.gpsSeconds, segp->end.gpsNanoSeconds/1000 );
    } else {
      fprintf( fp, "%10d.%09d %10d.%09d",
	       segp->start.gpsSeconds, segp->start.gpsNanoSeconds,
	       segp->end.gpsSeconds, segp->end.gpsNanoSeconds );
    }

    /* Include duration if specified; otherwise, just print the newline */
    if ( includeDuration ) {
      duration = (REAL8)(segp->end.gpsSeconds - segp->start.gpsSeconds) +
	1.e-9 * (REAL8)(segp->end.gpsNanoSeconds - segp->start.gpsNanoSeconds);
      if ( seglist->dplaces == 0 ) {
	fprintf( fp, " %9d\n", (INT4) duration );
      } else if ( seglist->dplaces <= 3 ) {
	fprintf( fp, " %13.3f\n", duration );
      } else if ( seglist->dplaces <= 6 ) {
	fprintf( fp, " %16.6f\n", duration );
      } else {
	fprintf( fp, " %19.9f\n", duration );
      }
    } else {
      fprintf( fp, "\n" );
    }

  }

  /*-- Close the file --*/
  LALFclose( fp );

  RETURN( status );
}
