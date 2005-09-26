/*----------------------------------- <lalVerbatim file="SegmentsHV">
 * Author: Peter Shawhan
 * Revision: $Id$
 *----------------------------------- </lalVerbatim> */

#if 0
------------------------------------- <lalLaTeX>
\section{Header \texttt{Segments.h}}
\label{s:Segments.h}

Provides data types and functions for manipulating lists of ``segments''
(GPS time intervals).

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Segments.h>
\end{verbatim}

\noindent This header defines data structures for segments and lists of
segments, as well as prototypes for functions that manipulate them.

A segment is a time interval with a start time and an end time.  The end time
must be equal to or later than the start time.  If the end time is equal to
the start time, then the segment represents a point in time.  If the end time
is later than the start time, then the segment represents a half-open time
interval, inclusive of its starting point and exclusive of its ending point.

All of the segment list manipulation functions are XLAL functions.
They handle error conditions by invoking the current XLAL error handler
and setting \texttt{xlalErrno} to a nonzero value.

\subsection*{Error conditions}

\begin{center}
\begin{tabular}{|cp{5.0in}|} \hline
\texttt{xlalErrno} & description \\ \hline
  \tt XLAL\_EFAULT & Null pointer passed for some argument \\
  \tt XLAL\_EINVAL & Attempted to use an uninitialized segment list structure \\
  \tt XLAL\_EDOM   & Pair of GPS times does not represent a valid segment \\
\hline
\end{tabular}
\end{center}

\subsection*{Structures}
\idx[Type]{LALSeg}
\idx[Type]{LALSegList}
\input{SegmentsHS}

\subsection*{Notes}

A \texttt{LALSegList} must be initialized before it is used.  Initialization
leaves it in an ``empty'' state, containing no segments.  Segments can then
be added to the list through an ``append'' operation.  The information about
each segment appended is copied to a memory location managed by the
\texttt{LALSegList} object.  In fact, the segments are stored in the form of
an array of \texttt{LALSeg} structures, with the \texttt{segs} field of
the segment list structure being the base address of the array.
This allows the segments to be accessed directly using a pointer as an
iterator, as in the following example code:
%
\begin{verbatim}
  LALSegList mylist;
  LALSeg *segp;
  ...
  /* (Append segments to the segment list 'mylist' here) */
  ...
  for ( segp=mylist.segs; segp<mylist.segs+mylist.length; segp++ ) {

    printf( "The end time of the segment is GPS %d.%09d\n",
            segp->end.gpsSeconds, segp->end.gpsNanoSeconds );

  }
\end{verbatim}
%
\ldots or by using an integer array index, as in the following example code:
%
\begin{verbatim}
  LALSegList mylist;
  LALSeg *segp;
  INT4 iseg;
  LIGOTimeGPS startgps;
  ...
  /* (Append segments to the segment list 'mylist' here) */
  ...
  for ( iseg=0; iseg<mylist.length; iseg++ ) {

    /* One way to access the segment... */
    startgps = mylist.segs[iseg].start;
    printf( "The start time of the segment is GPS %d.%09d\n",
            startgps.gpsSeconds, startgps.gpsNanoSeconds );

    /* Another way to access the segment... */
    segp = mylist.segs + iseg;
    printf( "The end time of the segment is GPS %d.%09d\n",
            segp->end.gpsSeconds, segp->end.gpsNanoSeconds );

  }
\end{verbatim}

Note that if the segment list is empty, then the \texttt{segs} field will
be NULL and the \texttt{length} field will be $0$.  So be careful not to
dereference the \texttt{segs} pointer unless you know that the length is
nonzero.

A segment list is considered ``sorted'' if the segments are in ascending
(or at least non-descending) order according to the comparison done by
the \texttt{XLALSegCmp} function.  A segment list is considered ``disjoint''
if no two segments in the list overlap, although they
may touch at an endpoint due to the half-open nature of the time intervals
represented by segments.  The \texttt{LALSegList} structure includes fields
which record whether the segment list is sorted and/or disjoint, and these
are used to search the segment list more efficiently when possible.  Note
that a segment list could in principle be disjoint but not sorted, but that
case is not of interest for the code; the \texttt{disjoint} field in the
structure specifically means that the list is sorted \emph{and} disjoint.


\vfill{\footnotesize\input{SegmentsHV}}
\newpage\input{SegmentsC}
------------------------------------- </lalLaTeX>
#endif

#ifndef _SEGMENTS_H
#define _SEGMENTS_H

#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( SEGMENTSH, "$Id$" );

/*------------------- Compile-time parameters -------------------*/

#define SEGMENTSH_ALLOCBLOCK 64  /**< Initial number of LALSeg spaces to
				  * allocate in memory at one time; this is
				  * intended to reduce the number of memory 
				  * reallocation calls needing to be made to
				  * build up a segment list.  For a large
				  * segment list, the reallocation size
				  * switches over to a multiplicative factor.
				  */

#define SEGMENTSH_INITMAGICVAL 729415386  /**< Distinctive value set in the
					   * 'initMagic' field to provide a
					   * check that the structure was
					   * properly initialized. */

/*------------------- Data structure definitions -------------------*/

/* <lalVerbatim file="SegmentsHS"> */
typedef struct
tagLALSeg
{
  LIGOTimeGPS start; /**< Beginning time of the segment */
  LIGOTimeGPS end;   /**< Ending time of the segment */
  INT4 id;           /**< Identifier (segment ID, array index, etc.) for user */
}
LALSeg;
/* </lalVerbatim> */

  
/* <lalVerbatim file="SegmentsHS"> */
typedef struct
tagLALSegList
{
  LALSeg *segs;      /**< Pointer to array of segments (LALSeg structures) */
  size_t arraySize;  /**< Size of array for which memory is allocated */
  UINT4 length;      /**< Number of segments in this segment list */
  UINT4 dplaces;     /**< Decimal places (0,3,6,9) to format GPS times */
  UINT4 sorted;      /**< Flag to indicate whether segment list is sorted */
  UINT4 disjoint;    /**< Flag to indicate whether segment list is disjoint */
  UINT4 initMagic;   /**< Internal value to help ensure list was initialized */
  LALSeg *lastFound; /**< Internal record of last segment found by a search */
}
LALSegList;
/* </lalVerbatim> */


/*----------------------- Function prototypes ----------------------*/

INT4
XLALSegSet( LALSeg *seg, const LIGOTimeGPS *start, const LIGOTimeGPS *end,
	    const INT4 id );

LALSeg *
XLALSegCreate( const LIGOTimeGPS *start, const LIGOTimeGPS *end,
	       const INT4 id );

int
XLALGPSInSeg( const LIGOTimeGPS *gps, const LALSeg *seg );

int
XLALSegCmp( const LALSeg *seg0, const LALSeg *seg1 );

INT4
XLALSegListInit( LALSegList *seglist );

INT4
XLALSegListClear( LALSegList *seglist );

INT4
XLALSegListAppend( LALSegList *seglist, const LALSeg *seg );

INT4
XLALSegListSort( LALSegList *seglist );

INT4
XLALSegListCoalesce( LALSegList *seglist );

LALSeg *
XLALSegListSearch( LALSegList *seglist, const LIGOTimeGPS *gps );

/*----------------------- Trailer stuff ----------------------------*/

#ifdef  __cplusplus
#pragma {
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
