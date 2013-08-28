/*
*  Copyright (C) 2007 Alexander Dietz, Jolien Creighton, Peter Shawhan
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

#ifndef _SEGMENTS_H
#define _SEGMENTS_H

#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup Segments_h
 * \author Peter Shawhan
 *
 * \brief Provides data types and functions for manipulating lists of ``segments'' (GPS time intervals).
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/Segments.h>
 * \endcode
 *
 * This header defines data structures for segments and lists of
 * segments, as well as prototypes for functions that manipulate them.
 *
 * A segment is a time interval with a start time and an end time.  The end time
 * must be equal to or later than the start time.  If the end time is equal to
 * the start time, then the segment represents a point in time.  If the end time
 * is later than the start time, then the segment represents a half-open time
 * interval, inclusive of its starting point and exclusive of its ending point.
 *
 * All of the segment list manipulation functions are XLAL functions.
 * They handle error conditions by invoking the current XLAL error handler
 * and setting \c xlalErrno to a nonzero value.
 *
 * \heading{Error conditions}
 *
 * <table><tr><th>xlalErrno</th><th>description</th></tr>
 * <tr><td>   XLAL_EFAULT</td><td>Null pointer passed for some argument</td></tr>
 * <tr><td>   XLAL_EINVAL</td><td>Attempted to use an uninitialized segment list structure</td></tr>
 * <tr><td>   XLAL_EDOM</td><td>Pair of GPS times does not represent a valid segment</td></tr>
 * </table>
 *
 * \heading{Notes}
 *
 * A \c LALSegList must be initialized before it is used.  Initialization
 * leaves it in an ``empty'' state, containing no segments. They also must be ''cleared''
 * after using \c XLALSegListClear(), and freed with \c LALFree() if it was dynamically allocated.
 * Segments can then be added to the list through an ``append'' operation.  The information about
 * each segment appended is copied to a memory location managed by the
 * \c LALSegList object.  In fact, the segments are stored in the form of
 * an array of \c LALSeg structures, with the \c segs field of
 * the segment list structure being the base address of the array.
 * This allows the segments to be accessed directly using a pointer as an
 * iterator, as in the following example code:
 *
 * \code
 * LALSegList mylist;
 * LALSeg *segp;
 * ...
 * /\* (Append segments to the segment list 'mylist' here) *\/
 * ...
 * for ( segp=mylist.segs; segp<mylist.segs+mylist.length; segp++ ) {
 *
 * printf( "The end time of the segment is GPS %d.%09d\n", segp->end.gpsSeconds, segp->end.gpsNanoSeconds );
 *
 * }
 * \endcode
 *
 * ... or by using an integer array index, as in the following example code:
 *
 * \code
 * LALSegList mylist;
 * LALSeg *segp;
 * INT4 iseg;
 * LIGOTimeGPS startgps;
 * ...
 * /\* (Append segments to the segment list 'mylist' here) *\/
 * ...
 * for ( iseg=0; iseg<mylist.length; iseg++ ) {
 *
 * /\* One way to access the segment... *\/
 * startgps = mylist.segs[iseg].start;
 * printf( "The start time of the segment is GPS %d.%09d\n", startgps.gpsSeconds, startgps.gpsNanoSeconds );
 *
 * /\* Another way to access the segment... *\/
 * segp = mylist.segs + iseg;
 * printf( "The end time of the segment is GPS %d.%09d\n", segp->end.gpsSeconds, segp->end.gpsNanoSeconds );
 *
 * }
 * \endcode
 *
 * Note that if the segment list is empty, then the \c segs field will
 * be NULL and the \c length field will be \f$0\f$.  So be careful not to
 * dereference the \c segs pointer unless you know that the length is
 * nonzero.
 *
 * A segment list is considered ``sorted'' if the segments are in ascending
 * (or at least non-descending) order according to the comparison done by
 * the \c XLALSegCmp() function.  A segment list is considered ``disjoint''
 * if no two segments in the list overlap, although they
 * may touch at an endpoint due to the half-open nature of the time intervals
 * represented by segments.  The \c LALSegList structure includes fields
 * which record whether the segment list is sorted and/or disjoint, and these
 * are used to search the segment list more efficiently when possible.  Note
 * that a segment list could in principle be disjoint but not sorted, but that
 * case is not of interest for the code; the \c disjoint field in the
 * structure specifically means that the list is sorted \e and disjoint.
 *
 * Also all segments in a segment list can be time-shifted using \c XLALSegListShift().
 *
 */
/*@{*/

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

/** Struct holding a single segment */
typedef struct
tagLALSeg
{
  LIGOTimeGPS start; /**< Beginning time of the segment */
  LIGOTimeGPS end;   /**< Ending time of the segment */
  INT4 id;           /**< Identifier (segment ID, array index, etc.) for user */
}
LALSeg;

/** Struct holding a segment list */
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

/*----------------------- Function prototypes ----------------------*/
INT4
XLALSegSet( LALSeg *seg, const LIGOTimeGPS *start, const LIGOTimeGPS *end,
	    const INT4 id );

LALSeg *
XLALSegCreate( const LIGOTimeGPS *start, const LIGOTimeGPS *end,
	       const INT4 id );

int
XLALGPSInSeg( const void *gps, const void *seg );

int
XLALSegCmp( const void *seg0, const void *seg1 );

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

INT4
XLALSegListShift( LALSegList *seglist, const LIGOTimeGPS *shift );

INT4
XLALSegListKeep(  LALSegList *seglist, const LIGOTimeGPS *start, const LIGOTimeGPS *end );

LALSeg *
XLALSegListGet( LALSegList *seglist, UINT4 indx );


int XLALSegListIsInitialized ( const LALSegList *seglist );
int XLALSegListInitSimpleSegments ( LALSegList *seglist, LIGOTimeGPS startTime, UINT4 Nseg, REAL8 Tseg );
char *XLALSegList2String ( const LALSegList *seglist );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */
