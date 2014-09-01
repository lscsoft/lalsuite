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

#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/Segments.h>
#include <lal/LALString.h>

/**
 * \addtogroup Segments_h
 *
 * ### Functions for handling segments ###
 *
 * The first few functions listed deal with \e segments:
 *
 * XLALSegSet(), XLALSegCreate(), XLALGPSInSeg(), XLALSegCmp()
 *
 * ### Functions for handling segment lists ###
 *
 * The rest of the functions listed deal with <em>segment lists</em>:
 *
 * XLALSegListInit(), XLALSegListClear(), XLALSegListAppend(), XLALSegListSort()
 * XLALSegListCoalesce(), XLALSegListSearch()
 *
 * ### Error codes and return values ###
 *
 * Each XLAL function listed above, if it fails invokes the current XLAL error
 * handler, sets \c xlalErrno to the appropriate XLAL error code, and
 * returns a particular value as noted below:
 *
 * <ul>
 * <li>
 * Functions which return an integer status code (XLALSegSet(),
 * XLALSegListInit(), XLALSegListClear(), XLALSegListAppend(), XLALSegListSort(),
 * XLALSegListCoalesce() return XLAL_SUCCESS if successful
 * or XLAL_FAILURE if an error occurs.</li>
 * <li>
 * XLALGPSInSeg() and XLALSegCmp() normally return a
 * comparison value (negative, 0, or positive).</li>
 * <li>
 * XLALSegCreate() normally returns a pointer to the created
 * segment.  If an error occurs, it returns NULL.</li>
 * <li>
 * XLALSegListSearch() returns a pointer to a segment in the list which
 * contains the time being searched for, or NULL if there is no such segment.
 * If more than one segment in the list contains the time, then this function
 * returns a pointer to \e one of the segments which contains it, not
 * necessarily the first such segment in the list.  (This is not an issue if
 * the list is ``disjoint'', which guarantees that it has no overlapping
 * segments.)  If no segment in the list contains the time, then this function
 * returns NULL; however, this is not really an error, use \c xlalErrno
 * to differentiate between failure and non-failure.
 * </li>
 * </ul>
 *
 */


/*---------------------------------------------------------------------------*/
/**
 * This function sets the start time, the end time, and the \a id of a segment.
 * The \a id can be any integer and is
 * solely for the use of the user, e.g. to store a segment ID code
 * or an index into some array containing additional information about the
 * segment.  XLALSegSet() checks to make sure the segment is valid,
 * i.e. the end time is later than or equal to the start time;
 * an error occurs if this condition is not true.
 */
int
XLALSegSet( LALSeg *seg, const LIGOTimeGPS *start, const LIGOTimeGPS *end,
            const INT4 id )
{

  /* Make sure a non-null pointer was passed for the segment to be set */
  if ( ! seg ) {
    XLALPrintError( "NULL LALSeg pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure non-null pointers were passed for the GPS times */
  if ( ! start || ! end ) {
    XLALPrintError( "NULL LIGOTimeGPS pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Check that segment end time is equal to or later than start time */
  if ( XLALGPSCmp(start,end) > 0 ) {
    XLALPrintError( "Inconsistent times passed to %s (%d.%09d > %d.%09d)\n",
        __func__, start->gpsSeconds, start->gpsNanoSeconds,
        end->gpsSeconds, end->gpsNanoSeconds );
    XLAL_ERROR( XLAL_EDOM );
  }

  /* Store the times and the id */
  seg->start = *start;
  seg->end = *end;
  seg->id = id;

  /* Return with success status code */
  return XLAL_SUCCESS;
}


/*---------------------------------------------------------------------------*/

/**
 * This function is similar to XLALSegSet()
 * except that it allocates memory for a new segment structure rather than setting
 * the fields of an existing segment structure.  It returns a pointer to the
 * new \a LALSeg structure.  When the structure is no longer needed, its
 * pointer should be passed to XLALFree().
 */
LALSeg *
XLALSegCreate( const LIGOTimeGPS *start, const LIGOTimeGPS *end,
               const INT4 id )
{
  LALSeg *segptr;

  /* Make sure non-null pointers were passed for the GPS times */
  if ( ! start || ! end ) {
    XLALPrintError( "NULL LIGOTimeGPS pointer passed to %s\n", __func__ );
    XLAL_ERROR_NULL( XLAL_EFAULT );
  }

  /* Check that segment end time is equal to or later than start time */
  if ( XLALGPSCmp(start,end) > 0 ) {
    XLALPrintError( "Inconsistent times passed to %s (%d.%09d > %d.%09d)\n",
        __func__, start->gpsSeconds, start->gpsNanoSeconds,
        end->gpsSeconds, end->gpsNanoSeconds );
    XLAL_ERROR_NULL( XLAL_EDOM );
  }

  /* Allocate memory for the new LALSeg object */
  segptr = (LALSeg *) LALMalloc( sizeof(LALSeg) );

  /* Check for a memory allocation failure */
  if ( ! segptr ) {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  /* Store the times and the id */
  segptr->start = *start;
  segptr->end = *end;
  segptr->id = id;

  /* Return a pointer to the newly-created segment */
  return segptr;
}


/*---------------------------------------------------------------------------*/

/**
 * This is designed to be usable as a comparison function for
 * <tt>bsearch()</tt> and therefore returns a negative value, 0, or a positive
 * value depending on whether the GPS time (the first argument) is before the
 * beginning of, within, or after the end of the segment (the second argument).
 * Note that a segment is a half-open interval, so the GPS time is considered
 * to be within the segment if it is equal to the start time of the segment
 * but not if it is equal to the end time of the segment.
 *
 * Returns a comparison value (negative, 0, or positive).
 */
int
XLALGPSInSeg( const void *pgps, const void *pseg )
{
  /* if time is < start of segment, return -1 */
  if ( XLALGPSCmp( (const LIGOTimeGPS *) pgps, &((const LALSeg *) pseg)->start ) < 0 )
    return -1;
  /* else if time is < end of segment, return 0 */
  if ( XLALGPSCmp( (const LIGOTimeGPS *) pgps, &((const LALSeg *) pseg)->end ) < 0 )
    return 0;
  /* time is >= end of segment, return +1 */
  return +1;
}


/*---------------------------------------------------------------------------*/

/**
 * This is designed to be usable as a comparison function for
 * <tt>qsort()</tt> and therefore returns a negative value, 0, or a positive value
 * depending on
 * whether the first argument is less than, equal to, or greater than the second.
 * The comparison is based on the start time of the segments, unless these are
 * equal, in which case the end times are compared.  Therefore, two segments
 * are considered equal only if their start \e and end times are identical.
 *
 * Returns a comparison value (negative, 0, or positive).
 */
int
XLALSegCmp( const void *pseg0, const void *pseg1 )
{
  const LALSeg *seg0 = pseg0;
  const LALSeg *seg1 = pseg1;
  int result;

  /* If segment start times are different, comparison is based on that */
  result = XLALGPSCmp( &seg0->start, &seg1->start );
  if ( !result ) {
    /* If we get here, then the segments have the same start time,
       so use the end times to do the comparison */
    result = XLALGPSCmp( &seg0->end, &seg1->end );
  }

  return result;
}


/*---------------------------------------------------------------------------*/

/**
 * This function must be called to initialize a segment
 * list structure before that structure can be used.
 */
int
XLALSegListInit( LALSegList *seglist )
{
  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "NULL LALSegList pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Initialize values for this segment list */
  seglist->segs = NULL;
  seglist->arraySize = 0;
  seglist->length = 0;

  /* So far, we do not need any decimal places to format GPS times */
  seglist->dplaces = 0;

  /* An empty list can be considered sorted and disjoint */
  seglist->sorted = 1;
  seglist->disjoint = 1;

  /* No search has been performed yet */
  seglist->lastFound = NULL;

  /* Store a distinctive integer value to serve as a check that initialization
     was performed */
  seglist->initMagic = SEGMENTSH_INITMAGICVAL;

  /* Return with success status code */
  return XLAL_SUCCESS;
}


/*---------------------------------------------------------------------------*/

/**
 * This function must be called when you are done
 * with a segment list, in order to free memory that was allocated to store
 * the segments in the list.  (Strictly speaking, this is only necessary if the
 * list contains one or more segments.)
 * The function leaves the segment list in a valid
 * state, but containing no segments.  After calling XLALSegListClear(),
 * it is OK to re-use the segment list structure (by using
 * XLALSegListAppend() to add segments to it, etc.).
 * You do not have to call XLALSegListInit() again before re-using it,
 * but it is OK to do so.  If you do call XLALSegListInit() again,
 * you must do so \e after calling XLALSegListClear() !
 */
int
XLALSegListClear( LALSegList *seglist )
{
  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "NULL LALSegList pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    XLALPrintError( "Passed unintialized LALSegList structure to %s\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* Clear the memory (if any) that has been allocated for segments */
  if ( seglist->segs ) {
    LALFree( seglist->segs );
  }

  /* Clear values related to the segment list */
  seglist->segs = NULL;
  seglist->arraySize = 0;
  seglist->length = 0;

  /* So far, we do not need any decimal places to format GPS times */
  seglist->dplaces = 0;

  /* An empty list can be considered sorted and disjoint */
  seglist->sorted = 1;
  seglist->disjoint = 1;

  /* No search has been performed yet */
  seglist->lastFound = NULL;

  /* Return with success status code */
  return XLAL_SUCCESS;
}


/*---------------------------------------------------------------------------*/

/**
 * This function appends a segment to a segment list.
 * It first checks to make sure the segment is valid,
 * i.e. the end time is later than or equal to the start time;
 * an error occurs if this condition is not true.
 * The input segment information is copied into an array of \c LALSeg
 * structures maintained ``internally'' by the segment list.  The function takes
 * care of extending this array when necessary.
 * It also checks whether the segment being appended preserves the ``sorted''
 * and/or ``disjoint'' properties of the segment list.  An empty segment list has
 * the ``sorted'' and ``disjoint'' properties to start with, and as long as
 * segments are appended in ascending time order and do not overlap, it retains
 * those properties.
 */
int
XLALSegListAppend( LALSegList *seglist, const LALSeg *seg )
{
  LALSeg *segptr;
  LALSeg *prev;
  size_t newSize;
  INT4 ns1, ns2;

  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "NULL LALSegList pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure a non-null pointer was passed for the segment */
  if ( ! seg ) {
    XLALPrintError( "NULL LALSeg pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    XLALPrintError( "Passed unintialized LALSegList structure to %s\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* Check that segment end time is equal to or later than start time */
  if ( XLALGPSCmp(&(seg->start),&(seg->end)) > 0 ) {
    XLALPrintError( "Invalid segment passed to %s (%d.%09d > %d.%09d)\n",
        __func__, seg->start.gpsSeconds, seg->start.gpsNanoSeconds,
        seg->end.gpsSeconds, seg->end.gpsNanoSeconds );
    XLAL_ERROR( XLAL_EDOM );
  }

  /* See whether we need to extend (or create) the segment array */
  if ( seglist->length == seglist->arraySize ) {

    if ( seglist->arraySize < 6*SEGMENTSH_ALLOCBLOCK ) {
      newSize = seglist->arraySize + SEGMENTSH_ALLOCBLOCK;
    } else {
      newSize = (seglist->arraySize * 6) / 5 ;
    }

    if ( seglist->arraySize ) {
      /* Extend the array */
      segptr = (LALSeg *) LALRealloc( seglist->segs, newSize*sizeof(LALSeg) );
    } else {
      /* Allocate the first block of memory for the array */
      segptr = (LALSeg *) LALMalloc( newSize*sizeof(LALSeg) );
    }

    /* Check for a memory allocation failure */
    if ( ! segptr ) {
      XLAL_ERROR( XLAL_ENOMEM );
    }

    /* If we get here, then the memory allocation succeeded */
    seglist->segs = segptr;
    seglist->arraySize = newSize;

  }

  /* Copy the segment information into the array */
  seglist->segs[seglist->length] = *seg ;
  seglist->length++;

  /* See whether more decimal places are needed to represent these times than
     were needed for segments already in the list.  Work with 0, 3, 6, or 9
     decimal places. */
  ns1 = seg->start.gpsNanoSeconds;
  ns2 = seg->end.gpsNanoSeconds;
  if ( seglist->dplaces < 9 ) {
    if ( ns1 % 1000 || ns2 % 1000 ) {
      /* 6 decimal places are not enough */
      seglist->dplaces = 9;
    } else if ( seglist->dplaces < 6 ) {
      if ( ns1 % 1000000 || ns2 % 1000000 ) {
        /* 3 decimal places are not enough */
        seglist->dplaces = 6;
      } else if ( seglist->dplaces < 3 ) {
        if ( ns1 || ns2 ) {
          /* At least one of the times does have a decimal part */
          seglist->dplaces = 3;
        }
      }
    }
  }

  /* See whether the "disjoint" and/or "sorted" properties still hold */
  if ( seglist->length > 1 ) {
    prev = seglist->segs + seglist->length - 2 ;

    if ( seglist->disjoint && XLALGPSCmp(&(prev->end),&(seg->start)) > 0 ) {
      /* Segments overlap, so segment list is no longer disjoint */
      seglist->disjoint = 0;
    }

    if ( seglist->sorted && XLALSegCmp(prev,seg) > 0 ) {
      /* Segments are no longer in ascending order */
      seglist->sorted = 0;
    }

  }

  /* Return with success status code */
  return XLAL_SUCCESS;
}


/*---------------------------------------------------------------------------*/

/**
 * This function sorts the segments in a segment list
 * into forward time order.  If the list is already sorted, then this function
 * returns promptly.
 */
int
XLALSegListSort( LALSegList *seglist )
{
  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "NULL LALSegList pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    XLALPrintError( "Passed unintialized LALSegList structure to %s\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* If segment list is known to be sorted already, just return */
  if ( seglist->sorted ) {
    return XLAL_SUCCESS;
  }

  /* Sort the segment list in place */
  /* This function call produces a compiler warning, because XLALSegCmp takes
     pointers to LALSeg objects, whereas qsort expects the comparison function
     passed to it to take void pointers as arguments.  Oh well. */
  qsort( (void *) seglist->segs, seglist->length, sizeof(LALSeg), XLALSegCmp );

  /* Now we can set the "sorted" flag */
  seglist->sorted = 1;

  /* Reset the 'lastFound' value, since the array has changed */
  seglist->lastFound = NULL;

  /* Return with success status code */
  return XLAL_SUCCESS;
}


/*---------------------------------------------------------------------------*/

/**
 * The function XLALSegListCoalesce() first sorts the segments in a
 * segment list (if not already sorted) and then joins together segments which
 * overlap or touch (i.e. share endpoints).
 * The result is a segment list which is sorted and is guaranteed
 * to not have any overlapping segments; thus it is ``disjoint''.
 * (Note, however, that a disjoint segment list is not necessarily coalesced,
 * since segments which touch at an endpoint are considered disjoint but will
 * be joined by XLALSegListCoalesce().)
 * If the list has the ``disjoint'' property to begin with,
 * then this function returns promptly.  Each segment in the output list is
 * assigned the \c id value taken from the first segment in the input list
 * (after it has been sorted) among those which were joined to make it.
 */
int
XLALSegListCoalesce( LALSegList *seglist )
{
  LALSeg *rp, *wp;   /* Read and write pointers for stepping through array */
  size_t newLength;
  LALSeg *segptr;

  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "NULL LALSegList pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    XLALPrintError( "Passed unintialized LALSegList structure to %s\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* If segment list is empty or has only one segment, just return */
  if ( seglist->length <= 1 ) {
    return XLAL_SUCCESS;
  }

  /* Make sure the segment list is sorted */
  XLALSegListSort( seglist );

  /* Step through the segment list with read and write pointers */
  wp = seglist->segs;
  newLength = 1;
  for ( rp = seglist->segs+1; rp < seglist->segs+seglist->length; rp++ ) {
    if ( XLALGPSCmp(&(wp->end),&(rp->start)) < 0 ) {
      /* Segment at read pointer does NOT overlap segment at write pointer */
      wp++;
      *wp = *rp;
      newLength++;
    } else {
      /* Segment at read pointer touches or overlaps segment at write pointer*/
      /* So extend the segment at the write pointer if necessary */
      if ( XLALGPSCmp(&(wp->end),&(rp->end)) < 0 ) {
        wp->end = rp->end;
      }
    }
  }

  /* Update the array length */
  seglist->length = (UINT4) newLength;

  /* Reduce the memory allocation for the array, if possible */
  if ( newLength < seglist->arraySize ) {
    /* Contract the array */
    segptr = (LALSeg *) LALRealloc( seglist->segs, newLength*sizeof(LALSeg) );

    /* Check for a memory allocation failure */
    if ( ! segptr ) {
      XLAL_ERROR( XLAL_ENOMEM );
    }

    /* If we get here, then the memory reallocation succeeded */
    seglist->segs = segptr;
    seglist->arraySize = newLength;
  }

  /* Now we can set the "disjoint" flag */
  seglist->disjoint = 1;

  /* Reset the 'lastFound' value, since the array has changed */
  seglist->lastFound = NULL;

  /* Return with success status code */
  return XLAL_SUCCESS;
}


/*---------------------------------------------------------------------------*/

/**
 * The function XLALSegListSearch() determines which segment in the
 * list, if any, contains the GPS time passed to this function.  It returns
 * a pointer to a segment containing the time, if there is one, otherwise
 * it returns NULL.  If more than one segment
 * in the list contains the time, then this function returns a pointer
 * to \e one of the segments which contains it, not necessarily the first
 * such segment in the list.  (This is not an issue if the list is ``disjoint'',
 * which guarantees that it has no overlapping segments.)
 * The following code shows how the XLALSegListSearch() function
 * might be used:
 * \code
 * LALSegList mylist;
 * LIGOTimeGPS tgps, startgps;
 * LALSeg *segp;
 * ...
 * /\* (Fill the segment list 'mylist' with segments here) *\/
 * /\* (Set the gps time, 'tgps', to search for) *\/
 * ...
 * segp = XLALSegListSearch( &mylist, &tgps );
 * if ( segp ) {
 *   startgps = segp->start;
 *   printf( "That time is within a segment which starts at GPS time %d.%09d\n",
 *   startgps.gpsSecconds, startgps.gpsNanoSeconds );
 * } else {
 *   printf( "That time is not within any segment in the list\n" );
 * }
 * \endcode
 * The search algorithm used by the XLALSegListSearch() function
 * depends on whether the segment list is ``sorted'' and/or ``disjoint''.  If
 * the segment list has both of these properties, then the function can use
 * a binary search to locate the segment containing the time, or to determine
 * that there is no such segment.  (Therefore, it is a good idea to pass
 * a segment list to XLALSegListCoalesce() before using it with
 * XLALSegListSearch(), unless segment list ordering or distinct
 * segments which touch/overlap are meaningful for what you are doing,
 * which is sometimes the case.)   Otherwise, it must use a linear search,
 * although a ``sorted'' list can still be searched slightly more efficiently than
 * an un-sorted list.  In all cases, XLALSegListSearch() first checks
 * whether the segment found by the last successful search contains the
 * specified time, and returns that promptly if so.
 *
 * \return a pointer to a segment in the list which
 * contains the time being searched for, or NULL if there is no such segment.
 * If more than one segment in the list contains the time, then this function
 * returns a pointer to \e one of the segments which contains it, not
 * necessarily the first such segment in the list.  (This is not an issue if
 * the list is ``disjoint'', which guarantees that it has no overlapping
 * segments.)  If no segment in the list contains the time, then this function
 * returns NULL; however, this is not really an error, use \c xlalErrno
 * to differentiate between failure and non-failure.
 */
LALSeg *
XLALSegListSearch( LALSegList *seglist, const LIGOTimeGPS *gps )
{
  int cmp;
  LALSeg *bstart = NULL;
  size_t bcount = 0;
  LALSeg *l1start=NULL, *l1end=NULL, *l2start=NULL, *l2end=NULL;
  LALSeg *segp;

  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "NULL LALSegList pointer passed to %s\n", __func__ );
    XLAL_ERROR_NULL( XLAL_EFAULT );
  }

  /* Make sure a non-null pointer was passed for the GPS time */
  if ( ! gps ) {
    XLALPrintError( "NULL LIGOTimeGPS pointer passed to %s\n", __func__ );
    XLAL_ERROR_NULL( XLAL_EFAULT );
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    XLALPrintError( "Passed unintialized LALSegList structure to %s\n", __func__ );
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  /* If the segment list is empty, simply return */
  if ( seglist->length == 0 ) {
    return NULL;
  }

  /* Do one or two quick checks based on the last segment found in a search,
     which might match this time too.  If we're not so lucky, then plan out
     the search strategy, depending on the segment list's 'disjoint' property.
     If the segment list is disjoint, then we will do a binary search;
     otherwise, we will do a linear search, possibly wrapping around from the
     end of the list to the beginning, in which case we set two pairs of
     search limits. */

  if ( seglist->disjoint ) {

    if ( ! seglist->lastFound ) {

      /* Do a binary search on the whole list */
      bstart = seglist->segs;
      bcount = seglist->length;

    } else {

      /* Check whether the time lies within the last segment found */
      cmp = XLALGPSInSeg( gps, seglist->lastFound );

      if ( cmp == 0 ) {
        return seglist->lastFound;

      } else if ( cmp < 0 ) {
        /* Do a binary search from the beginning of the list up to here */
        bstart = seglist->segs;
        bcount = ( seglist->lastFound - seglist->segs );
        /* If there are no earlier segments, just return */
        if ( bcount == 0 ) {
          return NULL;
        }

      } else {
        /* The time is later than the last segment found */

        /* If there are no later segments, just return */
        if ( seglist->lastFound == seglist->segs + seglist->length - 1 ) {
          return NULL;
        }

        /* Check whether the time lies within the next later segment */
        cmp = XLALGPSInSeg( gps, seglist->lastFound+1 );
        if ( cmp == 0 ) {
          /* We found a match, so update lastFound and return */
          seglist->lastFound++;
          return seglist->lastFound;
        } else if ( cmp < 0 ) {
          /* The time lies between the two segments */
          return NULL;
        } else {
          /* The time is later than this segment, so search to end of list */
          bstart = seglist->lastFound+2;
          bcount = ( seglist->segs - seglist->lastFound
                     + seglist->length-2 );
          /* If there are no later segments, just return */
          if ( bcount == 0 ) {
            return NULL;
          }
        }

      }

    }

  } else {

    /* The search limits are chosen the same for sorted vs. unsorted lists,
       but the actual search will behave differently for given limits */

    if ( ! seglist->lastFound ) {

      /* Do a single linear search on the whole list */
      l1start = seglist->segs;
      l1end = seglist->segs + seglist->length;
      l2start = l2end = NULL;

    } else {

      /* Check whether the time lies within the last segment found */
      cmp = XLALGPSInSeg( gps, seglist->lastFound );

      if ( cmp == 0 ) {
        return seglist->lastFound;

      } else if ( cmp < 0 ) {
        /* Do a linear search, starting from the beginning of the list.
           If the list is sorted, then we already know that we only have to
           search up to the lastFound segment, but the actual search is
           guaranteed to bail out at that point, so we can go ahead and set
           the search parameters as if the whole list would be searched. */
        l1start = seglist->segs;
        l1end = seglist->segs + seglist->length;
        l2start = l2end = NULL;

      } else {
        /* First search from here to the end of the list, then wrap around
           and search from the beginning of the list up to here */
        l1start = seglist->lastFound + 1;
        l1end = seglist->segs + seglist->length;
        l2start = seglist->segs;
        l2end = seglist->lastFound;

      }

    }

  }

  /* If we get here, then we have to do one or more searches */

  if ( seglist->disjoint ) {

    /* Do a binary search */
    /* This function call produces a compiler warning, because XLALGPSInSeg
       takes a pointer to a LIGOTimeGPS and a pointer to a LALSeg, whereas
       bsearch expects the comparison function passed to it to take void
       pointers as arguments.  Oh well. */
    segp = bsearch( (const void *) gps, (const void *) bstart,
                    bcount, sizeof(LALSeg), XLALGPSInSeg );
    /* If we found a match, update lastFound and return the pointer.
       Otherwise, return NULL. */
    if ( segp ) {
      seglist->lastFound = segp;
      return segp;
    } else {
      return NULL;
    }

  } else if ( seglist->sorted ) {

    /* Do a linear search, but bail out if we reach a segment whose start time
       is later than the target time */
    for ( segp=l1start; segp<l1end; segp++ ) {
      cmp = XLALGPSInSeg( gps, segp );
      if ( cmp == 0 ) {
        seglist->lastFound = segp;
        return segp;
      } else if ( cmp < 0 ) {
        /* This segment is beyond the target time, so we can bail out */
        break;
      }
    }
    /* Set the lastFound pointer to the last segment which is not beyond the
       target time, or the last segment in the list */
    seglist->lastFound = segp - 1;

    /* If necessary, do a second linear search, starting from the beginning
       of the list and going up to lastFound.  In this case, we know that the
       comparison value will never be negative. */
    for ( segp=l2start; segp<l2end; segp++ ) {
      if ( XLALGPSInSeg(gps,segp) == 0 ) {
        seglist->lastFound = segp;
        return segp;
      }
    }

    /* If we get here, then we didn't find a match.
       Continue to the end of the function. */

  } else {

    /* Do a linear search, looking for a match */
    for ( segp=l1start; segp<l1end; segp++ ) {
      if ( XLALGPSInSeg(gps,segp) == 0 ) {
        seglist->lastFound = segp;
        return segp;
      }
    }

    /* If necessary, do a second linear search, starting from the beginning
       of the list and going up to lastFound. */
    for ( segp=l2start; segp<l2end; segp++ ) {
      if ( XLALGPSInSeg(gps,segp) == 0 ) {
        seglist->lastFound = segp;
        return segp;
      }
    }

    /* If we get here, then we didn't find a match.
       Continue to the end of the function. */

  }

  /* If we get here, then we didn't find a match, so return NULL */
  return NULL;
}



/*---------------------------------------------------------------------------*/
/**
 * UNDOCUMENTED
 */
int
XLALSegListShift(  LALSegList *seglist, const LIGOTimeGPS *shift )
{
  unsigned i;

  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "NULL LALSegList pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure a non-null pointer was passed for the GPS time */
  if ( ! shift ) {
    XLALPrintError( "NULL LIGOTimeGPS pointer passed to %s\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    XLALPrintError( "Passed unintialized LALSegList structure to %s\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* time shift each segment in this list */
  for (i=0; i<seglist->length; i++) {
    XLALGPSAddGPS( &seglist->segs[i].start, shift);
    XLALGPSAddGPS( &seglist->segs[i].end, shift);
  }

  /* done */
  return 0;
}


/**
 * UNDOCUMENTED
 */
int
XLALSegListKeep(  LALSegList *seglist, const LIGOTimeGPS *start, const LIGOTimeGPS *end )
{
  LALSegList workspace;
  unsigned i;
  INT8 startNS, endNS;

  /* Make sure a non-null pointer was passed for the segment list */
  if ( ! seglist ) {
    XLALPrintError( "%s(): NULL LALSegList\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure a non-null pointer was passed for the GPS time */
  if ( !start || !end ) {
    XLALPrintError( "%s(): NULL boundaries\n", __func__ );
    XLAL_ERROR( XLAL_EFAULT );
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    XLALPrintError( "%s(): unintialized LALSegList\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* init the temporary list */
  XLALSegListInit( &workspace );

  /* convert the start and stop times to nanoseconds */
  startNS = XLALGPSToINT8NS( start );
  endNS = XLALGPSToINT8NS( end );
  if ( startNS >= endNS ) {
    XLALPrintError( "%s(): zero-length or improper interval\n", __func__ );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* loop over all segments in this list */
  for (i=0; i<seglist->length; i++) {
    LALSeg newSeg;
    INT8 segStart = XLALGPSToINT8NS( &seglist->segs[i].start );
    INT8 segEnd = XLALGPSToINT8NS( &seglist->segs[i].end );

    /* if segment lies entirely outside interval, discard */
    if ( segEnd <= startNS || segStart >= endNS )
      continue;

    /* intersect */
    if ( segStart < startNS )
      segStart = startNS;
    if ( segEnd > endNS )
      segEnd = endNS;

    /* store in new list */
    XLALINT8NSToGPS( &newSeg.start, segStart );
    XLALINT8NSToGPS( &newSeg.end, segEnd );
    newSeg.id = seglist->segs[i].id;
    XLALSegListAppend( &workspace, &newSeg );
  }

  /* clear the old list */
  XLALSegListClear( seglist );

  /* put the new segments into the old lsit */
  for ( i = 0; i < workspace.length; i++ )
    XLALSegListAppend( seglist, &workspace.segs[i] );

  /* clear the work space */
  XLALSegListClear( &workspace );

  /* done */
  return 0;
}

/**
 * Simple method to check whether a LALSegList is in an initialized state.
 *
 * Avoid the user having to deal with LALSegList internal 'magic'.
 */
int
XLALSegListIsInitialized ( const LALSegList *seglist )
{
  XLAL_CHECK ( seglist != NULL, XLAL_EINVAL, "Invalid NULL input 'seglist'\n");

  return (seglist->initMagic == SEGMENTSH_INITMAGICVAL);

}  /* XLALSegListIsInitialized() */

/**
 * (Re-)Initialize a segment list with Nseg 'simple' segments of length 'Tseg',
 * starting at t0 = startTime, ie
 * { [t0, t0+Tseg), [t0+Tseg, t0+2*Tseg), ... [t0+(N-1)*Tseg, t0+N*Tseg) }
 *
 * Note: accepts un-initialized segments list, as well as existing segment-lists.
 * The latter will be properly cleared before re-use.
 *
 * The 'Id' field of segment k is set to 'k'
 */
int
XLALSegListInitSimpleSegments ( LALSegList *seglist, LIGOTimeGPS startTime, UINT4 Nseg, REAL8 Tseg )
{
  XLAL_CHECK ( seglist != NULL, XLAL_EINVAL, "Invalid NULL input 'seglist'\n" );
  XLAL_CHECK ( Tseg > 0, XLAL_EDOM, "Invalid non-positive input 'Tseg=%g'\n", Tseg );

  if ( XLALSegListIsInitialized ( seglist ) )
    XLALSegListClear ( seglist );
  else
    XLALSegListInit ( seglist );

  for ( UINT4 k = 0; k < Nseg; k ++ )
    {
      LALSeg seg_k;

      LIGOTimeGPS start_k = startTime;
      XLALGPSAdd( &start_k, k * Tseg );		// t0_k = t0 + k * Tseg
      LIGOTimeGPS end_k = startTime;
      XLALGPSAdd( &end_k, (k+1) * Tseg );	// t1_k = t0 + (k+1) * Tseg

      XLAL_CHECK ( XLALSegSet ( &seg_k, &start_k, &end_k, k ) == XLAL_SUCCESS, XLAL_EFUNC );

      XLAL_CHECK ( XLALSegListAppend ( seglist, &seg_k ) == XLAL_SUCCESS, XLAL_EFUNC );

    } // for k < Nseg

  return XLAL_SUCCESS;

} /* XLALSegListInitSimpleSegments() */


/**
 * Output an (octave) formatting of 'seglist' as a string
 */
char *
XLALSegList2String ( const LALSegList *seglist )
{
  XLAL_CHECK_NULL ( seglist != NULL, XLAL_EINVAL, "Invalid NULL input 'seglist'\n" );
  XLAL_CHECK_NULL ( XLALSegListIsInitialized ( seglist ), XLAL_EINVAL, "Got invalid un-initialized seglist\n" );

  char *ret = NULL;
  XLAL_CHECK_NULL ( (ret = XLALStringAppend ( ret, "[ " )) != NULL, XLAL_ENOMEM, "Failed to ret=XLALStringAppend()\n" );

  char segfmt[64];
  sprintf ( segfmt, "%%.%df, %%.%df, %%d; ", seglist->dplaces, seglist->dplaces );	// seglist tells us output precision for GPS times to use

  UINT4 Nseg = seglist->length;
  for ( UINT4 k = 0; k < Nseg; k ++ )
    {
      char seg_buf[512];
      REAL8 t0_k = XLALGPSGetREAL8 ( &(seglist->segs[k].start) );
      REAL8 t1_k = XLALGPSGetREAL8 ( &(seglist->segs[k].end) );
      sprintf ( seg_buf, segfmt, t0_k, t1_k, seglist->segs[k].id );

      XLAL_CHECK_NULL ( (ret = XLALStringAppend ( ret, seg_buf ) ) != NULL, XLAL_ENOMEM, "Failed to ret=XLALStringAppend() for segment %d/%d\n", k+1, Nseg );

    } // for k < Nseg

  XLAL_CHECK_NULL ( (ret = XLALStringAppend ( ret, "]" ) ) != NULL, XLAL_ENOMEM, "Failed to ret=XLALStringAppend()" );

  return ret;

} /* XLALSegList2String() */

/*---------------------------------------------------------------------------*/
/**
 * Get a copy of the segment at indx in the internal array. If the segment is
 * beyond the bounds of the array, return NULL.
 */
LALSeg *
XLALSegListGet( LALSegList *seglist, UINT4 indx )
{
        if( indx >= seglist->length ){
                return NULL;
        }
        LALSeg* tmp = LALMalloc(sizeof(LALSeg));
        memcpy(tmp, &seglist->segs[indx], sizeof(LALSeg));
        return tmp;

}  /* XLALSegListGet() */
