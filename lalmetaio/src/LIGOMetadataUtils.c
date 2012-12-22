/*
*  Copyright (C) 2007 Drew Keppel, Duncan Brown, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Thomas Cokelaer
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
 * File Name: LIGOMetadataUtils.c
 *
 * Author: Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**

\author Brown, D. A.
\file
\ingroup lalmetaio
\brief General routines for manipulating LIGO metadatabase tables.

\heading{Description}

The function <tt>LALPlaygroundInSearchSummary()</tt> determines the
ammount of time in the search summary table \c ssTable that overlaps
with playground data. The time between \c in_start_time and
\c in_end_time that overlaps with playground is returned in
\c inPlayTime and the time between \c out_start_time and
\c out_end_time that overlaps with playground is returned in
\c outPlayTime.

<tt>LALCompareSearchSummaryByInTime()</tt> is a function to compare the in
times in two search summary tables.  It returns 1 if the
\c in_start_time of the first table is after the
\c in_start_time of the second and -1 if it is before.  If the two
\c in_start_times are identical, the test is repeated on the
\c in_end_times.  If these are also equal, the comparison returns 0.
<tt>LALCompareSearchSummaryByOutTime()</tt> operates in a similar manner, but
uses the out, rather than in, times.

<tt>LALTimeSortSearchSummary()</tt> will time sort a linked list of search
summary tables.  You can sort on in our out start time depending which
\c comparfunc is specified.

<tt>LALIfoScanSearchSummary()</tt> steps through a linked list of search
summary tables and returns a pointer \c output to a linked list of those
tables whos \c ifos field matches the string \c ifos.


<tt>LALIfoScanSummValue()</tt>, <tt>LALCompareSummValueByTime()</tt> and
<tt>LALTimeSortSummValue()</tt> performs the same functions as described
above.  The only difference being that they act on summ value tables.

<tt>LALCheckOutTimeFromSearchSummary()</tt> verifies that all times
between the specified \c startTime and \c endTime have been
searched precisely once for the given \c ifo.

Finally, <tt>LALDistanceScanSummValueTable()</tt> scan a summ value table
 searching for a trigger belonging to a given ifo and englobing a give GPS
 time.


\heading{Algorithm}

None.

\heading{Uses}

LALCalloc, LALMalloc, LALFree.

\heading{Notes}
%% Any relevant notes.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>

int XLALCountProcessTable(ProcessTable *head)

{
	int length;

	/* count the number of events in the list */
	for(length = 0; head; head = head->next)
		length++;

	return(length);
}


int XLALCountProcessParamsTable(ProcessParamsTable *head)

{
	int length;

	/* count the number of events in the list */
	for(length = 0; head; head = head->next)
		length++;

	return(length);
}


int XLALCountMultiInspiralTable(MultiInspiralTable *head)

{
	int length;
	/* count the number of events in the list */
	for(length = 0; head; head = head->next)
		length++;

	return(length);
}



int
XLALIFONumber(
    const char *ifo
    )

{
  switch( ifo[0] )
  {
    case 'G':
      return LAL_IFO_G1;
      break;

    case 'H':
      if ( !strcmp( ifo, "H1" ) )
      {
        return LAL_IFO_H1;
      }
      else if (!strcmp( ifo, "H2" ) )
      {
        return LAL_IFO_H2;
      }
      else
      {
        /* Invalid Hanford Detector */
        return LAL_UNKNOWN_IFO ;
      }
      break;

    case 'L':
      return LAL_IFO_L1;
      break;

    case 'T':
      return LAL_IFO_T1;
      break;

    case 'V':
      return LAL_IFO_V1;
      break;

    default:
      /* Invalid Detector Site */
      return LAL_UNKNOWN_IFO ;
  }
}


void
XLALReturnIFO(
    char                *ifo,
    InterferometerNumber IFONumber
    )

{
  switch( IFONumber )
  {
    case LAL_IFO_G1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "G1");
      break;

    case LAL_IFO_H1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "H1");
      break;

    case LAL_IFO_H2:
      snprintf( ifo, LIGOMETA_IFO_MAX, "H2");
      break;

    case LAL_IFO_L1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "L1");
      break;

    case LAL_IFO_T1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "T1");
      break;

    case LAL_IFO_V1:
      snprintf( ifo, LIGOMETA_IFO_MAX, "V1");
      break;

    default:
      /* Invalid Detector Site */
      snprintf( ifo, LIGOMETA_IFO_MAX, " ");
  }
}



void
XLALReturnDetector(
    LALDetector           *det,
    InterferometerNumber   IFONumber
    )

{
  switch( IFONumber )
  {
    case LAL_IFO_G1:
      *det = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      break;

    case LAL_IFO_H1:
      *det = lalCachedDetectors[LALDetectorIndexLHODIFF];
      break;

    case LAL_IFO_H2:
      *det = lalCachedDetectors[LALDetectorIndexLHODIFF];
      break;

    case LAL_IFO_L1:
      *det = lalCachedDetectors[LALDetectorIndexLLODIFF];
      break;

    case LAL_IFO_T1:
      *det = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      break;

    case LAL_IFO_V1:
      *det = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      break;

    default:
      /* Invalid Detector Site */
      memset(det, 0, sizeof(LALDetector) );
  }
}



static INT8 PlaygroundOverlap( INT8 seg_end, INT8 seg_length )
{
  const INT8 play_length = LAL_INT8_C(600000000000);
  const INT8 S2_start = LAL_INT8_C(729273613000000000);
  const INT8 mod_play = LAL_INT8_C(6370000000000);
  INT8 end_mod_play;

  /* if the end precedes the start of S2, then there is no playground, since
   * we only define playground from S2 onwards */
  if ( seg_end < S2_start )
  {
    return LAL_INT8_C(0);
  }

  /* handle a segment that contains two or more playground */
  /* segments by recursively bisecting it                  */
  if ( seg_length >= (mod_play - play_length) )
  {
    INT8 low_len, high_len, low_play, high_play;
    low_len = high_len = seg_length / LAL_INT8_C(2);
    if ( seg_length % LAL_INT8_C(2) ) ++low_len;

    low_play = PlaygroundOverlap( seg_end - high_len, low_len );
    high_play = PlaygroundOverlap( seg_end, high_len );

    return low_play + high_play;
  }

  end_mod_play = ((seg_end - S2_start) % mod_play );

  /* if no overlap with playground, return zero */
  if ( end_mod_play >= play_length + seg_length )
  {
    return LAL_INT8_C(0);
  }
  else
  {
    if ( seg_length > play_length )
    {
      if ( (seg_length < end_mod_play) &&
          (end_mod_play < play_length + seg_length ) )
      {
        return play_length + seg_length - end_mod_play;
      }
      else if ( (play_length <= end_mod_play) &&
          (end_mod_play <= seg_length ) )
      {
        return play_length;
      }
      else if ( end_mod_play < play_length )
      {
        return end_mod_play;
      }
    }
    else if ( seg_length <= play_length )
    {
      if ( (play_length < end_mod_play) &&
          (end_mod_play < seg_length + play_length ) )
      {
        return play_length + seg_length - end_mod_play;
      }
      else if ( (seg_length <= end_mod_play) &&
          (end_mod_play <= play_length ) )
      {
        return seg_length;
      }
      else if ( end_mod_play < seg_length )
      {
        return end_mod_play;
      }
    }
    else
    {
      LALPrintError( "Error determining playground overlap\n" );
      return LAL_INT8_C(-1);
    }
  }
  /* JC: HOPEFULLY NEVER GET HERE */
  return LAL_INT8_C(-1);
}


void
LALPlaygroundInSearchSummary (
    LALStatus          *status,
    SearchSummaryTable *ssTable,
    LIGOTimeGPS        *inPlayTime,
    LIGOTimeGPS        *outPlayTime
    )

{
  INT4 playCheck = 0;
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  playCheck = XLALPlaygroundInSearchSummary ( ssTable, inPlayTime,
      outPlayTime );

  if ( playCheck < 0 )
  {
    ABORT( status, LIGOMETADATAUTILSH_ETIME, LIGOMETADATAUTILSH_MSGETIME );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


int
XLALPlaygroundInSearchSummary (
    SearchSummaryTable *ssTable,
    LIGOTimeGPS        *inPlayTime,
    LIGOTimeGPS        *outPlayTime
    )

{
  INT8 startNS, endNS, lengthNS, playNS;

  startNS = XLALGPSToINT8NS( &(ssTable->in_start_time) );
  endNS = XLALGPSToINT8NS( &(ssTable->in_end_time) );
  lengthNS = endNS - startNS;
  playNS = PlaygroundOverlap( endNS, lengthNS );
  if ( playNS < 0 )
  {
    XLAL_ERROR(XLAL_EIO);
  }
  inPlayTime = XLALINT8NSToGPS( inPlayTime, playNS );

  startNS = XLALGPSToINT8NS( &(ssTable->out_start_time) );
  endNS = XLALGPSToINT8NS( &(ssTable->out_end_time) );
  lengthNS = endNS - startNS;

  playNS = PlaygroundOverlap( endNS, lengthNS );
  if ( playNS < 0 )
  {
    XLAL_ERROR(XLAL_EIO);
  }
  outPlayTime = XLALINT8NSToGPS( outPlayTime, playNS );

  return( 0 );
}



int
LALCompareSearchSummaryByInTime (
    const void *a,
    const void *b
    )

{
  LALStatus     status;

  const SearchSummaryTable *aPtr = *((const SearchSummaryTable * const *)a);
  const SearchSummaryTable *bPtr = *((const SearchSummaryTable * const *)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the in start times */
  memset( &status, 0, sizeof(LALStatus) );
  ta = XLALGPSToINT8NS( &(aPtr->in_start_time) );
  tb = XLALGPSToINT8NS( &(bPtr->in_start_time) );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    /* determine the in end times */
    memset( &status, 0, sizeof(LALStatus) );
    ta = XLALGPSToINT8NS( &( aPtr->in_end_time) );
    tb = XLALGPSToINT8NS( &( bPtr->in_end_time) );

    if ( ta > tb )
    {
      return 1;
    }
    else if ( ta < tb )
    {
      return -1;
    }
    else
    {
      return 0;
    }
  }
}



int
LALCompareSearchSummaryByOutTime (
    const void *a,
    const void *b
    )

{
  LALStatus     status;

  const SearchSummaryTable *aPtr = *((const SearchSummaryTable * const *)a);
  const SearchSummaryTable *bPtr = *((const SearchSummaryTable * const *)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the out start times */
  memset( &status, 0, sizeof(LALStatus) );
  ta = XLALGPSToINT8NS( &(aPtr->out_start_time) );
  tb = XLALGPSToINT8NS( &(bPtr->out_start_time) );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    /* determine the out end times */
    memset( &status, 0, sizeof(LALStatus) );
    ta = XLALGPSToINT8NS( &(aPtr->out_end_time) );
    tb = XLALGPSToINT8NS( &(bPtr->out_end_time) );

    if ( ta > tb )
    {
      return 1;
    }
    else if ( ta < tb )
    {
      return -1;
    }
    else
    {
      return 0;
    }
  }
}


int
XLALTimeSortSearchSummary(
    SearchSummaryTable  **summHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INT4                  i;
  INT4                  numSumms = 0;
  SearchSummaryTable    *thisSearchSumm = NULL;
  SearchSummaryTable   **summHandle = NULL;

  if ( !summHead )
  {
    XLAL_ERROR(XLAL_EIO);
  }


  /* count the number of summs in the linked list */
  for ( thisSearchSumm = *summHead; thisSearchSumm;
      thisSearchSumm = thisSearchSumm->next )
  {
    ++numSumms;
  }
  if ( ! numSumms )
  {
    return( 0 );
  }

  /* allocate memory for an array of ptrs to sort and populate array */
  summHandle = (SearchSummaryTable **)
    LALCalloc( numSumms, sizeof(SearchSummaryTable *) );
  for ( i = 0, thisSearchSumm = *summHead; i < numSumms;
      ++i, thisSearchSumm = thisSearchSumm->next )
  {
    summHandle[i] = thisSearchSumm;
  }

  /* qsort the array using the specified function */
  qsort( summHandle, numSumms, sizeof(summHandle[0]), comparfunc );

  /* re-link the linked list in the right order */
  thisSearchSumm = *summHead = summHandle[0];
  for ( i = 1; i < numSumms; ++i, thisSearchSumm = thisSearchSumm->next )
  {
    thisSearchSumm->next = summHandle[i];
  }
  thisSearchSumm->next = NULL;

  /* free the internal memory */
  LALFree( summHandle );

  return( 0 );
}




void
LALTimeSortSearchSummary (
    LALStatus            *status,
    SearchSummaryTable  **summHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( summHead, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  XLALTimeSortSearchSummary( summHead, comparfunc );

  DETATCHSTATUSPTR (status);
  RETURN( status );
}



SearchSummaryTable *
XLALIfoScanSearchSummary(
    SearchSummaryTable         *input,
    CHAR                       *ifos
    )

{
  SearchSummaryTable    *output = NULL;
  SearchSummaryTable    *thisSearchSumm = NULL;
  SearchSummaryTable    *keptSumm = NULL;


  if ( !input )
  {
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* Scan through a linked list of search_summary tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  for( thisSearchSumm = input; thisSearchSumm;
      thisSearchSumm = thisSearchSumm->next )
  {

    if ( !strcmp(thisSearchSumm->ifos, ifos) )
    {
      /* IFOs match so write this entry to the output table */
      if ( ! output  )
      {
        output = keptSumm = (SearchSummaryTable *)
          LALMalloc( sizeof(SearchSummaryTable) );
      }
      else
      {
        keptSumm = keptSumm->next = (SearchSummaryTable *)
          LALMalloc( sizeof(SearchSummaryTable) );
      }
      if ( !keptSumm )
      {
        while ( output )
        {
          thisSearchSumm = output;
          output = (output)->next;
          LALFree( thisSearchSumm );
        }
        XLAL_ERROR_NULL(XLAL_ENOMEM);
      }
      memcpy(keptSumm, thisSearchSumm, sizeof(SearchSummaryTable));
      keptSumm->next = NULL;
    }
  }
  return( output);
}




void
LALIfoScanSearchSummary(
    LALStatus                  *status,
    SearchSummaryTable        **output,
    SearchSummaryTable         *input,
    CHAR                       *ifos
    )

{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  *output = XLALIfoScanSearchSummary( input, ifos );

  DETATCHSTATUSPTR (status);
  RETURN (status);

}


void
LALDistanceScanSummValueTable (
    LALStatus            *status,
    SummValueTable       *summValueList,
    LIGOTimeGPS          gps,
    const CHAR           *ifo,
    REAL4                *distance)

{
  SummValueTable    *thisSummValue = NULL;
  /*INT4 test=0;*/
  INT8 ta=0, tb=0, tc=0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* initialize diatnce to zero */
  *distance = 0;

  /* convert the input GPS time into INT8 */
  ta = XLALGPSToINT8NS( &(gps) );

  /* scan the summ value table */
  for( thisSummValue = summValueList; thisSummValue;
      thisSummValue = thisSummValue->next )
  {

    /* if this summvar is actually the distance */
    if ( ! strncmp( thisSummValue->name, "inspiral_effective_distance",
		    LIGOMETA_SUMMVALUE_NAME_MAX ) )
      {
	/* if this is the requested ifo */
	if ( !strcmp(thisSummValue->ifo, ifo) )
	  {
	    /* IFOs match so now let us check if this entry coincides
	       with the requested GPS time */

	    tb = XLALGPSToINT8NS( &(thisSummValue->start_time) );
	    tc = XLALGPSToINT8NS( &(thisSummValue->end_time) );
	    if ( ta >= tb && ta<=tc )
	      {
		*distance = thisSummValue->value;
		break;
	      }
	  }
      }
  }

  if ( *distance == 0 )
    {
      ABORT ( status, LIGOMETADATAUTILSH_EDIST, LIGOMETADATAUTILSH_MSGEDIST );
    }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALCheckOutTimeFromSearchSummary (
    LALStatus            *status,
    SearchSummaryTable   *summList,
    CHAR                 *ifo,
    LIGOTimeGPS          *startTime,
    LIGOTimeGPS          *endTime
    )

{
  SearchSummaryTable   *thisIFOSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;
  INT8  startTimeNS = 0;
  INT8  endTimeNS = 0;
  INT8  unsearchedStartNS = 0;
  INT8  outStartNS = 0;
  INT8  outEndNS = 0;


  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* check that the data has been searched once
     and only once for the given IFO */

  /* first, create a list of search summary tables applicable to this IFO */
  LALIfoScanSearchSummary( status->statusPtr,  &thisIFOSummList, summList,
      ifo );
  CHECKSTATUSPTR( status );

  /* now, time sort the output list */
  LALTimeSortSearchSummary ( status->statusPtr,  &thisIFOSummList,
      *LALCompareSearchSummaryByOutTime );
  CHECKSTATUSPTR( status );

  /* calculate requested start and end time in NS */
  startTimeNS = XLALGPSToINT8NS( startTime );
  endTimeNS = XLALGPSToINT8NS( endTime );

  unsearchedStartNS = startTimeNS;

  /* check that all times are searched */
  for ( thisSearchSumm = thisIFOSummList; thisSearchSumm;
      thisSearchSumm = thisSearchSumm->next )
  {
    outStartNS = XLALGPSToINT8NS( &(thisSearchSumm->out_start_time) );

    if ( outStartNS < startTimeNS )
    {
      /* file starts before requested start time */
      outEndNS = XLALGPSToINT8NS( &(thisSearchSumm->out_end_time) );

      if ( outEndNS > startTimeNS )
      {
        /* file is partially in requested times, update unsearchedStart */
        unsearchedStartNS = outEndNS;
      }
    }
    else if ( outStartNS == unsearchedStartNS )
    {
      /* this file starts at the beginning of the unsearched data */
      /* calculate the end time and set unsearched start to this */
      outEndNS = XLALGPSToINT8NS( &(thisSearchSumm->out_end_time) );

      unsearchedStartNS = outEndNS;
    }
    else if ( outStartNS > unsearchedStartNS )
    {
      /* there is a gap in the searched data between unsearchedStart
         and outStart */
      ABORT( status, LIGOMETADATAUTILSH_ESGAP, LIGOMETADATAUTILSH_MSGESGAP );
    }
    else if ( outStartNS < unsearchedStartNS )
    {
      /* there is a region of data which was searched twice */
      ABORT( status, LIGOMETADATAUTILSH_ESDUB, LIGOMETADATAUTILSH_MSGESDUB );
    }
  }

  /* check that we got to the end of the requested time */
  if ( unsearchedStartNS < endTimeNS )
  {
    ABORT( status, LIGOMETADATAUTILSH_ESGAP, LIGOMETADATAUTILSH_MSGESGAP );
  }

  /* free memory allocated in LALIfoScanSearchSummary */
  while ( thisIFOSummList )
  {
    thisSearchSumm = thisIFOSummList;
    thisIFOSummList = thisIFOSummList->next;
    LALFree( thisSearchSumm );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

}




void
LALIfoScanSummValue(
    LALStatus                  *status,
    SummValueTable            **output,
    SummValueTable             *input,
    CHAR                       *ifo
    )

{
  SummValueTable    *thisSummValue = NULL;
  SummValueTable    *keptSumm = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* check that output is null and input non-null */
  ASSERT( !(*output), status,
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  ASSERT( input, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of search_summary tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  for( thisSummValue = input; thisSummValue;
      thisSummValue = thisSummValue->next )
  {

    if ( !strcmp(thisSummValue->ifo, ifo) )
    {
      /* IFOs match so write this entry to the output table */
      if ( ! *output  )
      {
        *output = keptSumm = (SummValueTable *)
          LALMalloc( sizeof(SummValueTable) );
      }
      else
      {
        keptSumm = keptSumm->next = (SummValueTable *)
          LALMalloc( sizeof(SummValueTable) );
      }
      memcpy(keptSumm, thisSummValue, sizeof(SummValueTable));
      keptSumm->next = NULL;
    }
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

}




int
LALCompareSummValueByTime (
    const void *a,
    const void *b
    )

{
  const SummValueTable *aPtr = *((const SummValueTable * const *)a);
  const SummValueTable *bPtr = *((const SummValueTable * const *)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the out start times */
  ta = XLALGPSToINT8NS( &(aPtr->start_time) );
  tb = XLALGPSToINT8NS( &(bPtr->start_time) );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    /* determine the out end times */
    ta = XLALGPSToINT8NS( &(aPtr->end_time) );
    tb = XLALGPSToINT8NS( &(bPtr->end_time) );

    if ( ta > tb )
    {
      return 1;
    }
    else if ( ta < tb )
    {
      return -1;
    }
    else
    {
      return 0;
    }
  }
}


int
XLALTimeSortSummValue(
    SummValueTable      **summHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INT4                  i;
  INT4                  numSumms = 0;
  SummValueTable    *thisSummValue = NULL;
  SummValueTable   **summHandle = NULL;

  if ( !summHead )
  {
    XLAL_ERROR(XLAL_EIO);
  }

  /* count the number of summs in the linked list */
  for ( thisSummValue = *summHead; thisSummValue;
      thisSummValue = thisSummValue->next )
  {
    ++numSumms;
  }
  if ( ! numSumms )
  {
    return(0);
  }

  /* allocate memory for an array of ptrs to sort and populate array */
  summHandle = (SummValueTable **)
    LALCalloc( numSumms, sizeof(SummValueTable *) );
  if ( !summHandle )
  {
    XLAL_ERROR(XLAL_ENOMEM);
  }

  for ( i = 0, thisSummValue = *summHead; i < numSumms;
      ++i, thisSummValue = thisSummValue->next )
  {
    summHandle[i] = thisSummValue;
  }

  /* qsort the array using the specified function */
  qsort( summHandle, numSumms, sizeof(summHandle[0]), comparfunc );

  /* re-link the linked list in the right order */
  thisSummValue = *summHead = summHandle[0];
  for ( i = 1; i < numSumms; ++i )
  {
    thisSummValue = thisSummValue->next = summHandle[i];
  }
  thisSummValue->next = NULL;

  /* free the internal memory */
  LALFree( summHandle );

  return(0);
}



void
LALTimeSortSummValue (
    LALStatus            *status,
    SummValueTable      **summHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( summHead, status,
       LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  XLALTimeSortSummValue( summHead, comparfunc );

  DETATCHSTATUSPTR (status);
  RETURN( status );
}



/**
 * Create a ProcessTable structure.
 */


ProcessTable *XLALCreateProcessTableRow(void)
{
  ProcessTable *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  memset(new->program, 0, sizeof(new->program));
  memset(new->version, 0, sizeof(new->version));
  memset(new->cvs_repository, 0, sizeof(new->cvs_repository));
  XLALGPSSet(&new->cvs_entry_time, 0, 0);
  memset(new->comment, 0, sizeof(new->comment));
  new->is_online = 0;
  memset(new->node, 0, sizeof(new->node));
  memset(new->username, 0, sizeof(new->username));
  XLALGPSSet(&new->start_time, 0, 0);
  XLALGPSSet(&new->end_time, 0, 0);
  new->jobid = 0;
  memset(new->domain, 0, sizeof(new->domain));
  new->unix_procid = 0;
  memset(new->ifos, 0, sizeof(new->ifos));
  new->process_id = -1;

  return new;
}


/**
 * Destroy a ProcessTable structure.
 */


void XLALDestroyProcessTableRow(ProcessTable *row)
{
  XLALFree(row);
}


/**
 * Destroy a ProcessTable linked list.
 */


void XLALDestroyProcessTable(ProcessTable *head)
{
  while(head)
  {
    ProcessTable *next = head->next;
    XLALDestroyProcessTableRow(head);
    head = next;
  }
}


/**
 * Return the next available process ID.
 */


long XLALProcessTableGetNextID(ProcessTable *head)
{
  long highest = -1;
  for(; head; head = head->next)
    if(head->process_id > highest)
      highest = head->process_id;
  return highest + 1;
}


/**
 * Create a ProcessParamsTable structure.
 */


ProcessParamsTable *XLALCreateProcessParamsTableRow(const ProcessTable *process)
{
  ProcessParamsTable *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  memset(new->program, 0, sizeof(new->program));
  if(process)
    new->process_id = process->process_id;
  else
    new->process_id = -1;
  memset(new->param, 0, sizeof(new->param));
  memset(new->type, 0, sizeof(new->type));
  memset(new->value, 0, sizeof(new->value));

  return new;
}


/**
 * Destroy a ProcessParamsTable structure.
 */


void XLALDestroyProcessParamsTableRow(ProcessParamsTable *row)
{
  XLALFree(row);
}


/**
 * Destroy a ProcessParamsTable linked list.
 */


void XLALDestroyProcessParamsTable(ProcessParamsTable *head)
{
  while(head)
  {
    ProcessParamsTable *next = head->next;
    XLALDestroyProcessParamsTableRow(head);
    head = next;
  }
}


/**
 * Create a TimeSlide structure.
 */


TimeSlide *XLALCreateTimeSlide(void)
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


void XLALDestroyTimeSlide(TimeSlide *row)
{
  XLALFree(row);
}


/**
 * Destroy a TimeSlide linked list.
 */


void XLALDestroyTimeSlideTable(TimeSlide *head)
{
  while(head)
  {
    TimeSlide *next = head->next;
    XLALDestroyTimeSlide(head);
    head = next;
  }
}


/**
 * Find and return the address of the first element in the linked list of
 * TimeSlide objects whose time_slide_id and instrument name equal the
 * values given.  TimeSlide elements whose instrument pointer is NULL are
 * skipped.  Returns NULL if no matching row is found.  There are two
 * versions, one for cost * TimeSlide rows and one for non-const (neither
 * modifies the TimeSlide rows, the two versions are identical, they are
 * provided to allow the const'ness to be "passed" through the function).
 */


const TimeSlide *XLALTimeSlideConstGetByIDAndInstrument(const TimeSlide *time_slide, long time_slide_id, const char *instrument)
{
	for(; time_slide && (time_slide->time_slide_id != time_slide_id || !time_slide->instrument || strcmp(time_slide->instrument, instrument)); time_slide = time_slide->next);
	return time_slide;
}


TimeSlide *XLALTimeSlideGetByIDAndInstrument(TimeSlide *time_slide, long time_slide_id, const char *instrument)
{
	for(; time_slide && (time_slide->time_slide_id != time_slide_id || !time_slide->instrument || strcmp(time_slide->instrument, instrument)); time_slide = time_slide->next);
	return time_slide;
}


/**
 * Create a SearchSummaryTable structure.
 */


SearchSummaryTable *XLALCreateSearchSummaryTableRow(const ProcessTable *process)
{
  SearchSummaryTable *new = XLALMalloc(sizeof(*new));

  if(!new)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  new->next = NULL;
  if(process)
    new->process_id = process->process_id;
  else
    new->process_id = -1;
  memset(new->comment, 0, sizeof(new->comment));
  XLALGPSSet(&new->in_start_time, 0, 0);
  XLALGPSSet(&new->in_end_time, 0, 0);
  XLALGPSSet(&new->out_start_time, 0, 0);
  XLALGPSSet(&new->out_end_time, 0, 0);
  new->nevents = -1;
  new->nnodes = -1;
  memset(new->ifos, 0, sizeof(new->ifos));

  return new;
}


/**
 * Destroy a SearchSummaryTable structure.
 */


void XLALDestroySearchSummaryTableRow(SearchSummaryTable *row)
{
  XLALFree(row);
}


/**
 * Destroy a SearchSummaryTable linked list.
 */


void XLALDestroySearchSummaryTable(SearchSummaryTable *head)
{
  while(head)
  {
    SearchSummaryTable *next = head->next;
    XLALDestroySearchSummaryTableRow(head);
    head = next;
  }
}
