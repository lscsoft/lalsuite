/*----------------------------------------------------------------------- 
 * 
 * File Name: LIGOMetadataUtils.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LIGOMetadataUtilsCV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>

NRCSID( LIGOMETADATAUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{LIGOMetadataUtils.c}}

\noindent General routines for manipulating LIGO metadatabase tables.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LIGOMetadataUtilsCP}
\idx{LALPlaygroundInSearchSummary()}

\subsubsection*{Description}

\noindent The function \texttt{LALPlaygroundInSearchSummary()} determines the
ammount of time in the search summary table \texttt{ssTable} that overlaps
with playground data. The time between \texttt{in\_start\_time} and
\texttt{in\_end\_time} that overlaps with playground is returned in
\texttt{inPlayTime} and the time between \texttt{out\_start\_time} and
\texttt{out\_end\_time} that overlaps with playground is returned in
\texttt{outPlayTime}.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{LIGOMetadataUtilsCV}}

</lalLaTeX>
#endif

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
  if ( seg_length >= mod_play )
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
}

/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void
LALPlaygroundInSearchSummary (
    LALStatus          *status,
    SearchSummaryTable *ssTable,
    LIGOTimeGPS        *inPlayTime,
    LIGOTimeGPS        *outPlayTime
    )
/* </lalVerbatim> */
{
  INT8 startNS, endNS, lengthNS, playNS;

  INITSTATUS( status, "LALPlaygroundInSearchSummary", LIGOMETADATAUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( ssTable, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ssTable, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ssTable, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  LALGPStoINT8( status->statusPtr, &startNS, &(ssTable->in_start_time) );
  CHECKSTATUSPTR( status );
  LALGPStoINT8( status->statusPtr, &endNS, &(ssTable->in_end_time) );
  CHECKSTATUSPTR( status );
  lengthNS = endNS - startNS;
  playNS = PlaygroundOverlap( endNS, lengthNS );
  if ( playNS < 0 )
  {
    ABORT( status, LIGOMETADATAUTILSH_ETIME, LIGOMETADATAUTILSH_MSGETIME );
  }
  LALINT8toGPS( status->statusPtr, inPlayTime, &playNS );
  CHECKSTATUSPTR( status );

  LALGPStoINT8( status->statusPtr, &startNS, &(ssTable->out_start_time) );
  CHECKSTATUSPTR( status );
  LALGPStoINT8( status->statusPtr, &endNS, &(ssTable->out_end_time) );
  CHECKSTATUSPTR( status );
  lengthNS = endNS - startNS;

  playNS = PlaygroundOverlap( endNS, lengthNS );
  if ( playNS < 0 )
  {
    ABORT( status, LIGOMETADATAUTILSH_ETIME, LIGOMETADATAUTILSH_MSGETIME );
  }
  LALINT8toGPS( status->statusPtr, outPlayTime, &playNS );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
int
LALCompareSearchSummaryByInTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;

  SearchSummaryTable *aPtr = *((SearchSummaryTable **)a);
  SearchSummaryTable *bPtr = *((SearchSummaryTable **)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the in start times */
  LALGPStoINT8( &status, &ta, &(aPtr->in_start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->in_start_time) );

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
    LALGPStoINT8( &status, &ta, &( aPtr->in_end_time) );
    LALGPStoINT8( &status, &tb, &( bPtr->in_end_time) );

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


/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
int
LALCompareSearchSummaryByOutTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;

  SearchSummaryTable *aPtr = *((SearchSummaryTable **)a);
  SearchSummaryTable *bPtr = *((SearchSummaryTable **)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the out start times */
  LALGPStoINT8( &status, &ta, &(aPtr->out_start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->out_start_time) );

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
    LALGPStoINT8( &status, &ta, &(aPtr->out_end_time) );
    LALGPStoINT8( &status, &tb, &(bPtr->out_end_time) );

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

/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void
LALTimeSortSearchSummary (
    LALStatus            *status,
    SearchSummaryTable  **summHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numSumms = 0;
  SearchSummaryTable    *thisSearchSumm = NULL;
  SearchSummaryTable   **summHandle = NULL;

  INITSTATUS( status, "LALTimeSortSearchSummary", LIGOMETADATAUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( summHead, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* count the number of summs in the linked list */
  for ( thisSearchSumm = *summHead; thisSearchSumm; 
      thisSearchSumm = thisSearchSumm->next )
  {
    ++numSumms;
  }
  if ( ! numSumms )
  {
    LALWarning( status, "No summs in list to sort" );
    RETURN( status );
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
  for ( i = 1; i < numSumms; ++i )
  {
    thisSearchSumm = thisSearchSumm->next = summHandle[i];
  }
  thisSearchSumm->next = NULL;

  /* free the internal memory */
  LALFree( summHandle );

  DETATCHSTATUSPTR (status);
  RETURN( status );
}


/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void
LALIfoScanSearchSummary(
    LALStatus                  *status,
    SearchSummaryTable        **output,
    SearchSummaryTable         *input,
    CHAR                       *ifos
    )
/* </lalVerbatim> */
{
  SearchSummaryTable    *thisSearchSumm = NULL;
  SearchSummaryTable    *keptSumm = NULL;

  INITSTATUS( status, "LALIfoScanSearchSummary", LIGOMETADATAUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that output is null and input non-null */
  ASSERT( !(*output), status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  ASSERT( input, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of search_summary tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  for( thisSearchSumm = input; thisSearchSumm; 
      thisSearchSumm = thisSearchSumm->next )
  {

    if ( !strcmp(thisSearchSumm->ifos, ifos) ) 
    {
      /* IFOs match so write this entry to the output table */
      if ( ! *output  )
      {
	*output = keptSumm = (SearchSummaryTable *) 
	  LALMalloc( sizeof(SearchSummaryTable) );
      }
      else
      {
	keptSumm = keptSumm->next = (SearchSummaryTable *) 
	  LALMalloc( sizeof(SearchSummaryTable) );
      }
      memcpy(keptSumm, thisSearchSumm, sizeof(SearchSummaryTable));
      keptSumm->next = NULL;
    }
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

}  

/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void
LALCheckOutTimeFromSearchSummary (
    LALStatus            *status,
    SearchSummaryTable   *summList,
    CHAR                 *ifo,
    LIGOTimeGPS          *startTime,
    LIGOTimeGPS          *endTime
    )
/* </lalVerbatim> */
{
  SearchSummaryTable   *thisIFOSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;
  INT8  startTimeNS = 0;
  INT8  endTimeNS = 0;
  INT8  unsearchedStartNS = 0;
  INT8  outStartNS = 0;
  INT8  outEndNS = 0;


  INITSTATUS( status, "LALCheckOutTimeSearchSummary", LIGOMETADATAUTILSC );
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
  LALGPStoINT8( status->statusPtr, &startTimeNS, startTime );
  CHECKSTATUSPTR( status );
  LALGPStoINT8( status->statusPtr, &endTimeNS, endTime );
  CHECKSTATUSPTR( status );

  unsearchedStartNS = startTimeNS;

  /* check that all times are searched */
  for ( thisSearchSumm = thisIFOSummList; thisSearchSumm; 
      thisSearchSumm = thisSearchSumm->next )
  {
    LALGPStoINT8( status->statusPtr, &outStartNS, 
	&(thisSearchSumm->out_start_time) );
    CHECKSTATUSPTR( status );

    if ( outStartNS < startTimeNS )    
    {
      /* file starts before requested start time */
      if ( outEndNS > startTimeNS )
      {
	/* file is partially in requested times, update unsearchedStart */
	LALGPStoINT8( status->statusPtr, &outEndNS, 
	    &(thisSearchSumm->out_end_time) );
	CHECKSTATUSPTR( status );

	unsearchedStartNS = outEndNS;
      }
    }
    else if ( outStartNS == unsearchedStartNS )
    {
      /* this file starts at the beginning of the unsearched data */
      /* calculate the end time and set unsearched start to this */
      LALGPStoINT8( status->statusPtr, &outEndNS, 
	  &(thisSearchSumm->out_end_time) );
      CHECKSTATUSPTR( status );

      unsearchedStartNS = outEndNS;
    }
    else if ( outStartNS > unsearchedStartNS )
    {
      /* there is a gap in the searched data between unsearchedStart
	 and outStart */
      ABORT( status, LIGOMETADATAUTILSH_ESSGAP, LIGOMETADATAUTILSH_MSGESSGAP );
    }
    else if ( outStartNS < unsearchedStartNS )    
    {
      /* there is a region of data which was searched twice */
      ABORT( status, LIGOMETADATAUTILSH_ESSDUB, LIGOMETADATAUTILSH_MSGESSDUB );
    }
  }

  /* check that we got to the end of the requested time */
  if ( unsearchedStartNS < outEndNS )
  {
    ABORT( status, LIGOMETADATAUTILSH_ESSGAP, LIGOMETADATAUTILSH_MSGESSGAP );
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



/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void
LALIfoScanSummValue(
    LALStatus                  *status,
    SummValueTable            **output,
    SummValueTable             *input,
    CHAR                       *ifo
    )
/* </lalVerbatim> */
{
  SummValueTable    *thisSummValue = NULL;
  SummValueTable    *keptSumm = NULL;

  INITSTATUS( status, "LALIfoScanSummValue", LIGOMETADATAUTILSC );
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



/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
int
LALCompareSummValueByTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;

  SummValueTable *aPtr = *((SummValueTable **)a);
  SummValueTable *bPtr = *((SummValueTable **)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the out start times */
  LALGPStoINT8( &status, &ta, &(aPtr->start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->start_time) );

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
    LALGPStoINT8( &status, &ta, &(aPtr->end_time) );
    LALGPStoINT8( &status, &tb, &(bPtr->end_time) );

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

/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void
LALTimeSortSummValue (
    LALStatus            *status,
    SummValueTable      **summHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numSumms = 0;
  SummValueTable    *thisSummValue = NULL;
  SummValueTable   **summHandle = NULL;

  INITSTATUS( status, "LALTimeSortSummValue", LIGOMETADATAUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( summHead, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* count the number of summs in the linked list */
  for ( thisSummValue = *summHead; thisSummValue; 
      thisSummValue = thisSummValue->next )
  {
    ++numSumms;
  }
  if ( ! numSumms )
  {
    LALWarning( status, "No summs in list to sort" );
    RETURN( status );
  }

  /* allocate memory for an array of ptrs to sort and populate array */
  summHandle = (SummValueTable **) 
    LALCalloc( numSumms, sizeof(SummValueTable *) );
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

  DETATCHSTATUSPTR (status);
  RETURN( status );
}

