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
\idx{LALCompareSearchSummaryByInTime()}
\idx{LALCompareSearchSummaryByOutTime ()}
\idx{LALTimeSortSearchSummary()}
\idx{LALIfoScanSearchSummary()}
\idx{LALIfoScanSummValue()}
\idx{LALCompareSummValueByTime()}
\idx{LALTimeSortSummValue()}
\idx{LALCheckOutTimeFromSearchSummary()}
  
\subsubsection*{Description}

The function \texttt{LALPlaygroundInSearchSummary()} determines the
ammount of time in the search summary table \texttt{ssTable} that overlaps
with playground data. The time between \texttt{in\_start\_time} and
\texttt{in\_end\_time} that overlaps with playground is returned in
\texttt{inPlayTime} and the time between \texttt{out\_start\_time} and
\texttt{out\_end\_time} that overlaps with playground is returned in
\texttt{outPlayTime}.

\texttt{LALCompareSearchSummaryByInTime()} is a function to compare the in
times in two search summary tables.  It returns 1 if the
\texttt{in\_start\_time} of the first table is after the
\texttt{in\_start\_time} of the second and -1 if it is before.  If the two
\texttt{in\_start\_time}s are identical, the test is repeated on the
\texttt{in\_end\_time}s.  If these are also equal, the comparison returns 0.
\texttt{LALCompareSearchSummaryByOutTime()} operates in a similar manner, but
uses the out, rather than in, times.

\texttt{LALTimeSortSearchSummary()} will time sort a linked list of search
summary tables.  You can sort on in our out start time depending which
\texttt{comparfunc} is specified.

\texttt{LALIfoScanSearchSummary()} steps through a linked list of search
summary tables and returns a pointer \texttt{output} to a linked list of those
tables whos \texttt{ifos} field matches the string \texttt{ifos}.


\texttt{LALIfoScanSummValue()}, \texttt{LALCompareSummValueByTime()} and
\texttt{LALTimeSortSummValue()} performs the same functions as described
above.  The only difference being that they act on summ value tables.

Finally, \texttt{LALCheckOutTimeFromSearchSummary()} verifies that all times
between the specified \texttt{startTime} and \texttt{endTime} have been
searched precisely once for the given \texttt{ifo}.

 
\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent LALGPStoINT8, LALCalloc, LALMalloc, LALFree.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{LIGOMetadataUtilsCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
int 
XLALIFONumber( 
    const char *ifo 
    )
/* </lalVerbatim> */
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

/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void 
XLALReturnIFO( 
    char                *ifo,
    InterferometerNumber IFONumber 
    )
/* </lalVerbatim> */
{
  switch( IFONumber )
  {
    case LAL_IFO_G1:
      LALSnprintf( ifo, LIGOMETA_IFO_MAX, "G1");
      break;

    case LAL_IFO_H1:
      LALSnprintf( ifo, LIGOMETA_IFO_MAX, "H1");
      break;

    case LAL_IFO_H2:
      LALSnprintf( ifo, LIGOMETA_IFO_MAX, "H2");
      break;

    case LAL_IFO_L1:
      LALSnprintf( ifo, LIGOMETA_IFO_MAX, "L1");
      break;

    case LAL_IFO_T1:
      LALSnprintf( ifo, LIGOMETA_IFO_MAX, "T1");
      break;

    case LAL_IFO_V1:
      LALSnprintf( ifo, LIGOMETA_IFO_MAX, "V1");
      break;

    default:
      /* Invalid Detector Site */
      LALSnprintf( ifo, LIGOMETA_IFO_MAX, "");
  }
}


/* <lalVerbatim file="LIGOMetadataUtilsCP"> */
void
XLALReturnDetector(
    LALDetector           *det,
    InterferometerNumber   IFONumber 
    )
/* </lalVerbatim> */
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
  /* JC: HOPEFULLY NEVER GET HERE */
  return LAL_INT8_C(-1);
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

  const SearchSummaryTable *aPtr = *((const SearchSummaryTable * const *)a);
  const SearchSummaryTable *bPtr = *((const SearchSummaryTable * const *)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the in start times */
  memset( &status, 0, sizeof(LALStatus) );
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
    memset( &status, 0, sizeof(LALStatus) );
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

  const SearchSummaryTable *aPtr = *((SearchSummaryTable * const *)a);
  const SearchSummaryTable *bPtr = *((SearchSummaryTable * const *)b);

  INT8 ta = 0;
  INT8 tb = 0;

  /* determine the out start times */
  memset( &status, 0, sizeof(LALStatus) );
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
    memset( &status, 0, sizeof(LALStatus) );
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
    DETATCHSTATUSPTR (status);
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
      LALGPStoINT8( status->statusPtr, &outEndNS, 
          &(thisSearchSumm->out_end_time) );
      CHECKSTATUSPTR( status );

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
      LALGPStoINT8( status->statusPtr, &outEndNS, 
          &(thisSearchSumm->out_end_time) );
      CHECKSTATUSPTR( status );

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

  const SummValueTable *aPtr = *((SummValueTable * const *)a);
  const SummValueTable *bPtr = *((SummValueTable * const *)b);

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
    DETATCHSTATUSPTR (status);
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

