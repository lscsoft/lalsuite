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
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>

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
  const INT8 play_length = 600000000000LL;
  const INT8 S2_start = 729273613000000000LL;
  const INT8 mod_play = 6370000000000LL;
  INT8 end_mod_play;

  /* handle the unhandled case */
  if ( seg_length >= mod_play )
  {
    LALPrintError( "Segment length is longer than mod_play: FIXME" );
    return -1LL;
  }
  
  end_mod_play = ((seg_end - S2_start) % mod_play );

  /* if no overlap with playground, return zero */
  if ( end_mod_play >= play_length + seg_length )
  {
    return 0LL;
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
    else if ( seg_length < play_length )
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
      LALPrintError( "Error determining playground overlap" );
      return -1LL;
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
