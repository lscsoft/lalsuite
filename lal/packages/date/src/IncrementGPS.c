/*
 * File Name: IncrementGPS.c
 *
 * Author: David Chin <dwchin@umich.edu> +1-734-730-1274
 *
 * Revision: $Id$
 *
 */

#if 0
<lalVerbatim file="IncrementGPSCV">
Author: David Chin
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{IncrementGPS.c}}
\label{ss:IncrementGPS.c}

Routines to perform arithmetic and comparisons on \texttt{LIGOTimeGPS} and
\texttt{LALTimeInterval} data types.

\subsection*{Prototypes}
\vspace{0.1in}
\input{IncrementGPSCP}
\idx{LALIncrementGPS()}
\idx{LALDecrementGPS()}
\idx{LALDeltaGPS()}
\idx{LALCompareGPS()}

\subsubsection*{Description}

This module contains a few utility routines to perform comparisons and 
arithmetic on \texttt{LIGOTimeGPS} GPS times. These routines do not convert 
the GPS times they operate on into a floating point representation.

\begin{itemize}
  \item \texttt{LALIncrementGPS()} increments a GPS time by a time interval
  \item \texttt{LALDecrementGPS()} decrements a GPS time by a time interval
  \item \texttt{LALDeltaGPS()} takes the difference between two GPS times
  \item \texttt{LALCompareGPS()} compares two GPS times, and returns a \texttt{LALGPSCompareResult} indicating if the first GPS time is earlier than, equal to,  or later than the second GPS time
\end{itemize}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{IncrementGPSCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/Date.h>

NRCSID( INCREMENTGPSC, "$Id$" );

/* module-scope constant */
static const INT4 oneBillion = 1000000000;

/* Increment a GPS time */
#pragma <lalVerbatim file="IncrementGPSCP">
void
LALIncrementGPS (LALStatus             *status,
                 LIGOTimeGPS           *pIncrementedGPS, /* output */
                 const LIGOTimeGPS     *pInitialGPS, /* input: GPS time */
                 const LALTimeInterval *pDeltaT) /* input: interval to increment by */
#pragma </lalVerbatim>
{
  INITSTATUS( status, "LALIncrementGPS", INCREMENTGPSC );
  ATTATCHSTATUSPTR(status);

  *pIncrementedGPS = *pInitialGPS;

  pIncrementedGPS->gpsSeconds += pDeltaT->seconds;

  if (pIncrementedGPS->gpsNanoSeconds + pDeltaT->nanoSeconds > oneBillion)
    {
      ++(pIncrementedGPS->gpsSeconds);
      pIncrementedGPS->gpsNanoSeconds += (pDeltaT->nanoSeconds - oneBillion);
    }
  else
    {
      pIncrementedGPS->gpsNanoSeconds += pDeltaT->nanoSeconds;
    }

  DETATCHSTATUSPTR(status);
  RETURN( status );
} /* END: LALIncrementGPS() */



/* Decrement a GPS time */
#pragma <lalVerbatim file="IncrementGPSCP">
void
LALDecrementGPS (LALStatus             *status,
                 LIGOTimeGPS           *pDecrementedGPS, /* output */
                 const LIGOTimeGPS     *pInitialGPS, /* input: GPS time */
                 const LALTimeInterval *pDeltaT) /* input: interval to decrement by */
#pragma </lalVerbatim>
{
  LIGOTimeGPS         deltaT_gps;
  LALGPSCompareResult comparison_result;
  
  INITSTATUS( status, "LALDecrementGPS", INCREMENTGPSC );
  ATTATCHSTATUSPTR(status);

  ASSERT(pDecrementedGPS, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pInitialGPS, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT(pDeltaT, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  *pDecrementedGPS = *pInitialGPS;

  /* Need to try to prevent going into negative GPS times, which
     are undefined; can do this easily by making the DeltaT into a GPS time
     and making sure that it is not later than the InitialGPS */
  deltaT_gps.gpsSeconds     = pDeltaT->seconds;
  deltaT_gps.gpsNanoSeconds = pDeltaT->nanoSeconds;

  TRY( LALCompareGPS(status->statusPtr, &comparison_result, &deltaT_gps,
                     pInitialGPS), status);

  ASSERT (comparison_result != LALGPS_LATER, status,
          DATEH_EDECRTIMETOOLARGE, DATEH_MSGEDECRTIMETOOLARGE);

  if (pInitialGPS->gpsNanoSeconds < pDeltaT->nanoSeconds)
    {
      pDecrementedGPS->gpsNanoSeconds = oneBillion
        + pInitialGPS->gpsNanoSeconds - pDeltaT->nanoSeconds;
      pDecrementedGPS->gpsSeconds = pInitialGPS->gpsSeconds - 1
        - pDeltaT->seconds;
    }
  else
    {
      pDecrementedGPS->gpsNanoSeconds = pInitialGPS->gpsNanoSeconds
        - pDeltaT->nanoSeconds;
      pDecrementedGPS->gpsSeconds = pInitialGPS->gpsSeconds
        - pDeltaT->seconds;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/* Return GPS1 - GPS2 */
#pragma <lalVerbatim file="IncrementGPSCP">
void
LALDeltaGPS (LALStatus         *status,
             LALTimeInterval   *pDeltaGPS, /* output: GPS1 - GPS2 */
             const LIGOTimeGPS *pGPS1, /* input: GPS1 */
             const LIGOTimeGPS *pGPS2) /* input: GPS2 */
#pragma </lalVerbatim>
{
  INITSTATUS( status, "LALDeltaGPS", INCREMENTGPSC );
  ATTATCHSTATUSPTR(status);

  ASSERT(pDeltaGPS, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pGPS1, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT(pGPS2, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /* 6 possibilities:

     1. GPS1.gpsSeconds == GPS2.gpsSeconds

        a. GPS1.gpsNanoSeconds >= GPS2.gpsNanoSeconds

           => GPS1 >= GPS2

        b. GPS1.gpsNanoSeconds < GPS2.gpsNanoSeconds

           => GPS1 < GPS2

     2. GPS1.gpsSeconds != GPS2.gpsSeconds

        a. GPS1.gpsSeconds > GPS2.gpsSeconds

           => GPS1 > GPS2

           i.  GPS1.gpsNanoSeconds >= GPS2.gpsNanoSeconds
           ii. GPS1.gpsNanoSeconds <  GPS2.gpsNanoSeconds

        b. GPS1.gpsSeconds < GPS2.gpsSeconds

           => GPS1 < GPS2

           i.  GPS1.gpsNanoSeconds > GPS2.gpsNanoSeconds
           ii. GPS1.gpsNanoSeconds <=  GPS2.gpsNanoSeconds
  */

  if (pGPS1->gpsSeconds == pGPS2->gpsSeconds) /* 1 */
    {
      if (pGPS1->gpsNanoSeconds >= pGPS2->gpsNanoSeconds) /* a */
        {
          pDeltaGPS->seconds = 0;
          pDeltaGPS->nanoSeconds = pGPS1->gpsNanoSeconds
            - pGPS2->gpsNanoSeconds;
        }
      else /* b */
        {
          pDeltaGPS->seconds = 0;
          pDeltaGPS->nanoSeconds = -(pGPS2->gpsNanoSeconds
                                     - pGPS1->gpsNanoSeconds);
        }
    }
  else /* 2 */
    {
      if (pGPS1->gpsSeconds > pGPS2->gpsSeconds) /* a */
        {
          if (pGPS1->gpsNanoSeconds >= pGPS2->gpsNanoSeconds) /* i */
            {
              pDeltaGPS->seconds = pGPS1->gpsSeconds - pGPS2->gpsSeconds;
              pDeltaGPS->nanoSeconds = pGPS1->gpsNanoSeconds
                - pGPS2->gpsNanoSeconds;
            }
          else /* ii */
            {
              pDeltaGPS->seconds = pGPS1->gpsSeconds - pGPS2->gpsSeconds
                - 1;
              pDeltaGPS->nanoSeconds = pGPS1->gpsNanoSeconds + oneBillion
                - pGPS2->gpsNanoSeconds;
            }
        }
      else /* b */
        {
          if (pGPS1->gpsNanoSeconds > pGPS2->gpsNanoSeconds) /* i */
            {
              pDeltaGPS->seconds = - (pGPS2->gpsSeconds - pGPS1->gpsSeconds
                                      - 1);
              pDeltaGPS->nanoSeconds = - (pGPS2->gpsNanoSeconds + oneBillion
                                          - pGPS1->gpsNanoSeconds);
            }
          else /* ii */
            {
              pDeltaGPS->seconds = - (pGPS2->gpsSeconds - pGPS1->gpsSeconds);
              pDeltaGPS->nanoSeconds = - (pGPS2->gpsNanoSeconds -
                                          pGPS1->gpsNanoSeconds);
            }
        }
    }

  DETATCHSTATUSPTR(status);
  RETURN( status );
} /* END: LALDeltaGPS() */




#pragma <lalVerbatim file="IncrementGPSCP">
void
LALCompareGPS(LALStatus           *status,
              LALGPSCompareResult *pResult, /* output: -1 => GPS1 < GPS2
                                                        0 => GPS1 = GPS2
                                                        1 => GPS1 > GPS2 */
              const LIGOTimeGPS   *pGPS1, /* input: GPS1 */
              const LIGOTimeGPS   *pGPS2) /* input: GPS2 */
#pragma </lalVerbatim>
{
  INITSTATUS( status, "LALCompareGPS", INCREMENTGPSC );
  ATTATCHSTATUSPTR(status);

  ASSERT(pResult, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pGPS1, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT(pGPS2, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*  Cases: (see DeltaGPS above, with some amendments) */

  if (pGPS1->gpsSeconds == pGPS2->gpsSeconds) /* 1 */
    {
      if (pGPS1->gpsNanoSeconds == pGPS2->gpsNanoSeconds)
          *pResult = 0;
      else if (pGPS1->gpsNanoSeconds > pGPS2->gpsNanoSeconds) /* a */
          *pResult = 1;
      else
          *pResult = -1;
    }
  else if (pGPS1->gpsSeconds > pGPS2->gpsSeconds) /* a */
      *pResult = 1;
  else
      *pResult = -1;

  DETATCHSTATUSPTR(status);
  RETURN( status );
} /* END: LALCompareGPS() */
