/*
 * File Name: IncrementGPS.c
 *
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
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
\idx{LALAddFloatToGPS()}
\idx{XLALAddFloatToGPS()}
\idx{LALDeltaGPS()}
\idx{LALDeltaFloatGPS()}
\idx{LALCompareGPS()}

\subsubsection*{Description}

This module contains a few utility routines to perform comparisons and 
arithmetic on \texttt{LIGOTimeGPS} GPS times. These routines do not convert 
the GPS times they operate on into a floating point representation.

\begin{itemize}
  \item \texttt{LALIncrementGPS()} increments a GPS time by a time interval
  \item \texttt{LALDecrementGPS()} decrements a GPS time by a time interval
  \item \texttt{LALAddFloatToGPS()} adds a REAL8 interval (in seconds) to a GPS time
  \item \texttt{XLALAddFloatToGPS()} adds a REAL8 interval (in seconds) to a GPS time
  \item \texttt{LALDeltaGPS()} returns the difference between two GPS times as a \texttt{LALTimeInterval}.
  \item \texttt{LALDeltaFloatGPS()} returns the difference between two GPS times in seconds as a \texttt{REAL8}.
  \item \texttt{LALCompareGPS()} compares two GPS times, and returns a \texttt{LALGPSCompareResult} indicating if the first GPS time is earlier than, equal to,  or later than the second GPS time
\end{itemize}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

In the \texttt{LALIncrementGPS()}, \texttt{LALDecrementGPS()} and \texttt{LALAddFloatToGPS()}
routines, it is legal to pass a pointer to the same \texttt{LIGOTimeGPS}
structure as input and output, \textit{e.g.}

\begin{verbatim}
  LALStatus status;
  LIGOTimeGPS gps;
  LALTimeInterval interval;

  ...

  gps.gpsSeconds = 10;
  gps.gpsNanoSeconds = 25;
  interval.seconds = 13;
  interval.nanoSeconds = 1000;

  LALIncrementGPS(&status, &gps, &gps, &interval);

\end{verbatim}


\vfill{\footnotesize\input{IncrementGPSCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <math.h>

NRCSID( INCREMENTGPSC, "$Id$" );

/* module-scope constant */
static const INT4 oneBillion = 1000000000;

/* Increment a GPS time */
/* <lalVerbatim file="IncrementGPSCP"> */
void
LALIncrementGPS (LALStatus             *status,
                 LIGOTimeGPS           *pIncrementedGPS, /* output */
                 const LIGOTimeGPS     *pInitialGPS, /* input: GPS time */
                 const LALTimeInterval *pDeltaT) /* input: interval to increment by */
/*  </lalVerbatim> */
{
  LIGOTimeGPS tmp_gps;

  INITSTATUS( status, "LALIncrementGPS", INCREMENTGPSC );
  ATTATCHSTATUSPTR(status);

  tmp_gps = *pInitialGPS;

  tmp_gps.gpsSeconds += pDeltaT->seconds;

  if (tmp_gps.gpsNanoSeconds + pDeltaT->nanoSeconds >= oneBillion)
    {
      ++(tmp_gps.gpsSeconds);
      tmp_gps.gpsNanoSeconds += (pDeltaT->nanoSeconds - oneBillion);
    }
  else
    {
      tmp_gps.gpsNanoSeconds += pDeltaT->nanoSeconds;
    }

   /* assign the computed values */
   *pIncrementedGPS = tmp_gps;

  DETATCHSTATUSPTR(status);
  RETURN( status );
} /* END: LALIncrementGPS() */



/* Decrement a GPS time */
/* <lalVerbatim file="IncrementGPSCP"> */
void
LALDecrementGPS (LALStatus             *status,
                 LIGOTimeGPS           *pDecrementedGPS, /* output */
                 const LIGOTimeGPS     *pInitialGPS, /* input: GPS time */
                 const LALTimeInterval *pDeltaT) /* input: interval to decrement by */
/* </lalVerbatim> */
{
  LIGOTimeGPS         deltaT_gps;
  LIGOTimeGPS         tmp_gps;    /* tmp so that we can use a call like:
                                     LALDecrementGPS(&stat, &gps, &gps, &interval)*/
  LALGPSCompareResult comparison_result;
  
  INITSTATUS( status, "LALDecrementGPS", INCREMENTGPSC );
  ATTATCHSTATUSPTR(status);

  ASSERT(pDecrementedGPS, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pInitialGPS, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT(pDeltaT, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  tmp_gps = *pInitialGPS;

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
      tmp_gps.gpsNanoSeconds = oneBillion
        + pInitialGPS->gpsNanoSeconds - pDeltaT->nanoSeconds;
      tmp_gps.gpsSeconds = pInitialGPS->gpsSeconds - 1
        - pDeltaT->seconds;
    }
  else
    {
      tmp_gps.gpsNanoSeconds = pInitialGPS->gpsNanoSeconds
        - pDeltaT->nanoSeconds;
      tmp_gps.gpsSeconds = pInitialGPS->gpsSeconds
        - pDeltaT->seconds;
    }

  /* assign the computed values */
  *pDecrementedGPS = tmp_gps;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/* Return GPS1 - GPS2 */
/* <lalVerbatim file="IncrementGPSCP"> */
void
LALDeltaGPS (LALStatus         *status,
             LALTimeInterval   *pDeltaGPS, /* output: GPS1 - GPS2 */
             const LIGOTimeGPS *pGPS1, /* input: GPS1 */
             const LIGOTimeGPS *pGPS2) /* input: GPS2 */
/* </lalVerbatim> */
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




/* <lalVerbatim file="IncrementGPSCP"> */
void
LALCompareGPS(LALStatus           *status,
              LALGPSCompareResult *pResult, /* output: -1 => GPS1 < GPS2
                                                        0 => GPS1 = GPS2
                                                        1 => GPS1 > GPS2 */
              const LIGOTimeGPS   *pGPS1, /* input: GPS1 */
              const LIGOTimeGPS   *pGPS2) /* input: GPS2 */
/* </lalVerbatim> */
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



/* Increment a GPS time by a float-interval */
/* <lalVerbatim file="IncrementGPSCP"> */
LIGOTimeGPS *
XLALAddFloatToGPS(
	LIGOTimeGPS *gps,
	REAL8 deltaT
)
/*  </lalVerbatim> */
{
  INT4 secs = deltaT;

  if(!gps)
	  return(NULL);

  gps->gpsSeconds += secs;
  gps->gpsNanoSeconds += floor((deltaT - secs) * oneBillion  + 0.5);

  if (gps->gpsNanoSeconds >= oneBillion)
    {
      gps->gpsNanoSeconds -= oneBillion;
      gps->gpsSeconds++;
    }
  else if (gps->gpsNanoSeconds < 0)
    {
      gps->gpsNanoSeconds += oneBillion;
      gps->gpsSeconds--;
    }

  return(gps);
} /* XLALAddFloatToGPS() */

/* Increment a GPS time by a float-interval */
/* <lalVerbatim file="IncrementGPSCP"> */
void
LALAddFloatToGPS(
    LALStatus          *status,
    LIGOTimeGPS        *outputGPS, /* outputGPS = startGPS + deltaT */
    const LIGOTimeGPS  *startGPS,  /* input: GPS time */
    REAL8 		deltaT     /* input: interval to increment by in sec */
    )
/*  </lalVerbatim> */
{
  INITSTATUS( status, "LALAddFloatToGPS", INCREMENTGPSC );

  *outputGPS = *startGPS;
  XLALAddFloatToGPS(outputGPS, deltaT);

  RETURN( status );

} /* LALAddFloatToGPS() */


/* Return GPS1 - GPS2 as a REAL8 !*/
/* <lalVerbatim file="IncrementGPSCP"> */
void
LALDeltaFloatGPS (LALStatus    *status,
		  REAL8 	*deltaT,	/* tGPS1 - tGPS2 */
		  const LIGOTimeGPS *tGPS1, 	/* input: tGPS1 */
		  const LIGOTimeGPS *tGPS2) 	/* input: tGPS2 */
/* </lalVerbatim> */
{

  INITSTATUS( status, "LALDeltaFloatGPS", INCREMENTGPSC );

  ASSERT(deltaT, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(tGPS1, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT(tGPS2, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  
  *deltaT = (REAL8)(tGPS1->gpsSeconds - tGPS2->gpsSeconds) + 1.0*(tGPS1->gpsNanoSeconds - tGPS2->gpsNanoSeconds)/oneBillion;

  RETURN( status );

} /* LALDeltaFloatGPS() */

