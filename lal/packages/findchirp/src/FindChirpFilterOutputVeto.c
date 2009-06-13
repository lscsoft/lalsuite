/*
*  Copyright (C) 2007 Andres C. Rodriguez, Chad Hanna, Drew Keppel, Duncan Brown, Jolien Creighton, Patrick Brady
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
 * File Name: FindChirpFilterOutputVeto.c
 *
 * Author: A.C. Rodriguez
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpFilterOutputVetoCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpFilterOutputVeto.c}}
\label{ss:FindChirpFilterOutputVeto.c}

Memory management functions for creating and destroying input data and
workspace memory for findchirp.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpFilterOutputVetoCP}
\idx{LALFindChirpFilterOutputVeto()}

\subsubsection*{Description}
The function \texttt{LALFindChirpFilterOutputVeto()} implements a signal based
veto, currently it is used primarily for testing. The function itself tests the
consistency of the triggers that survive the bns and macho inspiral search pipeline
by monitoring the behavior of the $r^2$ time series which is calculated for each
segment of data (256 seconds).

\subsubsection*{Thresholds for Searches}
Two thresholds are currently employed in the binary neutron star and primordial
black hole inspiral searches: a signal to noise ratio threshold
($\rho^*$(\emph{t$_j$}), ($^*$ denoting the threshold), and a threshold on the
consistency of the template chirp waveform with the data ($r^{2*}$). At a given instant
in time, \emph{t$_j$}. $r^2$(\emph{t$_j$}) is defined as:
\\
\begin{equation}
{r^2(t_j)} = \frac{\chi^2(t_j)}{p}
\end{equation}
\\
where:
\\
\\
p = number of $\chi^2$ bins
\\
\\
The search code calculates
$\rho$(\emph{t$_j$}) and $r^2$(\emph{t$_j$}) for a given segment of data and looks
for:
\\
\begin{equation}
\rho (t_j) > \rho^*(t_j)
\end{equation}
and
\begin{equation}
r^2(t_j) < r^{2*}(t_j) * (1 + \frac{\rho^2(t_j)  \delta^2}{p})
\end{equation}
\\
where:
\\
$^*$ = threshold used in the search\\
$\rho$ = signal to noise ratio\\
$\delta$ = mismatch between your data and template waveform\\
p = number of $\chi^2$ bins
\\
\\
If both these criteria are met at a given \emph{t$_j$}, an inspiral "trigger" is
recorded.

\subsubsection*{Algorithm}
The algorithm inputs the the vector \texttt{chisqVec} (which is actually $r^2$)
for the whole data segment and searches a time window (\texttt{rsqvetoWindow})
prior to the inferred coalescence time of the trigger up to the trigger time
and counts the number of time samples above a given $r^{2**}$ threshold
(\texttt{rsqvetoThresh}) different than the search pipeline employs. Note as well
that the threshold we impose does not get multiplied by the factor:
(1 + {$\rho^2$(t$_j$)$\delta^2$$/p$). The outputted value from this test is
stored in the \texttt{rsqveto\_duration} field in the \texttt{sngl\_inspiral}
xml table. Future implementation of this function will have it take the
calculated value and decide whether or not to store the trigger
for future analysis.

\subsubsection*{Uses}

\subsubsection*{Notes}
The same test described here could also be employed for monitoring the
behavior of the signal to noise time series, $\rho(\emph{t$_j$})$, about a
trigger, therefore the inclusion of \texttt{qVec} and \texttt{qNorm} as
input to the function for future work.

\vfill{\footnotesize\input{FindChirpFilterOutputVetoCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindChirp.h>
/*#include <lal/FindChirpFilterOutputVeto.h>*/

double rint(double x);

NRCSID (FINDCHIRPFILTEROUTPUTVETOC, "$Id$");

/* <lalVerbatim file="FindChirpFilterOutputVetoCP"> */
void
LALFindChirpFilterOutputVeto(
    LALStatus                          *status,
    SnglInspiralTable                 **eventList,
    FindChirpFilterInput               *input,
    FindChirpFilterParams              *fcParams
    )
/* </lalVerbatim> */
{

  UINT4                 x;                  /* for loops */
  UINT4		        eventIndex;         /* store the event as an index */
  UINT4		        windowIndex;        /* r^2 veto window index */
  UINT4                 timeaboversqThresh; /* time above the veto threshold */
  REAL4                 rsqvetoWindow;      /* the r^2 veto window */
  REAL4		        rsqvetoThresh;      /* the r^2 veto threshold */
  SnglInspiralTable    *event = NULL;
  SnglInspiralTable    *lastevent = NULL;
  FindChirpSegment     *segment;
  REAL4Vector          *chisqVec;
  REAL8                 deltaT;
  COMPLEX8Vector       *qVec;
  FindChirpFilterOutputVetoParams *params;

  event = *eventList;

  INITSTATUS( status, "LALFindChirpFilterFilterOutputVeto",
      FINDCHIRPFILTEROUTPUTVETOC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* store the values needed to compute the veto */
  segment = input->segment;
  chisqVec = fcParams->chisqVec;
  deltaT = fcParams->deltaT;
  qVec = fcParams->qVec;
  params = fcParams->filterOutputVetoParams;

  /* check that the filterOutputVeto parameter structure exist */
  ASSERT( params->rsqvetoWindow, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);

  /* check that the workspace vectors exist */
  ASSERT( chisqVec->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);

  /* check that deltaT is reasonable */
  ASSERT( deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );

  /* make sure that the segment structure exists */
  ASSERT( segment, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the segment structure contains some input */
  ASSERT( segment->data->epoch.gpsSeconds, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* If an event exists, begin the r^2 veto calculations */
  while ( event )
  {
    /* convert the time of the event into an index */
    eventIndex = (UINT4) floor( ( (REAL8) ( event->end_time.gpsSeconds -
            segment->data->epoch.gpsSeconds ) +
          (REAL8) ( event->end_time.gpsNanoSeconds + 1 ) * 1.0e-9 ) / deltaT );

    /* convert the window duration into an index */
    windowIndex =
      (UINT4) floor( (REAL8) ( ( params->rsqvetoWindow ) / deltaT ) );

    /* Initialize output to zero */
    timeaboversqThresh = 0;

    /* Count the number of time samples above the given r^2 veto threshold */
    for( x = eventIndex - windowIndex; x <= eventIndex; ++x )
    {
      if ( chisqVec->data[x] >= params->rsqvetoThresh )
        ++timeaboversqThresh;
    }

    /* Convert the output to seconds and record the computed values in the
       sngl_inspiral event */
    event->rsqveto_duration = (REAL4) timeaboversqThresh * deltaT;

    /* do the rsq veto if only rsqvetoTimeThresh and rsqvetoMaxSNR are given */
    if ( ( params->rsqvetoTimeThresh > 0 ) && ( params->rsqvetoMaxSNR > 0 ) &&
        ( ( params->rsqvetoCoeff < 0 ) || ( params->rsqvetoPow < 0 ) ) )
    {
      if ( ( event->snr < params->rsqvetoMaxSNR ) &&
          ( event->rsqveto_duration > params->rsqvetoTimeThresh ) )
      {
        /* reassign eventList if vetoing first event */
        if ( event == *eventList )
        {
          SnglInspiralTable    *tmpevent = event;
          *eventList = event->next;
          free( tmpevent );
          event = *eventList;
        }
        else
        {
          SnglInspiralTable    *tmpevent = event;
          lastevent->next = event->next;
          free( tmpevent );
          event = lastevent->next;
        }
      }
      else
      {
        lastevent = event;
        event = event->next;
      }
    }
    /* do the rsq veto if all rsqveto parameters are specified */
    else if ( ( params->rsqvetoTimeThresh > 0 ) &&
        ( params->rsqvetoMaxSNR > 0 ) && ( params->rsqvetoCoeff > 0 ) &&
        ( params->rsqvetoPow > 0 ) )
    {
      if ( ( ( event->snr < params->rsqvetoMaxSNR ) &&
            ( event->rsqveto_duration > params->rsqvetoTimeThresh ) ) ||
          ( event->rsqveto_duration >
            ( params->rsqvetoCoeff * pow( event->snr, params->rsqvetoPow ) ) ) )
      {
        /* reassign eventList if vetoing first event */
        if ( event == *eventList )
        {
          SnglInspiralTable *tmpevent = event;
          *eventList = event->next;
          free( tmpevent );
          event = *eventList;
        }
        else
        {
          SnglInspiralTable *tmpevent = event;
          lastevent->next = event->next;
          free( tmpevent );
          event = lastevent->next;
        }
      }
      else
      {
        lastevent = event;
        event = event->next;
      }
    }
    else
    {
      lastevent = event;
      event = event->next;
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
