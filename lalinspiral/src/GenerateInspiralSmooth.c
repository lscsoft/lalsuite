/*
*  Copyright (C) 2007 Stephen Fairhurst
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

/************************** <lalVerbatim file="GenerateInspiralSmoothCV">
Author: Fairhurst, S.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GenerateInspiralSmooth.c}}
\label{ss:GenerateInspiralSmooth.c}

Smooths the end of an inspiral waveform by adding an exponential
ringdown at the end.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateInspiralSmoothCP}
\idx{LALGenerateInspiralSmooth()}

\subsubsection*{Description}

This function creates a smooth ending to an inspiral waveform from
GeneratePPNInspiral.  It works by reading in \verb@**output@ and the
damping factor \verb@*qfactor@.

It is assumed that the \verb@**output@ structure is the output from
\verb@LAL_GeneratePPNInspiral@ so that it contains amplitude, frequency
and phase information in \verb@(*output)->a@, \verb@(*output)->f@ and
\verb@(*output)->phi@ respectively.  These data is then extended by
keeping the frequency fixed and exponentially damping the amplitude with
damping factor \verb@qfactor@.

Note:  The length of the injection stored in \verb@**waveform@ will be
correct.  However, the length and time of the inspiral are not updated
in the \verb@PPNParamStruc@ \verb@params@.  Therefore,
they will still contain the actual end time of the inspiral part of the
waveform.

\subsubsection*{Algorithm}

The function reads in $f_{\mathrm{final}}$ and
$(a_{+,\times})_{\mathrm{final}}$ then it populates additional data
entries by:
\begin{eqnarray}
  f &=& f_{\mathrm{final}} \\
  a_{+,\times} &=& (a_{+,\times})_{\mathrm{final}} \,
    \exp( - \pi \, f_{\mathrm{final}} \, t / \mathrm{qfactor}) \\
  \phi &=& \phi_{\mathrm{final}} + (f_{\mathrm{final}}) \, t \, .
\end{eqnarray}
Here, $t$ is the elapsed time after the end of the inspiral.  The waveform
ends when its amplitude has been suppressed by a factor of $\exp(-10)$.

\subsubsection*{Uses}
\begin{verbatim}
LALRealloc()                  LALFree()
LALSResizeVector()            LALDResizeVector()
LALSDestroyVectorSequence()   LALSDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GenerateInspiralSmoothCV}}

******************************************************* </lalLaTeX> */

/*********************************************************************
 * PREAMBLE                                                          *
 *********************************************************************/

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/FindRoot.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID( GENERATEINSPIRALSMOOTHC, "$Id$" );

#define BUFFSIZE 16384     /* Number of timesteps buffered */

/*********************************************************************
 * MAIN FUNCTION                                                     *
 *********************************************************************/

/* <lalVerbatim file="GenerateInspiralSmoothCP"> */
void
LALGenerateInspiralSmooth(  LALStatus     *stat,
			    CoherentGW    **output,
			    PPNParamStruc *params,
			    REAL4	  *qfactor )
{ /* </lalVerbatim> */

  REAL4 *a, *f; /* pointers to generated amplitude and frequency data */
  REAL8 *phi;   /* pointer to generated phase data */
  REAL4 *p; /* temporary pointer */

  INT4 dataIn;
  INT4 dataLength;
  REAL4 phase;
  REAL4 ampPlus;
  REAL4 ampCross;
  REAL4 finalFreq;
  REAL4 dampFac;
  REAL4 phaseFac;
  REAL4 Qfac;
  INT4 nSmooth,n;
  CoherentGW *waveform;

  INITSTATUS( stat, "LALGenerateInspiralSmooth", GENERATEINSPIRALSMOOTHC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter, qfactor and output structures exist. */
  ASSERT( params, stat, GENERATEPPNINSPIRALH_ENUL,
      GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( qfactor, stat, GENERATEPPNINSPIRALH_ENUL,
      GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( output, stat, GENERATEPPNINSPIRALH_ENUL,
      GENERATEPPNINSPIRALH_MSGENUL );

  waveform = *output;
  /* Make sure there is data for the amplitude, frequency and phase */
  ASSERT( (waveform->a->data->data), stat, GENERATEPPNINSPIRALH_ENUL,
      GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( (waveform->f->data->data), stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( (waveform->phi->data->data), stat, GENERATEPPNINSPIRALH_ENUL,
	  GENERATEPPNINSPIRALH_MSGENUL );


  /* read in the number of points already present */
  dataIn = waveform->a->data->length;

  /* record the final frequency, amplitude and phase */
  ampPlus = waveform->a->data->data[2*dataIn - 2];
  ampCross = waveform->a->data->data[2*dataIn - 1];
  finalFreq = waveform->f->data->data[dataIn - 1];
  phase = waveform->phi->data->data[dataIn - 1];
  Qfac = *qfactor;

  /* We will damp the amplitude according to A = A_{0} exp(-pi f t / Q )
     therefore, at each timestep, we need to multiply the previous
     amplitude by exp(-pi f dt / Q).  Compute this factor once right now */

  phaseFac = LAL_TWOPI * finalFreq * params->deltaT;
  dampFac = exp( - phaseFac / (2 * Qfac));

  /* calculate number of additional points necessary to damp by exp(10). */
  nSmooth = ceil((20 * Qfac)/( phaseFac));
  dataLength = dataIn + nSmooth;

  /* Reallocate the waveform data fields. */
  {
    waveform->a->data->length = dataLength;
    p = (REAL4 *) LALRealloc( ( waveform->a->data->data ),
		    2*dataLength*sizeof(REAL4) );
    if ( !p )
    {
      LALFree( waveform->a );
      waveform->a = NULL;
      LALFree( waveform->f );
      waveform->f = NULL;
      LALFree( waveform->phi );
      waveform->phi = NULL;
      ABORT( stat, GENERATEPPNINSPIRALH_EMEM, GENERATEPPNINSPIRALH_MSGEMEM );
    }
    waveform->a->data->data = p;

    LALSResizeVector( stat->statusPtr, &(waveform->f->data), dataLength );
    BEGINFAIL( stat )
    {
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &(waveform->a->data) ),
	   stat );
      LALFree( waveform->a );
      waveform->a = NULL;
      LALFree( waveform->f );
      waveform->f = NULL;
      LALFree( waveform->phi );
      waveform->phi = NULL;
    }
    ENDFAIL( stat );

    LALDResizeVector( stat->statusPtr, &( waveform->phi->data ), dataLength );
    BEGINFAIL( stat )
    {
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &(waveform->a->data) ),
	   stat );
      TRY( LALSDestroyVector( stat->statusPtr, &(waveform->f->data) ),
	   stat );
      LALFree( waveform->a );
      waveform->a = NULL;
      LALFree( waveform->f );
      waveform->f = NULL;
      LALFree( waveform->phi );
      waveform->phi = NULL;
    }
    ENDFAIL( stat );
  }



  a = &(waveform->a->data->data[2*dataIn]);
  phi = &(waveform->phi->data->data[dataIn]);
  f = &(waveform->f->data->data[dataIn]);

  for (n = 1; n <= nSmooth; n++)
  {
      /* Set frequency equal to final frequency */
      *(f++) = finalFreq;

      /* Compute the amplitude. */
      ampPlus = *(a++) = ampPlus * dampFac;
      ampCross = *(a++) = ampCross * dampFac;

      /* Compute the phase. */
      *(phi++) = phase + n * phaseFac;
  }

  /* Everything's been stored and cleaned up, so there's nothing left
     to do but quit! */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
