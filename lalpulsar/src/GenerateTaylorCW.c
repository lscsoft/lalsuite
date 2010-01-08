/*
*  Copyright (C) 2007 Teviet Creighton
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

/************************** <lalVerbatim file="GenerateTaylorCWCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GenerateTaylorCW.c}}
\label{ss:GenerateTaylorCW.c}

Computes a Taylor-parametrized continuous waveform.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateTaylorCWCP}
\idx{LALGenerateTaylorCW()}

\subsubsection*{Description}

This function computes a quaiperiodic waveform using the parameters in
\verb@*params@, storing the result in \verb@*output@.

In the \verb@*params@ structure, the routine uses all the ``input''
fields specified in \verb@GenerateTaylorCW.h@, and sets all of the
``output'' fields.  If \verb@params->f=NULL@, a precisely periodic
(monochromatic) waveform is generated.

In the \verb@*output@ structure, the field \verb@output->h@ is
ignored, but all other pointer fields must be set to \verb@NULL@.  The
function will create and allocate space for \verb@output->a@,
\verb@output->f@, and \verb@output->phi@ as necessary.  The
\verb@output->shift@ field will remain set to \verb@NULL@.

\subsubsection*{Algorithm}

This function is a fairly straightforward calculation of
Eqs.~\ref{eq:taylorcw-freq} and~\ref{eq:taylorcw-phi} in
\verb@GenerateTaylorCW.h@.  There are no real tricks involved, except
to note that the phase $\phi$ and the time elapsed $t-t_0$ are
computed and stored as \verb@REAL8@s in order to provide waveforms
that are accurate to small fractions of a cycle over many years.

Since the waveform does not include any effects such as precession,
the amplitudes $A_+$, $A_\times$ and the shift angle $\Phi$, as
defined in \verb@SimulateCoherentGW.h@, are constant.  This is dealt
with by setting \verb@output->a@ to be a
\verb@REAL4TimeVectorSequence@ of two identical vectors
$(A_+,A_\times)$ spanning the requested time of the waveform, under
the assumption that any routine using this output data (such as the
routines in \verb@SimulateCoherentGW.h@) will interpolate these two
endpoints to get the instantaneous values of $A_+$, $A_\times$.  The
field \verb@output->shift@ is left as \verb@NULL@, so the constant
value of $\Phi$ must be subsumed into the polarization angle $\psi$.

The fields \verb@output->f@ and \verb@output->phi@ are created and
filled at the requested sampling interval \verb@params->deltaT@; it is
up to the calling routine to ensure that this sampling interval is
reasonable.  As a guideline, we want to be able to determine the
instantaneous wave phase accurately to within a fraction of a cycle.
For functions that compute the phase by linear interpolation of
\verb@output->phi@, this means sampling on timescales $\Delta
t\lessim\dot{f}^{-1/2}\sim\max\{\sqrt{kf_0f_kT^{k-1}}\}$, where $T$ is
the duration of the waveform.  More precisely, the largest deviation
from linear phase evolution will typically be on the order of
$\Delta\phi\approx(1/2)\ddot{\phi}(\Delta t/2)^2\approx(\pi/4)\Delta
f\Delta t$, where $\Delta f$ is the frequency shift over the timestep.
So if we want our interpolated phase to agree with the true phase to
within, say, $\pi/2$ radians, then we would like to have
$$
\Delta f \Delta t \lessim 2 \;.
$$
This routine provides a check by setting the output parameter field
\verb@params->dfdt@ equal to the maximum value of $\Delta f\Delta t$
encountered during the integration.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                   LALFree()
LALSCreateVectorSequence()    LALSDestroyVectorSequence()
LALSCreateVector()            LALSDestroyVector()
LALDCreateVector()            LALDDestroyVector()
snprintf()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GenerateTaylorCWCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateTaylorCW.h>

NRCSID( GENERATETAYLORCWC, "$Id$" );

/* <lalVerbatim file="GenerateTaylorCWCP"> */
void
LALGenerateTaylorCW( LALStatus          *stat,
		     CoherentGW         *output,
		     TaylorCWParamStruc *params )
{ /* </lalVerbatim> */
  UINT4 n, i;          /* number of and index over samples */
  UINT4 nSpin = 0, j;  /* number of and index over spindown terms */
  REAL8 t, tPow, dt;   /* time, interval, and t raised to a power */
  REAL8 f0, phi0;      /* initial phase and frequency */
  REAL8 twopif0;       /* 2*pi*f0 */
  REAL8 f, fPrev;      /* current and previous values of frequency */
  REAL4 df = 0.0;      /* maximum difference between f and fPrev */
  REAL8 phi;           /* current value of phase */
  REAL8 *fSpin = NULL; /* pointer to Taylor coefficients */
  REAL4 *fData;        /* pointer to frequency data */
  REAL8 *phiData;      /* pointer to phase data */

  INITSTATUS( stat, "LALGenerateTaylorCW", GENERATETAYLORCWC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATETAYLORCWH_ENUL,
	  GENERATETAYLORCWH_MSGENUL );
  ASSERT( output, stat, GENERATETAYLORCWH_ENUL,
	  GENERATETAYLORCWH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->a ), stat, GENERATETAYLORCWH_EOUT,
	  GENERATETAYLORCWH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATETAYLORCWH_EOUT,
	  GENERATETAYLORCWH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATETAYLORCWH_EOUT,
	  GENERATETAYLORCWH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATETAYLORCWH_EOUT,
	  GENERATETAYLORCWH_MSGEOUT );

  /* If Taylor coeficients are specified, make sure they exist. */
  if ( params->f ) {
    ASSERT( params->f->data, stat, GENERATETAYLORCWH_ENUL,
	    GENERATETAYLORCWH_MSGENUL );
    nSpin = params->f->length;
    fSpin = params->f->data;
  }

  /* Set up some other constants, to avoid repeated dereferencing. */
  n = params->length;
  dt = params->deltaT;
  f0 = fPrev = params->f0;
  twopif0 = f0*LAL_TWOPI;
  phi0 = params->phi0;

  /* Allocate output structures. */
  if ( ( output->a = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    ABORT( stat, GENERATETAYLORCWH_EMEM, GENERATETAYLORCWH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( output->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATETAYLORCWH_EMEM, GENERATETAYLORCWH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( output->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATETAYLORCWH_EMEM, GENERATETAYLORCWH_MSGEMEM );
  }
  memset( output->phi, 0, sizeof(REAL8TimeSeries) );

  /* Set output structure metadata fields. */
  output->position = params->position;
  output->psi = params->psi;
  output->a->epoch = output->f->epoch = output->phi->epoch
    = params->epoch;
  output->a->deltaT = n*params->deltaT;
  output->f->deltaT = output->phi->deltaT = params->deltaT;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  snprintf( output->a->name, LALNameLength, "CW amplitudes" );
  snprintf( output->f->name, LALNameLength, "CW frequency" );
  snprintf( output->phi->name, LALNameLength, "CW phase" );

  /* Allocate phase and frequency arrays. */
  LALSCreateVector( stat->statusPtr, &( output->f->data ), n );
  BEGINFAIL( stat ) {
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );
  LALDCreateVector( stat->statusPtr, &( output->phi->data ), n );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	 stat );
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );

  /* Allocate and fill amplitude array. */
  {
    CreateVectorSequenceIn in; /* input to create output->a */
    in.length = 2;
    in.vectorLength = 2;
    LALSCreateVectorSequence( stat->statusPtr, &(output->a->data), &in );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	   stat );
      TRY( LALDDestroyVector( stat->statusPtr, &( output->phi->data ) ),
	   stat );
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
    output->a->data->data[0] = output->a->data->data[2] = params->aPlus;
    output->a->data->data[1] = output->a->data->data[3] = params->aCross;
  }

  /* Fill frequency and phase arrays. */
  fData = output->f->data->data;
  phiData = output->phi->data->data;
  for ( i = 0; i < n; i++ ) {
    t = tPow = i*dt;
    f = 1.0;
    phi = t;
    j = 0;
    while ( j < nSpin ) {
      f += fSpin[j]*tPow;
      phi += fSpin[j]*( tPow*=t )/( j + 2.0 );
      j++;
    }
    f *= f0;
    phi *= twopif0;
    if ( fabs( f - fPrev ) > df )
      df = fabs( f - fPrev );
    *(fData++) = fPrev = f;
    *(phiData++) = phi + phi0;
  }

  /* Set output field and return. */
  params->dfdt = df*dt;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
