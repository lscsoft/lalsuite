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

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/PulsarSimulateCoherentGW.h>
#include <lal/GenerateTaylorCW.h>

/**
\author Creighton, T. D.

\brief Computes a Taylor-parametrized continuous waveform.

\heading{Description}

This function computes a quaiperiodic waveform using the parameters in
<tt>*params</tt>, storing the result in <tt>*output</tt>.

In the <tt>*params</tt> structure, the routine uses all the "input"
fields specified in ::TaylorCWParamStruc, and sets all of the
"output" fields.  If <tt>params->f=NULL</tt>, a precisely periodic
(monochromatic) waveform is generated.

In the <tt>*output</tt> structure, the field <tt>output->h</tt> is
ignored, but all other pointer fields must be set to \c NULL.  The
function will create and allocate space for <tt>output->a</tt>,
<tt>output->f</tt>, and <tt>output->phi</tt> as necessary.  The
<tt>output->shift</tt> field will remain set to \c NULL.

\heading{Algorithm}

This function is a fairly straightforward calculation of
Eqs.\eqref{eq_taylorcw-freq} and\eqref{eq_taylorcw-phi} in
\ref GenerateTaylorCW_h.  There are no real tricks involved, except
to note that the phase \f$\phi\f$ and the time elapsed \f$t-t_0\f$ are
computed and stored as \c REAL8s in order to provide waveforms
that are accurate to small fractions of a cycle over many years.

Since the waveform does not include any effects such as precession,
the amplitudes \f$A_+\f$, \f$A_\times\f$ and the shift angle \f$\Phi\f$, as
defined in \ref PulsarSimulateCoherentGW_h, are constant.  This is dealt
with by setting <tt>output->a</tt> to be a
\c REAL4TimeVectorSequence of two identical vectors
\f$(A_+,A_\times)\f$ spanning the requested time of the waveform, under
the assumption that any routine using this output data (such as the
routines in \ref PulsarSimulateCoherentGW_h) will interpolate these two
endpoints to get the instantaneous values of \f$A_+\f$, \f$A_\times\f$.  The
field <tt>output->shift</tt> is left as \c NULL, so the constant
value of \f$\Phi\f$ must be subsumed into the polarization angle \f$\psi\f$.

The fields <tt>output->f</tt> and <tt>output->phi</tt> are created and
filled at the requested sampling interval <tt>params->deltaT</tt>; it is
up to the calling routine to ensure that this sampling interval is
reasonable.  As a guideline, we want to be able to determine the
instantaneous wave phase accurately to within a fraction of a cycle.
For functions that compute the phase by linear interpolation of
<tt>output->phi</tt>, this means sampling on timescales \f$\Delta
t\lesssim\dot{f}^{-1/2}\sim\max\{\sqrt{kf_0f_kT^{k-1}}\}\f$, where \f$T\f$ is
the duration of the waveform.  More precisely, the largest deviation
from linear phase evolution will typically be on the order of
\f$\Delta\phi\approx(1/2)\ddot{\phi}(\Delta t/2)^2\approx(\pi/4)\Delta
f\Delta t\f$, where \f$\Delta f\f$ is the frequency shift over the timestep.
So if we want our interpolated phase to agree with the true phase to
within, say, \f$\pi/2\f$ radians, then we would like to have
\f[
\Delta f \Delta t \lesssim 2 \;.
\f]
This routine provides a check by setting the output parameter field
<tt>params->dfdt</tt> equal to the maximum value of \f$\Delta f\Delta t\f$
encountered during the integration.

*/
void
LALGenerateTaylorCW( LALStatus          *stat,
		     PulsarCoherentGW         *output,
		     TaylorCWParamStruc *params )
{
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

  INITSTATUS(stat);
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
