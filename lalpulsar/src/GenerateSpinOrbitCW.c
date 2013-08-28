/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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
#include <lal/GenerateSpinOrbitCW.h>

static LALStatus empty_LALStatus;


/**
 * FIXME: Temporary XLAL-wapper to LAL-function LALGenerateSpinOrbitCW()
 *
 * NOTE: This violates the current version of the XLAL-spec, but is unavoidable at this time,
 * as LALGenerateSpinOrbitCW() hasn't been properly XLALified yet, and doing this would be beyond
 * the scope of this patch.
 * However, doing it here in this way is better than calling LALGenerateSpinOrbitCW() from various
 * far-flung XLAL-functions, as in this way the "violation" is localized in one place, and serves
 * as a reminder for future XLAL-ification at the same time.
 */
int
XLALGenerateSpinOrbitCW ( PulsarCoherentGW *sourceSignal,		///< [out] output signal
                          SpinOrbitCWParamStruc *sourceParams	///< [in] input parameters
                          )
{
  XLAL_CHECK ( sourceSignal != NULL, XLAL_EINVAL );
  XLAL_CHECK ( sourceParams != NULL, XLAL_EINVAL );

  LALStatus status = empty_LALStatus;

  LALGenerateSpinOrbitCW ( &status, sourceSignal, sourceParams );

  XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALGenerateSpinOrbitCW() failed with code=%d, msg='%s'\n", status.statusCode, status.statusDescription );

  return XLAL_SUCCESS;

} // XLALGenerateSpinOrbitCW()



/* First, define a function to compute C(a,b) = (a!)/[(b!)*(a-b)!] */
static UINT4
choose( UINT4 a, UINT4 b );
static UINT4
choose( UINT4 a, UINT4 b )
{
  UINT4 numer = 1;
  UINT4 denom = 1;
  UINT4 myindex = b + 1;
  while ( --myindex ) {
    numer *= a - b + myindex;
    denom *= myindex;
  }
  return numer/denom;
}


/**
 * \author Creighton, T. D.
 *
 * \brief Computes a spindown- and Doppler-modulated continuous waveform.
 *
 * This function computes a quasiperiodic waveform using the spindown and
 * orbital parameters in <tt>*params</tt>, storing the result in
 * <tt>*output</tt>.
 *
 * In the <tt>*params</tt> structure, the routine uses all the "input"
 * fields specified in \ref GenerateSpinOrbitCW_h, and sets all of the
 * "output" fields.  If <tt>params-\>f</tt>=\c NULL, no spindown
 * modulation is performed.  If <tt>params-\>rPeriNorm</tt>=0, no Doppler
 * modulation is performed.
 *
 * In the <tt>*output</tt> structure, the field <tt>output-\>h</tt> is
 * ignored, but all other pointer fields must be set to \c NULL.  The
 * function will create and allocate space for <tt>output-\>a</tt>,
 * <tt>output-\>f</tt>, and <tt>output-\>phi</tt> as necessary.  The
 * <tt>output-\>shift</tt> field will remain set to \c NULL.
 *
 * \heading{Algorithm}
 *
 * This routine calls <tt>LALGenerateCircularSpinOrbitCW()</tt>,
 * <tt>LALGenerateCircularSpinOrbitCW()</tt>,
 * <tt>LALGenerateCircularSpinOrbitCW()</tt>, or
 * <tt>LALGenerateCircularSpinOrbitCW()</tt>, depending on the value of
 * <tt>params-\>oneMinusEcc</tt>.  See the other modules under
 * \ref GenerateSpinOrbitCW_h for descriptions of these routines'
 * algorithms.
 *
 * If <tt>params-\>rPeriNorm</tt>=0, this routine will call
 * <tt>LALGenerateTaylorCW()</tt> to generate the waveform.  It creates a
 * \c TaylorCWParamStruc from the values in <tt>*params</tt>, adjusting
 * the values of <tt>params-\>phi0</tt>, <tt>params-\>f0</tt>, and
 * <tt>params-\>f</tt> from the reference time <tt>params-\>spinEpoch</tt> to
 * the time <tt>params-\>epoch</tt>, as follows: Let \f$\Delta
 * t=t^{(2)}-t^{(1)}\f$ be the time between the old epoch \f$t^{(1)}\f$ and the
 * new one \f$t^{(2)}\f$.  Then the phase, base frequency, and spindown
 * parameters for the new epoch are:
 * \f{eqnarray}{
 * \phi_0^{(2)} & = & \phi_0^{(1)} + 2\pi f_0^{(1)}t \left( 1 +
 * \sum_{k=1}^N \frac{1}{k+1}f_k^{(1)} \Delta t^k \right)
 * \nonumber\\
 * f_0^{(2)} & = & f_0^{(1)} \left( 1 +
 * \sum_{k=1}^N f_k^{(1)} \Delta t^k \right)
 * \nonumber\\
 * f_k^{(2)} & = & \frac{f_0^{(1)}}{f_0^{(2)}} \left( f_k^{(1)} +
 * \sum_{j=k+1}{N} \binom{j}{k} f_j^{(1)}\Delta t^{j-k} \right)
 * \nonumber
 * \f}
 * The phase function \f$\phi(t)=\phi_0^{(i)}+2\pi
 * f_0^{(i)}\left[t-t^{(i)}+\sum_{k=1}^N
 * \frac{f_k^{(i)}}{k+1}\left(t-t^{(i)}\right)^{k+1}\right]\f$ then has the
 * same functional dependence on \f$t\f$ for either \f$i=1\f$ or~2.
 */
void
LALGenerateSpinOrbitCW( LALStatus             *stat,
			PulsarCoherentGW            *output,
			SpinOrbitCWParamStruc *params )
{

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structure exists (output structure will be
     tested by subroutine). */
  ASSERT( params, stat, GENERATESPINORBITCWH_ENUL,
	  GENERATESPINORBITCWH_MSGENUL );

  /* If there is no orbital motion, use LALGenerateTaylorCW() to
     compute the waveform. */
  if ( params->rPeriNorm == 0.0 ) {
    TaylorCWParamStruc taylorParams; /* subroutine parameters */
    REAL8 t;                         /* time shift */
    memset( &taylorParams, 0, sizeof(TaylorCWParamStruc) );
    taylorParams.position = params->position;
    taylorParams.psi = params->psi;
    taylorParams.epoch = params->epoch;
    taylorParams.deltaT = params->deltaT;
    taylorParams.length = params->length;
    taylorParams.aPlus = params->aPlus;
    taylorParams.aCross = params->aCross;
    taylorParams.phi0 = params->phi0;
    taylorParams.f0 = params->f0;
    t = (REAL8)( params->epoch.gpsSeconds -
		 params->spinEpoch.gpsSeconds );
    t += ( 1.0e-9 )*(REAL8)( params->epoch.gpsNanoSeconds -
			     params->spinEpoch.gpsNanoSeconds );

    /* Adjust epochs. */
    if ( params->f ) {
      UINT4 length = params->f->length; /* number of coefficients */
      UINT4 i, j;       /* indecies over coefficients */
      REAL8 tN = 1.0;   /* t raised to various powers */
      REAL8 fFac = 1.0; /* fractional change in frequency */
      REAL8 tFac = 1.0; /* time integral of fFac */
      REAL8 *f1Data = params->f->data; /* pointer to coeficients */
      REAL8 *f2Data;    /* pointer to corrected coeficients */

      TRY( LALDCreateVector( stat->statusPtr, &(taylorParams.f),
			     length ), stat );
      f2Data = taylorParams.f->data;
      memcpy( f2Data, f1Data, length*sizeof(REAL8) );
      for ( i = 0; i < length; i++ ) {
        REAL8 tM = 1.0; /* t raised to various powers */
        fFac += f1Data[i]*( tN *= t );
        tFac += f1Data[i]*tN/( i + 2.0 );
        for ( j = i + 1; j < length; j++ )
          f2Data[i] += choose( j + 1, i + 1 )*f1Data[j]*( tM *= t );
      }
      taylorParams.phi0 += LAL_TWOPI*taylorParams.f0*t*tFac;
      taylorParams.f0 *= fFac;
      for ( i = 0; i < length; i++ )
        f2Data[i] /= fFac;
    } else
      taylorParams.phi0 += LAL_TWOPI*taylorParams.f0*t;

    /* Generate waveform. */
    LALGenerateTaylorCW( stat->statusPtr, output, &taylorParams );
    BEGINFAIL( stat ) {
      if ( taylorParams.f ) {
	TRY( LALDDestroyVector( stat->statusPtr, &(taylorParams.f) ),
	     stat );
      }
    } ENDFAIL( stat );
    if ( taylorParams.f ) {
      TRY( LALDDestroyVector( stat->statusPtr, &(taylorParams.f) ),
	   stat );
    }
    params->dfdt = taylorParams.dfdt;
  }

  /* If there is orbital motion, call the appropriate subroutine. */

  /* e < 0 is out of range. */
  else if ( params->oneMinusEcc > 1.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EECC,
	   GENERATESPINORBITCWH_MSGEECC );
  }

  /* 0 <= e < 1 is elliptic. */
  else if ( params->oneMinusEcc > 0.0 ) {
    TRY( LALGenerateEllipticSpinOrbitCW( stat->statusPtr, output,
					 params ), stat );
  }

  /* e = 1 is parabolic. */
  else if ( params->oneMinusEcc == 0.0 ) {
    TRY( LALGenerateParabolicSpinOrbitCW( stat->statusPtr, output,
					  params ), stat );
  }

  /* e > 1 is hyperbolic. */
  else {
    TRY( LALGenerateHyperbolicSpinOrbitCW( stat->statusPtr, output,
					   params ), stat );
  }

  /* That is all. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
