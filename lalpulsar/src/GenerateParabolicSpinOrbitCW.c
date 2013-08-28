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
#include <lal/GenerateSpinOrbitCW.h>

/**
 * \author Creighton, T. D.
 *
 * \brief Computes a continuous waveform with frequency drift and Doppler
 * modulation from a parabolic orbital trajectory.
 *
 * This function computes a quaiperiodic waveform using the spindown and
 * orbital parameters in <tt>*params</tt>, storing the result in
 * <tt>*output</tt>.
 *
 * In the <tt>*params</tt> structure, the routine uses all the "input"
 * fields specified in \ref GenerateSpinOrbitCW_h, and sets all of the
 * "output" fields.  If <tt>params-\>f</tt>=\c NULL, no spindown
 * modulation is performed.  If <tt>params-\>oneMinusEcc</tt>\f$\neq0\f$, or if
 * <tt>params-\>rPeriNorm</tt>\f$\times\f$<tt>params-\>angularSpeed</tt>\f$\geq1\f$
 * (faster-than-light speed at periapsis), an error is returned.
 *
 * In the <tt>*output</tt> structure, the field <tt>output-\>h</tt> is
 * ignored, but all other pointer fields must be set to \c NULL.  The
 * function will create and allocate space for <tt>output-\>a</tt>,
 * <tt>output-\>f</tt>, and <tt>output-\>phi</tt> as necessary.  The
 * <tt>output-\>shift</tt> field will remain set to \c NULL.
 *
 * ### Algorithm ###
 *
 * For parabolic orbits, we combine Eqs.\eqref{eq_spinorbit-tr},
 * \eqref{eq_spinorbit-t}, and\eqref{eq_spinorbit-upsilon} to get \f$t_r\f$
 * directly as a function of \f$E\f$:
 * \anchor eq_cubic-e \f{equation}{
 * \tag{eq_cubic-e}
 * t_r = t_p + \frac{r_p\sin i}{c} \left[ \cos\omega +
 * \left(\frac{1}{v_p} + \cos\omega\right)E -
 * \frac{\sin\omega}{4}E^2 + \frac{1}{12v_p}E^3\right] \;,
 * \f}
 * where \f$v_p=r_p\dot{\upsilon}_p\sin i/c\f$ is a normalized velocity at
 * periapsis.  Following the prescription for the general analytic
 * solution to the real cubic equation, we substitute
 * \f$E=x+3v_p\sin\omega\f$ to obtain:
 * \anchor eq_cubic-x \f{equation}{
 * \tag{eq_cubic-x}
 * x^3 + px = q \;,
 * \f}
 * where:
 * \anchor eq_cubic-p \anchor eq_cubic-q \f{eqnarray}{
 * \tag{eq_cubic-p}
 * p & = & 12 + 12v_p\cos\omega - 3v_p^2\sin^2\omega \;, \\
 * \tag{eq_cubic-q}
 * q & = & 12v_p^2\sin\omega\cos\omega - 24v_p\sin\omega +
 * 2v_p^3\sin^3\omega + 12\dot{\upsilon}_p(t_r-t_p) \;.
 * \f}
 * We note that \f$p>0\f$ is guaranteed as long as \f$v_p<1\f$, so the right-hand
 * side of Eq.\eqref{eq_cubic-x} is monotonic in \f$x\f$ and has exactly one
 * root.  However, \f$p\rightarrow0\f$ in the limit \f$v_p\rightarrow1\f$ and
 * \f$\omega=\pi\f$.  This may cause some loss of precision in subsequent
 * calculations.  But \f$v_p\sim1\f$ means that our solution will be
 * inaccurate anyway because we ignore significant relativistic effects.
 *
 * Since \f$p>0\f$, we can substitute \f$x=y\sqrt{3/4p}\f$ to obtain:
 * \f{equation}{
 * 4y^3 + 3y = \frac{q}{2}\left(\frac{3}{p}\right)^{3/2} \equiv C \;.
 * \f}
 * Using the triple-angle hyperbolic identity
 * \f$\sinh(3\theta)=4\sinh^3\theta+3\sinh\theta\f$, we have
 * \f$y=\sinh\left(\frac{1}{3}\sinh^{-1}C\right)\f$.  The solution to the
 * original cubic equation is then:
 * \f{equation}{
 * E = 3v_p\sin\omega + 2\sqrt{\frac{p}{3}}
 * \sinh\left(\frac{1}{3}\sinh^{-1}C\right) \;.
 * \f}
 * To ease the calculation of \f$E\f$, we precompute the constant part
 * \f$E_0=3v_p\sin\omega\f$ and the coefficient \f$\Delta E=2\sqrt{p/3}\f$.
 * Similarly for \f$C\f$, we precompute a constant piece \f$C_0\f$ evaluated at
 * the epoch of the output time series, and a stepsize coefficient
 * \f$\Delta C=6(p/3)^{3/2}\dot{\upsilon}_p\Delta t\f$, where \f$\Delta t\f$ is
 * the step size in the (output) time series in \f$t_r\f$.  Thus at any
 * timestep \f$i\f$, we obtain \f$C\f$ and hence \f$E\f$ via:
 * \f{eqnarray}{
 * C & = & C_0 + i\Delta C \;, \nonumber\\
 * E & = & E_0 + \Delta E\times\left\{\begin{array}{l@{\qquad}c}
 * \sinh\left[\frac{1}{3}\ln\left(
 * C + \sqrt{C^2+1} \right) \right]\;, & C\geq0 \;,\\
 * \\
 * \sinh\left[-\frac{1}{3}\ln\left(
 * -C + \sqrt{C^2+1} \right) \right]\;, & C\leq0 \;,\\
 * \end{array}\right. \nonumber
 * \f}
 * where we have explicitly written \f$\sinh^{-1}\f$ in terms of functions in
 * \c math.h.  Once \f$E\f$ is found, we can compute
 * \f$t=E(12+E^2)/(12\dot{\upsilon}_p)\f$ (where again \f$1/12\dot{\upsilon}_p\f$
 * can be precomputed), and hence \f$f\f$ and \f$\phi\f$ via
 * Eqs.\eqref{eq_taylorcw-freq} and\eqref{eq_taylorcw-phi}.  The
 * frequency \f$f\f$ must then be divided by the Doppler factor:
 * \f[
 * 1 + \frac{\dot{R}}{c} = 1 + \frac{v_p}{4+E^2}\left(
 * 4\cos\omega - 2E\sin\omega \right)
 * \f]
 * (where once again \f$4\cos\omega\f$ and \f$2\sin\omega\f$ can be precomputed).
 *
 * This routine does not account for relativistic timing variations, and
 * issues warnings or errors based on the criterea of
 * Eq.\eqref{eq_relativistic-orbit} in GenerateEllipticSpinOrbitCW().
 * The routine will also warn if
 * it seems likely that \c REAL8 precision may not be sufficient to
 * track the orbit accurately.  We estimate that numerical errors could
 * cause the number of computed wave cycles to vary by
 * \f[
 * \Delta N \lesssim f_0 T\epsilon\left[
 * \sim6+\ln\left(|C|+\sqrt{|C|^2+1}\right)\right] \;,
 * \f]
 * where \f$|C|\f$ is the maximum magnitude of the variable \f$C\f$ over the
 * course of the computation, \f$f_0T\f$ is the approximate total number of
 * wave cycles over the computation, and \f$\epsilon\approx2\times10^{-16}\f$
 * is the fractional precision of \c REAL8 arithmetic.  If this
 * estimate exceeds 0.01 cycles, a warning is issued.
 */
void
LALGenerateParabolicSpinOrbitCW( LALStatus             *stat,
				 PulsarCoherentGW            *output,
				 SpinOrbitCWParamStruc *params )
{
  UINT4 n, i;              /* number of and index over samples */
  UINT4 nSpin = 0, j;      /* number of and index over spindown terms */
  REAL8 t, dt, tPow;       /* time, interval, and t raised to a power */
  REAL8 phi0, f0, twopif0; /* initial phase, frequency, and 2*pi*f0 */
  REAL8 f, fPrev;      /* current and previous values of frequency */
  REAL4 df = 0.0;      /* maximum difference between f and fPrev */
  REAL8 phi;           /* current value of phase */
  REAL8 vp;            /* projected speed at periapsis */
  REAL8 argument;      /* argument of periapsis */
  REAL8 fourCosOmega;  /* four times the cosine of argument */
  REAL8 twoSinOmega;   /* two times the sine of argument */
  REAL8 vpCosOmega;    /* vp times cosine of argument */
  REAL8 vpSinOmega;    /* vp times sine of argument */
  REAL8 vpSinOmega2;   /* vpSinOmega squared */
  REAL8 vDot6;         /* 6 times angular speed at periapsis */
  REAL8 oneBy12vDot;   /* one over (12 times angular speed) */
  REAL8 pBy3;          /* constant sqrt(p/3) in cubic equation */
  REAL8 p32;           /* constant (p/3)^1.5 in cubic equation */
  REAL8 c, c0, dc;     /* C variable, offset, and step increment */
  REAL8 e, e2, e0;     /* E variable, E^2, and constant piece of E */
  REAL8 de;            /* coefficient of sinh() piece of E */
  REAL8 tpOff;         /* orbit epoch - time series epoch (s) */
  REAL8 spinOff;       /* spin epoch - orbit epoch (s) */
  REAL8 *fSpin = NULL; /* pointer to Taylor coefficients */
  REAL4 *fData;        /* pointer to frequency data */
  REAL8 *phiData;      /* pointer to phase data */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATESPINORBITCWH_ENUL,
	  GENERATESPINORBITCWH_MSGENUL );
  ASSERT( output, stat, GENERATESPINORBITCWH_ENUL,
	  GENERATESPINORBITCWH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->a ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );

  /* If Taylor coeficients are specified, make sure they exist. */
  if ( params->f ) {
    ASSERT( params->f->data, stat, GENERATESPINORBITCWH_ENUL,
	    GENERATESPINORBITCWH_MSGENUL );
    nSpin = params->f->length;
    fSpin = params->f->data;
  }

  /* Set up some constants (to avoid repeated calculation or
     dereferencing), and make sure they have acceptable values. */
  vp = params->rPeriNorm*params->angularSpeed;
  vDot6 = 6.0*params->angularSpeed;
  n = params->length;
  dt = params->deltaT;
  f0 = fPrev = params->f0;
  if ( params->oneMinusEcc != 0.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EECC,
	   GENERATESPINORBITCWH_MSGEECC );
  }
  if ( vp >= 1.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EFTL,
	   GENERATESPINORBITCWH_MSGEFTL );
  }
  if ( vp <= 0.0 || dt <= 0.0 || f0 <= 0.0 || vDot6 <= 0.0 ||
       n == 0 ) {
    ABORT( stat, GENERATESPINORBITCWH_ESGN,
	   GENERATESPINORBITCWH_MSGESGN );
  }
#ifndef NDEBUG
  if ( lalDebugLevel & LALWARNING ) {
    if ( f0*n*dt*vp*vp > 0.5 )
      LALWarning( stat, "Orbit may have significant relativistic"
		  " effects that are not included" );
  }
#endif

  /* Compute offset between time series epoch and periapsis, and
     betweem periapsis and spindown reference epoch. */
  tpOff = (REAL8)( params->orbitEpoch.gpsSeconds -
		   params->epoch.gpsSeconds );
  tpOff += 1.0e-9 * (REAL8)( params->orbitEpoch.gpsNanoSeconds -
			     params->epoch.gpsNanoSeconds );
  spinOff = (REAL8)( params->orbitEpoch.gpsSeconds -
		     params->spinEpoch.gpsSeconds );
  spinOff += 1.0e-9 * (REAL8)( params->orbitEpoch.gpsNanoSeconds -
			       params->spinEpoch.gpsNanoSeconds );

  /* Set up some other constants. */
  twopif0 = f0*LAL_TWOPI;
  phi0 = params->phi0;
  argument = params->omega;
  oneBy12vDot = 0.5/vDot6;
  fourCosOmega = 4.0*cos( argument );
  twoSinOmega = 2.0*sin( argument );
  vpCosOmega = 0.25*vp*fourCosOmega;
  vpSinOmega = 0.5*vp*twoSinOmega;
  vpSinOmega2 = vpSinOmega*vpSinOmega;
  pBy3 = sqrt( 4.0*( 1.0 + vpCosOmega ) - vpSinOmega2 );
  p32 = 1.0/( pBy3*pBy3*pBy3 );
  c0 = p32*( vpSinOmega*( 6.0*vpCosOmega - 12.0 + vpSinOmega2 ) -
	     tpOff*vDot6 );
  dc = p32*vDot6*dt;
  e0 = 3.0*vpSinOmega;
  de = 2.0*pBy3;

  /* Check whether REAL8 precision is good enough. */
#ifndef NDEBUG
  if ( lalDebugLevel & LALWARNING ) {
    REAL8 x = fabs( c0 + n*dc ); /* a temporary computation variable */
    if ( x < fabs( c0 ) )
      x = fabs( c0 );
    x = 6.0 + log( x + sqrt( x*x + 1.0 ) );
    if ( LAL_REAL8_EPS*f0*dt*n*x > 0.01 )
      LALWarning( stat, "REAL8 arithmetic may not have sufficient"
		  " precision for this orbit" );
  }
#endif

  /* Allocate output structures. */
  if ( ( output->a = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    ABORT( stat, GENERATESPINORBITCWH_EMEM,
	   GENERATESPINORBITCWH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( output->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATESPINORBITCWH_EMEM,
	   GENERATESPINORBITCWH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( output->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATESPINORBITCWH_EMEM,
	   GENERATESPINORBITCWH_MSGEMEM );
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

    /* Compute emission time. */
    c = c0 + dc*i;
    if ( c > 0 )
      e = e0 + de*sinh( log( c + sqrt( c*c + 1.0 ) )/3.0 );
    else
      e = e0 + de*sinh( -log( -c + sqrt( c*c + 1.0 ) )/3.0 );
    e2 = e*e;
    phi = t = tPow = oneBy12vDot*e*( 12.0 + e2 );

    /* Compute source emission phase and frequency. */
    f = 1.0;
    for ( j = 0; j < nSpin; j++ ) {
      f += fSpin[j]*tPow;
      phi += fSpin[j]*( tPow*=t )/( j + 2.0 );
    }

    /* Appy frequency Doppler shift. */
    f *= f0 / ( 1.0 + vp*( fourCosOmega - e*twoSinOmega )
		/( 4.0 + e2 ) );
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
