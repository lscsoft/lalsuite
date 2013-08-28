/*
*  Copyright (C) 2007 Reinhard Prix, Teviet Creighton
*  Copyright (C) 2012 Teviet Creighton, Evan Goetz
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
#include <gsl/gsl_roots.h>


static REAL8 gsl_E_solver(REAL8 e, void *p);

struct E_solver_params {
   REAL8 a, b, x;
};

static REAL8 gsl_E_solver(REAL8 e, void *p) {
   struct E_solver_params *params = (struct E_solver_params*)p;
   return e + params->a*sin(e) + params->b*(cos(e) - 1.0) - params->x;
}

/**
 * \author Creighton, T. D.
 *
 * \brief Computes a continuous waveform with frequency drift and Doppler
 * modulation from an elliptical orbital trajectory.
 *
 * This function computes a quaiperiodic waveform using the spindown and
 * orbital parameters in <tt>*params</tt>, storing the result in
 * <tt>*output</tt>.
 *
 * In the <tt>*params</tt> structure, the routine uses all the "input"
 * fields specified in \ref GenerateSpinOrbitCW_h, and sets all of the
 * "output" fields.  If <tt>params-\>f</tt>=\c NULL, no spindown
 * modulation is performed.  If <tt>params-\>oneMinusEcc</tt>\f$\notin(0,1]\f$
 * (an open orbit), or if
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
 * For elliptical orbits, we combine Eqs.\eqref{eq_spinorbit-tr},
 * \eqref{eq_spinorbit-t}, and\eqref{eq_spinorbit-upsilon} to get \f$t_r\f$
 * directly as a function of the eccentric anomaly \f$E\f$:
 * \anchor eq_tr-e1 \f{eqnarray}{
 * t_r = t_p & + & \left(\frac{r_p \sin i}{c}\right)\sin\omega \nonumber\\
 * \tag{eq_tr-e1}
 * & + & \left(\frac{P}{2\pi}\right) \left( E +
 * \left[v_p(1-e)\cos\omega - e\right]\sin E
 * + \left[v_p\sqrt{\frac{1-e}{1+e}}\sin\omega\right]
 * [\cos E - 1]\right) \;,
 * \f}
 * where \f$v_p=r_p\dot{\upsilon}_p\sin i/c\f$ is a normalized velocity at
 * periapsis and \f$P=2\pi\sqrt{(1+e)/(1-e)^3}/\dot{\upsilon}_p\f$ is the
 * period of the orbit.  For simplicity we write this as:
 * \anchor eq_tr-e2 \f{equation}{
 * \tag{eq_tr-e2}
 * t_r = T_p + \frac{1}{n}\left( E + A\sin E + B[\cos E - 1] \right) \;,
 * \f}
 *
 * \image html  inject_eanomaly.png "Fig. [fig_binary_orbit_ell]: Function to be inverted to find eccentric anomaly"
 * \image latex inject_eanomaly.pdf "Function to be inverted to find eccentric anomaly" width=0.23\textwidth
 *
 * where \f$T_p\f$ is the \e observed time of periapsis passage and
 * \f$n=2\pi/P\f$ is the mean angular speed around the orbit.  Thus the key
 * numerical procedure in this routine is to invert the expression
 * \f$x=E+A\sin E+B(\cos E - 1)\f$ to get \f$E(x)\f$.  We note that
 * \f$E(x+2n\pi)=E(x)+2n\pi\f$, so we only need to solve this expression in
 * the interval \f$[0,2\pi)\f$, sketched to the right.
 *
 * We further note that \f$A^2+B^2<1\f$, although it approaches 1 when
 * \f$e\rightarrow1\f$, or when \f$v_p\rightarrow1\f$ and either \f$e=0\f$ or
 * \f$\omega=\pi\f$.  Except in this limit, Newton-Raphson methods will
 * converge rapidly for any initial guess.  In this limit, though, the
 * slope \f$dx/dE\f$ approaches zero at the point of inflection, and an
 * initial guess or iteration landing near this point will send the next
 * iteration off to unacceptably large or small values.  However, by
 * restricting all initial guesses and iterations to the domain
 * \f$E\in[0,2\pi)\f$, one will always end up on a trajectory branch that
 * will converge uniformly.  This should converge faster than the more
 * generically robust technique of bisection. (Note: the danger with Newton's method
 * has been found to be unstable for certain binary orbital parameters. So if
 * Newton's method fails to converge, a bisection algorithm is employed.)
 *
 * In this algorithm, we start the computation with an arbitrary initial
 * guess of \f$E=0\f$, and refine it until the we get agreement to within
 * 0.01 parts in part in \f$N_\mathrm{cyc}\f$ (where \f$N_\mathrm{cyc}\f$ is the
 * larger of the number of wave cycles in an orbital period, or the
 * number of wave cycles in the entire waveform being generated), or one
 * part in \f$10^{15}\f$ (an order of magnitude off the best precision
 * possible with \c REAL8 numbers).  The latter case indicates that
 * \c REAL8 precision may fail to give accurate phasing, and one
 * should consider modeling the orbit as a set of Taylor frequency
 * coefficients \'{a} la <tt>LALGenerateTaylorCW()</tt>.  On subsequent
 * timesteps, we use the previous timestep as an initial guess, which is
 * good so long as the timesteps are much smaller than an orbital period.
 * This sequence of guesses will have to readjust itself once every orbit
 * (as \f$E\f$ jumps from \f$2\pi\f$ down to 0), but this is relatively
 * infrequent; we don't bother trying to smooth this out because the
 * additional tests would probably slow down the algorithm overall.
 *
 * Once a value of \f$E\f$ is found for a given timestep in the output
 * series, we compute the system time \f$t\f$ via Eq.\eqref{eq_spinorbit-t},
 * and use it to determine the wave phase and (non-Doppler-shifted)
 * frequency via Eqs.\eqref{eq_taylorcw-freq}
 * and\eqref{eq_taylorcw-phi}.  The Doppler shift on the frequency is
 * then computed using Eqs.\eqref{eq_spinorbit-upsilon}
 * and\eqref{eq_orbit-rdot}.  We use \f$\upsilon\f$ as an intermediate in
 * the Doppler shift calculations, since expressing \f$\dot{R}\f$ directly in
 * terms of \f$E\f$ results in expression of the form \f$(1-e)/(1-e\cos E)\f$,
 * which are difficult to simplify and face precision losses when
 * \f$E\sim0\f$ and \f$e\rightarrow1\f$.  By contrast, solving for \f$\upsilon\f$ is
 * numerically stable provided that the system <tt>atan2()</tt> function is
 * well-designed.
 *
 * The routine does not account for variations in special relativistic or
 * gravitational time dilation due to the elliptical orbit, nor does it
 * deal with other gravitational effects such as Shapiro delay.  To a
 * very rough approximation, the amount of phase error induced by
 * gravitational redshift goes something like \f$\Delta\phi\sim
 * fT(v/c)^2\Delta(r_p/r)\f$, where \f$f\f$ is the typical wave frequency, \f$T\f$
 * is either the length of data or the orbital period (whichever is
 * \e smaller), \f$v\f$ is the \e true (unprojected) speed at
 * periapsis, and \f$\Delta(r_p/r)\f$ is the total range swept out by the
 * quantity \f$r_p/r\f$ over the course of the observation.  Other
 * relativistic effects such as special relativistic time dilation are
 * comparable in magnitude.  We make a crude estimate of when this is
 * significant by noting that \f$v/c\gtrsim v_p\f$ but
 * \f$\Delta(r_p/r)\lesssim 2e/(1+e)\f$; we take these approximations as
 * equalities and require that \f$\Delta\phi\lesssim\pi\f$, giving:
 * \anchor eq_relativistic-orbit \f{equation}{
 * \tag{eq_relativistic-orbit}
 * f_0Tv_p^2\frac{4e}{1+e}\lesssim1 \;.
 * \f}
 * When this critereon is violated, a warning is generated.  Furthermore,
 * as noted earlier, when \f$v_p\geq1\f$ the routine will return an error, as
 * faster-than-light speeds can cause the emission and reception times to
 * be non-monotonic functions of one another.
 */
void
LALGenerateEllipticSpinOrbitCW( LALStatus             *stat,
				PulsarCoherentGW            *output,
				SpinOrbitCWParamStruc *params )
{
  UINT4 n, i;              /* number of and index over samples */
  UINT4 nSpin = 0, j;      /* number of and index over spindown terms */
  REAL8 t, dt, tPow;       /* time, interval, and t raised to a power */
  REAL8 phi0, f0, twopif0; /* initial phase, frequency, and 2*pi*f0 */
  REAL8 f, fPrev;          /* current and previous values of frequency */
  REAL4 df = 0.0;          /* maximum difference between f and fPrev */
  REAL8 phi;               /* current value of phase */
  REAL8 p, vDotAvg;        /* orbital period, and 2*pi/(period) */
  REAL8 vp;                /* projected speed at periapsis */
  REAL8 upsilon, argument; /* true anomaly, and argument of periapsis */
  REAL8 eCosOmega;         /* eccentricity * cosine of argument */
  REAL8 tPeriObs;          /* time of observed periapsis */
  REAL8 spinOff;           /* spin epoch - orbit epoch */
  REAL8 x;                 /* observed mean anomaly */
  REAL8 dx, dxMax;         /* current and target errors in x */
  REAL8 a, b;              /* constants in equation for x */
  REAL8 ecc;               /* orbital eccentricity */
  REAL8 oneMinusEcc, onePlusEcc; /* 1 - ecc and 1 + ecc */
  REAL8 e = 0.0;                 /* eccentric anomaly */
  REAL8 de = 0.0;                /* eccentric anomaly step */
  REAL8 sine = 0.0, cose = 0.0;  /* sine of e, and cosine of e minus 1 */
  REAL8 *fSpin = NULL;           /* pointer to Taylor coefficients */
  REAL4 *fData;                  /* pointer to frequency data */
  REAL8 *phiData;                /* pointer to phase data */

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
  oneMinusEcc = params->oneMinusEcc;
  ecc = 1.0 - oneMinusEcc;
  onePlusEcc = 1.0 + ecc;
  if ( ecc < 0.0 || oneMinusEcc <= 0.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EECC,
	   GENERATESPINORBITCWH_MSGEECC );
  }
  vp = params->rPeriNorm*params->angularSpeed;
  n = params->length;
  dt = params->deltaT;
  f0 = fPrev = params->f0;
  vDotAvg = params->angularSpeed
    *sqrt( oneMinusEcc*oneMinusEcc*oneMinusEcc/onePlusEcc );
  if ( vp >= 1.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EFTL,
	   GENERATESPINORBITCWH_MSGEFTL );
  }
  if ( vp <= 0.0 || dt <= 0.0 || f0 <= 0.0 || vDotAvg <= 0.0 ||
       n == 0 ) {
    ABORT( stat, GENERATESPINORBITCWH_ESGN,
	   GENERATESPINORBITCWH_MSGESGN );
  }

  /* Set up some other constants. */
  twopif0 = f0*LAL_TWOPI;
  phi0 = params->phi0;
  argument = params->omega;
  p = LAL_TWOPI/vDotAvg;
  a = vp*oneMinusEcc*cos( argument ) + oneMinusEcc - 1.0;
  b = vp*sqrt( oneMinusEcc/( onePlusEcc ) )*sin( argument );
  eCosOmega = ecc*cos( argument );
  if ( n*dt > p )
    dxMax = 0.01/( f0*n*dt );
  else
    dxMax = 0.01/( f0*p );
  if ( dxMax < 1.0e-15 ) {
    dxMax = 1.0e-15;
#ifndef NDEBUG
    LALWarning( stat, "REAL8 arithmetic may not have sufficient"
		" precision for this orbit" );
#endif
  }
#ifndef NDEBUG
  if ( lalDebugLevel & LALWARNING ) {
    REAL8 tau = n*dt;
    if ( tau > p )
      tau = p;
    if ( f0*tau*vp*vp*ecc/onePlusEcc > 0.25 )
      LALWarning( stat, "Orbit may have significant relativistic"
		  " effects that are not included" );
  }
#endif

  /* Compute offset between time series epoch and observed periapsis,
     and betweem true periapsis and spindown reference epoch. */
  tPeriObs = (REAL8)( params->orbitEpoch.gpsSeconds -
		      params->epoch.gpsSeconds );
  tPeriObs += 1.0e-9 * (REAL8)( params->orbitEpoch.gpsNanoSeconds -
				params->epoch.gpsNanoSeconds );
  tPeriObs += params->rPeriNorm*sin( params->omega );
  spinOff = (REAL8)( params->orbitEpoch.gpsSeconds -
		     params->spinEpoch.gpsSeconds );
  spinOff += 1.0e-9 * (REAL8)( params->orbitEpoch.gpsNanoSeconds -
			       params->spinEpoch.gpsNanoSeconds );

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
    INT4 nOrb; /* number of orbits since the specified orbit epoch */

    /* First, find x in the range [0,2*pi]. */
    x = vDotAvg*( i*dt - tPeriObs );
    nOrb = (INT4)( x/LAL_TWOPI );
    if ( x < 0.0 )
      nOrb -= 1;
    x -= LAL_TWOPI*nOrb;

    /* Newton-Raphson iteration to find E(x). Maximum of 100 iterations. */
    INT4 maxiter = 100, iter = 0;
    while ( iter<maxiter && fabs( dx = e + a*sine + b*cose - x ) > dxMax ) {
      iter++;
      //Make a check on the step-size so we don't step too far
      de = dx/( 1.0 + a*cose + a - b*sine );
      if ( de > LAL_PI )
        de = LAL_PI;
      else if ( de < -LAL_PI )
        de = -LAL_PI;
      e -= de;

      if ( e < 0.0 )
        e = 0.0;
      else if ( e > LAL_TWOPI )
        e = LAL_TWOPI;
      sine = sin( e );
      cose = cos( e ) - 1.0;
    }
    /* Bisection algorithm from GSL if Newton's method (above) fails to converge. */
    if (iter==maxiter && fabs( dx = e + a*sine + b*cose - x ) > dxMax ) {
       //Initialize solver
       const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
       gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
       REAL8 e_lo = 0.0, e_hi = LAL_TWOPI;
       gsl_function F;
       struct E_solver_params pars = {a, b, x};
       F.function = &gsl_E_solver;
       F.params = &pars;

       if (gsl_root_fsolver_set(s, &F, e_lo, e_hi) != 0) {
          LALFree( output->a );   output->a = NULL;
          LALFree( output->f );   output->f = NULL;
          LALFree( output->phi ); output->phi = NULL;
          ABORT( stat, -1, "GSL failed to set initial points" );
       }

       INT4 keepgoing = 1;
       INT4 success = 0;
       INT4 root_status = keepgoing;
       e = 0.0;
       iter = 0;
       while (root_status==keepgoing && iter<maxiter) {
          iter++;
          root_status = gsl_root_fsolver_iterate(s);
          if (root_status!=keepgoing && root_status!=success) {
             LALFree( output->a );   output->a = NULL;
             LALFree( output->f );   output->f = NULL;
             LALFree( output->phi ); output->phi = NULL;
             ABORT( stat, -1, "gsl_root_fsolver_iterate() failed" );
          }
          e = gsl_root_fsolver_root(s);
          sine = sin(e);
          cose = cos(e) - 1.0;
          if (fabs( dx = e + a*sine + b*cose - x ) > dxMax) root_status = keepgoing;
          else root_status = success;
       }

       if (root_status!=success) {
          LALFree( output->a );   output->a = NULL;
          LALFree( output->f );   output->f = NULL;
          LALFree( output->phi ); output->phi = NULL;
          gsl_root_fsolver_free(s);
          ABORT( stat, -1, "Could not converge using bisection algorithm" );
       }

       gsl_root_fsolver_free(s);
    }

    /* Compute source emission time, phase, and frequency. */
    phi = t = tPow =
      ( e + LAL_TWOPI*nOrb - ecc*sine )/vDotAvg + spinOff;
    f = 1.0;
    for ( j = 0; j < nSpin; j++ ) {
      f += fSpin[j]*tPow;
      phi += fSpin[j]*( tPow*=t )/( j + 2.0 );
    }

    /* Appy frequency Doppler shift. */
    upsilon = 2.0*atan2( sqrt( -onePlusEcc*cose ), sqrt( oneMinusEcc*( cose + 2.0 ) ) );
    f *= f0 / ( 1.0 + vp*( cos( argument + upsilon ) + eCosOmega ) /onePlusEcc );
    phi *= twopif0;
    if ( (i > 0) && (fabs( f - fPrev ) > df) )
      df = fabs( f - fPrev );
    *(fData++) = fPrev = f;
    *(phiData++) = phi + phi0;

  } /* for i < n */

  /* Set output field and return. */
  params->dfdt = df*dt;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
