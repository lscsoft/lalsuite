/*
 *  Copyright (C) 2007, 2009 Chris Messenger
 *  Copyright (C) 2006 John T. Whelan, Badri Krishnan
 *  Copyright (C) 2005, 2006, 2007, 2010, 2014 Reinhard Prix
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

/*---------- INCLUDES ----------*/
#include <math.h>

#include <lal/SSBtimes.h>
#include <lal/AVFactories.h>

/* GSL includes */
#include <lal/LALGSL.h>
#include <gsl/gsl_roots.h>

/*---------- local DEFINES ----------*/
#define EA_ACC          1E-9                    /* the timing accuracy of LALGetBinaryTimes in seconds */

/*----- Macros ----- */

/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

/*---------- Global variables ----------*/

/*---------- internal prototypes ----------*/

static double gsl_E_solver ( double E, void *p );

struct E_solver_params {
  double A, B, x0;
};

/*==================== FUNCTION DEFINITIONS ====================*/

/** Compute extra time-delays for a CW source in a (Keplerian) binary orbital system.
 *
 *
  \f[
  \newcommand{\Det}{{\mathrm{det}}}
  \newcommand{\tDet}{t_\Det}
  \newcommand{\tSSB}{t_{\mathrm{SSB}}}
  \newcommand{\tEM}{t_{\mathrm{em}}}
  \newcommand{\tRef}{t_{\mathrm{ref}}}
  \newcommand{\tPeri}{t_{\mathrm{p}}}
  \newcommand{\TPeri}{T_{\mathrm{p}}}
  \newcommand{\argp}{\omega}
  \newcommand{\sini}{\sin i}
  \newcommand{\fdot}{\dot{f}}
  \f]

 * For given binary-orbital parameters in \a Doppler, compute the source-frame emission times \f$\tEM(\tDet)\f$
 * of CW wavefronts arriving at the detector frame at times \f$\tDet\f$.
 * The input to this function should contain the wavefront arrival times \f$\tSSB(\tDet)\f$ in the solar-system barycenter (SSB) frame,
 * the additional time differences due to the binary orbital motion are then added by this function.
 *
 * NOTE: the output vector \a tSSBOut can be passed either as
 * - unallocated (in this case it must be <tt>(*tSSBOut)==NULL</tt>): it gets allocated here, or
 * - an allocated vector of the same length as the input vector \a tSSBIn
 *   (this is useful in order to minimize unnecessary memory allocs+frees on repeated calls using the same number of timesteps).
 *
 * NOTE2: it is OK to pass an identical in- and output-vector, i.e. <tt>(*tSSBOut) = tSSBIn</tt>: in this case the
 * binary-orbital time offsets will be added, modifying the input.
 *
 * NOTE3: (FIXME!) the SSBtimes structure naming is very misleading: 'SSB' here really refers to the *source emission time*,
 * not generally the solar-system-barycenter time, so this naming is only correct in the case of isolated NSs.
 * The additional time-delay computed by this function is between the inertial-SSB time and the source emission time, so
 * this is quite confusing right now and the corresponding variable-names and documentation should be fixed ASAP.
 *

 ## CW Timing Generalities ##

Notation: (all times are referring to a global Newtonian time axis):
\f{align*}{
& \tDet \ldots \text{wave-front arrival time at detector}\\
& \tSSB \ldots \text{arrival time at SSB}\\
& \tEM  \ldots \text{wave-front emission time at the CW source}\\
& \tRef \ldots \text{reference time at which} \{\phi_0, f, \fdot, \ldots\} \text{are defined}
\f}

The intrinsic CW phase in the source frame can be written as
\f{equation}{
\label{eq:sourcePhase}
\begin{split}
\Phi(\tEM) &= \phi_0 + 2 \pi \left[ f \, \Delta\tau + \frac{1}{2} \fdot \, \Delta\tau^2 + \ldots \right]\,,\\
\Delta\tau &\equiv \tEM - \tRef\,.
\end{split}
\f}

In order to relate this to the CW phase \f$\Phi_\Det\f$ arriving at the detector, we need the timing relation \f$\tEM(\tDet)\f$,
such that \f$\Phi_\Det(\tDet) = \Phi\left(\tEM(\tDet)\right)\f$.
We also compute the explicit values of the time-derivative, i.e. \f$\frac{d\tEM}{d \tDet}(\tDet)\f$, which are required
by XLALComputeFaFb(), for example.

XLALGetSSBtimes() computes the timing relation between detector time \f$\tDet\f$ and SSB time \f$\tSSB\f$, namely
\f{align}{
\Delta\tau_0(\tDet) &\equiv \tSSB(\tDet) - \tRef\,,\\
\dot{\tau}_0(\tDet) &\equiv \frac{d\tSSB}{d \tDet}(\tDet)\,
\f}
which is passed to this function as an input in \a tSSBIn. We therefore need to compute the extra time delay due to the binary-orbital motion, namely
\f{align}{
\Delta\tau(\tDet) &\equiv \tEM(\tDet) - \tRef = [\tEM(\tDet) - \tSSB(\tDet)] + [\tSSB(\tDet) - \tRef] \notag\\
                  &= [\tEM - \tSSB](\tDet) + \Delta\tau_0(\tDet) \,, \label{eq:7a}\\
\dot{\tau}(\tDet) &\equiv \frac{d\tEM}{d \tDet}(\tDet) = \frac{d\tEM}{d\tSSB}\left(\tSSB(\tDet)\right) \, \frac{d\tSSB}{d\tDet}(\tDet) \notag\\
                  &= \left[ \frac{d\tEM}{d\tSSB}\left(\tSSB(\tDet)\right) \right]\,\dot{\tau}_0(\tDet)\,, \label{eq:7b}
\f}

The relation between \f$\tEM\f$ and \f$\tSSB\f$ contains an unknown offset due to the distance of the binary system from the
SSB, which can be either constant, or changing at a constant or variable rate (eg due to accelerations in globular
clusters). However, we can absorb this unknown effect into the definition of the pulsar parameters
\f$\{\phi_0,f,\fdot,\ldots\}\f$, which are either unknown anyway, or affected by the exact same re-definitions if they are
known from photon astronomy. We therefore ignore this effect and pretend that there is no time delay between the SSB and
the binary-system barycenter (BSB), ie effectively \f$\tSSB = t_{\mathrm{BSB}}\f$.

### (Newtonian) Binary NS Timing equations ###

The extra time-delay from the binary orbital motion can be written as
\f{equation}{
  \label{eq:1}
  \tSSB(\tEM) = \tEM + R(\tEM)\,,
\f}
where \f$R\f$ is the radial distance from the BSB to the emitting NS along the line of sight, and here and in the following we are using
units where \f$c=1\f$. The sign convention is such that \f$R>0\f$ means that the NS is further away than the BSB, when \f$R<0\f$ it is closer.

\figure{inject_binary,eps,0.47,Binary orbit orientation parameters [T. Creighton]}

In terms of orbital parameters, the radial distance \f$R\f$ can be expressed as (see \figref{inject_binary}).
\f{equation}{
  \label{eq:2}
  R = r\,\sini\,\sin(\argp + \upsilon)\,,
\f}
where \f$r\f$ is the distance of the NS from the BSB (ie the focus of the ellips), \f$i\f$ is the inclination angle
between the orbital plane and the sky, \f$\argp\f$ is the argument of periapse, and \f$\upsilon\f$ is the <em>true anomaly</em> (ie
the angle from the periapse to the current NS location around the BSB.

\figure{Eccentric_and_true_anomaly,png,0.2,Definition of true anomaly 'v' and eccentric anomaly 'E' to describe an ellipse [Wikipedia]}

Using elementary trigonometry (cf. \figref{Eccentric_and_true_anomaly} and https://en.wikipedia.org/wiki/Eccentric_anomaly), one can see that the elliptical
orbit can be described in terms of the true anomaly \f$\upsilon\f$ as
\f{equation}{
  \label{eq:3}
  r(\upsilon) = \frac{a\,(1- e^2)}{1 +  e\cos\upsilon}\,,
\f}
and in terms of the <em>eccentric anomaly</em> \f$E\f$ as
\f{equation}{
  \label{eq:4}
  r(E) = a \, \left( 1 -  e\,\cos E\right)\,.
\f}
From \eqref{eq:3} and \eqref{eq:4} we easily obtain the relations
\f{equation}{
  \begin{split}
    \cos\upsilon &= \frac{\cos E -  e}{1 -  e\cos E}\,,\\
    \sin\upsilon &= \sin E \frac{\sqrt{1 -  e^2}}{1 -  e\cos E}\,,
  \end{split}
    \label{eq:9}
\f}
where in the second equation we have used the fact that \f$\sin E\f$ and \f$\sin\upsilon\f$ always have the same sign, as can be
seen from \figref{Eccentric_and_true_anomaly}.

The (Keplerian) motion of the NS on this elliptical orbit is described by Kepler's equation:
\f{equation}{
  \label{eq:6}
  \tEM - \tPeri = \frac{P}{2\pi}\left( E -  e\,\sin E\right)\,,
\f}
where we fixed the (discrete) gauge as \f$E=0\f$ at \f$\tEM=\tPeri\f$.

#### Algorithm to solve for \f$\tEM(\tSSB)\f$ ####

Following Teviet's strategy explained in LALGenerateEllipticSpinOrbitCW() we proceed by first expressing \f$\tSSB(E)\f$,
solving it numerically for \f$E(\tSSB)\f$, from which we obtain by substitution \f$\tEM(\tSSB)\f$ and
\f$d\tEM/d\tSSB\f$ required for \eqref{eq:7a},\eqref{eq:7b}.

We start by writing \eqref{eq:1} using \eqref{eq:6} as
\f{equation}{
  \label{eq:8}
  \tSSB(E) = \tPeri + \frac{P}{2\pi}\left( E -  e\,\sin E\right) + R(E)\,,
\f}
and using \eqref{eq:2}, \eqref{eq:3} with \eqref{eq:9}, we obtain
\f{equation}{
  \label{eq:R_E}
  R(E) = a\sini\left[ \sin\argp ( \cos E - e) + \cos\argp\,\sin E \sqrt{1 -  e^2} \right]\,,
\f}
Before solving \eqref{eq:8} for \f$E(\tSSB)\f$, it is useful to rewrite it a bit further:
First we note that at periapse (i.e. \f$E=0\f$),
\f{equation}{
\TPeri \equiv \tSSB(E=0) = \tPeri + a\sini \,\sin\argp\,(1-e)\,, \label{eq:defTPeri}
\f}
which corresponds to the <em>observed</em> (SSB) time of periapse passage.

Furthermore, as pointed out in LALGenerateEllipticSpinOrbitCW(), the (Newtonian) binary NS timing is periodic, namely
\f$\tSSB(E + m\,2\pi) = \tSSB(E) + m\,P\f$ for integer \f$m\f$, and we can therefore simplify the numerical solution of \eqref{eq:8} for \f$E(\tSSB)\f$ by
restricting the solution to the interval \f$E\in[0,2\pi)\f$ by mapping \f$\tSSB\f$ into \f$[\TPeri,\,\TPeri + P)\f$, via
\f{equation}{
x_0 \equiv \frac{2\pi}{P} (\tSSB - \TPeri) \mod 2\pi\,,
\f}
where we add \f$2\pi\f$ if \f$x_0 < 0\f$ to ensure that \f$x_0 \in [0,\,2\pi)\f$.
We can therefore rewrite \eqref{eq:8} as
\f{align}{
x_0 &= E + A\,\sin E + B \,( \cos E -1 )\,, \label{eq:E_tSSB}\\
A &\equiv \frac{2\pi}{P}\, a\sini \,\cos\argp\,\sqrt{1-e^2} - e\,, \label{eq:defA}\\
B &\equiv \frac{2\pi}{P}\, a\sini \,\sin\argp\,, \label{eq:defB}
\f}
which (after some substitutions) can be seen to agree with Teviet's equation found in LALGenerateEllipticSpinOrbitCW().

This is solved numerically for \f$E(\tSSB) \in [0,\,2\pi)\f$, and plugging this back into \eqref{eq:R_E} and \eqref{eq:1}
we obtain \f$\tEM(E)\f$, and from this the required timing relation \eqref{eq:7a}.

#### Computing the derivate \f$d\tEM/d\tSSB\f$ ####

We start from \eqref{eq:1} and write,
\f{equation}{
  \label{eq:12}
  \frac{d\tEM}{d\tSSB} = \left[\frac{d\tSSB}{d\tEM}\right]^{-1} = \left[1 + \frac{dR}{d E}\,\frac{d E}{d\tEM}\right]^{-1}\,.
\f}
From \eqref{eq:6} we obtain
\f{equation}{
  \label{eq:11}
  \frac{d E}{d\tEM} = \frac{2\pi}{P} \frac{1}{1 -  e\cos E}\,,
\f}
and \eqref{eq:R_E} yields
\f{equation}{
\label{eq:13}
\begin{split}
  \frac{dR}{d E} &= a\sini\left[ - \sin E \, \sin\argp + \cos E\,\cos\argp\,\sqrt{1- e^2}\right]\\
  &= \frac{P}{2\pi}\left( (A + e)\,\cos E - B\,\sin E \right)\,,
\end{split}
\f}
which results in the derivative of \eqref{eq:12} to be expressible as
\f{equation}{
\label{eq:dtEM_dtSSB}
\frac{d\tEM}{d\tSSB} = \frac{1 - e\cos E}{1 + A\,\cos E - B \,\sin E}\,.
\f}
*/
int
XLALAddBinaryTimes ( SSBtimes **tSSBOut,			//!< [out] reference-time offsets in emission frame: \f$\tEM(\tDet)-\tRef\f$ and \f$d\tEM/d\tDet\f$
                     const SSBtimes *tSSBIn,			//!< [in] reference-time offsets in SSB frame: \f$\tSSB(\tDet)-\tRef\f$ and \f$d\tSSB/d\tDet\f$
                     const PulsarDopplerParams *Doppler		//!< [in] pulsar Doppler parameters, includes binary orbit parameters */
                     )
{
  XLAL_CHECK ( tSSBIn != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSB'\n" );
  XLAL_CHECK ( tSSBIn->DeltaT != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSBIn->DeltaT'\n" );
  XLAL_CHECK ( tSSBIn->Tdot != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSBIn->Tdot'\n" );

  XLAL_CHECK ( Doppler != NULL, XLAL_EINVAL, "Invalid NULL input 'Doppler'\n");
  XLAL_CHECK ( Doppler->asini >= 0, XLAL_EINVAL );

  UINT4 numSteps = tSSBIn->DeltaT->length;		/* number of timesteps */
  XLAL_CHECK (tSSBIn->Tdot->length == numSteps, XLAL_EINVAL,
                   "Length tSSBIn->DeltaT = %d, while tSSBIn->Tdot = %d\n", numSteps, tSSBIn->Tdot->length );

  SSBtimes *binaryTimes;
  // ----- prepare output timeseries: either allocate or re-use existing
  if ( (*tSSBOut) == NULL )	// creating new output vector
    {
      XLAL_CHECK ( (binaryTimes = XLALDuplicateSSBtimes ( tSSBIn )) != NULL, XLAL_EFUNC );
    }
  else if ( (*tSSBOut) == tSSBIn )	// input==output vector
    {
      binaryTimes = (*tSSBOut);
    }
  else // input vector given, but not identical to output vector
    {
      binaryTimes = (*tSSBOut);
      // need to do a few more sanity checks
      XLAL_CHECK ( binaryTimes->DeltaT->length == numSteps, XLAL_EINVAL,
                   "Length (*tSSBOut)->DeltaT = %d, while tSSBIn->DeltaT = %d\n", binaryTimes->DeltaT->length, numSteps );
      XLAL_CHECK ( binaryTimes->Tdot->length == numSteps, XLAL_EINVAL,
                        "Length tSSBOut->Tdot = %d, while tSSBIn->Tdot = %d\n", binaryTimes->Tdot->length, numSteps );
      // ... and copy the vector contents from the input SSB vector
      binaryTimes->refTime = tSSBIn->refTime;
      memcpy ( binaryTimes->DeltaT->data, tSSBIn->DeltaT->data, numSteps * sizeof(binaryTimes->DeltaT->data[0]) );
      memcpy ( binaryTimes->Tdot->data,   tSSBIn->Tdot->data,   numSteps * sizeof(binaryTimes->Tdot->data[0]) );
    } // re-using input vector

  /* ----- convenience variables */
  REAL8 Porb 	= Doppler->period;		/* binary orbital period */
  REAL8 e 	= Doppler->ecc;			/* the eccentricity */
  REAL8 sqrtome2= sqrt(1.0 - e*e);
  REAL8 asini 	= Doppler->asini;		/* the projected orbital semimajor axis */
  REAL8 argp	= Doppler->argp;
  REAL8 sinw 	= sin ( argp );		/* the sin and cos of the argument of periapsis */
  REAL8 cosw 	= cos ( argp );
  REAL8 n    	= LAL_TWOPI / Porb;

  REAL8 refTimeREAL8 = XLALGPSGetREAL8 ( &tSSBIn->refTime );

  // compute time-independent coefficients for tSSB(E) equation
  REAL8 A = n * asini * cosw * sqrtome2 - e;  	// see Eq.(eq:defA)
  REAL8 B = n * asini * sinw;			// see Eq.(eq:defB)
  REAL8 tp = XLALGPSGetREAL8 ( &(Doppler->tp) );
  REAL8 Tp = tp + asini * sinw * (1 - e);	// see Eq.(eq:defTPeri)
  REAL8 acc = LAL_TWOPI * EA_ACC / Porb;   // root-finding accuracy, EA_ACC represents the required timing precision in seconds

  /* loop over the SFTs i */
  for ( UINT4 i = 0; i < numSteps; i++ )
    {
      REAL8 tSSB_i    = refTimeREAL8 + tSSBIn->DeltaT->data[i];	// SSB time for the current SFT midpoint
      REAL8 fracOrb_i = fmod ( tSSB_i - Tp, Porb ) / Porb; 	// fractional orbit
      if ( fracOrb_i < 0 ) {
        fracOrb_i += 1;	// enforce fracOrb to be within [0, 1)
      }
      REAL8 x0 = fracOrb_i * LAL_TWOPI;
      REAL8 E_i;              // eccentric anomaly at emission of the wavefront arriving in SSB at tSSB
      { // ---------- use GSL for the root-finding
        const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
        gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        REAL8 E_lo = 0, E_hi = LAL_TWOPI;	// gauge-choice mod (2pi)
        gsl_function F;
        struct E_solver_params pars = {A, B, x0};
        F.function = &gsl_E_solver;
        F.params = &pars;

        XLAL_CHECK ( gsl_root_fsolver_set(s, &F, E_lo, E_hi) == 0, XLAL_EFAILED );

        XLALPrintInfo ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "abstol", "err(est)");
        int max_iter = 100;
        int iter = 0;
        int status;
        do
          {
            iter++;
            status = gsl_root_fsolver_iterate(s);
            XLAL_CHECK ( (status == GSL_SUCCESS) || (status == GSL_CONTINUE), XLAL_EFAILED );
            E_i = gsl_root_fsolver_root(s);
            E_lo = gsl_root_fsolver_x_lower (s);
            E_hi = gsl_root_fsolver_x_upper (s);
            status = gsl_root_test_interval ( E_lo, E_hi, acc, 0 );

            if (status == GSL_SUCCESS) { XLALPrintInfo ("Converged:\n"); }
            XLALPrintInfo ("%5d [%.7f, %.7f] %.7f %+10.7g %10.7g\n", iter, E_lo, E_hi, E_i, acc, E_hi - E_lo);

          } while ( (status == GSL_CONTINUE) && (iter < max_iter) );

        XLAL_CHECK ( status == GSL_SUCCESS, XLAL_EMAXITER, "Eccentric anomaly: failed to converge within %d iterations\n", max_iter );
        gsl_root_fsolver_free(s);
      } // gsl-root finding block

      // use this value of E(tSSB) to compute the additional binary time delay
      REAL8 sinE = sin(E_i);
      REAL8 cosE = cos(E_i);

      REAL8 R 	    	= asini * ( sinw * ( cosE - e ) + cosw * sinE * sqrtome2 );	// see Eq.(eq:R_E)
      REAL8 dtEM_dtSSB 	= ( 1.0 - e * cosE ) / ( 1.0 + A * cosE - B * sinE );		// see Eq.(eq:dt_EMdtSSB)

      binaryTimes->DeltaT->data[i] -= R;
      binaryTimes->Tdot->data[i]   *= dtEM_dtSSB;

    } /* for i < numSteps */

  // pass back output SSB timings
  (*tSSBOut) = binaryTimes;

  return XLAL_SUCCESS;

} /* XLALAddBinaryTimes() */

/** Function implementing \eqref{eq:E_tSSB} to be solved via numerical root-finder for \f$E(\tSSB)\f$
 */
static double
gsl_E_solver ( REAL8 E, void *par )
{
  struct E_solver_params *params = (struct E_solver_params*) par;
  double A  = params->A;
  double B  = params->B;
  double x0 = params->x0;

  double diff = - x0 + ( E + A * sin(E) + B * ( cos(E) - 1.0 ) );

  return diff;
} // gsl_E_solver()


/**
 * Multi-IFO version of XLALAddBinaryTimes().

 * For given binary-orbital parameters in \a Doppler, compute the source-frame emission times \f$\tEM(\tDet)\f$
 * of CW wavefronts arriving at the detector frame at times \f$\tDet\f$.
 * The input to this function should contain the wavefront arrival times \f$\tSSB(\tDet)\f$ in the solar-system barycenter (SSB) frame,
 * the additional time differences due to the binary orbital motion are then added by this function.
 *
 * NOTE: the output vector \a multiSSBOut can be passed either as
 * - unallocated (in this case it must be <tt>(*tSSBOut)==NULL</tt>): it gets allocated here, or
 * - an allocated vector of the same length as the input vector \a multiSSBIn
 *   (this is useful in order to minimize unnecessary memory allocs+frees on repeated calls using the same number of timesteps).
 *
 * NOTE2: it is OK to pass an identical in- and output-vector, i.e. <tt>(*tSSBOut) = tSSBIn</tt>: in this case the
 * binary-orbital time offsets will be added, modifying the input.
 *
 * NOTE3: (FIXME!) the SSBtimes structure naming is very misleading: 'SSB' here really refers to the *source emission time*,
 * not generally the solar-system-barycenter time, so this naming is only correct in the case of isolated NSs.
 * The additional time-delay computed by this function is between the inertial-SSB time and the source emission time, so
 * this is quite confusing right now and the corresponding variable-names and documentation should be fixed ASAP.
 */
int
XLALAddMultiBinaryTimes ( MultiSSBtimes **multiSSBOut,		/**< [out] output SSB times */
                          const MultiSSBtimes *multiSSBIn,	/**< [in] SSB-timings for all input detector-state series */
                          const PulsarDopplerParams *Doppler	/**< [in] pulsar Doppler parameters, includes binary orbit parameters */
                          )
{
  /* check input */
  XLAL_CHECK ( multiSSBIn != NULL, XLAL_EINVAL, "Invalid NULL input 'multiSSB'\n");
  XLAL_CHECK ( Doppler != NULL, XLAL_EINVAL, "Invalid NULL input 'Doppler'\n");
  XLAL_CHECK ( Doppler->asini >= 0, XLAL_EINVAL );

  UINT4 numDetectors = multiSSBIn->length;

  MultiSSBtimes *multiBinaryTimes;
  // ----- prepare output timeseries: either allocate or re-use existing
  if ( (*multiSSBOut) == NULL )	// creating new output vector
    {
      XLAL_CHECK ( (multiBinaryTimes = XLALDuplicateMultiSSBtimes ( multiSSBIn )) != NULL, XLAL_EFUNC );
    }
  else // input vector given
    {
      multiBinaryTimes = (*multiSSBOut);
      XLAL_CHECK ( multiBinaryTimes->length == numDetectors, XLAL_EINVAL,
                   "Inconsistent length (*multiSSBOut)->length = %d, while multiSSBIn->length = %d\n", (*multiSSBOut)->length, numDetectors );
      // we'll leave all other sanity-checks to XLALAddBinaryTimes() calls
    }

  // ----- simply loop over detectors for XLALAddBinaryTimes()
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      int ret = XLALAddBinaryTimes ( &(multiBinaryTimes->data[X]), multiSSBIn->data[X], Doppler );
      XLAL_CHECK( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALAddBinaryTimes() failed for X=%d\n", X );
    } /* for X < numDet */

  // pass back result-vector
  (*multiSSBOut) = multiBinaryTimes;

  return XLAL_SUCCESS;

} /* XLALAddMultiBinaryTimes() */


/** Duplicate (ie allocate + copy) an input SSBtimes structure.
 * This can be useful for creating a copy before adding binary-orbital corrections in XLALAddBinaryTimes()
 */
SSBtimes *
XLALDuplicateSSBtimes ( const SSBtimes *tSSB )
{
  XLAL_CHECK_NULL ( tSSB != NULL, XLAL_EINVAL, "Invalid NULL input 'tSSB'\n" );

  UINT4 len;
  SSBtimes *ret;
  ret = XLALCalloc ( 1, len = sizeof (*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( 1, %d )\n", len );

  ret->refTime = tSSB->refTime;

  if ( tSSB->DeltaT )
    {
      len = tSSB->DeltaT->length;
      ret->DeltaT = XLALCreateREAL8Vector ( len );
      XLAL_CHECK_NULL ( ret->DeltaT != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed\n", len );
      memcpy ( ret->DeltaT->data, tSSB->DeltaT->data, len * sizeof(ret->DeltaT->data[0]) );
    }

  if ( tSSB->Tdot )
    {
      len = tSSB->Tdot->length;
      ret->Tdot = XLALCreateREAL8Vector ( len );
      XLAL_CHECK_NULL ( ret->Tdot != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed\n", len );
      memcpy ( ret->Tdot->data, tSSB->Tdot->data, len * sizeof(ret->Tdot->data[0]) );
    }

  return ret;

} /* XLALDuplicateSSBtimes() */


/** Duplicate (ie allocate + copy) an input MultiSSBtimes structure.
 */
MultiSSBtimes *
XLALDuplicateMultiSSBtimes ( const MultiSSBtimes *multiSSB )
{
  XLAL_CHECK_NULL ( multiSSB != NULL, XLAL_EINVAL, "Invalid NULL input 'multiSSB'\n");

  UINT4 len;
  MultiSSBtimes *ret;
  ret = XLALCalloc ( 1, len=sizeof(*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( 1, %d )\n", len );

  UINT4 numDetectors = multiSSB->length;
  ret->length = numDetectors;
  ret->data = XLALCalloc ( numDetectors, len = sizeof(ret->data[0]) );
  XLAL_CHECK_NULL ( ret->data != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( %d, %d )\n", numDetectors, len );

  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      ret->data[X] = XLALDuplicateSSBtimes ( multiSSB->data[X] );
      XLAL_CHECK_NULL( ret->data[X] != NULL, XLAL_EFUNC, "XLALDuplicateSSBtimes() failed for detector X=%d\n", X );
    } // for X < numDetectors

  return ret;

} /* XLALDuplicateMultiSSBtimes() */

/** For a given DetectorStateSeries, calculate the time-differences
 *  \f$\Delta T_\alpha\equiv T(t_\alpha) - T_0\f$, and their
 *  derivatives \f$\dot{T}_\alpha \equiv d T / d t (t_\alpha)\f$.
 *
 *  \note The return-vector is allocated here
 *
 */
SSBtimes *
XLALGetSSBtimes ( const DetectorStateSeries *DetectorStates,	/**< [in] detector-states at timestamps t_i */
                  SkyPosition pos,				/**< source sky-location */
                  LIGOTimeGPS refTime,				/**< SSB reference-time T_0 of pulsar-parameters */
                  SSBprecision precision			/**< relativistic or Newtonian SSB transformation? */
                  )
{
  XLAL_CHECK_NULL ( DetectorStates != NULL, XLAL_EINVAL, "Invalid NULL input 'DetectorStates'\n" );
  XLAL_CHECK_NULL ( precision < SSBPREC_LAST, XLAL_EDOM, "Invalid value precision=%d, allowed are [0, %d]\n", precision, SSBPREC_LAST -1 );
  XLAL_CHECK_NULL ( pos.system == COORDINATESYSTEM_EQUATORIAL, XLAL_EDOM, "Only equatorial coordinate system (=%d) allowed, got %d\n", COORDINATESYSTEM_EQUATORIAL, pos.system );

  UINT4 numSteps = DetectorStates->length;		/* number of timestamps */

  // prepare output SSBtimes struct
  int len;
  SSBtimes *ret = XLALCalloc ( 1, len = sizeof(*ret) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,%d)\n", len );
  ret->DeltaT = XLALCreateREAL8Vector ( numSteps );
  XLAL_CHECK_NULL ( ret->DeltaT != NULL, XLAL_EFUNC, "ret->DeltaT = XLALCreateREAL8Vector(%d) failed\n", numSteps );
  ret->Tdot = XLALCreateREAL8Vector ( numSteps );
  XLAL_CHECK_NULL ( ret->Tdot != NULL, XLAL_EFUNC, "ret->Tdot = XLALCreateREAL8Vector(%d) failed\n", numSteps );

  /* convenience variables */
  REAL8 alpha = pos.longitude;
  REAL8 delta = pos.latitude;
  REAL8 refTimeREAL8 = XLALGPSGetREAL8 ( &refTime );
  BarycenterBuffer *bBuffer = NULL;

  BarycenterInput XLAL_INIT_DECL(baryinput);

  /*----- now calculate the SSB transformation in the precision required */
  switch (precision)
    {
      REAL8 vn[3];		/* unit-vector pointing to source in Cart. coord. */

    case SSBPREC_NEWTONIAN:	/* use simple vr.vn to calculate time-delay */

      /*----- get the cartesian source unit-vector */
      vn[0] = cos(alpha) * cos(delta);
      vn[1] = sin(alpha) * cos(delta);
      vn[2] = sin(delta);

      for (UINT4 i = 0; i < numSteps; i++ )
	{
	  LIGOTimeGPS *ti = &(DetectorStates->data[i].tGPS);
	  /* DeltaT_alpha */
	  ret->DeltaT->data[i]  = XLALGPSGetREAL8 ( ti );
	  ret->DeltaT->data[i] += SCALAR(vn, DetectorStates->data[i].rDetector);
	  ret->DeltaT->data[i] -= refTimeREAL8;

	  /* Tdot_alpha */
	  ret->Tdot->data[i] = 1.0 + SCALAR(vn, DetectorStates->data[i].vDetector);

	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTIC:	/* use LALBarycenter() to get SSB-times and derivative */

      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = alpha;
      baryinput.delta = delta;
      baryinput.dInv = 0;

      for (UINT4 i = 0; i < numSteps; i++ )
	{
	  EmissionTime emit;
	  DetectorState *state = &(DetectorStates->data[i]);

	  baryinput.tgps = state->tGPS;

          if ( XLALBarycenter ( &emit, &baryinput, &(state->earthState) ) != XLAL_SUCCESS )
            XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALBarycenter() failed with xlalErrno = %d\n", xlalErrno );

	  ret->DeltaT->data[i] = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
	  ret->Tdot->data[i] = emit.tDot;

	} /* for i < numSteps */

      break;

    case SSBPREC_RELATIVISTICOPT:	/* use optimized version XLALBarycenterOpt() */

      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = alpha;
      baryinput.delta = delta;
      baryinput.dInv = 0;

      for ( UINT4 i = 0; i < numSteps; i++ )
        {
          EmissionTime emit;
          DetectorState *state = &(DetectorStates->data[i]);
          baryinput.tgps = state->tGPS;

          if ( XLALBarycenterOpt ( &emit, &baryinput, &(state->earthState), &bBuffer ) != XLAL_SUCCESS )
            XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALBarycenterOpt() failed with xlalErrno = %d\n", xlalErrno );

          ret->DeltaT->data[i] = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
          ret->Tdot->data[i] = emit.tDot;

        } /* for i < numSteps */
      // free buffer memory
      XLALFree ( bBuffer );

      break;

    default:
      XLAL_ERROR_NULL (XLAL_EFAILED, "\n?? Something went wrong.. this should never be called!\n\n" );
      break;
    } /* switch precision */

  /* finally: store the reference-time used into the output-structure */
  ret->refTime = refTime;

  return ret;

} /* XLALGetSSBtimes() */

/** Multi-IFO version of LALGetSSBtimes().
 * Get all SSB-timings for all input detector-series.
 *
 * NOTE: this functions *allocates* the output-vector,
 * use XLALDestroyMultiSSBtimes() to free this.
 */
MultiSSBtimes *
XLALGetMultiSSBtimes ( const MultiDetectorStateSeries *multiDetStates, /**< [in] detector-states at timestamps t_i */
                       SkyPosition skypos,		/**< source sky-position [in equatorial coords!] */
                       LIGOTimeGPS refTime,		/**< SSB reference-time T_0 for SSB-timing */
                       SSBprecision precision		/**< use relativistic or Newtonian SSB timing?  */
                       )
{
  /* check input */
  XLAL_CHECK_NULL ( multiDetStates != NULL, XLAL_EINVAL, "Invalid NULL input 'multiDetStates'\n");
  XLAL_CHECK_NULL ( multiDetStates->length > 0, XLAL_EINVAL, "Invalid zero-length 'multiDetStates'\n");

  UINT4 numDetectors = multiDetStates->length;

  // prepare return struct
  int len;
  MultiSSBtimes *ret = XLALCalloc ( 1, len = sizeof( *ret ) );
  XLAL_CHECK_NULL ( ret != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,%d)\n", len );
  ret->length = numDetectors;
  ret->data = XLALCalloc ( numDetectors, len=sizeof ( *ret->data ) );
  XLAL_CHECK_NULL ( ret->data != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(%d,%d)\n", numDetectors, len );

  // loop over detectors
  for ( UINT4 X = 0; X < numDetectors; X ++ )
    {
      ret->data[X] = XLALGetSSBtimes ( multiDetStates->data[X], skypos, refTime, precision );
      XLAL_CHECK_NULL ( ret->data[X] != NULL, XLAL_EFUNC, "ret->data[%d] = XLALGetSSBtimes() failed with xlalErrno = %d\n", X, xlalErrno );

    } /* for X < numDet */

  return ret;

} /* XLALGetMultiSSBtimes() */

/** Find the earliest timestamp in a multi-SSB data structure
 *
*/
int XLALEarliestMultiSSBtime ( LIGOTimeGPS *out,              /**< output earliest GPS time */
                               const MultiSSBtimes *multiSSB,      /**< input multi SSB SFT-midpoint timestamps */
                               const REAL8 Tsft                    /**< the length of an SFT */
                               )
{
  UINT4 i,j;
  LIGOTimeGPS t;
  REAL8 delta;

  /* check sanity of input */
  if ( !multiSSB || (multiSSB->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  if ( !multiSSB->data[0] || (multiSSB->data[0]->DeltaT->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }


  /* initialise the earliest and latest sample value */
  out->gpsSeconds = multiSSB->data[0]->refTime.gpsSeconds;
  out->gpsNanoSeconds = multiSSB->data[0]->refTime.gpsNanoSeconds;
  delta = multiSSB->data[0]->DeltaT->data[0] - 0.5*Tsft*multiSSB->data[0]->Tdot->data[0];
  if ( (XLALGPSAdd(out,delta)) == NULL) {
    XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  /* loop over detectors */
  for (i=0;i<multiSSB->length;i++) {

    /* loop over all SSB times and find the earliest SSB SFT start time */
    for (j=0;j<multiSSB->data[i]->DeltaT->length;j++) {

      /* reset the reference time */
      t.gpsSeconds = multiSSB->data[i]->refTime.gpsSeconds;
      t.gpsNanoSeconds = multiSSB->data[i]->refTime.gpsNanoSeconds;

      /* compute SSB time - we approximate the SFT start time in the SSB as t_mid_SSB - 0.5*Tsft*dt_SSB/dt_det */
      delta = multiSSB->data[i]->DeltaT->data[j] - 0.5*Tsft*multiSSB->data[i]->Tdot->data[j];
      if ( (XLALGPSAdd(&t,delta)) == NULL) {
        XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
        XLAL_ERROR (XLAL_ENOMEM);
      }

      /* compare it to the existing earliest */
      if ( (XLALGPSCmp(out,&t) == 1 ) ) {
        out->gpsSeconds = t.gpsSeconds;
        out->gpsNanoSeconds = t.gpsNanoSeconds;
      }

    }

  }

  /* success */
  return XLAL_SUCCESS;

} /* XLALEarliestMultiSSBtime() */

/** Find the latest timestamp in a multi-SSB data structure
 *
*/
int XLALLatestMultiSSBtime ( LIGOTimeGPS *out,                   /**< output latest GPS time */
                             const MultiSSBtimes *multiSSB,      /**< input multi SSB SFT-midpoint timestamps */
                             const REAL8 Tsft                    /**< the length of an SFT */
                             )
{
  UINT4 i,j;
  LIGOTimeGPS t;
  REAL8 delta;

  /* check sanity of input */
  if ( !multiSSB || (multiSSB->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  if ( !multiSSB->data[0] || (multiSSB->data[0]->DeltaT->length == 0) )
    {
      XLALPrintError ("%s: empty multiSSBtimes input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }


  /* initialise the earliest and latest sample value */
  out->gpsSeconds = multiSSB->data[0]->refTime.gpsSeconds;
  out->gpsNanoSeconds = multiSSB->data[0]->refTime.gpsNanoSeconds;
  delta = multiSSB->data[0]->DeltaT->data[0] + 0.5*Tsft*multiSSB->data[0]->Tdot->data[0];
  if ( (XLALGPSAdd(out,delta)) == NULL) {
    XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
    XLAL_ERROR (XLAL_ENOMEM);
  }
  /* loop over detectors */
  for (i=0;i<multiSSB->length;i++) {

    /* loop over all SSB times and find the latest SSB SFT start time */
    for (j=0;j<multiSSB->data[i]->DeltaT->length;j++) {

      /* reset the reference time */
      t.gpsSeconds = multiSSB->data[i]->refTime.gpsSeconds;
      t.gpsNanoSeconds = multiSSB->data[i]->refTime.gpsNanoSeconds;

      /* compute SSB time - we approximate the SFT end time in the SSB as t_mid_SSB + 0.5*Tsft*dt_SSB/dt_det */
      delta = multiSSB->data[i]->DeltaT->data[j] + 0.5*Tsft*multiSSB->data[i]->Tdot->data[j];
      if ( (XLALGPSAdd(&t,delta)) == NULL) {
        XLALPrintError ("%s: XLALGPSAdd() failed!  errno = %d!\n", __func__, xlalErrno );
        XLAL_ERROR (XLAL_ENOMEM);
      }

      /* compare it to the existing earliest */
      if ( (XLALGPSCmp(out,&t) == -1 ) ) {
        out->gpsSeconds = t.gpsSeconds;
        out->gpsNanoSeconds = t.gpsNanoSeconds;
      }

    }

  }

  /* success */
  return XLAL_SUCCESS;

} /* XLALLatestMultiSSBtime() */

/* ===== Object creation/destruction functions ===== */

/** Destroy a MultiSSBtimes structure.
 * Note, this is "NULL-robust" in the sense that it will not crash
 * on NULL-entries anywhere in this struct, so it can be used
 * for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiSSBtimes ( MultiSSBtimes *multiSSB )
{
  UINT4 X;
  SSBtimes *tmp;

  if ( ! multiSSB )
    return;

  if ( multiSSB->data )
    {
      for ( X=0; X < multiSSB->length; X ++ )
	{
	  if ( (tmp = multiSSB->data[X]) != NULL )
	    {
	      if ( tmp->DeltaT )
		XLALDestroyREAL8Vector ( tmp->DeltaT );
	      if ( tmp->Tdot )
		XLALDestroyREAL8Vector ( tmp->Tdot );
	      LALFree ( tmp );
	    } /* if multiSSB->data[X] */
	} /* for X < numDetectors */
      LALFree ( multiSSB->data );
    }
  LALFree ( multiSSB );

  return;

} /* XLALDestroyMultiSSBtimes() */
