/*
*  Copyright (C) 2007 Philip Charlton, Duncan Brown, Jolien Creighton, David McKechan, Stephen Fairhurst, Teviet Creighton, Thomas Cokelaer, John Whelan
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

#ifndef _GENERATEPPNINSPIRAL_H
#define _GENERATEPPNINSPIRAL_H

#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>
#include <lal/Random.h>
#include <lal/LALSimInspiral.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 \addtogroup GeneratePPNInspiral_h
 \author Creighton, T. D.

 \brief Provides routines to generate restricted parametrized
post\f${}^{5/2}\f$-Newtonian inspiral waveforms.

\heading{Synopsis}
\code
#include <lal/GeneratePPNInspiral.h>
\endcode

\heading{Description}

This header covers routines to generate a "restricted" parametrized
post\f${}^{5/2}\f$-Newtonian binary inspiral waveform in the time domain.
That is, the calculation of the wave phase is accurate to
post\f${}^{5/2}\f$-Newtonian order (including corrections up to order
\f$v^5/c^5\f$, where \f$v\f$ is the orbital speed), but the wave amplitudes
are accurate only to leading (post\f${}^0\f$-Newtonian) order.
Furthermore, at each order the post\f${}^{n/2}\f$-Newtonian correction can
be turned on, off, or set to an unphysical value, by adjusting a
parameter \f$p_n\f$.

The post-Newtonian expansion implicitly assumes an \e adiabatic
inspiral, where one can represent the waveform by an "instantaneous"
amplitude and frequency that vary over timescales longer than one wave
period.  The \e orbital frequency of the system to
post\f${}^{5/2}\f$-Newtonian order is given in Eqs.\ 6.4.1 and 6.9.1
of [\ref GRASP2000]; here we work entirely in terms of the
<em>gravitational-wave</em> frequency, which is twice the orbital
frequency:
\anchor eq_ppn_freq \f{eqnarray}{
f(t) & = & \frac{M_\odot}{8\pi T_\odot m_\mathrm{tot}}\left\{
	p_0\Theta^{-3/8}+
	p_1\Theta^{-1/2}+
	p_2\left(\frac{743}{2688}+\frac{11}{32}\eta\right)\Theta^{-5/8}-
	p_3\frac{3\pi}{10}\Theta^{-3/4} \right. \nonumber \\
\tag{eq_ppn_freq}
& & \left.+ p_4\left(\frac{1855099}{14450688}+\frac{56975}{258048}\eta+
		\frac{371}{2048}\eta^2\right)\Theta^{-7/8}-
	p_5\left(\frac{7729}{21504}+\frac{3}{256}\eta\right)\pi\Theta^{-1}
	\right\} \; ,
\f}
where \f$M_\odot\f$ is the mass of the Sun,
\f$T_\odot=GM_\odot/c^3=4.925491\times10^{-6}\f$s is the "geometrized"
solar mass in time units, \f$m_\mathrm{tot}=m_1+m_2\f$ is the total mass
of the binary, \f$\eta=m_1m_2/m_\mathrm{tot}^2\f$ is the (symmetric) mass
ratio parameter, and \f$\Theta\f$ is a dimensionless time parameter:
\anchor eq_ppn_theta \f{equation}{
\tag{eq_ppn_theta}
\Theta(t) = \frac{\eta M_\odot}{5T_\odot m_\mathrm{tot}}(t_c-t) \; .
\f}
Here \f$t_c\f$ is the time of coalescence of the two masses in the
point-mass approximation.  The post-Newtonian parameters \f$p_k\f$ are
defined such that in a normal (physical) post\f${}^{n/2}\f$-Newtonian
expansion, one sets \f$p_1=0\f$ and \f$p_{k>n}=0\f$, and \f$p_k=1\f$ for all other
\f$k\f$.  However, changing this convention can be used to model in an
approximate way things such as spin, eccentricity, or non-GR theories
of gravity.  We also note that while most terms are normalized to
their normal post-Newtonian values, the normalization on the \f$p_1\f$
term is completely arbitrary, since it is zero in a normal
post-Newtonian expansion.

The wave phase as a function of time can be computed analytically from
Eq.\eqref{eq_ppn_freq} as \f$\phi_\mathrm{orb}=2\pi\int f\,dt\f$:
\anchor eq_ppn_phi \f{eqnarray}{
\phi(t) & = & \phi_c - \frac{2}{\eta}\left\{
	p_0\Theta^{5/8}+
	p_1\frac{5}{4}\Theta^{1/2}+
	p_2\left(\frac{3715}{8064}+\frac{55}{96}\eta\right)\Theta^{3/8}-
	p_3\frac{3\pi}{4}\Theta^{1/4} \right. \nonumber \\
\tag{eq_ppn_phi}
& & \left.+ p_4\left(\frac{9275495}{14450688}+\frac{284875}{258048}\eta+
		\frac{1855}{2048}\eta^2\right)\Theta^{1/8}-
	p_5\left(\frac{38645}{172032}+\frac{15}{2048}\eta\right)\pi
		\log\left(\frac{\Theta}{\Theta_0}\right)\right\} \; .
\f}
Here \f$\Theta_0\f$ is an arbitrary constant; changing it is equivalent to
changing \f$\phi_c\f$.  We note that the post\f${}^{5/2}\f$-Newtonian term
introduces a late-time divergence in phase which renders meaningless
the interpretation of \f$\phi_c\f$ as "phase at coalescence"; in our
convention we define \f$\phi_c\f$ to correspond to the case \f$\Theta_0=1\f$.

We refer the interested reader to Sec.\ 6.6 of\ [\ref GRASP2000]
for a discussion of how propagation effects shift the phase of the
waveform relative to the orbital phase.  To summarize, though: A
changing propagation delay does introduce a time-dependent phase shift
in the waveform, but the dependence on \f$t\f$ is weak except at very late
times; although it looks like a post\f${}^{3/2}\f$-Newtonian phase
correction, it can in fact be represented as a post\f${}^{3}\f$-Newtonian
phase correction combined with a post\f${}^{3/2}\f$-Newtonian amplitude
correction.  Since we are concerned with \e restricted
post\f${}^{5/2}\f$-Newtonian waveforms, which model the amplitude only to
leading (post\f${}^0\f$-Newtonian) order, we can ignore these propagation
effects.

To leading order, then, the amplitude of the + and \f$\times\f$
polarizations of the wave are given by Eqs.\ 6.6.1--6.6.4
of\ [\ref GRASP2000] as:
\anchor eq_ppn_aplus \anchor eq_ppn_across \f{eqnarray}{
\tag{eq_ppn_aplus}
A_+(t) & = & -\frac{2T_\odot c}{D}(1+\cos^2 i)
	\left(\frac{\eta m_\mathrm{tot}}{M_\odot}\right)
	\left[\frac{\pi T_\odot	m_\mathrm{tot}f(t)}{M_\odot}
	\right]^{2/3} \; , \\
\tag{eq_ppn_across}
A_\times(t) & = & -\frac{2T_\odot c}{D}(2\cos i)
	\left(\frac{\eta m_\mathrm{tot}}{M_\odot}\right)
	\left[\frac{\pi T_\odot	m_\mathrm{tot}f(t)}{M_\odot}
	\right]^{2/3} \; ,
\f}
where \f$D\f$ is the distance to the source and \f$i\f$ is the inclination of
the axis of the source to the line of sight.  The normal polarization
convention in\ [\ref Will_C_1996] is used, where the reference
\f$x\f$-coordinate axis for the + and \f$\times\f$ polarization tensors is the
ascending node of the rotational plane as it crosses the plane
transverse to the propagation direction.  This convention implies that
the + and \f$\times\f$ waveforms are elliptically polarized as follows:
\anchor eq_ppn_hplus \anchor eq_ppn_hcross \f{eqnarray}{
\tag{eq_ppn_hplus}
h_+(t) & = & A_+(t)\cos\phi(t) \; , \\
\tag{eq_ppn_hcross}
h_\times(t) & = & A_\times(t)\sin\phi(t) \; .
\f}

*/
/*@{*/

/** \name Error Codes */
/*@{*/
#define GENERATEPPNINSPIRALH_ENUL  1	/**< Unexpected null pointer in arguments */
#define GENERATEPPNINSPIRALH_EOUT  2	/**< output field a, f, phi, or shift already exists */
#define GENERATEPPNINSPIRALH_ETBAD 3	/**< Bad sampling interval */
#define GENERATEPPNINSPIRALH_EFBAD 4	/**< Bad starting frequency; could not get valid start time */
#define GENERATEPPNINSPIRALH_EPBAD 5	/**< Bad post-Newtonian parameters */
#define GENERATEPPNINSPIRALH_EMBAD 6	/**< Bad masses */
#define GENERATEPPNINSPIRALH_EDBAD 7	/**< Bad distance */
#define GENERATEPPNINSPIRALH_EMEM  8	/**< Out of memory */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATEPPNINSPIRALH_MSGENUL  "Unexpected null pointer in arguments"
#define GENERATEPPNINSPIRALH_MSGEOUT  "output field a, f, phi, or shift already exists"
#define GENERATEPPNINSPIRALH_MSGETBAD "Bad sampling interval"
#define GENERATEPPNINSPIRALH_MSGEFBAD "Bad starting frequency; could not get valid start time"
#define GENERATEPPNINSPIRALH_MSGEPBAD "Bad post-Newtonian parameters"
#define GENERATEPPNINSPIRALH_MSGEMBAD "Bad masses"
#define GENERATEPPNINSPIRALH_MSGEDBAD "Bad distance"
#define GENERATEPPNINSPIRALH_MSGEMEM  "Out of memory"
/** \endcond */

/** \name More Termination conditions

In addition to the error conditions above, there are a number of ways
that the signal generation routine can terminate gracefully while
still returning a valid waveform.  In many cases one \e wants to
continue generating a waveform "until things fall apart"; the
following codes, returned in the \c PPNParamStruc below, allow the
waveform generator to report exactly \e how things fell apart.

For the sake of LAL namespace conventions, these termination codes are
<tt>\#define</tt>d and autodocumented exactly like error codes.
*/
/*@{*/
#define GENERATEPPNINSPIRALH_EFSTOP     0	/**< Reached requested termination frequency */
#define GENERATEPPNINSPIRALH_ELENGTH    1	/**< Reached maximum length, or end of provided time series vector */
#define GENERATEPPNINSPIRALH_EFNOTMON   2	/**< Frequency no longer increasing monotonically */
#define GENERATEPPNINSPIRALH_EPNFAIL    3	/**< Evolution dominated by higher-order PN terms */
#define GENERATEPPNINSPIRALH_ERTOOSMALL 4	/**< Orbital radius too small for PN approximation */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATEPPNINSPIRALH_MSGEFSTOP     "Reached requested termination frequency"
#define GENERATEPPNINSPIRALH_MSGELENGTH    "Reached maximum length, or end of provided time series vector"
#define GENERATEPPNINSPIRALH_MSGEFNOTMON   "Frequency no longer increasing monotonically"
#define GENERATEPPNINSPIRALH_MSGEPNFAIL    "Evolution dominated by higher-order PN terms"
#define GENERATEPPNINSPIRALH_MSGERTOOSMALL "Orbital radius too small for PN approximation"
/** \endcond */

/*
 * FIXME: SWIG with -Werror won't let const CHAR *termDescription pass,
 *        as it leaves a lot of potential for memory leaks. I choose to
 *        make this an opaque struct.
 */
#ifndef SWIG
/** This structure stores the parameters for constructing a restricted
 * post-Newtonian waveform.  It is divided into three parts: parameters
 * passed along to the output structure but not used by waveform
 * generator, parameters used as input to the waveform generator, and
 * parameters set by the generator to evaluate its success.
 */
typedef struct tagPPNParamStruc {
  /** \name Passed parameters. */
  /*@{*/
  SkyPosition position; /**< location of source on sky */
  REAL4 psi;            /**< polarization angle (radians) */
  LIGOTimeGPS epoch;    /**< start time of output time series */
  /*@}*/

  /**\name Input parameters. */
  /*@{*/
  REAL4 mTot;       	/**< The total mass \f$m_\mathrm{tot}=m_1+m_2\f$ of the binary system, in solar masses */
  REAL4 eta;        	/**< The mass ratio \f$\eta=m_1m_2/m_\mathrm{tot}^2\f$ of the binary system;  Physically this
                         * parameter must lie in the range \f$\eta\in(0,1/4]\f$; values outside of
                         * this range may be permitted in order to represent "nonphysical"
                         * post-Newtonian expansions
                         */
  REAL4 d;          	/**< The distance to the system, in metres */
  REAL4 inc;        	/**< The inclination of the system to the line of sight, in radians */
  REAL4 phi;        	/**< The phase at coalescence \f$\phi_c\f$ (or arbitrary reference phase for a post\f${}^{5/2}\f$-Newtonian
                         * approximation), in radians
                         */
  REAL8 deltaT;     	/**< The requested sampling interval of the waveform, in s */
  REAL4 fStartIn;   	/**< The requested starting frequency of the waveform, in Hz */
  REAL4 fStopIn;    	/**< The requested termination frequency of
                         * the waveform, in Hz;  If set to 0, the waveform will be generated
                         * until a termination condition (above) is met;  If set to a negative
                         * number, the generator will use its absolute value as the terminating
                         * frequency, but will ignore post-Newtonian breakdown; it will terminate
                         * only at the requested frequency \f$-\mathtt{fStopIn}\f$, a local maximum
                         * frequency, or the central singularity
                         */
  UINT4 lengthIn;   	/**< The maximum number of samples in the generated waveform;  If zero, the waveforms can be arbitrarily long */
  REAL4Vector *ppn; 	/**< The parameters \f$p_n\f$ selecting the type of post-Newtonian expansion;  If \c ppn=\c NULL, a "normal" (physical) expansion is assumed */
  INT4 ampOrder;    	/**< PN amplitude selection 0-5 */
  /*@}*/

  /** \name Output parameters. */
  /*@{*/
  REAL8 tc;         	/**< The time \f$t_c-t\f$ from the start of the waveform to coalescence (in the point-mass approximation), in s */
  REAL4 dfdt;       	/**< The maximum value of \f$\Delta f\Delta t\f$ encountered over any timestep \f$\Delta t\f$ used in generating the waveform */
  REAL4 fStart;     	/**< The actual starting frequency of the waveform, in Hz (normally close but not identical to \c fStartIn) */
  REAL4 fStop;      	/**< The frequency at the termination of the waveform, in Hz */
  UINT4 length;     	/**< The length of the generated waveform */
  INT4 termCode;    	/**< The termination condition (above) that stopped computation of the waveform */
  const CHAR *termDescription; /**< The termination code description (above) */
  /*@}*/
} PPNParamStruc;
#else  /* SWIG */
typedef struct tagPPNParamStruc PPNParamStruc;
#endif  /* SWIG */

/** This structure stores the position and mass parameters of a galactic inspiral event.
 */
typedef struct tagGalacticInspiralParamStruc {
  REAL4 rho;    		/**< The distance of the binary system from the Galactic axis, in kpc */
  REAL4 z;	      		/**< The distance of the system from the Galactic plane, in kpc */
  REAL4 lGal;   		/**< The Galactocentric Galactic longitude of the system (ie the Galactic longitude of the direction <em>from
                                 * the Galactic centre</em> through the system), in radians; See\ \ref SkyCoordinates_h for the definition of this quantity
                                 */
  REAL4 m1, m2; 		/**< The masses of the binary components, in solar masses */
  LIGOTimeGPS geocentEndTime; 	/**< The geocentric end time of the inspiral event */
} GalacticInspiralParamStruc;

/** UNDOCUMENTED */
typedef struct tagAmpSwitchStruc {
	UINT4 q0, q1, q2, q3, q4, q5;
} AmpSwitchStruc;



/* ---------- Function prototypes. ----------  */
void
LALGeneratePPNInspiral( LALStatus     *,
			CoherentGW    *output,
			PPNParamStruc *params );




void
LALGeneratePPNAmpCorInspiral( LALStatus     *,
			CoherentGW    *output,
			PPNParamStruc *params );






void
LALGetInspiralParams( LALStatus                  *,
		      PPNParamStruc              *output,
		      GalacticInspiralParamStruc *input,
		      RandomParams               *params );




void
LALGenerateInspiralSmooth( LALStatus            *,
		      	   CoherentGW		**output,
			   PPNParamStruc	*params,
			   REAL4		*qfactor);


/*@}*/ /* end:GeneratePPNInspiral_h */


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATEPPNINSPIRAL_H */

