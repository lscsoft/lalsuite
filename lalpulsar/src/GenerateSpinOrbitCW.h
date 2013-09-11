/*
*  Copyright (C) 2007 Reinhard Prix, Teviet Creighton
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

#ifndef _GENERATESPINORBITCW_H
#define _GENERATESPINORBITCW_H

#include <lal/LALStdlib.h>
#include <lal/PulsarSimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup GenerateSpinOrbitCW_h
 * \author Creighton, T. D.
 *
 * \brief Provides routines to generate continuous waveforms with spindown and orbital modulation.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/GenerateSpinOrbitCW.h>
 * \endcode
 *
 * This header covers routines to generate continuous quasiperiodic
 * waveforms with a smoothly-varying intrinsic frequency modulated by
 * orbital motions around a binary companion.  The intrinsic frequency is
 * modeled by Taylor series coefficients as in \ref GenerateTaylorCW.h,
 * and the orbital modulation is described by a reduced set of orbital
 * parameters.  Note that the routines do \e not account for spin
 * precession, accretion processes, or other complicating factors; they
 * simply Doppler-modulate a polynomial frequency function.
 *
 * The frequency and phase of the wave in the source's rest frame are
 * given by \eqref{eq_taylorcw-freq} and \eqref{eq_taylorcw-phi} of
 * \ref GenerateTaylorCW_h, where \f$t\f$ is the proper time in this rest
 * frame.  The frequency and phase of the wave fronts crossing a
 * reference point in an inertial frame (e.g.\ the Solar system
 * barycentre) are simply \f$f[t(t_r)]\f$ and \f$\phi[t(t_r)]\f$, where
 * \f{equation}{
 * \label{eq_spinorbit-tr}
 * t_r = t + R(t)/c
 * \f}
 * is the (retarded) time measured at the inertial reference point a
 * distance \f$r\f$ from the source.
 *
 * The generation of the waveform thus consists of computing the radial
 * component \f$R(t)\f$ of the orbital motion of the source in the binary
 * centre-of-mass frame, inverting \eqref{eq_spinorbit-tr} to find
 * the "emission time" \f$t\f$ for a given "detector time" \f$t_r\f$, and
 * plugging this into the Taylor expansions to generate the instantaneous
 * frequency and phase.  The received frequency is also multiplied by the
 * instantaneous Doppler shift \f$[1+\dot{R}(t)/c]^{-1}\f$ at the time of
 * emission.
 *
 * Since we do not include precession effects, the polarization state of
 * the wave is constant: we simply specify the polarization amplitudes
 * \f$A_+\f$, \f$A_\times\f$ and the polarization phase \f$\psi\f$ based on the
 * (constant) orientation of the source's \e rotation.  The following
 * discussion defines a set of parameters for the source's orbital
 * \e revolution, which we regard as completely independent from its
 * rotation.
 *
 * ### Orbital motion ###
 *
 * \figure{inject_binary,eps,0.47,Binary orbit orientation parameters}
 *
 * \figref{inject_binary} illustrates the notation conventions
 * defining a binary orbit.  We define a radial axis \f$R\f$ directed
 * \e from the observer (Earth) \e to the source, as shown.  The
 * horizontal plane is thus the plane of the sky, and the direction
 * marked \f$N\f$ is the direction along a meridian towards the North
 * celestial pole.  The tilted plane is the plane of the binary orbit,
 * and the axis labeled \f$z\f$ is the normal to this plane directed such
 * that the orbit is right-handed about this axis.  The <em>ascending
 * node</em> of the orbit, denoted by \f$\Omega\f$, is the direction
 * defined by \f$\hat{\mathbf{\mathit{R}}}\times\hat{\mathbf{\mathit{z}}}\f$.
 * The binary orbit itself is shown as an off-centred ellipse, with the
 * barycentre at one of its foci; the wave-emitting source is also shown.
 *
 * The <em>inclination angle</em> \f$i\f$ is the angle between the sky and
 * orbital planes.  The <em>longitude of the ascending node</em> \f$\Omega\f$
 * is the angle in the plane of the sky from the North direction to the
 * ascending node, measured right-handed about
 * \f$\hat{\mathbf{\mathit{R}}}\f$.  The <em>argument of the periapsis</em>
 * \f$\omega\f$ is the angle in the orbital plane from the ascending node to
 * the direction of periapsis (point where the source is closest to the
 * system barycentre), and the <em>true anomaly</em> \f$\upsilon(t)\f$ of the
 * source is the angle from the periapsis to the current location of the
 * source; both angles are measured right-handed about
 * \f$\hat{\mathbf{\mathit{z}}}\f$ (i.e.\ prograde).  The <em>periapsis
 * separation</em> \f$r_p\f$ is the distance from the periapsis to the
 * barycentre, and we denote the \e eccentricity of the orbital
 * ellipse as \f$e\f$, so that the separation between the source and the
 * barycentre at any time is \f$r=r_p(1+e)/(1+e\cos\upsilon)\f$.
 *
 * In this convention, \f$i\in[0,\pi]\f$ and \f$\Omega\in[0,2\pi)\f$.  Another
 * convention common in astronomy is to restrict \f$\Omega\f$ to the range
 * \f$[0,\pi)\f$, refering to whichever node (ascending or descending) lies
 * in this range.  The argument of the periapsis \f$\omega\f$ is then also
 * measured from this node.  In this case the range of \f$i\f$ must be
 * extended to \f$(-\pi,\pi]\f$; it is negative if the reference node is
 * descending, and positive if it is ascending.  The formulae that follow
 * are the same in either convention, though, since one can verify that
 * adding \f$\pi\f$ to \f$\Omega\f$ and \f$\omega\f$ is equivalent to reversing the
 * sign on \f$i\f$.
 *
 * Some spherical trigonometry gives us \f$R=r\sin(\omega+\upsilon)\sin i\f$.
 * We can differentiate \f$R\f$ with respect to \f$t\f$, and apply Keplers
 * second law
 * \f$r^2\dot{\upsilon}=r_p^2\dot{\upsilon}_p=\mathrm{constant}\f$, where
 * \f$\dot{\upsilon}_p\f$ is the angular speed at periapsis, to get:
 * \f{eqnarray}{
 * \label{eq_orbit-r}
 * R & = & R_0 + \frac{(1+e) r_p\sin i}{1+e\cos\upsilon}
 * \sin(\omega+\upsilon) \;,\\
 * \label{eq_orbit-rdot}
 * \dot{R} & = & \dot{R}_0 + \frac{\dot{\upsilon}_p r_p\sin i}{1+e}
 * \left[ \cos(\omega+\upsilon) + e\cos\omega \right] \;.
 * \f}
 * Without loss of generality, we will henceforth drop the offsets \f$R_0\f$
 * and (constant) \f$\dot{R}_0\f$ from these equations.  This means that we
 * ignore the overall propagation delay between the \f$R=R_0\f$ plane and the
 * observer, and incorporate any (constant) Doppler shifts due to net
 * centre-of-mass motions into the values of \f$f\f$ and \f$\dot{\upsilon}_p\f$.
 * The resulting times and parameter values are referred to as being in
 * the \e barycentric frame.  The only time delays and Doppler shifts
 * that we explicitly treat are those arising from the motion of the
 * source relative to the \f$R=R_0\f$ sky plane passing through the system
 * barycentre.
 *
 * All we need now to determine the orbital motion is an equation for
 * \f$\upsilon(t)\f$.  Many basic astronomy textbooks give exact but
 * transcendental expressions relating \f$\upsilon\f$ and \f$t\f$ for elliptical
 * orbits with \f$0\leq e<1\f$, and/or series expansions of \f$\upsilon(t)\f$ for
 * \f$e\ll1\f$.  However, for a generic binary system we cannot guarantee
 * that \f$e\ll1\f$, and for now we would like to retain the possibility of
 * modeling open orbits with \f$e\geq1\f$.  For now we will simply present
 * the exact formulae, and discuss the numerical solution methods in the
 * modules under this header.
 *
 * Let \f$t_p\f$ be the time of a periapsis passage (preferably a recent one
 * in the case of closed orbits).  We express both \f$t\f$ and \f$\upsilon\f$ in
 * terms of an intermediate variable \f$E\f$ (called the <em>eccentric
 * anomaly</em> for elliptic orbits, unnamed for open orbits).  The formulae
 * are:
 * \f{equation}{
 * \label{eq_spinorbit-t}
 * t - t_p = \left\{ \begin{array}{l@{\qquad}c}
 * \frac{1}{\dot{\upsilon}_p} \sqrt{\frac{1+e}{(1-e)^3}}
 * \left( E - e\sin E \right) & 0 \leq e < 1 \\ & \\
 * \frac{1}{\dot{\upsilon}_p} E
 * \left( 1 + \frac{E^2}{12} \right) & e = 1 \\ & \\
 * \frac{1}{\dot{\upsilon}_p} \sqrt{\frac{e+1}{(e-1)^3}}
 * \left( e\sinh E - E \right) & e > 1
 * \end{array} \right.
 * \f}
 * \f{equation}{
 * \label{eq_spinorbit-upsilon}
 * \begin{array}{c} \tan\left(\frac{\upsilon}{2}\right) \end{array}
 * = \left\{ \begin{array}{l@{\qquad}c}
 * \sqrt{\frac{1+e}{1-e}}\tan\left(\frac{E}{2}\right)
 * & 0 \leq e < 1 \\ & \\
 * \frac{E}{2} & e = 1 \\ & \\
 * \sqrt{\frac{e+1}{e-1}}\tanh\left(\frac{E}{2}\right) & e > 1
 * \end{array} \right.
 * \f}
 *
 * Thus to solve for \f$\upsilon(t)\f$ one typically inverts the equation for
 * \f$t-t_p\f$ numerically or by series expansion, finds the corresponding
 * \f$E\f$, and then plugs this into the expression for \f$\upsilon\f$.  However,
 * in our case we would then need to do another numerical inversion to
 * find the retarded time \f$t_r\f$ from \eqref{eq_spinorbit-tr}.  A more
 * efficient approach is thus to take an initial guess for \f$E\f$, compute
 * both \f$t\f$, \f$\upsilon\f$, and hence \f$t_r\f$, and then refine directly on
 * \f$E\f$.
 *
 * ### Other notation conventions ###
 *
 * Since we may deal with highly eccentric or open orbits, we will
 * specify these orbits with parameters that are definable for all
 * classes of orbit.  Thus we specify the size of the orbit with the
 * periapsis separation \f$r_p\f$ rather than the semimajor axis \f$a\f$, and the
 * speed of the orbit with the angular speed at periapsis
 * \f$\dot{\upsilon}_p\f$ rather than with the period \f$P\f$.  These parameters
 * are related by:
 * \f{eqnarray}{
 * \label{eq_spinorbit-a}
 * a & = & \frac{r_p}{1-e} \;,\\
 * \label{eq_spinorbit-p}
 * P & = & \frac{2\pi}{\dot{\upsilon}_p} \sqrt{\frac{1+e}{(1-e)^3}} \;.
 * \f}
 * Furthermore, for improved numerical precision when dealing with
 * near-parabolic orbits, we specify the value of \f$1-e\f$ rather than the
 * value of \f$e\f$.  We note that \f$1-e\f$ has a maximum value of \f$1\f$ for a
 * circular orbit, positive for closed elliptical orbits, zero for
 * parabolic orbits, and negative (unbounded) for hyperbolic orbits.
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define GENERATESPINORBITCWH_ENUL 1	/**< Unexpected null pointer in arguments */
#define GENERATESPINORBITCWH_EOUT 2	/**< Output field a, f, phi, or shift already exists */
#define GENERATESPINORBITCWH_EMEM 3	/**< Out of memory */
#define GENERATESPINORBITCWH_EECC 4	/**< Eccentricity out of range */
#define GENERATESPINORBITCWH_EFTL 5	/**< Periapsis motion is faster than light */
#define GENERATESPINORBITCWH_ESGN 6	/**< Sign error: positive parameter expected */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATESPINORBITCWH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATESPINORBITCWH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATESPINORBITCWH_MSGEMEM "Out of memory"
#define GENERATESPINORBITCWH_MSGEECC "Eccentricity out of range"
#define GENERATESPINORBITCWH_MSGEFTL "Periapsis motion is faster than light"
#define GENERATESPINORBITCWH_MSGESGN "Sign error: positive parameter expected"
/** \endcond */

/**
 * This structure stores the parameters for constructing a gravitational
 * waveform with both a Taylor-polynomial intrinsic frequency and phase,
 * and a binary-orbit modulation.  As with the ::PPNParamStruc type
 * in \ref GeneratePPNInspiral_h, we divide the fields into passed
 * fields (which are supplied to the final PulsarCoherentGW structure
 * but not used in any calculations), input fields (that are used by the
 * waveform generator), and output fields (that are set by the waveform
 * generator).
 */
typedef struct tagSpinOrbitCWParamStruc {
  /** \name Passed parameters. */
  /*@{*/
  SkyPosition position;   /**< The location of the source on the sky, normally in equatorial coordinates */
  REAL4 psi;              /**< The polarization angle of the source, in radians */
  /*@}*/

  /** \name Input parameters. */
  /*@{*/
  LIGOTimeGPS epoch;      /**< The start time of the output series */
  LIGOTimeGPS spinEpoch;  /**< A reference time \f$t_\mathrm{ref}\f$ (in the barycentric frame) at which the rotational properties of the source are specified */
  LIGOTimeGPS orbitEpoch; /**< A time \f$t_\mathrm{peri}\f$ (in the barycentric frame) at which the source passes through periapsis.
                           * Note that this is the proper or "true" time of passage; the \e observed periapsis passage occurs at
                           * time \f$t_\mathrm{peri}+r(t_\mathrm{peri})/c\f$ */
  REAL8 deltaT;           /**< The requested sampling interval of the waveform, in s */
  UINT4 length;           /**< The number of samples in the generated waveform */
  REAL4 aPlus, aCross;    /**< The polarization amplitudes \f$A_+\f$, \f$A_\times\f$, in dimensionless strain units */
  REAL8 phi0;             /**< The phase of the wave emitted at time \f$t_\mathrm{ref}\f$, in radians */
  REAL8 f0;               /**< The frequency of the wave emitted at time \f$t_\mathrm{ref}\f$ (and incorporating any Doppler shift due to \f$\dot{R}_0\f$), in Hz */
  REAL8Vector *f;         /**< The spin-normalized Taylor parameters \f$f_k\f$, as defined in \eqref{eq_taylorcw-freq} of \ref GenerateTaylorCW_h.
                           * If \c f=\c NULL, the (proper) spin of the source is assumed to be constant */
  REAL8 omega;            /**< The argument of the periapsis, \f$\omega\f$, in radians */
  REAL8 rPeriNorm;        /**< The projected, speed-of-light-normalized periapsis separation of the orbit, \f$(r_p/c)\sin i\f$, in s */
  REAL8 oneMinusEcc;      /**< The value of \f$1-e\f$ */
  REAL8 angularSpeed;     /**< The angular speed at periapsis, \f$\dot{\upsilon}_p\f$, in Hz */
  /*@}*/

  /** \name Output parameters. */
  /*@{*/
  REAL4 dfdt;             /**< The maximum value of \f$\Delta f\Delta t\f$ encountered over any timestep \f$\Delta t\f$ used in generating the waveform */
  /*@}*/
} SpinOrbitCWParamStruc;


/* ---------- Function prototypes. ---------- */
int XLALGenerateSpinOrbitCW ( PulsarCoherentGW *sourceSignal, SpinOrbitCWParamStruc *sourceParams );

void
LALGenerateSpinOrbitCW( LALStatus             *,
			PulsarCoherentGW            *output,
			SpinOrbitCWParamStruc *params );

void
LALGenerateEllipticSpinOrbitCW( LALStatus             *,
				PulsarCoherentGW            *output,
				SpinOrbitCWParamStruc *params );

void
LALGenerateParabolicSpinOrbitCW( LALStatus             *,
				 PulsarCoherentGW            *output,
				 SpinOrbitCWParamStruc *params );

void
LALGenerateHyperbolicSpinOrbitCW( LALStatus             *,
				  PulsarCoherentGW            *output,
				  SpinOrbitCWParamStruc *params );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATESPINORBITCW_H */
