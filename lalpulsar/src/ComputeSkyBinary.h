/*
*  Copyright (C) 2007 Chris Messenger, Jolien Creighton
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

/**
\author Messenger, C.J., Berukoff, S.J., Papa, M.A.
\file
\ingroup pulsarTODO

\brief Computes phase coefficients necessary for a correct demodulation for a source in a binary system.

\heading{Synopsis}
\code
#include <lal/ComputeSkyBinary.h>
\endcode

The methods employed here follow very closely those used within <tt>ComputeSky()</tt>.
Note that at present this code simply corrects for the Doppler modulation present in a polynomial
frequency function for signals from sources in elliptical orbits.  It does not account for general
relativistic effects.

At the risk of repeating existing documentation, but in the interests of clarity much of the
following can also be found in the <tt>ComputeSky()</tt> documentation.  Recall that a demodulated
Fourier Transform (DeFT) is given by

\anchor eqB_demod_FT \f{equation}{
\hat{x}_b({\vec{\lambda}})=
\sum_{\alpha =0}^{M-1}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}\left[\frac{1}{N}\sum_{j=0}^{N-1}e^{-2\pi i(\Phi_{\alpha jb}(\vec{\lambda})-\frac{jk}{N})}\right]
\tag{eqB_demod_FT}
\f}

The index \f$b\f$ defines the DeFT frequency bin, the index \f$\alpha\f$ loops through
the SFTs that build the DeFT, \f$k\f$ runs on all the SFT frequency bins, and \f$j\f$
is a time index that runs on each SFT.  As shown in
\ref LALDemod.h, the next step in the development of the demodulation
technique involves Taylor expanding the phase model about the temporal
midpoint of each short segment of data, while retaining only first order
terms.  At this point it is neccessary to clearly define some quantities.
Times as defined at the chosen detector are denoted by \f$t\f$, times defined at the
solar system barycenter (SSB) are denoted by \f$T\f$, and the retarded time measured at an inertial
reference point (chosen as the SSB) at a distance from the source are denote by \f$t^{\prime}\f$.

The Taylor expansion of \f$\Phi (t)\f$ about the temporal midpoint
\f$t_{\alpha,1/2}\f$ is

\anchor eqB_taylor2 \f{equation}{
\Phi_{\alpha}(t) = \Phi(t_{\alpha,1/2})+\left[t-t_{\alpha,1/2}\right]\frac{d\Phi}{dt}(t_{\alpha,1/2})\tag{eqB_taylor2} \\
\f}

For each value of \f$\alpha\f$, this expression consists of either constant or linear terms in time.
With the particular time discretization chosen in this code, \f$t=t_{0}+(N\alpha+j)\ T_{obs}/NM\f$, we have

\anchor eqB_time \f{equation}{
\tag{eqB_time}
\left[t-t_{\alpha,1/2}\right]=\frac{\ T_{obs}}{M}\left(\frac{j}{N}-\frac{1}{2}\right)=\mathcal{T}_{s}\left(\frac{j}{N}-\frac{1}{2}\right),
\f}

where \f$\mathcal{T}_{s}\f$ is the short time baseline of the \f$M\f$ short FTs.  On
the other hand, the phase can also be expressed as a function of SSB time \f$T\f$
(i.e. the time at the solar system barycenter).  We will assume the source to
be at rest in this reference frame.  If we now adopt the notation \f$\Delta
t^{\prime}_{\alpha}\equiv\left[t^{\prime}(T(t_{\alpha,1/2}))-
t^{\prime}(T(t_{0}))\right]\f$ and \f$\dot{t^{\prime}}_{\alpha}\equiv dt^{\prime}/dt({\alpha,1/2})\f$,
the phase terms described in Eq.\eqref{eqB_taylor2} become (neglecting constants)

\anchor eqB_phi \anchor eqB_dphi \f{eqnarray}{
\Phi(t_{\alpha,1/2})   & = & f_{0}(\Delta t^{\prime}_{\alpha})+\frac{1}{2}f_{1}(\Delta t^{\prime}_{\alpha})^{2}
+\frac{1}{3}f_{2}(\Delta t^{\prime}_{\alpha})^{3}+\frac{1}{4}f_{3}(\Delta t^{\prime}_{\alpha})^{4}+\frac{1}{5}f_{4}(\Delta t^{\prime}_{\alpha})^{5}
+\frac{1}{6}f_{5}(\Delta t^{\prime}_{\alpha})^{6}, \nonumber\tag{eqB_phi} \\
                                         &   & \\
\frac{d\Phi}{dt}(t_{\alpha,1/2})         & = & \dot{t^{\prime}}_{\alpha}\left(f_{0}+ f_{1}(\Delta t^{\prime}_{\alpha})
+f_{2}(\Delta t^{\prime}_{\alpha})^{2}+f_{3}(\Delta t^{\prime}_{\alpha})^{3}
+f_{4}(\Delta t^{\prime}_{\alpha})^{4}+f_{5}(\Delta t^{\prime}_{\alpha})^{5}\right). \tag{eqB_dphi}
\f}

Note that the polynomial phase function is expressed as a function of the retarded time \f$t^{\prime}\f$ and subsequently
the intrinsic frequency and its derivitives (\f$f_{i}\f$) are defined at the chosen inertial reference frame (SSB).

In order to calculate, for each value of \f$\alpha\f$, the quantities \f$\dot{t^{\prime}}_{\alpha}\f$ and
\f$\Delta t^{\prime}_{\alpha}\f$, we must now look at the binary system in more detail.  At present
we reference \ref GenerateSpinOrbitCW.h, where the definition of all the following orbital variables
can be found.  For a given set of orbital input parameters we obtain the eccentric anomoly \f$E\f$ by
numerically solving

\anchor eqB_bintime \f{eqnarray}{\tag{eqB_bintime}
  T(t_{\alpha})-T_{p}&=&\frac{P}{2\pi}\left[E+\left(p\sin{E}+q\left(\cos{E}-1\right)\right)\right].
\f}

where the quantities \f$p\f$ and \f$q\f$, dependent only on the orbital parameters
of the source system, are given by

\anchor eqB_pq \f{eqnarray}{
\tag{eqB_pq}
p&=&\frac{2\pi}{P}\frac{a}{c}\sin{i}\sqrt{1-e^{2}}\cos{\omega}-e \nonumber \\
q&=&\frac{2\pi}{P}\frac{a}{c}\sin{i}\sin{\omega}.
\f}

\f$T(t_{\alpha})\f$ is returned via a call to \c LALBarycenter and \f$a\sin{i},P,T_{p},\omega,e\f$
are the projected semi-major axis (projected along the line of sight), the orbital period, the
time of observed periapse passage as measured
in the SSB, the argument of periapse, and the orbital eccentricity respectively.  Having defined \f$E\f$
(where the source is in it`s orbit) at a given detector time \f$t_{\alpha}\f$ we can calculate the derivitive
of the retarded source time \f$\prime{t}\f$ with respect to the SSB time \f$T\f$.  This is given by

\anchor eqB_binderiv \f{equation}{\tag{eqB_binderiv}
 \frac{dt^{\prime}}{dT}=\frac{[1-e\cos{E}]}{\left[1+pcos{E}-q\sin{E}\right]}.
\f}

The quantity \f$\dot{t^{\prime}}_{\alpha}\f$ can now be expressed as

\anchor eqB_binderivtwo \f{equation}{\tag{eqB_binderivtwo}
\dot{t^{\prime}}_{\alpha}=\frac{dT}{dt}\frac{dt^{\prime}}{dT},
\f}

where \f$dT/dt\f$ is returned via a call to <tt>LALBarycenter()</tt>.

We can now rewrite Eq.\eqref{eqB_taylor2} and by grouping together the terms in \f$j\f$ (linear
in \f$t\f$) in order to save computations, we have

\anchor eqB_phasecalc \f{equation}{
\Phi_{\alpha}(t)=\sum_{s=0}^{n_{spin}}f_{s}A_{s\alpha}+\frac{j}{N}\sum_{s=0}^{n_{spin}}f_{s}B_{s\alpha},
\tag{eqB_phasecalc}
\f}

where \f$n_{spin}\f$ is the maximum order of spindown parameter.

Thus, for a given sky position and set of orbital parameters, the quantities \f$\dot{t^{\prime}}_{\alpha}\f$ and
\f$\Delta t^{\prime}_{\alpha}\f$ are calculated only once, just as in <tt>ComputeSky()</tt>.  The analytical constants
defined in Eq.\eqref{eqB_phasecalc} now become

\f{equation}{
A_{s \alpha}=\frac{1}{s+1}\Delta (t^{\prime}_{\alpha})^{s+1}-\frac{1}{2}\mathcal{T}_{SFT}\dot{t^{\prime}}_{\alpha}\Delta (t^{\prime}_{\alpha})^{s}
\f}
\f{equation}{
B_{s \alpha}=\mathcal{T}_{SFT}\dot{t^{\prime}}_{\alpha}\Delta (t^{\prime}_{\alpha})^{s}.
\f}

It is these constants that form the input to the function <tt>LALDemod()</tt>.

*/

#ifndef _COMPUTESKYBINARY_H
#define _COMPUTESKYBINARY_H

#include <lal/LALStdlib.h>
#include <lal/LALBarycenter.h>
#include <lal/Date.h>

#ifdef __cplusplus
extern "C" {
#endif

/**\name Error Codes */ /*@{*/
#define COMPUTESKYBINARYH_ENULL 1
#define COMPUTESKYBINARYH_ENNUL 2
#define COMPUTESKYBINARYH_ERANG 3
#define COMPUTESKYBINARYH_ENEGA 4
#define COMPUTESKYBINARYH_MSGENULL "Null Pointer"
#define COMPUTESKYBINARYH_MSGENNUL "Non-Null Pointer"
#define COMPUTESKYBINARYH_MSGERANG "Input parameter out of range"
#define COMPUTESKYBINARYH_MSGENEGA "Bad Negative Value"
/*@}*/

/** The following quantity represents the required accuracy for the timing of the binary orbit in seconds */
#define ACC 1e-9

  /** This structure contains the parameters for the LALComputeSkyBinary() routine */
typedef struct
tagCSBParams
{
  INT8		spinDwnOrder;	/**< The maximal number of spindown parameters per spindown parameter set */
  INT8		mObsSFT;	/**< The number of SFTs in the observation time */
  REAL8		tSFT;		/**< The timescale of one SFT */
  LIGOTimeGPS	*tGPS;		/**< An array containing the GPS times of the first datum from each SFT */
  REAL8		*skyPos; 	/**< The array containing the sky patch coordinates */
  REAL8		SemiMajorAxis;  /**< The projected speed-of-light-normalised semi-major axis of the orbit (in seconds) */
  REAL8		OrbitalPeriod;  /**< The orbital period (in seconds) */
  REAL8         OrbitalEccentricity; /**<  The orbital eccentricity */
  REAL8         ArgPeriapse;    /**< The argument of the periapse (in radians) */
  LIGOTimeGPS   TperiapseSSB;   /**< A time of observed periapse passage as defined in the SSB */
  BarycenterInput *baryinput;	/**< TO BE DOCUMENTED */
  EmissionTime  *emit;		/**< TO BE DOCUMENTED */
  EarthState    *earth;		/**< TO BE DOCUMENTED */
  const EphemerisData *edat;	/**< ephemeris data */
}
CSBParams;


void LALComputeSkyBinary (LALStatus *status,
			REAL8 		*skyConst,
			INT8 		iSkyCoh,
			CSBParams 	*params);


#ifdef __cplusplus
}
#endif

#endif /* _COMPUTESKYBINARY_H */
