/*
*  Copyright (C) 2007 Jolien Creighton, Steve Berukoff
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
 * \author Berukoff, S.J., Papa, M.A.
 * \file
 * \ingroup pulsarTODO
 *
 * \brief Computes phase coefficients necessary for a correct demodulation.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/ComputeSky.h>
 * \endcode
 *
 * This is a short summary of the analytical calculations which form the basis for the code in this routine.
 *
 * Recall that a demodulated Fourier Transform (DeFT) is given by
 * \f{equation}{
 * \label{eq_e4a}
 * \hat{x}_b({\vec{\lambda}})=
 * \sum_{\alpha =0}^{M-1}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}\left[\frac{1}{N}\sum_{j=0}^{N-1}e^{-2\pi i(\Phi_{\alpha jb}(\vec{\lambda})-\frac{jk}{N})}\right]
 * \f}
 * The index \f$b\f$ defines the DeFT frequency bin, the index \f$\alpha\f$ loops through
 * the SFTs that build the DeFT, \f$k\f$ runs on all the SFT frequency bins, and \f$j\f$
 * is a time index that runs on each SFT.  The next step in the development of the demodulation
 * technique involves Taylor expanding the phase model about the temporal
 * midpoint of each short segment of data, while retaining only first order
 * terms.  The Taylor expansion of \f$\Phi (t)\f$ about the temporal midpoint
 * \f$t_{\alpha,1/2}\f$ is
 * \f{equation}{
 * \label{eq_taylor2}
 * \Phi_{\alpha}(t) = \Phi(t_{\alpha,1/2})+\left[t-t_{\alpha,1/2}\right]\frac{d\Phi}{dt}(t_{\alpha,1/2}) \\
 * \f}
 * For each value of \f$\alpha\f$, this expression consist of either constant or linear terms in time.  With the particular time discretization chosen in this code, \f$t=t_{0}+(N\alpha+j)\ T_{obs}/NM\f$, we have
 * \f{equation}{
 * \label{eq_time}
 * \left[t-t_{\alpha,1/2}\right]=\frac{\ T_{obs}}{M}\left(\frac{j}{N}-\frac{1}{2}\right)=\mathcal{T}_{s}\left(\frac{j}{N}-\frac{1}{2}\right),
 * \f}
 * where \f$\mathcal{T}_{s}\f$ is the short time baseline of the \f$M\f$ short FTs.  On
 * the other hand, the phase can also be expressed as a function of SSB time \f$T\f$
 * (i.e. the time at the solar system barycenter).  We will assume the source to
 * be at rest in this reference frame.  Now, if one adopts the notation \f$\Delta
 * T_{\alpha}\equiv\left[T(t_{\alpha,1/2})-
 * T(t_{0})\right]\f$ and \f$\dot{T}_{\alpha}\equiv
 * dT/dt(t_{\alpha,1/2})\f$
 * the phase terms in the above equation are (neglecting constants)
 * \f{eqnarray}{
 * \label{eq_phi}
 * \Phi(t_{\alpha,1/2})                     & = & f_{0}\Delta T_{\alpha}+\frac{1}{2}f_{1}\Delta T_{\alpha}^{2}
 * +\frac{1}{3}f_{2}\Delta T_{\alpha}^{3}+\frac{1}{4}f_{3}\Delta T_{\alpha}^{4}+\frac{1}{5}f_{4}\Delta T_{\alpha}^{5}
 * +\frac{1}{6}f_{5}\Delta T_{\alpha}^{6} \\
 * &   & \\
 * \label{eq_dphi}
 * \frac{d\Phi}{dt}(t_{\alpha,1/2})         & = & \dot{T}_{\alpha}\left(f_{0}+ f_{1}\Delta T_{\alpha}
 * +f_{2}\Delta T_{\alpha}^{2}+f_{3}\Delta T_{\alpha}^{3}
 * +f_{4}\Delta T_{\alpha}^{4}+f_{5}\Delta T_{\alpha}^{5}\right).
 * \f}
 * These constants, for each value of \f$\alpha\f$, require \f$\dot{T}_{\alpha}\f$ and
 * \f$\Delta T_{\alpha}\f$, which are calculated by a suitable timing routine.  For
 * this demodulation package, this timing routine is provided by <tt>tdb()</tt>.
 * Thus, for a given sky position, the timing routine will be called once for
 * each short time chunk, each call returning a specific  \f$\dot{T}_{\alpha}\f$ and
 * \f$\Delta T_{\alpha}\f$.  By substituting \eqref{eq_time}, \eqref{eq_phi} and
 * \eqref{eq_dphi} in \eqref{eq_taylor2} and grouping together the terms in \f$j\f$ (linear
 * in \f$t\f$) in order to save computations, we have
 * \f{equation}{
 * \label{eq_phasecalc}
 * \Phi_{\alpha}(t)=\sum_{s=0}^{n_{spin}}f_{s}A_{s\alpha}+\frac{j}{N}\sum_{s=0}^{n_{spin}}f_{s}B_{s\alpha},
 * \f}
 * where \f$n_{spin}\f$ is the maximum order of spindown parameter.  Rather than
 * store the values of \f$\dot{T}_{\alpha}\f$ and \f$\Delta T_{\alpha}\f$ for each value
 * of \f$\alpha\f$, it is more efficient to calculate the constants \f$A_{s\alpha}\f$ and
 * \f$B_{s\alpha}\f$ only once, and then use these values for every spindown
 * parameter set used when searching in a given sky position.  Analytical
 * formulae for these constants are easily derived:
 * \f{equation}{
 * A_{s \alpha}=\frac{1}{s+1}\Delta T_{\alpha}^{s+1}-\frac{1}{2}\mathcal{T}_{SFT}\dot{T}_{\alpha}\Delta T_{\alpha}^{s}
 * \f}
 * \f{equation}{
 * B_{s \alpha}=\mathcal{T}_{SFT}\dot{T}_{\alpha}\Delta T_{\alpha}^{s}
 * \f}
 *
 */




#ifndef _COMPUTESKY_H
#define _COMPUTESKY_H

#include <lal/LALStdlib.h>
#include <lal/LALBarycenter.h>

#ifdef __cplusplus
extern "C" {
#endif

/**\name Error Codes */ /*@{*/
#define COMPUTESKYH_ENULL 1
#define COMPUTESKYH_ENNUL 2
#define COMPUTESKYH_ENEGA 4
#define COMPUTESKYH_MSGENULL "Null Pointer"
#define COMPUTESKYH_MSGENNUL "Non-Null Pointer"
#define COMPUTESKYH_MSGENEGA "Bad Negative Value"
/*@}*/

/**
 * This structure contains the parameters for the LALComputeSky() routine.
 */
typedef struct
tagCSParams
{
  INT8			spinDwnOrder;	/**< The maximal number of spindown parameters per spindown parameter set */
  INT8			mObsSFT;	/**< The number of SFTs in the observation time */
  REAL8			tSFT;		/**< The timescale of one SFT */
  LIGOTimeGPS		*tGPS;		/**< An array containing the GPS times of the first datum from each SFT */
  REAL8 		*skyPos; 	/**< The array containing the sky patch coordinates */
  BarycenterInput 	*baryinput;	/**< A switch which turns modulation on/off */
  EmissionTime 		*emit;		/**< TO BE DOCUMENTED */
  EarthState 		*earth;		/**< TO BE DOCUMENTED */
  const EphemerisData 	*edat;		/**< ephemeris data */
}
CSParams;


void LALComputeSky (LALStatus *status,
			REAL8 		*skyConst,
			INT8 		iSkyCoh,
			CSParams 	*params);


#ifdef __cplusplus
}
#endif

#endif /* _COMPUTESKY_H */
