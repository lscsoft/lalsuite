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

#ifndef _GENERATETAYLORCW_H
#define _GENERATETAYLORCW_H

#include <lal/LALStdlib.h>
#include <lal/PulsarSimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup GenerateTaylorCW_h
 * \author Creighton, T. D.
 *
 * \brief Provides routines to generate Taylor-parameterized continuous waveforms.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/GenerateTaylorCW.h>
 * \endcode
 *
 * This header covers routines to generate continuous quasiperiodic
 * waveforms whose frequency varies slowly and smoothly with time.  For
 * such sources the frequency function is normally described by its
 * Taylor "spindown" (or spin-up) coefficients.  This type of waveform
 * may be typical of objects such as neutron stars that are gradually
 * shedding angular momentum, or are accelerating in the gravitational
 * potential of a star cluster.  The Taylor expansion is likely
 * \e not suitable for long-term modelling of the frequency of waves
 * from glitching neutron stars, neutron stars in close binary orbits, or
 * neutron stars that are accreting or shedding angular momentum in a
 * stochastic manner.
 *
 * The frequency and phase of such slowly-varying quasiperiodic sources
 * are given by their Taylor series:
 * \anchor eq_taylorcw-freq \anchor eq_taylorcw-phi \f{eqnarray}{
 * \tag{eq_taylorcw-freq}
 * f(t)    & = & f_0 \left[ 1 + \sum_{k=1}^n f_k(t-t_0)^k \right] \;, \\
 * \tag{eq_taylorcw-phi}
 * \phi(t) & = & \phi_0 + 2\pi f_0 \left[ (t-t_0) +
 * \sum_{k=1}^n \frac{f_k}{k+1}(t-t_0)^{k+1} \right] \;,
 * \f}
 * where \f$f_k\f$ are the spin-normalized Taylor coefficients.  If the
 * source's spin is varying over some timescale \f$\tau\f$, one typically
 * expects that \f$f_k\sim\tau^{-k}\f$.  Note that in this and later
 * discussions, \f$f\f$ and \f$\phi\f$ refer to the frequency and phase of the
 * gravitational wave, which are typically some constant multiple of
 * (often twice) the frequency and phase of the rotating source.
 *
 * The \c PulsarCoherentGW structure allows for a very general
 * description of waveforms with modulations in the amplitudes or
 * relative phases of the wave polarizations, as described in
 * \ref PulsarSimulateCoherentGW_h.  However, in this simplest model of
 * quasiperiodic waveforms, we neglect such phenomena as precession that
 * would produce these effects.  Thus for any given source one can choose
 * a polarization basis (described by some polarization angle \f$\psi\f$) in
 * which the wave has a constant elliptical polarization of the form:
 * \anchor eq_taylorcw-hplus \anchor eq_taylorcw-hcross \f{eqnarray}{
 * \tag{eq_taylorcw-hplus}
 * h_+(t)      & = & A_+      \cos\phi(t) \;, \\
 * \tag{eq_taylorcw-hcross}
 * h_\times(t) & = & A_\times \sin\phi(t) \;.
 * \f}
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define GENERATETAYLORCWH_ENUL 1	/**< Unexpected null pointer in arguments */
#define GENERATETAYLORCWH_EOUT 2	/**< Output field a, f, phi, or shift already exists */
#define GENERATETAYLORCWH_EMEM 3	/**< Out of memory */
/*@}*/

/** \cond DONT_DOXYGEN */
#define GENERATETAYLORCWH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATETAYLORCWH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATETAYLORCWH_MSGEMEM "Out of memory"
/** \endcond */

/**
 * This structure stores the parameters for constructing a gravitational
 * waveform with a Taylor-polynomial frequency and phase.  As with the
 * \c PPNParamStruc type in \ref GeneratePPNInspiral_h, we divide
 * the fields into passed fields (which are supplied to the final
 * \c PulsarCoherentGW structure but not used in any calculations), input
 * fields (that are used by the waveform generator), and output fields
 * (that are set by the waveform generator).
 */
typedef struct tagTaylorCWParamStruc {
  /** \name Passed parameters. */
  /*@{*/
  SkyPosition position; /**< The location of the source on the sky, normally in equatorial coordinates */
  REAL4 psi;            /**< The polarization angle of the source, in radians */
  LIGOTimeGPS epoch;    /**< The start time \f$t_0\f$ of the output series */
  /*@}*/

  /** \name Input parameters. */
  /*@{*/
  REAL8 deltaT;         /**< The requested sampling interval of the waveform, in s */
  UINT4 length;         /**< The number of samples in the generated waveform */
  REAL4 aPlus, aCross;  /**<  The polarization amplitudes \f$A_+\f$, \f$A_\times\f$, in dimensionless strain units */
  REAL8 phi0;           /**< The wave phase at time \f$t_0\f$, in radians */
  REAL8 f0;             /**< The wave frequency at time \f$t_0\f$, in Hz */
  REAL8Vector *f;       /**< The spin-normalized Taylor parameters \f$f_k\f$, as defined in Eq.\eqref{eq_taylorcw-freq}; If \c f=\c NULL, a monochromatic wave is generated */
  /*@}*/

  /** \name Output parameters. */
  /*@{*/
  REAL4 dfdt;           /**< The maximum value of \f$\Delta f\Delta t\f$ encountered over any timestep \f$\Delta t\f$ used in generating the waveform */
  /*@}*/
} TaylorCWParamStruc;


/* ---------- Function prototypes. ---------- */
void
LALGenerateTaylorCW( LALStatus          *,
		     PulsarCoherentGW         *output,
		     TaylorCWParamStruc *params );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATETAYLORCW_H */
