/*
 * Copyright (C) 2008 J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, E. Ochsner, C. Pankow, 2104 A. Klein
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

#ifndef _LALSIMINSPIRAL_H
#define _LALSIMINSPIRAL_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimSphHarmSeries.h>
#include <lal/LALSimInspiralTestGRParams.h>
#include <lal/LALSimInspiralWaveformFlags.h>
#include <gsl/gsl_matrix.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @addtogroup lalsimulation_inspiral
 * @details
 *
 * Various types of inspiral approximants are supported for producing
 * waveforms in the time- or frequency-domains.  The high-level routines
 * for generating simulated inspiral waveforms are
 * XLALSimInspiralChooseTDWaveform() (for time-domain waveforms) and
 * XLALSimInspiralChooseFDWaveform() (for frequency-domain waveforms).
 * The following examples show their basic usage.
 *
 * Generate a time-domain waveform:
 * @code
 * #include <lal/LALSimInspiral.h>
 * #include <lal/LALConstants.h>
 * ...
 * double phiRef = 0;             // gravitational wave phase at end
 * double deltaT = 1.0/16384.0;   // series sampling interval
 * double m1 = 1.4 * LAL_MSUN_SI; // mass of body 1
 * double m2 = 1.4 * LAL_MSUN_SI; // mass of body 2
 * double S1x = 0.0;              // x-component of dimensionless spin of body 1
 * double S1y = 0.0;              // y-component of dimensionless spin of body 1
 * double S1z = 0.0;              // z-component of dimensionless spin of body 1
 * double S2x = 0.0;              // x-component of dimensionless spin of body 2
 * double S2y = 0.0;              // y-component of dimensionless spin of body 2
 * double S2z = 0.0;              // z-component of dimensionless spin of body 2
 * double f_min = 40.0;           // start frequency of inspiral
 * double f_ref = 0.0;            // reference frequency: 0 means waveform end
 * double r = 1e6 * LAL_PC_SI;    // distance
 * double i = 0.0;                // inclination
 * double lambda1 = 0.0;          // dimensionless tidal parameter of body 1
 * double lambda2 = 0.0;          // dimensionless tidal parameter of body 2
 * LALSimInspiralWaveformFlags *waveFlags = NULL;  // no extra flags
 * LALSimInspiralTestGRParam *nonGRparams = NULL;  // no non-GR parameters
 * int amplitudeO = -1;           // amplitude pN order: -1 means include all
 * int phaseO = -1;               // phase pN order: -1 means include all
 * Approximant approximant = TaylorT2;  // post-Newtonian approximant
 * REAL8TimeSeries *hplus = NULL;  // plus polarization to be returned
 * REAL8TimeSeries *hcross = NULL; // cross polarization to be returned
 * ...
 * XLALSimInspiralChooseTDWaveform(&hplus, &hcross, phiRef, deltaT, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_ref, r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
 * @endcode
 *
 * Generate a frequency-domain waveform:
 * @code
 * #include <lal/LALSimInspiral.h>
 * #include <lal/LALConstants.h>
 * ...
 * double phiRef = 0;             // gravitational wave phase at end
 * double deltaF = 1.0/64.0;      // series sampling interval
 * double m1 = 1.4 * LAL_MSUN_SI; // mass of body 1
 * double m2 = 1.4 * LAL_MSUN_SI; // mass of body 2
 * double S1x = 0.0;              // x-component of dimensionless spin of body 1
 * double S1y = 0.0;              // y-component of dimensionless spin of body 1
 * double S1z = 0.0;              // z-component of dimensionless spin of body 1
 * double S2x = 0.0;              // x-component of dimensionless spin of body 2
 * double S2y = 0.0;              // y-component of dimensionless spin of body 2
 * double S2z = 0.0;              // z-component of dimensionless spin of body 2
 * double f_min = 40.0;           // start frequency of inspiral
 * double f_max = 0.0;            // end frequency of inspiral: 0 means use default
 * double f_ref = 0.0;            // reference frequency: 0 means waveform end
 * double r = 1e6 * LAL_PC_SI;    // distance
 * double i = 0.0;                // inclination
 * double lambda1 = 0.0;          // dimensionless tidal parameter of body 1
 * double lambda2 = 0.0;          // dimensionless tidal parameter of body 2
 * LALSimInspiralWaveformFlags *waveFlags = NULL;  // no extra flags
 * LALSimInspiralTestGRParam *nonGRparams = NULL;  // no non-GR parameters
 * int amplitudeO = -1;           // amplitude pN order: -1 means include all
 * int phaseO = -1;               // phase pN order: -1 means include all
 * Approximant approximant = TaylorF2;  // post-Newtonian approximant
 * COMPLEX16FrequencySeries *hptilde = NULL;  // plus polarization to be returned
 * COMPLEX16FrequencySeries *hctilde = NULL; // cross polarization to be returned
 * ...
 * XLALSimInspiralChooseFDWaveform(&hptilde, &hctilde, phiRef, deltaF, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max, f_ref, r, i, lambda1, lambda2, waveFlags, nonGRparams, amplitudeO, phaseO, approximant);
 * @endcode
 *
 * ### Coordinate Systems
 *
 * The diagram below illustrates how the source frame (x,y,z) of the binary is
 * related to the wave frame (X,Y,Z) in which the gravitational waveform is
 * defined.
 *
 * @anchor lalsiminspiral_orbitelements
 * @image html lalsiminspiral_orbitelements.svg "Orbital Elements"
 *
 * The origin of the coordinate systems is the instantaneous center-of-mass
 * of the binary system.  The orbiting body shown in the diagram is body 1.
 *
 * The binary's instantaneous orbital angular momentum @b L at the
 * reference gravitational wave frequency @p f_ref defines the z-axis of the
 * binary system, while the direction to the periapsis defines the x-axis of
 * the binary system.  The x-y-plane is therefore the orbital plane, at least
 * at the moment the binary system is at the reference gravitational wave
 * frequency.
 *
 * The spin components for body 1, (@p S1x,@p S1y, @p S1z), and for body 2,
 * (@p S2x,@p S2y, @p S2z), are defined in the source-frame.  Therefore,
 * when the spins are aligned with the orbital angular momentum,
 * @p S1x = @p S1y = @p S2x = @p S2y = 0.
 *
 * The wave frame is defined by the Z-axis, which points toward the Earth,
 * and some reference direction, defining the X-axis.  The X-Y-plane is
 * therefore the plane of the sky.
 *
 * The plus- and cross-polarizations of the gravitational waveform are defined
 * in this wave frame.  Specifically, if \f$ h^{ij} \f$ is computed in the
 * source frame, then
 * \f[ h_+ = \frac12 ( \hat{p}_i \hat{p}_j - \hat{q}_i \hat{q}_j ) h^{ij} \f]
 * and
 * \f[ h_\times = \frac12 ( \hat{p}_i \hat{q}_j + \hat{q}_i \hat{p}_j ) h^{ij} \f]
 * where \f$ \hat{p}_i \f$ are the components of the unit vector pointing
 * along the X-axis and \f$ \hat{q}_i \f$ are the components of the unit
 * vector pointing along the Y-axis.
 *
 * The orbital elements are:
 *
 *  * Inclination (&iota;).  The angle between the Z-axis of the wave frame
 *    and the z-axis of the source frame.
 *  * Longitude of ascending node (&Omega;).  The angle on the plane of the
 *    sky from the X-axis of the reference direction in the wave frame to the
 *    ascending node @htmlonly &#x260A; @endhtmlonly.
 *    @note This angle is entirely degenerate with the polarization angle &psi;.
 *  * Argument of pariapsis (&omega;).  The angle on the orbital plane from
 *    the ascending node @htmlonly &#x260A; @endhtmlonly to the x-axis in the
 *    source frame.
 *  * True anomaly (&phi;).  The angle along the orbital plane from the
 *    periapsis to the present position of the orbiting body (body 1).
 *    The reference phase @p phiRef is @e twice the true anomaly of body 1
 *    at the moment when the system reaches the gravitational wave frequency
 *    @p f_ref which is @e twice the orbital frequency.
 *
 * @attention
 * At present, eccentric orbits are not fully supported, and the x-axis
 * of the source frame is defined to be the ascending node
 * @htmlonly &#x260A; @endhtmlonly.  Therefore, &omega;=0 by definition.
 *
 * @attention
 * In the present implementation, the reference direction in the wave frame,
 * i.e., the X-axis, is defined to be the ascending node
 * @htmlonly &#x260A; @endhtmlonly.  Therefore, &Omega;=0 by definition.  At
 * present, then, the X-axis and the x-axis coincide.
 *
 * @sa
 * The coordinate systems used here follow those of
 * > Clifford M. Will and Alan G. Wiseman
 * > "Gravitational radiation from compact binary systems: Gravitational
 * > waveforms and energy loss to second post-Newtonian order"
 * > Phys. Rev. D @b 54, 4813 (1996)
 * > http://dx.doi.org/10.1103/PhysRevD.54.4813
 *
 */

/**
 * @defgroup LALSimInspiral_h Header LALSimInspiral.h
 * @ingroup lalsimulation_inspiral
 *
 * @brief Routines for generating binary inspiral gravitational waveforms.
 *
 * @{
 * @defgroup LALSimInspiral_c                      Module LALSimInspiral.c
 * @defgroup LALSimInspiralPNMode_c                Module LALSimInspiralPNMode.c
 * @defgroup LALSimInspiralTaylorXX_c              Module LALSimInspiralTaylorXX.c
 * @defgroup LALSimInspiralSpinTaylor_c            Module LALSimInspiralSpinTaylor.c
 * @defgroup LALSimInspiralEccentricTD_c           Module LALSimInspiralEccentricTD.c
 * @defgroup LALSimInspiralEccentricityFD_c        Module LALSimInspiralEccentricityFD.c
 * @defgroup LALSimInspiralSpinDominatedWaveform_c Module LALSimInspiralSpinDominatedWaveform.c
 * @defgroup LALSimInspiralTaylorF2ReducedSpin_c   Module LALSimInspiralTaylorF2ReducedSpin.c
 * @defgroup LALSimInspiralHGimri_c                Module LALSimInspiralHGimri.c
 * @defgroup LALSimInspiralWaveformFlags_c         Module LALSimInspiralWaveformFlags.c
 * @defgroup LALSimInspiralTestGRParams_c          Module LALSimInspiralTestGRParams.c
 * @defgroup LALSimInspiralWaveformTaper_c         Module LALSimInspiralWaveformTaper.c
 * @}
 *
 * @addtogroup LALSimInspiral_h
 * @{
 */

#define LAL_PN_MODE_L_MAX 3
/* (2x) Highest available PN order - UPDATE IF NEW ORDERS ADDED!!*/
#define LAL_MAX_PN_ORDER 8

/**
 * Enum that specifies the PN approximant to be used in computing the waveform.
 */
typedef enum {
   TaylorT1, 		/**< Time domain Taylor approximant in which the energy and flux are both kept
                         * as Taylor expansions and a first order ordinary differential equation is solved
                         * or the GW phase as a function of \f$t\f$; Outputs a time-domain wave.
                         * @remarks Implemented in lalsimulation (time domain).
                         */
   TaylorT2,		/**< Time domain Taylor approximant in which the phase evolution \f$\varphi(t)\f$ is
                         * obtained by iteratively solving post-Newtonian expansions \f$\varphi(v)\f$ and \f$t(v)\f$;
                         * Outputs a time-domain wave.
                         * @remarks Implemented in lalsimulation (time domain).
                         */
   TaylorT3,		/**< Time domain Taylor approximant in which phase is explicitly given as a function
                         * of time; outputs a time-domain wave.
                         * @remarks Implemented in lalsimulation (time domain).
                         */
   TaylorF1,		/**< The stationary phase approximation that correctly represents, in the Fourier domain,
                         * the waveform given by \c TaylorT1 approximant (see \cite dis2000 for details);
                         * Outputs a frequency-domain wave.
                         * @attention Not implemented in lalsimulation. */
   EccentricFD,         /**< Frequency domain waveform in the SPA to describe low eccentricity systems.
                         * @remarks Implemented in lalsimulation (frequency domain). */
   TaylorF2,		/**< The standard stationary phase approximation; Outputs a frequency-domain wave.
                         * @remarks Implemented in lalsimulation (frequency domain). */
   TaylorR2F4,		/**< A frequency domain model closely related to TaylorT4
                         * @attention Not implemented in lalsimulation. */
   TaylorF2RedSpin,		/**< TaylorF2 waveforms for non-precessing spins, defined in terms of a single (reduced-spin) parameter [Ajith_2011ec].
                                 * @remarks Implemented in lalsimulation (frequency domain). */
   TaylorF2RedSpinTidal,	/**< TaylorF2 waveforms for non-precessing spins, defined in terms of a single (reduced-spin) parameter [Ajith_2011ec] plus tidal terms (http://arxiv.org/abs/1101.1673).
                                 * @remarks Implemented in lalsimulation (frequency domain).  */
   PadeT1,		/**< Time-domain P-approximant; Outputs a time-domain wave.
                         * @attention Not implemented in lalsimulation. */
   PadeF1,		/**< Frequency-domain P-approximant (not yet implemented).
                         * @attention Not implemented in lalsimulation. */
   EOB,			/**< Effective one-body waveform; Outputs a time-domain wave.
                         * @attention Not implemented in lalsimulation. */
   BCV,			/**< Detection template family of Buonanno, Chen and Vallisneri \cite BCV03; Outputs a frequency-domain wave.
                         * @attention Not implemented in lalsimulation. */
   BCVSpin,		/**< Detection template family of Buonanno, Chen and Vallisneri including  spin effects \cite BCV03b; Outputs a frequency-domain wave.
                         * @attention Not implemented in lalsimulation. */
   SpinTaylorT1,	/**< Spinning case T1 models.
                         * @remarks Implemented in lalsimulation (time domain). */
   SpinTaylorT2,	/**< Spinning case T2 models.
                         * @remarks Implemented in lalsimulation (time domain). */
   SpinTaylorT3,	/**< Spinning case T3 models
                         * @attention Not implemented in lalsimulation. */
   SpinTaylorT4,	/**< Spinning case T4 models (lalsimulation's equivalent of SpinTaylorFrameless).
                         * @remarks Implemented in lalsimulation (time domain). */
   SpinTaylorT5,       /**< Spinning case T5. Ref. Sec III of P. Ajith, Phys Rev D (2011)
                         * @attention Not implemented in lalsimulation.  */
   SpinTaylorF2,	/**< Spinning case F2 models (single spin only).
                         * @remarks Implemented in lalsimulation (frequency domain). */
   SpinTaylorFrameless,	/**< Spinning case PN models (replace SpinTaylor by removing the coordinate singularity)
                         * @attention Not implemented in lalsimulation. */
   SpinTaylor,		/**< Spinning case PN models (should replace SpinTaylorT3 in the future)
                         * @attention Not implemented in lalsimulation. */
   PhenSpinTaylor,      /**< Inspiral part of the PhenSpinTaylorRD.
                         * @remarks Implemented in lalsimulation (time domain). */
   PhenSpinTaylorRD,	/**< Phenomenological waveforms, interpolating between a T4 spin-inspiral and the ringdown.
                         * @remarks Implemented in lalsimulation (time domain). */
   SpinQuadTaylor,	/**< Spinning case PN models with quadrupole-monopole and self-spin interaction.
                         * @attention Not implemented in lalsimulation. */
   FindChirpSP,		/**< The stationary phase templates implemented by FindChirpSPTemplate in the findchirp package (equivalent to TaylorF2 at twoPN order).
                         * @attention Not implemented in lalsimulation. */
   FindChirpPTF,	/**< UNDOCUMENTED
                         * @attention Not implemented in lalsimulation. */
   GeneratePPN,		/**< The time domain templates generated by LALGeneratePPNInspiral() in the inject package (equivalent to TaylorT3 at twoPN order).
                         * @attention Not implemented in lalsimulation. */
   BCVC,		/**< UNDOCUMENTED
                         * @attention Not implemented in lalsimulation. */
   FrameFile,		/**< The waveform contains arbitrary data read from a frame file.
                         * @attention Not implemented in lalsimulation. */
   AmpCorPPN,		/**< UNDOCUMENTED
                         * @attention Not implemented in lalsimulation. */
   NumRel,		/**< UNDOCUMENTED
                         * @attention Not implemented in lalsimulation. */
   NumRelNinja2,	/**< The waveform contains REAL8 data generated by lalapps_fr_ninja from a file in the format described in arXiv:0709.0093v3
                         * @attention Not implemented in lalsimulation. */
   Eccentricity,	/**< UNDOCUMENTED
                         * @attention Not implemented in lalsimulation. */
   EOBNR,		/**< UNDOCUMENTED
                         * @attention Not implemented in lalsimulation. */
   EOBNRv2,		/**< UNDOCUMENTED
                         * @remarks Implemented in lalsimulation (time domain). */
   EOBNRv2HM,		/**< UNDOCUMENTED
                         * @remarks Implemented in lalsimulation (time domain). */
   EOBNRv2_ROM,       /**< Frequency domain reduced order model of model EOBNRv2HM, no spin neither higher modes.
                         * @attention Not implemented in lalsimulation. */
   EOBNRv2HM_ROM,       /**< Frequency domain reduced order model of model EOBNRv2HM, no spin but with higher modes.
                         * @attention Not implemented in lalsimulation. */
   SEOBNRv1,		/**< Spin-aligned EOBNR model
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv2,		/**< Spin-aligned EOBNR model v2
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv2_opt,	/**< Optimized Spin-aligned EOBNR model v2
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv3,		/**< Spin precessing EOBNR model v3
                         * @todo Fix implementation in lalsimulation (time domain). */
   SEOBNRv1_ROM_EffectiveSpin, /**< Single-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv1 See [Purrer:2014fza]
                                * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv1_ROM_DoubleSpin, /**< Double-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv1 See [Purrer:2014fza]
                             * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv2_ROM_EffectiveSpin, /**< Single-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv2
                                * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv2_ROM_DoubleSpin, /**< Double-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv2
                             * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv2_ROM_DoubleSpin_HI, /**< High resolution low-mass double-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv2
                                * @remarks Implemented in lalsimulation (frequency domain). */
   HGimri,		/**< Time domain inspiral-merger-ringdown waveform for quasi-circular intermediate mass-ratio inspirals [Huerta & Gair arXiv:1009.1985]
                         * @remarks Implemented in lalsimulation (time domain). */
   IMRPhenomA,		/**< Time domain (non-spinning) inspiral-merger-ringdown waveforms generated from the inverse FFT of IMRPhenomFA.
                         * @remarks Implemented in lalsimulation (time domain and frequency domain). */
   IMRPhenomB,		/**< Time domain (non-precessing spins) inspiral-merger-ringdown waveforms generated from the inverse FFT of IMRPhenomFB.
                         * @remarks Implemented in lalsimulation (time domain and frequency domain). */
   IMRPhenomFA,		/**< Frequency domain (non-spinning) inspiral-merger-ringdown templates of Ajith et al [Ajith_2007kx] with phenomenological coefficients defined in the Table I of [Ajith_2007xh]
                         * @attention Not implemented in lalsimulation.*/
   IMRPhenomFB,		/**< Frequency domain (non-precessing spins) inspiral-merger-ringdown templates of Ajith et al [Ajith_2009bn]
                         * @attention Not implemented in lalsimulation. */
   IMRPhenomC,		/**< Frequency domain (non-precessing spins) inspiral-merger-ringdown templates of Santamaria et al [Santamaria:2010yb] with phenomenological coefficients defined in the Table II of [Santamaria:2010yb].
                         * @remarks Implemented in lalsimulation (time domain and frequency domain). */
   IMRPhenomD,		/**< Frequency domain (non-precessing spins) inspiral-merger-ringdown templates of Husa et al, arXiv:1508.07250 and Khan et al, arXiv:1508.07253 with phenomenological coefficients defined in the Table ...
                         * @remarks Implemented in lalsimulation (frequency domain). */
   IMRPhenomP,		/**< Frequency domain (generic spins) inspiral-merger-ringdown templates of Hannam et al., arXiv:1308.3271 [gr-qc]. Based on IMRPhenomC.
                         * @remarks Implemented in lalsimulation (frequency domain).  */
   IMRPhenomPv2,		/**< Frequency domain (generic spins) inspiral-merger-ringdown templates of Hannam et al., arXiv:1308.3271 [gr-qc]. Based on IMRPhenomD, arXiv:1508.07250 and arXiv:1508.07253.
                         * @remarks Implemented in lalsimulation (frequency domain).  */
   IMRPhenomFC,		/**< Frequency domain (non-precessing spins) inspiral-merger-ringdown templates of Santamaria et al [Santamaria:2010yb] with phenomenological coefficients defined in the Table II of [Santamaria:2010yb]
                         * @attention Not implemented in lalsimulation.*/
   TaylorEt,		/**< UNDOCUMENTED
                         * @remarks Implemented in lalsimulation (time domain). */
   TaylorT4,		/**< UNDOCUMENTED
                         * @remarks Implemented in lalsimulation (time domain). */
   EccentricTD,		/**< Time domain Taylor T4 approximant including orbital eccentricity effects
                         * @remarks Implemented in lalsimulation (time domain). */
   TaylorN,		/**< UNDOCUMENTED
                         * @attention Not implemented in lalsimulation. */
   SpinTaylorT4Fourier, /**< Frequency domain (generic spins) inspiral only waveforms based on TaylorT4, arXiv: 1408.5158
                         * @remarks Implemented in lalsimulation (frequency domain). */
   SpinTaylorT2Fourier, /**< Frequency domain (generic spins) inspiral only waveforms based on TaylorT2, arXiv: 1408.5158
                         * @remarks Implemented in lalsimulation (frequency domain). */
   SpinDominatedWf,     /**< Time domain, inspiral only, 1 spin, precessing waveform, Tapai et al, arXiv: 1209.1722
                         * @remarks Implemented in lalsimulation (time domain). */
   NR_hdf5,              /**< Time domain, NR waveform from HDF file. From INSERT LINKS HERE */
   NumApproximants	/**< Number of elements in enum, useful for checking bounds */
 } Approximant;

/** Enum of various frequency functions */
typedef enum {
    fSchwarzISCO, /**< Schwarzschild ISCO */
    fIMRPhenomAFinal, /**< Final frequency of IMRPhenomA */
    fIMRPhenomBFinal, /**< Final of IMRPhenomB */
    fIMRPhenomCFinal, /**< Final of IMRPhenomC */
    fIMRPhenomDPeak, /**< Frequency of the peak amplitude in IMRPhenomD */
    fEOBNRv2RD, /**< Ringdown frequency of EOBNRv2 */
    fEOBNRv2HMRD, /**< Ringdown frequency of highest harmonic in EOBNRv2HM */
    fSEOBNRv1Peak, /**< Frequency of the peak amplitude in SEOBNRv1 */
    fSEOBNRv1RD, /**< Dominant ringdown frequency in SEOBNRv1 */
    fSEOBNRv2Peak, /**< Frequency of the peak amplitude in SEOBNRv2 */
    fSEOBNRv2RD, /**< Dominant ringdown frequency in SEOBNRv2 */
    NumFreqFunctions /**< Number of elements in the enum */
 } FrequencyFunction;

/** Enum of possible values to use for post-Newtonian order. */
typedef enum {
  LAL_PNORDER_NEWTONIAN,	/**< Newtonain (leading) order */
  LAL_PNORDER_HALF,		/**< 0.5PN <==> O(v) */
  LAL_PNORDER_ONE,		/**< 1PN <==> O(v^2) */
  LAL_PNORDER_ONE_POINT_FIVE,	/**< 1.5PN <==> O(v^3) */
  LAL_PNORDER_TWO,		/**< 2PN <==> O(v^4) */
  LAL_PNORDER_TWO_POINT_FIVE,	/**< 2.5PN <==> O(v^5) */
  LAL_PNORDER_THREE,		/**< 3PN <==> O(v^6) */
  LAL_PNORDER_THREE_POINT_FIVE,	/**< 3.5PN <==> O(v^7)  */
  LAL_PNORDER_PSEUDO_FOUR,	/**< pseudo-4PN tuning coefficients included, true 4PN terms currently unknown */
  LAL_PNORDER_NUM_ORDER		/**< Number of elements in enum, useful for checking bounds */
 } LALPNOrder;

/** Enumeration to specify the tapering method to apply to the waveform */
typedef enum
{
  LAL_SIM_INSPIRAL_TAPER_NONE,		/**< No tapering */
  LAL_SIM_INSPIRAL_TAPER_START,		/**< Taper the start of the waveform */
  LAL_SIM_INSPIRAL_TAPER_END,		/**< Taper the end of the waveform */
  LAL_SIM_INSPIRAL_TAPER_STARTEND,	/**< Taper the start and the end of the waveform */
  LAL_SIM_INSPIRAL_TAPER_NUM_OPTS	/**< Number of elements in enum, useful for checking bounds */
}  LALSimInspiralApplyTaper;

/** Enumeration to specify time or frequency domain */
typedef enum {
  LAL_SIM_DOMAIN_TIME,
  LAL_SIM_DOMAIN_FREQUENCY
 } LALSimulationDomain;

typedef enum {
   LAL_SIM_INSPIRAL_SPINLESS, /** These approximants cannot include spin terms */
   LAL_SIM_INSPIRAL_SINGLESPIN, /** These approximants support a signle spin (by default that is the object 1)*/
   LAL_SIM_INSPIRAL_ALIGNEDSPIN, /** These approximants can include spins aligned with L_N */
   LAL_SIM_INSPIRAL_PRECESSINGSPIN, /** These approximant support fully precessing spins */
   LAL_SIM_INSPIRAL_NUMSPINSUPPORT	/**< Number of elements in enum, useful for checking bounds */
 } SpinSupport;

typedef enum {
  LAL_SIM_INSPIRAL_NO_TESTGR_PARAMS,   /** These approximants cannot accept testGR params as input params */
  LAL_SIM_INSPIRAL_TESTGR_PARAMS,      /** These approximants accept testGR params as input params */
  LAL_SIM_INSPIRAL_NUM_TESTGR_ACCEPT  /**< Number of elements in enum, useful for checking bounds */
 } TestGRaccept;


/**
 * Structure for passing around PN phasing coefficients.
 * For use with the TaylorF2 waveform.
 */
#define PN_PHASING_SERIES_MAX_ORDER 12
typedef struct tagPNPhasingSeries
{
    REAL8 v[PN_PHASING_SERIES_MAX_ORDER+1];
    REAL8 vlogv[PN_PHASING_SERIES_MAX_ORDER+1];
    REAL8 vlogvsq[PN_PHASING_SERIES_MAX_ORDER+1];
}
PNPhasingSeries;

/** @} */

/* general waveform switching generation routines  */
int XLALSimInspiralChooseTDWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant);
int XLALSimInspiralChooseFDWaveform(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_max, REAL8 f_ref, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant);
int XLALSimInspiralTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 z, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant);
SphHarmTimeSeries* XLALSimInspiralTDModesFromPolarizations(REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 z, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant);
int XLALSimInspiralFD(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hcross, REAL8 phiRef, REAL8 deltaF, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_max, REAL8 f_ref, REAL8 r, REAL8 z, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant);
int XLALSimInspiralChooseWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, Approximant approximant); /* DEPRECATED */

/* general waveform switching mode generation routines */
SphHarmTimeSeries *XLALSimInspiralChooseTDModes(REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, int lmax, Approximant approximant);
SphHarmTimeSeries *XLALSimInspiralModesTD(REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, int lmax, Approximant approximant);
COMPLEX16TimeSeries *XLALSimInspiralChooseTDMode(REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, int l, int m, Approximant approximant);

/* routines for generating inspiral waveforms from orbital data */
int XLALSimInspiralPNPolarizationWaveforms(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, REAL8 i, int ampO);
int XLALSimInspiralPNPolarizationWaveformsFromModes(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *v, REAL8TimeSeries *phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, REAL8 i, int O);
int XLALSimInspiralPolarizationsFromSphHarmTimeSeries(REAL8TimeSeries **hp, REAL8TimeSeries **hc, SphHarmTimeSeries *hlms, REAL8 iota, REAL8 psi);
int XLALSimInspiralPNPolarizationWaveformsEccentric(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *V, REAL8TimeSeries *Ecc, REAL8TimeSeries *U, REAL8TimeSeries *Phi, REAL8 m1, REAL8 m2, REAL8 r, REAL8 i, int ampO, int ph_O);
int XLALSimInspiralPrecessingPolarizationWaveforms(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8TimeSeries *LNhatx, REAL8TimeSeries *LNhaty, REAL8TimeSeries *LNhatz, REAL8TimeSeries *E1x, REAL8TimeSeries *E1y, REAL8TimeSeries *E1z, REAL8 m1, REAL8 m2, REAL8 r, REAL8 v0, INT4 ampO);
int XLALSimInspiralPrecessingPolarizationWaveformHarmonic(COMPLEX16 *hplus, COMPLEX16 *hcross, REAL8 v, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhx, REAL8 lnhy, REAL8 lnhz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 dm, REAL8 eta, REAL8 v0, INT4 n, INT4 ampO);

/* approximant, order, axis, and modes handling routines */
int XLALSimInspiralImplementedTDApproximants(Approximant approximant);
int XLALSimInspiralImplementedFDApproximants(Approximant approximant);
int XLALSimInspiralDecomposeWaveformString(int *approximant, int *order, int *axis, const char *waveform);
int XLALSimInspiralGetApproximantFromString(const char *waveform);
int XLALSimInspiralGetPNOrderFromString(const char *waveform);
int XLALSimInspiralGetFrameAxisFromString(const char *waveform);
int XLALSimInspiralGetTaperFromString(const char *string);
int XLALSimInspiralGetHigherModesFromString(const char *string);
int XLALSimInspiralGetSpinSupportFromApproximant(Approximant approx);
int XLALSimInspiralApproximantAcceptTestGRParams(Approximant approx);
const char * XLALSimInspiralGetStringFromApproximant(Approximant approximant);
const char * XLALSimInspiralGetStringFromPNOrder(LALPNOrder order);
const char * XLALSimInspiralGetStringFromTaper(LALSimInspiralApplyTaper taper);
const char * XLALSimInspiralGetStringFromFrameAxis(LALSimInspiralFrameAxis axis);
const char * XLALSimInspiralGetStringFromModesChoice(LALSimInspiralModesChoice modes);

int XLALGetApproximantFromString(const char *waveform); /* DEPRECATED */
int XLALGetOrderFromString(const char *waveform); /* DEPRECATED */
int XLALGetFrameAxisFromString(const char *waveform); /* DEPRECATED */
int XLALGetTaperFromString(const char *string); /* DEPRECATED */
int XLALGetHigherModesFromString(const char *string); /* DEPRECATED */
const char * XLALGetStringFromApproximant(Approximant approximant); /* DEPRECATED */

/* routines for finding information about waveform durations or frequencies */
REAL8 XLALSimInspiralChirpTimeBound(REAL8 fstart, REAL8 m1, REAL8 m2, REAL8 s1, REAL8 s2);
REAL8 XLALSimInspiralMergeTimeBound(REAL8 m1, REAL8 m2);
REAL8 XLALSimInspiralRingdownTimeBound(REAL8 M, REAL8 s);
REAL8 XLALSimInspiralFinalBlackHoleSpinBound(REAL8 S1z, REAL8 S2z);
REAL8 XLALSimInspiralChirpStartFrequencyBound(REAL8 tchirp, REAL8 m1, REAL8 m2);
double XLALSimInspiralGetFrequency(REAL8 m1, REAL8 m2, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, FrequencyFunction freqFunc);
double XLALSimInspiralGetFinalFreq(REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, Approximant approximant);
REAL8 XLALSimInspiralTaylorLength(REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, int O);

/* routines for conditioning waveforms */
int XLALSimInspiralTDConditionStage1(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, REAL8 textra, REAL8 f_min);
int XLALSimInspiralTDConditionStage2(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, REAL8 f_min, REAL8 f_max);

/* routines for transforming initial conditions of precessing waveforms */
int XLALSimInspiralTransformPrecessingNewInitialConditions(REAL8 *incl, REAL8 *S1x, REAL8 *S1y, REAL8 *S1z, REAL8 *S2x, REAL8 *S2y, REAL8 *S2z, REAL8 thetaJN, REAL8 phiJL, REAL8 theta1, REAL8 theta2, REAL8 phi12, REAL8 chi1, REAL8 chi2, REAL8 m1, REAL8 m2, REAL8 fRef);
int XLALSimInspiralTransformPrecessingObsoleteInitialConditions(REAL8 *incl, REAL8 *S1x, REAL8 *S1y, REAL8 *S1z, REAL8 *S2x, REAL8 *S2y, REAL8 *S2z, REAL8 thetaJN, REAL8 phiJL, REAL8 theta1, REAL8 theta2, REAL8 phi12, REAL8 chi1, REAL8 chi2, REAL8 m1, REAL8 m2, REAL8 fRef);

/* routines for generating PN modes based on orbital data */
/* in module LALSimInspiralPNMode.c */

COMPLEX16TimeSeries *XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(REAL8TimeSeries *v, REAL8TimeSeries *phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O, int l, int m);

COMPLEX16TimeSeries *XLALSimInspiralPNMode22(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode21(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode20(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);

COMPLEX16TimeSeries *XLALSimInspiralPNMode33(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode32(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode31(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode30(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);

COMPLEX16TimeSeries *XLALSimInspiralPNMode44(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode43(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode42(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode41(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode40(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);

COMPLEX16TimeSeries *XLALSimInspiralPNMode55(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode54(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode53(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode52(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode51(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode50(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);

COMPLEX16TimeSeries *XLALSimInspiralPNMode66(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode65(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode64(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode63(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode62(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode61(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);
COMPLEX16TimeSeries *XLALSimInspiralPNMode60(REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O);


/* TaylorT1 functions */
/* in module LALSimInspiralTaylorT1.c */

int XLALSimInspiralTaylorT1PNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT1PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT1PNModes(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT1PNMode(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);
int XLALSimInspiralTaylorT1PN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT1PNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);


/* TaylorT2 functions */
/* in module LALSimInspiralTaylorT2.c */

int XLALSimInspiralTaylorT2PNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT2PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT2PNModes(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT2PNMode(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);
int XLALSimInspiralTaylorT2PN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT2PNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);


/* TaylorT3 functions */
/* in module LALSimInspiralTaylorT3.c */

int XLALSimInspiralTaylorT3PNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT3PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT3PNModes(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT3PNMode(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);

int XLALSimInspiralTaylorT3PN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT3PNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);


/* TaylorT4 functions */
/* in module LALSimInspiralTaylorT4.c */

int XLALSimInspiralTaylorT4PNEvolveOrbit(REAL8TimeSeries **v, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT4PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT4PNModes(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT4PNMode(REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);
int XLALSimInspiralTaylorT4PN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT4PNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);


/* TaylorEt functions */
/* in module LALSimInspiralTaylorEt.c */

int XLALSimInspiralTaylorEtPNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **phi, REAL8 phic, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, int O);
int XLALSimInspiralTaylorEtPNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phic, REAL8 x0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 r, REAL8 i, int amplitudeO, int phaseO);
int XLALSimInspiralTaylorEtPN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phic, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 r, REAL8 i, int O);
int XLALSimInspiralTaylorEtPNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phic, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 r, REAL8 i, int O);


/* HGimri functions */
/* in module LALSimInspiralHGimri.c */

int XLALHGimriGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 r, REAL8 i, REAL8 S1z);


/* TaylorF2 functions */
/* in module LALSimInspiralTaylorF2.c */

int XLALSimInspiralTaylorF2AlignedPhasing(PNPhasingSeries **pfa, const REAL8 m1, const REAL8 m2, const REAL8 chi1, const REAL8 chi2, const REAL8 qm_def1, const REAL8 qm_def2, const LALSimInspiralSpinOrder spinO, const LALSimInspiralTestGRParam *nonGRparams);
int XLALSimInspiralTaylorF2Core(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, const REAL8 phi_ref, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 S1z, const REAL8 S2z, const REAL8 f_ref, const REAL8 shft, const REAL8 r, const REAL8 quadparam1, const REAL8 quadparam2, const REAL8 lambda1, const REAL8 lambda2, const LALSimInspiralSpinOrder spinO, const LALSimInspiralTidalOrder tideO, const INT4 phaseO, const INT4 amplitudeO, const LALSimInspiralTestGRParam *nonGRparams);
int XLALSimInspiralTaylorF2(COMPLEX16FrequencySeries **htilde, const REAL8 phi_ref, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 S1z, const REAL8 S2z, const REAL8 fStart, const REAL8 fEnd, const REAL8 f_ref, const REAL8 r, const REAL8 quadparam1, const REAL8 quadparam2, const REAL8 lambda1, const REAL8 lambda2, const LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, const INT4 phaseO, const INT4 amplitudeO, const LALSimInspiralTestGRParam *nonGRparams);


/* SpinTaylor precessing waveform functions */
/* in module LALSimInspiralSpinTaylor.c */

/* Struct containing all of the non-dynamical coefficients needed
 * to evolve a TaylorTx spinning, precessing binary and produce a waveform.
 * This struct is passed to the static Derivatives and StoppingTest functions.*/
typedef struct tagXLALSimInspiralSpinTaylorTxCoeffs
{
  REAL8 M; ///< total mass in solar mass units
  REAL8 Mchirp; ///< chirp mass in solar mass units
  REAL8 eta; ///< symmetric mass ratio
  REAL8 m1M; ///< m1 / M
  REAL8 m2M; ///< m2 / M
  REAL8 wdotnewt; ///< leading order coefficient of wdot = \f$\dot{\omega}\f$
  REAL8 wdotcoeff[LAL_MAX_PN_ORDER]; ///< coeffs. of PN corrections to wdot
  REAL8 wdotlogcoeff; ///< coefficient of log term in wdot
  REAL8 wdot3S1O, wdot3S2O; ///< non-dynamical 1.5PN SO corrections
  REAL8 wdot4S1S2, wdot4S1OS2O; ///< non-dynamical 2PN SS corrections
  REAL8 wdot4S1S1,wdot4S2S2; ///< non-dynamical self S^2 2PN correction
  REAL8 wdot4S1OS1O, wdot4S2OS2O; ///< non-dynamical self SO^2 2PN correction
  REAL8 wdot4QMS1S1; ///< non-dynamical S1^2 2PN quadrupole-monopole correct
  REAL8 wdot4QMS1OS1O; ///< non-dynamical (S1.L)^2 2PN quadrupole-monopole co
  REAL8 wdot4QMS2S2; ///< non-dynamical S2^2 2PN quadrupole-monopole correct
  REAL8 wdot4QMS2OS2O; ///< non-dynamical (S2.L)^2 2PN quadrupole-monopole c
  REAL8 wdot5S1O, wdot5S2O; ///< non-dynamical 2.5PN SO corrections
  REAL8 wdot6S1O, wdot6S2O; ///< non-dynamical 3PN SO corrections
  REAL8 wdot6S1S2, wdot6S1OS2O; ///< non-dynamical 3PN S1-S2 corrections
  REAL8 wdot6S1S1, wdot6S1OS1O; ///< non-dynamical 3PN Spin^2 corrections
  REAL8 wdot6S2S2, wdot6S2OS2O; ///< non-dynamical 3PN Spin^2 corrections
  REAL8 wdot6QMS1S1, wdot6QMS1OS1O; ///< non-dynamical 3PN quadrupole-monopole S1^2 corrections
  REAL8 wdot6QMS2S2, wdot6QMS2OS2O; ///< non-dynamical 3PN quadrupole-monopole S2^2 corrections
  REAL8 wdot7S1O, wdot7S2O; ///< non-dynamical 3.5PN SO corrections
  REAL8 wdottidal10;	///< leading order tidal correction
  REAL8 wdottidal12;	///< next to leading order tidal correction
  REAL8 Ecoeff[LAL_MAX_PN_ORDER]; ///< coeffs. of PN corrections to energy
  REAL8 E3S1O, E3S2O; ///< non-dynamical 1.5PN SO corrections
  REAL8 E4S1S2,E4S1OS2O; ///< non-dynamical 2PN SS correction
  REAL8 E4QMS1S1; ///< non-dynamical S1^2 2PN quadrupole-monopole correction
  REAL8 E4QMS1OS1O;///< non-dynamical (S1.L)^2 2PN quadrupole-monopole correction
  REAL8 E4QMS2S2; ///< non-dynamical S2^2 2PN quadrupole-monopole correction
  REAL8 E4QMS2OS2O;///< non-dynamical (S2.L)^2 2PN quadrupole-monopole correction
  REAL8 E5S1O, E5S2O; ///< non-dynamical 2.5PN SO corrections
  REAL8 E6S1S2;  ///< non-dynamical 3PN S1-S2 correction
  REAL8 E6S1OS2O; ///< non-dynamical 3PN S1.LN S2.LN correction
  REAL8 E6S1S1, E6S1OS1O; ///< non-dynamical 3PN slef-spin^2 corrections
  REAL8 E6S2S2, E6S2OS2O; ///< non-dynamical 3PN self-spin^2 corrections
  REAL8 E6QMS1S1, E6QMS1OS1O; ///< non-dynamical 3PN quadrupole-monopole spin^2 corrections
  REAL8 E6QMS2S2, E6QMS2OS2O; ///< non-dynamical 3PN quadrupole-monopole spin^2 corrections
  REAL8 E7S1O, E7S2O; ///< non-dynamical 3.5PN SO corrections
  REAL8 Etidal10; ///< leading order 5PN tidal correction to energy
  REAL8 Etidal12; ///< next to leading order 6PN tidal correction to energy
  REAL8 dEdvnewt;
  REAL8 Fcoeff[LAL_MAX_PN_ORDER];///<FluxCoeff
  REAL8 Fnewt; ///<newtonian term in Flux
  REAL8 Flogcoeff; ///<log coeff in flux
  REAL8 F3S1O;  ///< Coefficient of S1.LN term
  REAL8 F3S2O;  ///< Coefficient of S2.LN term
  REAL8 F4S1S2; ///< Coefficient of S1.S2 term
  REAL8 F4S1OS2O;///< Coefficient of S1.LN S2.LN term
  REAL8 F4S1S1;  ///< Coefficient of S1.S1 term
  REAL8 F4S1OS1O;///< Coefficient of (S1.LN)^2 term
  REAL8 F4S2S2;  ///< Coefficient of S1.S2 term
  REAL8 F4S2OS2O; ///< Coefficient of (S2.LN)^2 term
  REAL8 F4QMS1S1; ///< Coefficient of S1.S1 term
  REAL8 F4QMS2S2; ///< Coefficient of S2.S2 term
  REAL8 F4QMS1OS1O; ///< Coefficient of quad-monop. (S1.LN)^2 term
  REAL8 F4QMS2OS2O; ///< Coefficient of quad-monop. (S2.LN)^2 term
  REAL8 F5S1O;  ///< Coefficient of (S1.LN)
  REAL8 F5S2O;  ///< Coefficient of (S1.LN) term
  REAL8 F6S1O, F6S2O; ///< Coefficient of (Si.LN) term
  REAL8 F6S1S2, F6S1OS2O; ///< Coefficients of S1.S2 and S1.LN S2.LN terms
  REAL8 F6S1S1, F6S1OS1O; ///< Coefficients of S1.S1 and (S1.LN)^2 terms
  REAL8 F6S2S2, F6S2OS2O; ///< Coefficients of S2.S2 and (S2.LN)^2 terms
  REAL8 F6QMS1S1, F6QMS2S2; ///< Coefficients of quad-monop. S1.S1 and S2.S2 terms
  REAL8 F6QMS1OS1O, F6QMS2OS2O; ///< Coefficients of quad-monop. (S1.LN)^2 and (S2.LN)^2 terms
  REAL8 F7S1O; ///< Coefficients of S1.LN term
  REAL8 F7S2O; ///< Coefficients of S2.LN term
  REAL8 Ftidal10;     ///< leading order 5PN tidal correction
  REAL8 Ftidal12;     ///< next-to-leading order 6PN tidal correction
  REAL8 Ldot3S1O, Ldot3S2O; ///< non-dynamical 1.5PN SO corrections
  REAL8 Ldot4S1S2; ///< non-dynamical 2PN coefficients of S1.LN S2xL and S2.LN S1xL
  REAL8 Ldot4QMS1; ///< non-dynamical quad-monop. 2PN coeff of S1.LN S1xL
  REAL8 Ldot4QMS2; ///< non-dynamical quad-monop. 2PN coeff of S2.LN S2xL
  REAL8 Ldot5S1O, Ldot5S2O; ///< non-dynamical 2.5PN SO corrections
  REAL8 Ldot6S1OS2, Ldot6S2OS1; ///< non-dynamical 3PN S1S2 corrections
  REAL8 Ldot6S1OS1, Ldot6S2OS2; ///< non-dynamical 3PN S^2 corrections
  REAL8 Ldot6QMS1O, Ldot6QMS2O; ///< non-dynamical 3PN quadrupole-monopole S^2 corrections
  REAL8 Ldot7S1, Ldot7S2; ///< non-dynamical 3.5PN SxL terms in Ldot
  REAL8 S1dot3; ///< coeff of LNxS1 term in S1dot
  REAL8 S2dot3; ///< coeff of LNxS2 term in S2dot
  REAL8 Sdot4S2;  ///< coeff of S2xS1 term in S1dot and of S1xS2 in S2dot
  REAL8 Sdot4S2O; ///< coeff of LN.S2 LNxS1 term in S1dot and of LN.S1 LNxS2 in S2dot
  REAL8 S1dot4QMS1O; ///< coeff of quad-monop. LN.S1 LNxS1 term in S1dot
  REAL8 S2dot4QMS2O; ///< coeff of quad-monop. LN.S2 LNxS2 term in S1dot
  REAL8 S1dot5S2; ///< coeff of LNxS1 term in S1dot
  REAL8 S2dot5S1; ///< coeff of LNxS2 term in S2dot
  REAL8 S1dot6S1O, S1dot6S2O;  ///< coeff of LN.Si LNxS1 term in S1dot
  REAL8 S1dot6S2; ///< coeff of S2xS1 term in S1dot
  REAL8 S1dot6QMS1O; ///< coeff of quad-monop. S1.LN LNxS1 term in S1dot
  REAL8 S2dot6S1O, S2dot6S2O; ///< coeff of LN.Si LNxS2 term in S2dot
  REAL8 S2dot6S1;    // Coefficient of S1 x S2 in S2dot
  REAL8 S2dot6QMS2O; //Coeff. of quad-monop. S2.LN LN X S2 term in S2dot
  REAL8 S1dot7S2;// Coefficient of S1 x S2 in S1dot
  REAL8 S2dot7S1;// Coefficient of S1 x S2 in S2dot
  REAL8 fStart; ///< starting GW frequency of integration
  REAL8 fEnd; ///< ending GW frequency of integration
  INT4 phaseO; ///< Twice PN order of GW-phase
  LALSimInspiralSpinOrder spinO; ///< Twice PN order of included spin effects
  LALSimInspiralTidalOrder tideO;///< Twice PN order of included tidal effects
  REAL8 prev_domega; ///< Previous value of domega/dt used in stopping test
} XLALSimInspiralSpinTaylorTxCoeffs;
int XLALSimInspiralSpinTaylorPNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **Phi, REAL8TimeSeries **S1x, REAL8TimeSeries **S1y, REAL8TimeSeries **S1z, REAL8TimeSeries **S2x, REAL8TimeSeries **S2y, REAL8TimeSeries **S2z, REAL8TimeSeries **LNhatx, REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz, REAL8TimeSeries **E1x, REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, Approximant approx);
int XLALSimInspiralSpinTaylorT1(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, int phaseO, int amplitudeO);
int XLALSimInspiralSpinTaylorT2(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, int phaseO, int amplitudeO);
int XLALSimInspiralSpinTaylorT4(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, int phaseO, int amplitudeO);
int XLALSimInspiralSpinTaylorT5(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 incAngle, int phaseO, int amplitudeO);
int XLALSimInspiralSpinTaylorT4PTFQVecs(REAL8TimeSeries **Q1, REAL8TimeSeries **Q2, REAL8TimeSeries **Q3, REAL8TimeSeries **Q4, REAL8TimeSeries **Q5, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 chi1, REAL8 kappa, REAL8 fStart, REAL8 lambda1, REAL8 lambda2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, int phaseO);
int XLALSimInspiralSpinTaylorT2Fourier(COMPLEX16FrequencySeries **hplus, COMPLEX16FrequencySeries **hcross, REAL8 fMin, REAL8 fMax, REAL8 deltaF, INT4 kMax, REAL8 phiRef, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, INT4 amplitudeO, INT4 phiRefAtEnd);
int XLALSimInspiralSpinTaylorT4Fourier(COMPLEX16FrequencySeries **hplus, COMPLEX16FrequencySeries **hcross, REAL8 fMin, REAL8 fMax, REAL8 deltaF, INT4 kMax, REAL8 phiRef, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, INT4 amplitudeO, INT4 phiRefAtEnd);
int XLALSimInspiralSpinTaylorF2(COMPLEX16FrequencySeries **hplus_out, COMPLEX16FrequencySeries **hcross_out, REAL8 phi_ref, REAL8 deltaF, REAL8 m1_SI, REAL8 m2_SI, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, const REAL8 fStart, const REAL8 fEnd, const REAL8 f_ref, const REAL8 r, LALSimInspiralTestGRParam *moreParams, const LALSimInspiralSpinOrder spinO, const INT4 phaseO, const INT4 amplitudeO);
int XLALSimInspiralPrecessingPTFQWaveforms(REAL8TimeSeries **Q1, REAL8TimeSeries **Q2, REAL8TimeSeries **Q3, REAL8TimeSeries **Q4, REAL8TimeSeries **Q5, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8TimeSeries *LNhatx, REAL8TimeSeries *LNhaty, REAL8TimeSeries *LNhatz, REAL8TimeSeries *E1x, REAL8TimeSeries *E1y, REAL8TimeSeries *E1z, REAL8 m1, REAL8 m2, REAL8 r);
int XLALSimInspiralInitialConditionsPrecessingApproxs(REAL8 *inc, REAL8 *S1x, REAL8 *S1y, REAL8 *S1z, REAL8 *S2x, REAL8 *S2y, REAL8 *S2z, const REAL8 inclIn, const REAL8 S1xIn, const REAL8 S1yIn, const REAL8 S1zIn, const REAL8 S2xIn, const REAL8 S2yIn, const REAL8 S2zIn, const REAL8 m1, const REAL8 m2, const REAL8 fRef, LALSimInspiralFrameAxis axisChoice);
INT4 XLALSimInspiralSpinDerivatives(REAL8 *dLNhx, REAL8 *dLNhy, REAL8 *dLNhz, REAL8 *dE1x, REAL8 *dE1y, REAL8 *dE1z, REAL8 *dS1x, REAL8 *dS1y, REAL8 *dS1z, REAL8 *dS2x, REAL8 *dS2y, REAL8 *dS2z, REAL8 *dphiExtra, const REAL8 v, const REAL8 LNhx, const REAL8 LNhy, const REAL8 LNhz, const REAL8 E1x, const REAL8 E1y, const REAL8 E1z, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 LNhdotS1, const REAL8 LNhdotS2, XLALSimInspiralSpinTaylorTxCoeffs *params);
INT4 XLALSimInspiralSpinTaylorT4Derivatives(REAL8 t, const REAL8 values[], REAL8 dvalues[], void *mparams);
INT4 XLALSimInspiralSpinTaylorT4Setup(XLALSimInspiralSpinTaylorTxCoeffs *params, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO);
INT4 XLALSimSpinTaylorEnergySpinDerivativeSetup(XLALSimInspiralSpinTaylorTxCoeffs *params, const REAL8 lambda1, const REAL8 lambda2, const REAL8 quadparam1, const REAL8 quadparam2);
INT4 XLALSimInspiralSetEnergyPNTerms(REAL8 *Espin3, REAL8 *Espin4, REAL8 *Espin5, REAL8 *Espin6, REAL8 *Espin7, XLALSimInspiralSpinTaylorTxCoeffs *params, const REAL8 LNhdotS1, const REAL8 LNhdotS2, const REAL8 S1sq, const REAL8 S2sq, const REAL8 S1dotS2);

/* time domain eccentric functions */
/* in module LALSimInspiralEccentricTD.c */

int XLALSimInspiralEccentricTDPNEvolveOrbit(REAL8TimeSeries **v, REAL8TimeSeries **et, REAL8TimeSeries **l, REAL8TimeSeries **lambda, REAL8TimeSeries **u, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 e_min, int O);
int XLALSimInspiralEccentricTDPNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 e_min, int amplitudeO, int phaseO);
int XLALSimInspiralEccentricTDPN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 e_min, int O);
int XLALSimInspiralEccentricTDPNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 e_min, int O);


/* frequency domain eccentric functions */
/* in module LALSimInspiralEccentricityFD.c */

int XLALSimInspiralEFD(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8 phiRef, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 fStart, const REAL8 fEnd, const REAL8 i, const REAL8 r, const REAL8 inclination_azimuth, const REAL8 e_min, int phaseO);


/* spin-dominated waveform functions */
/* in module LALSimInspiralSpinDominatedWaveform.c */

int XLALSimInspiralSpinDominatedWaveformInterfaceTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 D, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, int phaseO, int amplitudeO, REAL8 phiRef);
int XLALSimInspiralSpinDominatedWaveformDriver(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 totalmass, REAL8 nu, REAL8 chi1, REAL8 D, REAL8 kappa1, REAL8 beta1, REAL8 theta, REAL8 fStart, REAL8 fRef, int phaseO, int amplitudeO, REAL8 deltaT, REAL8 phiRef, REAL8 phin0);


/* TaylorF2 Reduced Spin routines */
/* in module LALSimInspiralTaylorF2ReducedSpin.c */

int XLALSimInspiralTaylorF2ReducedSpin(COMPLEX16FrequencySeries **htilde, const REAL8 phic, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 fStart, const REAL8 fEnd, const REAL8 r, const INT4 phaseO, const INT4 ampO);
int XLALSimInspiralTaylorF2ReducedSpinTidal(COMPLEX16FrequencySeries **htilde, const REAL8 phic, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const REAL8 lam1, const REAL8 lam2, const REAL8 fStart, const REAL8 fEnd, const REAL8 r, const INT4 phaseO, const INT4 ampO);
REAL8 XLALSimInspiralTaylorF2ReducedSpinChirpTime(const REAL8 fStart, const
REAL8 m1_SI, const REAL8 m2_SI, const REAL8 chi, const INT4 O);
REAL8 XLALSimInspiralTaylorF2ReducedSpinComputeChi(const REAL8 m1, const REAL8 m2, const REAL8 s1z, const REAL8 s2z);
int XLALSimInspiralTaylorF2RedSpinMetricMChirpEtaChi(REAL8 *gamma00, REAL8 *gamma01, REAL8 *gamma02, REAL8 *gamma11, REAL8 *gamma12, REAL8 *gamma22, const REAL8 mc, const REAL8 eta, const REAL8 chi, const REAL8 fLow, const REAL8FrequencySeries *Sh);
gsl_matrix *XLALSimInspiralTaylorF2RedSpinFisherMatrixChirpTimes(const REAL8 theta0, const REAL8 theta3, const REAL8 theta3s, const REAL8 fLow, const REAL8 df, REAL8Vector *momI_0, REAL8Vector *momI_2, REAL8Vector *momI_3, REAL8Vector *momI_4, REAL8Vector *momI_5, REAL8Vector *momI_6, REAL8Vector *momI_7, REAL8Vector *momI_8, REAL8Vector *momI_9, REAL8Vector *momI_10, REAL8Vector *momI_11, REAL8Vector *momI_12, REAL8Vector *momI_13, REAL8Vector *momI_14, REAL8Vector *momI_15, REAL8Vector *momI_16, REAL8Vector *momJ_5, REAL8Vector *momJ_6, REAL8Vector *momJ_7, REAL8Vector *momJ_8, REAL8Vector *momJ_9, REAL8Vector *momJ_10, REAL8Vector *momJ_11, REAL8Vector *momJ_12, REAL8Vector *momJ_13, REAL8Vector *momJ_14, REAL8Vector *momK_10, REAL8Vector *momK_11, REAL8Vector *momK_12);
int XLALSimInspiralTaylorF2RedSpinMetricChirpTimes(REAL8 *gamma00, REAL8 *gamma01, REAL8 *gamma02, REAL8 *gamma11, REAL8 *gamma12, REAL8 *gamma22, const REAL8 theta0, const REAL8 theta3, const REAL8 theta3s, const REAL8 fLow, const REAL8 df, REAL8Vector *momI_0, REAL8Vector *momI_2, REAL8Vector *momI_3, REAL8Vector *momI_4, REAL8Vector *momI_5, REAL8Vector *momI_6, REAL8Vector *momI_7, REAL8Vector *momI_8, REAL8Vector *momI_9, REAL8Vector *momI_10, REAL8Vector *momI_11, REAL8Vector *momI_12, REAL8Vector *momI_13, REAL8Vector *momI_14, REAL8Vector *momI_15, REAL8Vector *momI_16, REAL8Vector *momJ_5, REAL8Vector *momJ_6, REAL8Vector *momJ_7, REAL8Vector *momJ_8, REAL8Vector *momJ_9, REAL8Vector *momJ_10, REAL8Vector *momJ_11, REAL8Vector *momJ_12, REAL8Vector *momJ_13, REAL8Vector *momJ_14, REAL8Vector *momK_10, REAL8Vector *momK_11, REAL8Vector *momK_12);
int XLALSimInspiralTaylorF2RedSpinComputeNoiseMoments(REAL8Vector *momI_0, REAL8Vector *momI_2, REAL8Vector *momI_3, REAL8Vector *momI_4, REAL8Vector *momI_5, REAL8Vector *momI_6, REAL8Vector *momI_7, REAL8Vector *momI_8, REAL8Vector *momI_9, REAL8Vector *momI_10, REAL8Vector *momI_11, REAL8Vector *momI_12, REAL8Vector *momI_13, REAL8Vector *momI_14, REAL8Vector *momI_15, REAL8Vector *momI_16, REAL8Vector *momJ_5, REAL8Vector *momJ_6, REAL8Vector *momJ_7, REAL8Vector *momJ_8, REAL8Vector *momJ_9, REAL8Vector *momJ_10, REAL8Vector *momJ_11, REAL8Vector *momJ_12, REAL8Vector *momJ_13, REAL8Vector *momJ_14, REAL8Vector *momK_10, REAL8Vector *momK_11, REAL8Vector *momK_12, REAL8Vector *Sh, REAL8 fLow, REAL8 df);
void XLALSimInspiralTaylorF2RedSpinChirpTimesFromMchirpEtaChi(double *theta0, double *theta3, double *theta3s, double mc, double eta, double chi, double fLow);
void XLALSimInspiralTaylorF2RedSpinMchirpEtaChiFromChirpTimes(double *mc, double *eta, double *chi, double theta0, double theta3, double theta3s, double fLow);

/* waveform tapering routines */
/* in module LALSimInspiralWaveformTaper.c */

int XLALSimInspiralREAL4WaveTaper(REAL4Vector *signalvec, LALSimInspiralApplyTaper bookends);
int XLALSimInspiralREAL8WaveTaper(REAL8Vector *signalvec, LALSimInspiralApplyTaper bookends);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRAL_H */
