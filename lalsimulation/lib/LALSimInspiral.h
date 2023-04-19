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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#ifndef _LALSIMINSPIRAL_H
#define _LALSIMINSPIRAL_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimSphHarmSeries.h>
#include <lal/LALSimInspiralTestGRParams.h>
#include <lal/LALSimInspiralWaveformFlags.h>
#include <lal/LALSimInspiralWaveformParams.h>

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
 * double m1 = 1.4 * LAL_MSUN_SI; // mass of body 1
 * double m2 = 1.4 * LAL_MSUN_SI; // mass of body 2
 * double S1x = 0.0;              // x-component of dimensionless spin of body 1
 * double S1y = 0.0;              // y-component of dimensionless spin of body 1
 * double S1z = 0.0;              // z-component of dimensionless spin of body 1
 * double S2x = 0.0;              // x-component of dimensionless spin of body 2
 * double S2y = 0.0;              // y-component of dimensionless spin of body 2
 * double S2z = 0.0;              // z-component of dimensionless spin of body 2
 * double r = 1e6 * LAL_PC_SI;    // distance
 * double inclination = 0.0;      // angle between L and view direction, \iota in @image
 * double phiRef = 0.0;           // orbital phase at reference, helf of main GW phase at reference
 * double longAscNodes=0.0;       // longitude of ascending nodes, degenerate with the polarization angle, related to Omega in @image by Omega=longAscNodes+pi/2
 * double eccentricity=0.0;       // eccentricity at reference epoch
 * double meanPerAno=0.0;         // mean anomaly at reference epoch, i.e. the ratio of time passed since last periastron passage to the time interval between two periastron passages, times 2pi.  Note: This is not a geometric angle that can be visualized in @image
 * double deltaT = 1.0/16384.0;   // series sampling interval
 * double f_min = 40.0;           // start frequency of inspiral
 * double f_ref = 0.0;            // reference frequency: 0 means waveform end
 * LALDict *LALpars = NULL;       // structure containing variable with default values
 * Approximant approximant = TaylorT2;  // post-Newtonian approximant
 * REAL8TimeSeries *hplus = NULL;  // plus polarization to be returned
 * REAL8TimeSeries *hcross = NULL; // cross polarization to be returned
 * ...
 * XLALSimInspiralChooseTDWaveform( &hplus, &hcross, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, r, i, phiRef, longAscNodes, eccentricity, meanPerAno, deltaT, f_min, f_ref, LALpars, approximant);
 * @endcode
 *
 * Generate a frequency-domain waveform:
 * @code
 * #include <lal/LALSimInspiral.h>
 * #include <lal/LALConstants.h>
 * ...
 * double m1 = 1.4 * LAL_MSUN_SI;     // mass of body 1
 * double m2 = 1.4 * LAL_MSUN_SI;     // mass of body 2
 * double S1x = 0.0;                  // x-component of dimensionless spin of body 1
 * double S1y = 0.0;                  // y-component of dimensionless spin of body 1
 * double S1z = 0.0;                  // z-component of dimensionless spin of body 1
 * double S2x = 0.0;                  // x-component of dimensionless spin of body 2
 * double S2y = 0.0;                  // y-component of dimensionless spin of body 2
 * double S2z = 0.0;                  // z-component of dimensionless spin of body 2
 * double distance = 1e6 * LAL_PC_SI; // distance
 * double inclination = 0.0;          // angle between L and view direction, \iota in @image
 * double phiRef = 0;                 // orbital phase at reference, half of main GW phase at reference
 * double longAscNodes=0.0;           // longitude of ascending nodes, degenerate with the polarization angle, related to Omega in @image by Omega=longAscNodes+pi/2
 * double eccentricity=0.0;           // eccentricity at reference epoch
 * double meanPerAno=0.0;             // mean anomaly at reference epoch, i.e. the ratio of time passed since last periastron passage to the time interval between two periastron passages, times 2pi.  Note: This is not a geometric angle that can be visualized in @image
 * double deltaF = 1.;                // frequency sampling interval
 * double f_min = 40.0;               // start frequency of inspiral
 * double f_max = 0.0;                // end frequency of inspiral: 0 means use default
 * double f_ref = 0.0;                // reference frequency: 0 means waveform end
 * LALDict *LALpars = NULL;       // structure containing variable with default values
 * Approximant approximant = TaylorF2;  // post-Newtonian approximant
 * COMPLEX16FrequencySeries *hptilde = NULL;  // plus polarization to be returned
 * COMPLEX16FrequencySeries *hctilde = NULL; // cross polarization to be returned
 * ...
 * XLALSimInspiralChooseFDWaveform(&hptilde, &hctilde, m1, m2, S1x, S1y, S1z, S2x, S2y, S2z, f_min, f_max, f_ref, r, i, phiRef, longAscNodes, eccentricity, meanPerAno, deltaF, f_min, f_ref, LALpars, approximant);
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
 * binary system, while the unit vector from body 2 to body 1 defines the x-axis of
 * the binary system.  The x-y-plane is therefore the orbital plane, at least
 * at the moment the binary system is at the reference gravitational wave
 * frequency.
 *
 * The spin components for body 1, (@p S1x,@p S1y, @p S1z), and for body 2,
 * (@p S2x,@p S2y, @p S2z), are defined in the source-frame.  Therefore,
 * when the spins are aligned with the orbital angular momentum,
 * @p S1x = @p S1y = @p S2x = @p S2y = 0.
 * @note
 * The spin components transverse to the orbital angular momentum @b L at the
 * reference gravitational wave frequency @p f_ref are given with respect to
 * the triad x-y-z, with x-axis parallel to the vector pointing from
 * body 2 to body 1.
 *
 * The wave frame is defined by the Z-axis, which points toward the Earth,
 * and some reference direction, defining the X-axis.  The X-Y-plane is
 * therefore the plane of the sky.
 *
 * The plus- and cross-polarizations of the gravitational waveform are defined
 * in this wave frame.  Specifically, if \f$ h^{ij} \f$ is computed in the
 * source frame, then
 * \f[ h_+ = \frac12 ( \hat{P}_i \hat{P}_j - \hat{Q}_i \hat{Q}_j ) h^{ij} \f]
 * and
 * \f[ h_\times = \frac12 ( \hat{P}_i \hat{Q}_j + \hat{Q}_i \hat{P}_j ) h^{ij} \f]
 * where \f$ \hat{P}_i \f$ are the components of the unit vector pointing
 * along the X-axis and \f$ \hat{Q}_i \f$ are the components of the unit
 * vector pointing along the Y-axis.
 *
 * The orbital elements are:
 *
 *  * Inclination (&iota;).  The angle between the Z-axis of the wave frame
 *    and the z-axis of the source frame.
 *  * Reference phase (&Phi;) : The angle on the plane of the orbit from the
 *    line of ascending nodes to the position of body 1
 *    (x axis in our convention).
 *    ascending node @htmlonly &#x260A; @endhtmlonly.
 *  * Longitude of ascending node (&Omega;).  The angle on the plane of the
 *    sky from the X-axis of the reference direction in the wave frame to the
 *    ascending node @htmlonly &#x260A; @endhtmlonly.
 *    @note This angle is entirely degenerate with the polarization angle &psi;.
 *    @attention
 *    In the present implementation, the Y-axis in the wave frame is defined to
 *    be the ascending node @htmlonly &#x260A; @endhtmlonly.
 *    Therefore, &Omega;=&pi; /2 by default with the consequences that
 *    the z axis lies in the X-Z plane, with positive projection over X.
 *    Value of &Omega; can be changed providing a non zero
 *    longAscNodes= &Omega; - &pi; /2, this corresponding to a rotation in the
 *    observation direction, i.e. a polarization rotation.
 *    Another consequence is that the Z axis lies in the plane spanned by z
 *    and the axis perpendicular both z and the line of ascending nodes
 *    (i.e. y at &Phi;=0) with positive projection over the latter.
 *  * True anomaly (&delta;). The angle along the orbital plane from the
 *    periapsis to the present position of the orbiting body 1
 *    (it only applies to eccentric orbits).
 *
 * @sa
 * The coordinate systems used here follows those of
 * > L. Blanchet, G. Faye, B. R. Iyer and S. Sinha,
 * > "The Third post-Newtonian gravitational wave polarisations and
 * > associated spherical harmonic modes for inspiralling compact binaries
 * > in quasi-circular orbits"
 * > Class. Quant. Grav. @b 25, (2008) 165003
 * > Erratum: [Class. Quant. Grav. @b 29, (2012) 239501,
 * > arXiv:0802.1249 [gr-qc].
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
 * @defgroup LALSimInspiralTaylorF2Ecc_c           Module LALSimInspiralTaylorF2Ecc.c
 * @defgroup LALSimInspiralSpinTaylor_c            Module LALSimInspiralSpinTaylor.c
 * @defgroup LALSimInspiralEccentricTD_c           Module LALSimInspiralEccentricTD.c
 * @defgroup LALSimInspiralEccentricityFD_c        Module LALSimInspiralEccentricityFD.c
 * @defgroup LALSimInspiralSpinDominatedWaveform_c Module LALSimInspiralSpinDominatedWaveform.c
 * @defgroup LALSimInspiralTaylorF2ReducedSpin_c   Module LALSimInspiralTaylorF2ReducedSpin.c
 * @defgroup LALSimInspiralHGimri_c                Module LALSimInspiralHGimri.c
 * @defgroup LALSimInspiralWaveformFlags_c         Module LALSimInspiralWaveformFlags.c
 * @defgroup LALSimInspiralTestGRParams_c          Module LALSimInspiralTestGRParams.c
 * @defgroup LALSimInspiralWaveformTaper_c         Module LALSimInspiralWaveformTaper.c
 * @defgroup LALSimInspiralNRSur4d2s_c             Module LALSimInspiralNRSur4d2s.c
 * @defgroup LALSimIMRNRHybSur3dq8_c               Module LALSimIMRNRHybSur3dq8.c
 * @}
 *
 * @addtogroup LALSimInspiral_h
 * @{
 */

#define LAL_PN_MODE_L_MAX 3
/* (2x) Highest available PN order - UPDATE IF NEW ORDERS ADDED!!*/
#define LAL_MAX_PN_ORDER 8
#define LAL_MAX_ECC_PN_ORDER 6
#define LAL_DEFAULT_F_ECC -1.0

/**
 * Enum that specifies the PN approximant to be used in computing the waveform.
 * Please add new approximants at the end of the list, so as to prevent the
 * values changing for existing approximants and breaking ABI compatibility.
 */
typedef enum tagApproximant {
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
   TaylorF2Ecc,		/**< The standard stationary phase approximation with eccentricity; Outputs a frequency-domain wave.
                         * @remarks Implemented in lalsimulation (frequency domain). */
   TaylorF2NLTides,     /**< The standard stationary phase approximation including a phenomenological model of nonlinear tidal effects; Outputs a frequency-domain wave.
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
   SpinTaylorT2,	/**< Spinning case T2 models
                         * @attention Not implemented in lalsimulation. */
   SpinTaylorT3,	/**< Spinning case T3 models
                         * @attention Not implemented in lalsimulation. */
   SpinTaylorT4,	/**< Spinning case T4 models (lalsimulation's equivalent of SpinTaylorFrameless).
                         * @remarks Implemented in lalsimulation (time domain). */
   SpinTaylorT5,        /**< Spinning case T5 models, which is a variant of the spinning version of the original TaylorT2 (see \cite Buonanno:2009zt) described in sec. III of \cite Ajith:2011ec. SpinTaylorT2 is NOT implemented in LALSimulation.
		         * @remarks Implemented in lalsimulation (time domain). */
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
   TEOBResum_ROM,         /**< Time domain reduced order model of EOB with tidal effects.
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv1,		/**< Spin-aligned EOBNR model
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv2,		/**< Spin-aligned EOBNR model v2
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv2_opt,	/**< Optimized Spin-aligned EOBNR model v2
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv3,		/**< Spin precessing EOBNR model v3
                         * @todo Fix implementation in lalsimulation (time domain). */
   SEOBNRv3_pert,        /**< Perturbed [m1 -> m1*(1+1e-15)] Spin precessing EOBNR model v3
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv3_opt,        /**< Optimized Spin precessing EOBNR model v3
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv3_opt_rk4,        /**< USE RK4 Optimized Spin precessing EOBNR model v3
                         * @todo Fix implementation in lalsimulation (time domain). */
   SEOBNRv4,		/**< Spin nonprecessing EOBNR model v4
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv4_opt,	/**< Optimized Spin-aligned EOBNR model v4
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv4P,		/**< Spin precessing EOBNR model based on SEOBNRv4
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv4PHM,		/**< Spin precessing EOBNR model based on SEOBNRv4HM
                         * @remarks Implemented in lalsimulation (time domain). */
   SEOBNRv2T,	/**< Tidal EOB model
                     * @remarks Implemented in lalsimulation (time domain). Parameter range: q=[1,3], Sz=[-0.5,0.5], Lambda2=[0,5000]. Initial conditions solver can fail when starting frequency is too low (rate of failure 0.3% at fmin=10Hz for M=3Msol). */
   SEOBNRv4T,	/**< Tidal EOB model
             * @remarks Implemented in lalsimulation (time domain). Parameter range: q=[1,3], Sz=[-0.5,0.5], Lambda2=[0,5000]. Initial conditions solver can fail when starting frequency is too low (rate of failure 0.3% at fmin=10Hz for M=3Msol). */
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
   Lackey_Tidal_2013_SEOBNRv2_ROM, /**< Frequency domain tidal model based on reduced order model of SEOBNRv2
                                * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv4_ROM, /**< Low-mass double-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv4
                                * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv4HM_ROM, /**< Low-mass double-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv4hm
		    * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv4_ROM_NRTidal, /**< Low-mass double-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv4 [Bohe et al, arXiv:1611.03703] with tidal phase corrections [Dietrich et al, arXiv:1706.02969
                                * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv4_ROM_NRTidalv2, /**< based on NRTidalv2; https://arxiv.org/abs/1905.06011.
                             * @remarks Implemented in lalsimulation (time domain and frequency domain). */
   SEOBNRv4_ROM_NRTidalv2_NSBH, /**< NSBH model based on SEOBNRv4_ROM_NRTidalv2
                                * @remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv4T_surrogate, /**< Double-spin frequency domain surrogate model of spin-aligned tidal EOBNR model SEOBNRv4T
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
   IMRPhenomD_NRTidal,   /**< Uses arxiv:1706.02969 to upgrad IMRPhenomD to a tidal approximant
                         * @remarks Implemented in lalsimulation (frequency domain). */
   IMRPhenomD_NRTidalv2, /**< NRTidalv2; https://arxiv.org/abs/1905.06011
                            * @remarks Implemented in lalsimulation (time domain and frequency domain).*/
   IMRPhenomNSBH,   /**< NSBH Tidal model.
                         * @remarks Implemented in lalsimulation (frequency domain). */
   IMRPhenomHM,     /**< Frequency domain with higher modes (non-precessing spins) inspiral-merger-ringdown templates, based on IMRPhenomD.
                   * @remarks Implemented in lalsimulation (frequency domain). Ref London et al, arXiv:1708.00404 */
   IMRPhenomP,		/**< Frequency domain (generic spins) inspiral-merger-ringdown templates of Hannam et al., arXiv:1308.3271 [gr-qc]. Based on IMRPhenomC.
                         * @remarks Implemented in lalsimulation (frequency domain).  */
   IMRPhenomPv2,		/**< Frequency domain (generic spins) inspiral-merger-ringdown templates of Hannam et al., arXiv:1308.3271 [gr-qc]. Based on IMRPhenomD, arXiv:1508.07250 and arXiv:1508.07253.
                         * @remarks Implemented in lalsimulation (frequency domain).  */
   IMRPhenomPv2_NRTidal, /**< Frequency domain tidal version of IMRPhenomPv2, using NRTidal framework from arXiv:1706.02969 */
   IMRPhenomPv2_NRTidalv2, /**< Frequency domain tidal version; based on https://arxiv.org/abs/1905.06011
                            * @remarks Implemented in lalsimulation (time domain and frequency domain).*/
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
   SpinTaylorT5Fourier, /**< Frequency domain (generic spins) inspiral only waveforms based on TaylorT5, \cite Klein:2014bua , (the paper refers to SpinTaylorT2, but it is actually SpinTaylorT5 which is being used.)
                         * @remarks Implemented in lalsimulation (frequency domain). */
   SpinDominatedWf,     /**< Time domain, inspiral only, 1 spin, precessing waveform, Tapai et al, arXiv: 1209.1722
                         * @remarks Implemented in lalsimulation (time domain). */
   NR_hdf5,              /**< Time domain, NR waveform from HDF file. From INSERT LINKS HERE */
   NRSur4d2s,
   NRSur7dq2,           /**< Time domain, fully precessing NR surrogate model with up to ell=4 modes, arxiv: 1705.07089 */
   NRSur7dq4,           /**< q=4 extension of NRSur7dq2, arxiv: 1905.09300 */
   SEOBNRv4HM,	/**< Spin nonprecessing EOBNR model v4 with higher modes, PhysRevD.98.084028 [arXiv:1803.10701]
                     * @remarks Implemented in lalsimulation (time domain). */
   NRHybSur3dq8,        /**< Time domain, aligned-spin, higher modes, hybridized. Paper arxiv:1812.07865 */
   IMRPhenomXAS, 		/**< Frequency domain, non-precessing phenomenological IMR waveform model ([arXiv:2001.11412]). */
   IMRPhenomXHM, 		/**< Frequency domain, non-precessing phenomenological IMR waveform model with subdominant modes ([arXiv:2001.10914 [gr-qc]]) and accelerated evaluation through adapted grids (arXiv:2001.10897 [gr-qc]) */
   IMRPhenomPv3,    /**< Frequency domain (generic spins) inspiral-merger-ringdown templates of Hannam et al., arXiv:1308.3271 [gr-qc]. Based on IMRPhenomD, arXiv:1508.07250 and arXiv:1508.07253. But updated the precession angles to use the ones in arXiv 1703.03967.
                         * @remarks Implemented in lalsimulation (frequency domain).  */
   IMRPhenomPv3HM,    /**< Frequency domain (generic spins) inspiral-merger-ringdown templates of Khan et al. PhysRevD.101.024056. Based on IMRPhenomHM arXiv:1708.00404. And the precession angles of IMRPhenomPv3 1809.10113 and arXiv 1703.03967.
                       * @remarks Implemented in lalsimulation (frequency domain).  */
   IMRPhenomXP, 		/**< Frequency domain, precessing phenomenological IMR waveform model. */
   IMRPhenomXPHM, 	/**< Frequency domain, precessing with subdominant modes phenomenological IMR waveform model. */
   TEOBResumS,          /**< Resummed Spin-aligned Tidal EOB
                         * @remarks Implemented in lalsimulation (time domain). */
   IMRPhenomT,      /** Time domain, non-precessing phenomenological IMR waveform model for the dominant (2,2) and (2,-2) modes ([arXiv: 20XY.ZZZZZ]). */
   IMRPhenomTHM,      /** Time domain, non-precessing phenomenological IMR waveform model with subdominant modes ([arXiv: 20XY.ZZZZZ]). */
   IMRPhenomTP,      /** Time domain, precessing phenomenological IMR waveform model for L=2 sector ([arXiv: 20XY.ZZZZZ]). */
   IMRPhenomTPHM,      /** Time domain, precessing phenomenological IMR waveform model with subdominant modes ([arXiv: 20XY.ZZZZZ]). */
   SEOBNRv5_ROM, /**< Low-mass double-spin frequency domain reduced order model of spin-aligned EOBNR model SEOBNRv5 remarks Implemented in lalsimulation (frequency domain). */
   SEOBNRv4HM_PA,    /** Spin non-precessing EOBNR model v4 with higher modes post-adiabatic dynamics (time domain), PhysRevD.104.124087 [arXiv:2105.06983] */
   pSEOBNRv4HM_PA,    /** Spin non-precessing EOBNR model v4 with higher modes post-adiabatic dynamics (time domain) and TGR ringdown effects, PhysRevD.104.124087 [arXiv:2105.06983] */
   IMRPhenomXAS_NRTidalv2,         /**< Tidal extension of IMRPhenomXAS based on [arXiv:1905.06011]. */
   IMRPhenomXP_NRTidalv2,         /**< Tidal extension of IMRPhenomXP based on [arXiv:1905.06011]. */
   IMRPhenomXO4a,    /**< Frequency domain, precessing with subdominant modes phenomenological IMR waveform model with NR-tuned precession angles. */
   ExternalPython, /** External Python model **/
   NumApproximants,	/**< Number of elements in enum, useful for checking bounds */
 } Approximant;

/** Enum of various frequency functions */
typedef enum tagFrequencyFunction {
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
    fSEOBNRv4Peak, /**< Frequency of the peak amplitude in SEOBNRv4 */
    fSEOBNRv4RD, /**< Dominant ringdown frequency in SEOBNRv4 */
    fTEOBResumSFinal, /**< Dominant ringdown frequency in TEOBResumS */
    fSEOBNRv5Peak, /**< Frequency of the peak amplitude in SEOBNRv5_ROM */
    fSEOBNRv5RD, /**< Dominant ringdown frequency in SEOBNRv5_ROM */
    NumFreqFunctions /**< Number of elements in the enum */
 } FrequencyFunction;

/** Enum of possible values to use for post-Newtonian order. */
typedef enum tagLALPNOrder {
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
typedef enum tagLALSimInspiralApplyTaper
{
  LAL_SIM_INSPIRAL_TAPER_NONE,		/**< No tapering */
  LAL_SIM_INSPIRAL_TAPER_START,		/**< Taper the start of the waveform */
  LAL_SIM_INSPIRAL_TAPER_END,		/**< Taper the end of the waveform */
  LAL_SIM_INSPIRAL_TAPER_STARTEND,	/**< Taper the start and the end of the waveform */
  LAL_SIM_INSPIRAL_TAPER_NUM_OPTS	/**< Number of elements in enum, useful for checking bounds */
}  LALSimInspiralApplyTaper;

/** Enumeration to specify time or frequency domain */
typedef enum tagLALSimulationDomain {
  LAL_SIM_DOMAIN_TIME,
  LAL_SIM_DOMAIN_FREQUENCY
 } LALSimulationDomain;

typedef enum tagSpinSupport {
   LAL_SIM_INSPIRAL_SPINLESS, /** These approximants cannot include spin terms */
   LAL_SIM_INSPIRAL_SINGLESPIN, /** These approximants support a signle spin (by default that is the object 1)*/
   LAL_SIM_INSPIRAL_ALIGNEDSPIN, /** These approximants can include spins aligned with L_N */
   LAL_SIM_INSPIRAL_PRECESSINGSPIN, /** These approximant support fully precessing spins */
   LAL_SIM_INSPIRAL_CASEBYCASE_SPINSUPPORT, /** This approximant (ExternalPython) has spin support determined by the external python module on a case-by-case basis **/
   LAL_SIM_INSPIRAL_NUMSPINSUPPORT	/**< Number of elements in enum, useful for checking bounds */
 } SpinSupport;

typedef enum tagSpinFreq {
  LAL_SIM_INSPIRAL_SPINS_F_REF,   /** These approximants are parameterized by the spins at f_ref */
  LAL_SIM_INSPIRAL_SPINS_FLOW,      /** These approximants are parameterized by the spins at flow */
  LAL_SIM_INSPIRAL_SPINS_NONPRECESSING, /** These approximants have nonprecessing spins */
  LAL_SIM_INSPIRAL_SPINS_CASEBYCASE, /** These approximants (NR waveforms) have spins parameterized at different frequencies on a case-by-case basis **/
  LAL_SIM_INSPIRAL_NUMSPINFREQ  /**< Number of elements in enum, useful for checking bounds */
} SpinFreq;

typedef enum tagAllowZeroMinFreq {
  LAL_SIM_INSPIRAL_ALLOW_ZERO_FMIN,   /** These approximants allow f_min=0, which means the full length of the available waveform is returned. */
  LAL_SIM_INSPIRAL_DISALLOW_ZERO_FMIN,   /** These approximants do not allow f_min=0. This is set as default. */
  LAL_SIM_INSPIRAL_NUMZEROFMIN  /**< Number of elements in enum, useful for checking bounds */
} AllowZeroMinFreq;

typedef enum tagTestGRaccept {
  LAL_SIM_INSPIRAL_NO_TESTGR_PARAMS,   /** These approximants cannot accept testGR params as input params */
  LAL_SIM_INSPIRAL_TESTGR_PARAMS,      /** These approximants accept testGR params as input params */
  LAL_SIM_INSPIRAL_CASEBYCASE_TESTGR_PARAMS, /** This approximant (ExternalPython) accept testGR parameters depending on the external python module loaded **/
  LAL_SIM_INSPIRAL_NUM_TESTGR_ACCEPT  /**< Number of elements in enum, useful for checking bounds */
 } TestGRaccept;


/**
 * Structure for passing around PN phasing coefficients.
 * For use with the TaylorF2 waveform.
 */
#define PN_PHASING_SERIES_MAX_ORDER 15
typedef struct tagPNPhasingSeries
{
    REAL8 v[PN_PHASING_SERIES_MAX_ORDER+1];
    REAL8 vlogv[PN_PHASING_SERIES_MAX_ORDER+1];
    REAL8 vlogvsq[PN_PHASING_SERIES_MAX_ORDER+1];
    REAL8 vneg[PN_PHASING_SERIES_MAX_ORDER+1];
}
PNPhasingSeries;

/** @} */

/* general waveform switching generation routines  */

int XLALSimInspiralChooseTDWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 m1, const REAL8 m2, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 meanPerAno, const REAL8 deltaT, const REAL8 f_min, REAL8 f_ref, LALDict *params, const Approximant approximant);
int XLALSimInspiralChooseTDWaveformOLD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 m1, const REAL8 m2, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 meanPerAno, const REAL8 deltaT, const REAL8 f_min, REAL8 f_ref, const REAL8 lambda1, const REAL8 lambda2, const REAL8 dQuadParam1, const REAL8 dQuadParam2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, const int phaseO, const Approximant approximant);
int XLALSimInspiralChooseFDWaveform(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8 m1, const REAL8 m2, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 meanPerAno,  const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, REAL8 f_ref, LALDict *LALpars, const Approximant approximant);
int XLALSimInspiralChooseFDWaveformOLD(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8 m1, const REAL8 m2, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 meanPerAno,  const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, REAL8 f_ref, const REAL8 lambda1, const REAL8 lambda2, const REAL8 dQuadParam1, const REAL8 dQuadParam2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, const Approximant approximant);
int XLALSimInspiralPolarizationsFromChooseFDModes(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8 m1, const REAL8 m2, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 meanPerAno,  const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, REAL8 f_ref, LALDict *LALpars, const Approximant approximant);
int XLALSimInspiralTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 distance, REAL8 inclination, REAL8 phiRef, REAL8 longAscNodes, REAL8 eccentricity, REAL8 meanPerAno, REAL8 deltaT, REAL8 f_min, REAL8 f_ref, LALDict *LALparams, Approximant approximant);
SphHarmTimeSeries * XLALSimInspiralTDModesFromPolarizations(REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 distance, REAL8 phiRef, REAL8 longAscNodes, REAL8 eccentricity, REAL8 meanPerAno, REAL8 deltaT, REAL8 f_min, REAL8 f_ref, LALDict *LALparams, Approximant approximant);
int XLALSimInspiralFD(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 distance, REAL8 inclination, REAL8 phiRef, REAL8 longAscNodes, REAL8 eccentricity, REAL8 meanPerAno, REAL8 deltaF, REAL8 f_min, REAL8 f_max, REAL8 f_ref, LALDict *LALparams, Approximant approximant);
int XLALSimInspiralChooseWaveform(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, const REAL8 m1, const REAL8 m2, const REAL8 s1x, const REAL8 s1y, const REAL8 s1z, const REAL8 s2x, const REAL8 s2y, const REAL8 s2z, const REAL8 inclination, const REAL8 phiRef, const REAL8 distance, const REAL8 longAscNodes, const REAL8 eccentricity, const REAL8 meanPerAno, const REAL8 deltaT, const REAL8 f_min, const REAL8 f_ref, LALDict *LALpars, const Approximant approximant);
/* DEPRECATED */

/* general waveform switching mode generation routines */
SphHarmTimeSeries *XLALSimInspiralChooseTDModes(REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 f_min, REAL8 f_ref, REAL8 r, LALDict* LALpars, int lmax, Approximant approximant);
SphHarmFrequencySeries *XLALSimInspiralChooseFDModes(REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 deltaF, REAL8 f_min, REAL8 f_max, REAL8 f_ref, REAL8 phiRef, REAL8 distance, REAL8 inclination, LALDict *LALpars, Approximant approximant);
SphHarmTimeSeries *XLALSimInspiralModesTD(REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, LALDict *LALpars, int lmax, Approximant approximant);
COMPLEX16TimeSeries *XLALSimInspiralChooseTDMode(REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_ref, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralWaveformFlags *waveFlags, LALSimInspiralTestGRParam *nonGRparams, int amplitudeO, int phaseO, int l, int m, Approximant approximant);

/* routines for generating inspiral waveforms from orbital data */
int XLALSimInspiralPNPolarizationWaveforms(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, REAL8 i, int ampO);
int XLALSimInspiralPNPolarizationWaveformsFromModes(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *v, REAL8TimeSeries *phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, REAL8 i, int O);
int XLALSimInspiralPolarizationsFromSphHarmTimeSeries(REAL8TimeSeries **hp, REAL8TimeSeries **hc, SphHarmTimeSeries *hlms, REAL8 iota, REAL8 psi);
int XLALSimInspiralPolarizationsFromSphHarmFrequencySeries(COMPLEX16FrequencySeries **hp, COMPLEX16FrequencySeries **hc, SphHarmFrequencySeries *hlms, REAL8 theta, REAL8 phi);
int XLALSimInspiralPNPolarizationWaveformsEccentric(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *V, REAL8TimeSeries *Ecc, REAL8TimeSeries *U, REAL8TimeSeries *Phi, REAL8 m1, REAL8 m2, REAL8 r, REAL8 i, int ampO, int ph_O);
int XLALSimInspiralPrecessingPolarizationWaveforms(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8TimeSeries *LNhatx, REAL8TimeSeries *LNhaty, REAL8TimeSeries *LNhatz, REAL8TimeSeries *E1x, REAL8TimeSeries *E1y, REAL8TimeSeries *E1z, REAL8 m1, REAL8 m2, REAL8 r, INT4 ampO);
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
int XLALSimInspiralGetSpinFreqFromApproximant(Approximant approx);
int XLALSimInspiralGetAllowZeroMinFreqFromApproximant(Approximant approx);
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
int XLALSimInspiralTransformPrecessingNewInitialConditions(REAL8 *incl, REAL8 *S1x, REAL8 *S1y, REAL8 *S1z, REAL8 *S2x, REAL8 *S2y, REAL8 *S2z, const REAL8 thetaJN, const REAL8 phiJL, const REAL8 theta1, const REAL8 theta2, const REAL8 phi12, const REAL8 chi1, const REAL8 chi2, const REAL8 m1, const REAL8 m2, const REAL8 fRef, REAL8 phiRef);
int XLALSimInspiralTransformPrecessingObsoleteInitialConditions(REAL8 *incl, REAL8 *S1x, REAL8 *S1y, REAL8 *S1z, REAL8 *S2x, REAL8 *S2y, REAL8 *S2z, REAL8 thetaJN, REAL8 phiJL, REAL8 theta1, REAL8 theta2, REAL8 phi12, REAL8 chi1, REAL8 chi2, REAL8 m1, REAL8 m2, REAL8 fRef);
int XLALSimInspiralTransformPrecessingWvf2PE( REAL8 *thetaJN, REAL8 *phiJL, REAL8 *theta1,  	REAL8 *theta2, REAL8 *phi12, REAL8 *chi1, REAL8 *chi2,const REAL8 incl, const REAL8 S1x,	const REAL8 S1y,const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 m1, const REAL8 m2,	const REAL8 fRef, const REAL8 phiRef);

/* routines for generating PN modes based on orbital data */
/* in module LALSimInspiralPNMode.c */

COMPLEX16TimeSeries *XLALCreateSimInspiralPNModeCOMPLEX16TimeSeries(REAL8TimeSeries *v, REAL8TimeSeries *phi, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 r, int O, int l, int m);
COMPLEX16TimeSeries *XLALCreateSimInspiralPNModeCOMPLEX16TimeSeriesLALConvention(REAL8TimeSeries *v, REAL8TimeSeries *phi, REAL8 m1, REAL8 m2, REAL8 r, int O, int l, int m);

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

INT4 XLALSimInspiralSpinPNMode2m(SphHarmTimeSeries **hlm, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *LNhx, REAL8TimeSeries *LNhy, REAL8TimeSeries *LNhz, REAL8TimeSeries *e1x, REAL8TimeSeries *e1y, REAL8TimeSeries *e1z, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8 m1, REAL8 m2, REAL8 distance, int ampO);
INT4 XLALSimInspiralSpinPNMode3m(SphHarmTimeSeries **hlm, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *LNhx, REAL8TimeSeries *LNhy, REAL8TimeSeries *LNhz, REAL8TimeSeries *e1x, REAL8TimeSeries *e1y, REAL8TimeSeries *e1z, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8 m1, REAL8 m2, REAL8 distance, int ampO);
INT4 XLALSimInspiralSpinPNMode4m(SphHarmTimeSeries **hlm, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *LNhx, REAL8TimeSeries *LNhy, REAL8TimeSeries *LNhz, REAL8TimeSeries *e1x, REAL8TimeSeries *e1y, REAL8TimeSeries *e1z, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8 m1, REAL8 m2, REAL8 distance, int ampO);


/* TaylorT1 functions */
/* in module LALSimInspiralTaylorT1.c */

int XLALSimInspiralTaylorT1PNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT1PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT1PNModes(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT1PNMode(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);
int XLALSimInspiralTaylorT1PN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT1PNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);


/* TaylorT2 functions */
/* in module LALSimInspiralTaylorT2.c */

int XLALSimInspiralTaylorT2PNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT2PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT2PNModes(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT2PNMode(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);
int XLALSimInspiralTaylorT2PN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT2PNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);


/* TaylorT3 functions */
/* in module LALSimInspiralTaylorT3.c */

int XLALSimInspiralTaylorT3PNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT3PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT3PNModes(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT3PNMode(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);

int XLALSimInspiralTaylorT3PN(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT3PNRestricted(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);


/* TaylorT4 functions */
/* in module LALSimInspiralTaylorT4.c */

int XLALSimInspiralTaylorT4PNEvolveOrbit(REAL8TimeSeries **v, REAL8TimeSeries **phi, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int O);
int XLALSimInspiralTaylorT4PNGenerator(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 i, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO);
SphHarmTimeSeries *XLALSimInspiralTaylorT4PNModes(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int lmax);
COMPLEX16TimeSeries *XLALSimInspiralTaylorT4PNMode(REAL8 v0, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 fRef, REAL8 r, REAL8 lambda1, REAL8 lambda2, LALSimInspiralTidalOrder tideO, int amplitudeO, int phaseO, int l, int m);
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
int XLALSimInspiralTaylorF2AlignedPhasing(PNPhasingSeries **pfa, const REAL8 m1, const REAL8 m2, const REAL8 chi1, const REAL8 chi2, LALDict *extraPars);
int XLALSimInspiralTaylorF2AlignedPhasingArray(REAL8Vector **phasingvals, REAL8Vector mass1, REAL8Vector mass2, REAL8Vector chi1, REAL8Vector chi2, REAL8Vector lambda1, REAL8Vector lambda2, REAL8Vector dquadmon1, REAL8Vector dquadmon2);
int XLALSimInspiralTaylorF2Core(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, const REAL8 phi_ref, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_ref, const REAL8 shft, const REAL8 r, LALDict *LALparams, PNPhasingSeries *pfaP);

int XLALSimInspiralTaylorF2(COMPLEX16FrequencySeries **htilde, const REAL8 phi_ref, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 S1z, const REAL8 S2z, const REAL8 fStart, const REAL8 fEnd, const REAL8 f_ref, const REAL8 r, LALDict *LALpars);

/* TaylorF2Ecc functions */
/* in module LALSimInspiralTaylorF2Ecc.c */
int XLALSimInspiralTaylorF2CoreEcc(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, const REAL8 phi_ref, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 f_ref, const REAL8 shft, const REAL8 r, const REAL8 eccentricity, LALDict *LALparams, PNPhasingSeries *pfaP);
int XLALSimInspiralTaylorF2Ecc(COMPLEX16FrequencySeries **htilde, const REAL8 phi_ref, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 S1z, const REAL8 S2z, const REAL8 fStart, const REAL8 fEnd, const REAL8 f_ref, const REAL8 r, const REAL8 eccentricity, LALDict *LALparams);

/* TaylorF2NLPhase functions */
/* in module LALSimInspiralTaylorF2NLTides.c */

int XLALSimInspiralTaylorF2AlignedPhasingNLTides(PNPhasingSeries **pfa, const REAL8 m1, const REAL8 m2, const REAL8 chi1, const REAL8 chi2, LALDict *extraPars);
int XLALSimInspiralTaylorF2NLPhase(REAL8Sequence *dphi, const REAL8Sequence *freqs, const REAL8 Anl1, const REAL8 n1, const REAL8 fo1, const REAL8 m1_SI, const REAL8 Anl2, const REAL8 n2, const REAL8 fo2, const REAL8 m2_SI);
int XLALSimInspiralTaylorF2CoreNLTides(COMPLEX16FrequencySeries **htilde, const REAL8Sequence *freqs, const REAL8 phi_ref, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 S1z, const REAL8 S2z, const REAL8 f_ref, const REAL8 shft, const REAL8 r, LALDict *LALparams);
int XLALSimInspiralTaylorF2NLTides(COMPLEX16FrequencySeries **htilde, const REAL8 phi_ref, const REAL8 deltaF, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 S1z, const REAL8 S2z, const REAL8 fStart, const REAL8 fEnd, const REAL8 f_ref, const REAL8 r, LALDict *LALpars);


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
  REAL8 wdot4S1S2Avg, wdot4S1OS2OAvg; ///< non-dynamical, averaged 2PN S1-S2 terms
  REAL8 wdot4S1S1Avg,wdot4S1OS1OAvg, wdot4S2S2Avg,wdot4S2OS2OAvg; ///< non-dynamical, averaged self Spin^2 2PN correction
  REAL8 wdot4QMS1S1Avg,wdot4QMS1OS1OAvg, wdot4QMS2S2Avg, wdot4QMS2OS2OAvg; ///< non-dynamical,averaged self Spin^2 2PN quadrupole-monopole corrections
  REAL8 wdot5S1O, wdot5S2O; ///< non-dynamical 2.5PN SO corrections
  REAL8 wdot6S1O, wdot6S2O; ///< non-dynamical, 3PN SO corrections
  REAL8 wdot6QMS1S1, wdot6QMS1nS1n, wdot6QMS1vS1v; ///< non-dynamical 3PN quadrupole-monopole (S_1)^2 corrections
  REAL8 wdot6QMS2S2, wdot6QMS2nS2n, wdot6QMS2vS2v; ///< non-dynamical 3PN quadrupole-monopole (S_2)^2 corrections
  REAL8 wdot6S1S2Avg, wdot6S1OS2OAvg, wdot6S1S1Avg, wdot6S1OS1OAvg, wdot6S2S2Avg, wdot6S2OS2OAvg; ///< non-dynamical, averaged self Spin^2 3PN correction
  REAL8 wdot6QMS1S1Avg, wdot6QMS1OS1OAvg, wdot6QMS2S2Avg, wdot6QMS2OS2OAvg; ///< non-dynamical, averaged self QM Spin^2 3PN correction
  REAL8 wdot7S1O, wdot7S2O; ///< non-dynamical 3.5PN SO corrections
  REAL8 wdottidal10;	///< leading order tidal correction
  REAL8 wdottidal12;	///< next to leading order tidal correction
  REAL8 Ecoeff[LAL_MAX_PN_ORDER]; ///< coeffs. of PN corrections to energy
  REAL8 E3S1O, E3S2O; ///< non-dynamical 1.5PN SO corrections
  REAL8 E4S1S2Avg,E4S1OS2OAvg; ///< non-dynamical, averaged 2PN S1-S2 corrections
  REAL8 E4QMS1S1Avg, E4QMS1OS1OAvg, E4QMS2S2Avg, E4QMS2OS2OAvg;///< non-dynamical, averaged (Spin)^2 2PN quadrupole-monopole correction
  REAL8 E5S1O, E5S2O; ///< non-dynamical 2.5PN SO corrections
  REAL8 E6S1S2Avg, E6S1OS2OAvg, E6S1S1Avg, E6S1OS1OAvg, E6S2S2Avg, E6S2OS2OAvg; ///< non-dynamical 3PN self-spin^2 averaged corrections
  REAL8 E6QMS1S1, E6QMS1nS1n, E6QMS1vS1v; ///< non-dynamical 3PN quadrupole-monopole spin^2 corrections
  REAL8 E6QMS2S2, E6QMS2nS2n, E6QMS2vS2v; ///< non-dynamical 3PN quadrupole-monopole spin^2 corrections
  REAL8 E6QMS1S1Avg, E6QMS1OS1OAvg, E6QMS2S2Avg, E6QMS2OS2OAvg; ///< non-dynamical 3PN quadrupole-monopole averaged spin^2 corrections
  REAL8 E7S1O, E7S2O; ///< non-dynamical 3.5PN SO corrections
  REAL8 Etidal10; ///< leading order 5PN tidal correction to energy
  REAL8 Etidal12; ///< next to leading order 6PN tidal correction to energy
  REAL8 dEdvnewt;
  REAL8 Fcoeff[LAL_MAX_PN_ORDER];///<FluxCoeff
  REAL8 Fnewt; ///<newtonian term in Flux
  REAL8 Flogcoeff; ///<log coeff in flux
  REAL8 F3S1O, F3S2O;  ///< Coefficient of S.LN terms
  REAL8 F4S1S2Avg, F4S1OS2OAvg;///< Coefficients of averaged S1.S2 terms
  REAL8 F4S1S1Avg,F4S1OS1OAvg, F4S2S2Avg,F4S2OS2OAvg;///< Coefficient of averaged Spin^2 terms
  REAL8 F4QMS1S1,F4QMS1nS1n,F4QMS1vS1v, F4QMS2S2,F4QMS2nS2n,F4QMS2vS2v; ///< Averaged coefficient of quad-monop. Spin^2 terms
  REAL8 F4QMS1S1Avg,F4QMS1OS1OAvg, F4QMS2S2Avg,F4QMS2OS2OAvg; ///< Averaged coefficient of quad-monop. Spin^2 terms
  REAL8 F5S1O;  ///< Coefficient of (S1.LN)
  REAL8 F5S2O;  ///< Coefficient of (S1.LN) term
  REAL8 F6S1O, F6S2O; ///< Coefficient of (Si.LN) term
  REAL8 F6S1S2Avg, F6S1OS2OAvg, F6S1S1Avg, F6S1OS1OAvg, F6S2S2Avg, F6S2OS2OAvg; ///< Coefficients of Si^2 term avged
  REAL8 F6QMS1S1Avg, F6QMS1OS1OAvg, F6QMS2S2Avg, F6QMS2OS2OAvg; ///< Coefficients of quad-monop. S1.S1 terms-avged
  REAL8 F7S1O; ///< Coefficients of S1.LN term
  REAL8 F7S2O; ///< Coefficients of S2.LN term
  REAL8 Ftidal10;     ///< leading order 5PN tidal correction
  REAL8 Ftidal12;     ///< next-to-leading order 6PN tidal correction
  REAL8 S1dot3; ///< coeff of LNxS1 term in S1dot
  REAL8 S2dot3; ///< coeff of LNxS2 term in S2dot
  REAL8 S1dot4S2Avg,S1dot4S2OAvg,S1dot4QMS1OAvg,S2dot4QMS2OAvg; ///< coeff of averaged S2xS1 and quad-monop. (L.Si) LxSi terms in Sidot
  REAL8 S1dot5; ///< coeff of LNxS1 term in S1dot
  REAL8 S2dot5; ///< coeff of LNxS2 term in S2dot
  REAL8 S1dot6S2Avg,S1dot6S2OAvg,S1dot6S1OAvg,S1dot6QMS1OAvg; // 6PN S1dot avged-coefficients
  REAL8 S2dot6S1Avg,S2dot6S1OAvg,S2dot6S2OAvg,S2dot6QMS2OAvg; // 6PN S2dot avged-coefficients
  REAL8 S1dot7S2;// Coefficient of S1 x S2 in S1dot
  REAL8 S2dot7S1;// Coefficient of S1 x S2 in S2dot
  REAL8 omegashiftS1,omegashiftS2;// non-dynamical coefficients of \omega shift wrt \dot\phi, see eq. (34) of https://dcc.ligo.org/LIGO-T1500554
  REAL8 fStart; ///< starting GW frequency of integration
  REAL8 fEnd; ///< ending GW frequency of integration
  INT4 phaseO; ///< Twice PN order of GW-phase
  LALSimInspiralSpinOrder spinO; ///< Twice PN order of included spin effects
  LALSimInspiralTidalOrder tideO;///< Twice PN order of included tidal effects
  REAL8 prev_domega; ///< Previous value of domega/dt used in stopping test
  INT4 lscorr; ///< Flag for including spin corrections to orb. ang. mom.
  INT4 phenomtp; ///< Flag for using spinO=7 and not spinO=6 teems with orbital-averaged quantities for phenomtphm approx
} XLALSimInspiralSpinTaylorTxCoeffs;

int XLALSimInspiralSpinTaylorPNEvolveOrbit(REAL8TimeSeries **V, REAL8TimeSeries **Phi, REAL8TimeSeries **S1x, REAL8TimeSeries **S1y, REAL8TimeSeries **S1z, REAL8TimeSeries **S2x, REAL8TimeSeries **S2y, REAL8TimeSeries **S2z, REAL8TimeSeries **LNhatx, REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz, REAL8TimeSeries **E1x, REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, INT4 lscorr, Approximant approx);
int XLALSimInspiralSpinTaylorPNEvolveOrbitOnlyFinal(REAL8TimeSeries **V, REAL8TimeSeries **Phi, REAL8TimeSeries **S1x, REAL8TimeSeries **S1y, REAL8TimeSeries **S1z, REAL8TimeSeries **S2x, REAL8TimeSeries **S2y, REAL8TimeSeries **S2z, REAL8TimeSeries **LNhatx, REAL8TimeSeries **LNhaty, REAL8TimeSeries **LNhatz, REAL8TimeSeries **E1x, REAL8TimeSeries **E1y, REAL8TimeSeries **E1z, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, INT4 lscorr, Approximant approx);
int XLALSimInspiralSpinTaylorT1(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, LALDict *LALparams);
int XLALSimInspiralSpinTaylorT4(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, LALDict *LALParams);
int XLALSimInspiralSpinTaylorT5(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, LALDict *LALparams);
int XLALSimInspiralSpinTaylorT5duplicate(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 incAngle, int phaseO, int amplitudeO);
int XLALSimInspiralSpinTaylorDriver(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8TimeSeries **Vout, REAL8TimeSeries **Phiout, REAL8TimeSeries **S1xout, REAL8TimeSeries **S1yout, REAL8TimeSeries **S1zout, REAL8TimeSeries **S2xout, REAL8TimeSeries **S2yout, REAL8TimeSeries **S2zout, REAL8TimeSeries **LNhxout, REAL8TimeSeries **LNhyout, REAL8TimeSeries **LNhzout, REAL8TimeSeries **E1xout, REAL8TimeSeries **E1yout, REAL8TimeSeries **E1zout, REAL8 phiRef, REAL8 deltaT, REAL8 m1_SI, REAL8 m2_SI, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z,	REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, LALDict *LALparams, Approximant approx);
int XLALSimInspiralSpinTaylorOrbitalDriver(REAL8TimeSeries **Vout, REAL8TimeSeries **Phiout, REAL8TimeSeries **S1xout, REAL8TimeSeries **S1yout, REAL8TimeSeries **S1zout, REAL8TimeSeries **S2xout, REAL8TimeSeries **S2yout, REAL8TimeSeries **S2zout, REAL8TimeSeries **LNhxout, REAL8TimeSeries **LNhyout, REAL8TimeSeries **LNhzout, REAL8TimeSeries **E1xout, REAL8TimeSeries **E1yout, REAL8TimeSeries **E1zout, REAL8 phiRef, REAL8 deltaT, REAL8 m1_SI, REAL8 m2_SI, REAL8 fStart, REAL8 fRef, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, LALDict *LALparams, Approximant approx);
int XLALSimInspiralSpinTaylorT1OLD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT,  REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams, int phaseO, int amplitudeO);
int XLALSimInspiralSpinTaylorT4OLD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT,  REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams, int phaseO, int amplitudeO);
int XLALSimInspiralSpinTaylorT4PTFQVecs(REAL8TimeSeries **Q1, REAL8TimeSeries **Q2, REAL8TimeSeries **Q3, REAL8TimeSeries **Q4, REAL8TimeSeries **Q5, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 chi1, REAL8 kappa, REAL8 fStart, REAL8 lambda1, REAL8 lambda2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, int phaseO);
int XLALSimInspiralSpinTaylorT5Fourier(COMPLEX16FrequencySeries **hplus, COMPLEX16FrequencySeries **hcross, REAL8 fMin, REAL8 fMax, REAL8 deltaF, INT4 kMax, REAL8 phiRef, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams, INT4 phaseO, INT4 amplitudeO, INT4 phiRefAtEnd);
int XLALSimInspiralSpinTaylorT4Fourier(COMPLEX16FrequencySeries **hplus, COMPLEX16FrequencySeries **hcross, REAL8 fMin, REAL8 fMax, REAL8 deltaF, INT4 kMax, REAL8 phiRef, REAL8 v0, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALDict *LALparams, INT4 phaseO, INT4 amplitudeO, INT4 phiRefAtEnd);
int XLALSimInspiralSpinTaylorF2(COMPLEX16FrequencySeries **hplus_out, COMPLEX16FrequencySeries **hcross_out, REAL8 phi_ref, REAL8 deltaF, REAL8 m1_SI, REAL8 m2_SI, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, const REAL8 fStart, const REAL8 fEnd, const REAL8 f_ref, const REAL8 r, LALDict *moreParams, INT4 phaseO, INT4 amplitudeO);
int XLALSimInspiralPrecessingPTFQWaveforms(REAL8TimeSeries **Q1, REAL8TimeSeries **Q2, REAL8TimeSeries **Q3, REAL8TimeSeries **Q4, REAL8TimeSeries **Q5, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8TimeSeries *LNhatx, REAL8TimeSeries *LNhaty, REAL8TimeSeries *LNhatz, REAL8TimeSeries *E1x, REAL8TimeSeries *E1y, REAL8TimeSeries *E1z, REAL8 m1, REAL8 m2, REAL8 r);
int XLALSimInspiralInitialConditionsPrecessingApproxs(REAL8 *inc, REAL8 *S1x, REAL8 *S1y, REAL8 *S1z, REAL8 *S2x, REAL8 *S2y, REAL8 *S2z, const REAL8 inclIn, const REAL8 S1xIn, const REAL8 S1yIn, const REAL8 S1zIn, const REAL8 S2xIn, const REAL8 S2yIn, const REAL8 S2zIn, const REAL8 m1, const REAL8 m2, const REAL8 fRef, const REAL8 phiRef, LALSimInspiralFrameAxis axisChoice);
INT4 XLALSimInspiralSpinDerivativesAvg(REAL8 *dLNhx, REAL8 *dLNhy, REAL8 *dLNhz, REAL8 *dE1x, REAL8 *dE1y, REAL8 *dE1z, REAL8 *dS1x, REAL8 *dS1y, REAL8 *dS1z, REAL8 *dS2x, REAL8 *dS2y, REAL8 *dS2z, const REAL8 v, const REAL8 LNhx, const REAL8 LNhy, const REAL8 LNhz, const REAL8 E1x, const REAL8 E1y, const REAL8 E1z, const REAL8 S1x, const REAL8 S1y, const REAL8 S1z, const REAL8 S2x, const REAL8 S2y, const REAL8 S2z, const REAL8 LNhdotS1, const REAL8 LNhdotS2, XLALSimInspiralSpinTaylorTxCoeffs *params);
INT4 XLALSimInspiralSpinTaylorT4DerivativesAvg(REAL8 t, const REAL8 values[], REAL8 dvalues[], void *mparams);
int XLALSimInspiralSpinTaylorT5Setup(XLALSimInspiralSpinTaylorTxCoeffs **params, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, INT4 lscorr);
int XLALSimInspiralSpinTaylorT1Setup(XLALSimInspiralSpinTaylorTxCoeffs **params, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, INT4 lscorr);
INT4 XLALSimInspiralSpinTaylorT4Setup(XLALSimInspiralSpinTaylorTxCoeffs **params, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO, INT4 lscorr, INT4 phenomtp);
INT4 XLALSimSpinTaylorEnergySpinDerivativeSetup(XLALSimInspiralSpinTaylorTxCoeffs **params, const REAL8 m1_SI, const REAL8 m2_SI, const REAL8 fStart, const REAL8 fEnd, const LALSimInspiralSpinOrder  spinO, const LALSimInspiralTidalOrder tideO, const INT4 phaseO, const REAL8 lambda1, const REAL8 lambda2, const REAL8 quadparam1, const REAL8 quadparam2, const INT4 lscorr, const INT4 phenomtp);
INT4 XLALSimInspiralSpinTaylorHlmModesFromOrbit(SphHarmTimeSeries **hlm, REAL8TimeSeries *V, REAL8TimeSeries *Phi, REAL8TimeSeries *LNhx, REAL8TimeSeries *LNhy, REAL8TimeSeries *LNhz, REAL8TimeSeries *e1x, REAL8TimeSeries *e1y, REAL8TimeSeries *e1z, REAL8TimeSeries *S1x, REAL8TimeSeries *S1y, REAL8TimeSeries *S1z, REAL8TimeSeries *S2x, REAL8TimeSeries *S2y, REAL8TimeSeries *S2z, REAL8 m1_SI, REAL8 m2_SI, REAL8 distance, int ampO, LALValue *modearray);
INT4 XLALSimInspiralSpinTaylorHlmModes(SphHarmTimeSeries **hlm, REAL8 phiRef, REAL8 dT, REAL8 m1_SI, REAL8 m2_SI, REAL8 fStart, REAL8 fRef, REAL8 dist_SI, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, int ampO, LALValue *modearray, LALDict *LALparams, Approximant approx);
INT4 XLALSimInspiralSetEnergyPNTermsAvg(REAL8 *Espin3, REAL8 *Espin4, REAL8 *Espin5, REAL8 *Espin6, REAL8 *Espin7, XLALSimInspiralSpinTaylorTxCoeffs *params, const REAL8 LNhdotS1, const REAL8 LNhdotS2, const REAL8 S1sq, const REAL8 S2sq, const REAL8 S1dotS2);

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

int XLALSimInspiralSpinDominatedWaveformInterfaceTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 deltaT, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 D, REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 incl, int phaseO, int amplitudeO, REAL8 phiRef);
int XLALSimInspiralSpinDominatedWaveformDriver(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 totalmass, REAL8 nu, REAL8 chi1, REAL8 D, REAL8 kappa1, REAL8 beta1, REAL8 theta, REAL8 fStart, REAL8 fRef, int phaseO, int amplitudeO, REAL8 deltaT, REAL8 phiRef, REAL8 phin0, REAL8 polarizationangle);


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

REAL8 XLALSimInspiralfLow2fStart(REAL8 fLow, INT4 ampOrder, INT4 approximant);

/* NRSur4d2s functions */
int XLALSimNRSur4d2s(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, REAL8 phiRef, REAL8 deltaF, REAL8 fLow, REAL8 fHigh, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z);
int XLALSimNRSur4d2sFrequencySequence(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, const REAL8Sequence *freqs, REAL8 phiRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z);

/* waveform tapering routines */
/* in module LALSimInspiralWaveformTaper.c */

int XLALSimInspiralREAL4WaveTaper(REAL4Vector *signalvec, LALSimInspiralApplyTaper bookends);
int XLALSimInspiralREAL8WaveTaper(REAL8Vector *signalvec, LALSimInspiralApplyTaper bookends);

/* in module LALSimInspiralTEOBResumROM.c */

int XLALSimInspiralTEOBResumROM(REAL8TimeSeries **hPlus, REAL8TimeSeries **hCross, REAL8 phiRef, REAL8 deltaT, REAL8 fLow, REAL8 fRef, REAL8 distance, REAL8 inclination, REAL8 m1SI, REAL8 m2SI, REAL8 lambda1, REAL8 lambda2);

int XLALSimInspiralSetQuadMonParamsFromLambdas(LALDict *LALpars);

/**
 * Evaluates the NRHybSur3dq8 surrogate model and combines different modes to
 * obtain the plus and cross polarizations.
 * In module LALSimIMRNRHybSur3dq8.c
 */
INT4 XLALSimIMRNRHybSur3dq8Polarizations(
    REAL8TimeSeries **hplus,        /**<Output: \f$h_+\f$ polarization. */
    REAL8TimeSeries **hcross,       /**<Output: \f$h_{\times}\f$ polarization.*/
    REAL8 phiRef,                   /**< azimuthal angle for Ylms */
    REAL8 inclination,              /**< Inclination angle. */
    REAL8 deltaT,                   /**< Sampling interval (s). */
    REAL8 m1,                       /**< Mass of Bh1 (kg). */
    REAL8 m2,                       /**< Mass of Bh2 (kg). */
    REAL8 distance,                 /**< Distance of source (m). */
    REAL8 fMin,                     /**< Start GW frequency (Hz). */
    REAL8 fRef,                     /**< Reference GW frequency (Hz). */
    REAL8 chi1z,                    /**< Dimensionless spin of Bh1. */
    REAL8 chi2z,                    /**< Dimensionless spin of Bh2. */
    LALDict* LALparams             /**< Dict with extra parameters */
);

/**
 * Evaluates the NRHybSur3dq8 surrogate model and returns all required modes.
 * In module LALSimIMRNRHybSur3dq8.c
 */
SphHarmTimeSeries *XLALSimIMRNRHybSur3dq8Modes(
    REAL8 deltaT,                   /**< Sampling interval (s). */
    REAL8 m1,                       /**< Mass of Bh1 (kg). */
    REAL8 m2,                       /**< Mass of Bh2 (kg). */
    REAL8 chi1z,                    /**< Dimensionless spin of Bh1. */
    REAL8 chi2z,                    /**< Dimensionless spin of Bh2. */
    REAL8 fMin,                     /**< Start GW frequency (Hz). */
    REAL8 fRef,                     /**< Reference GW frequency (Hz). */
    REAL8 distance,                 /**< Distance of source (m). */
    LALDict* LALparams             /**< Dict with extra parameters */
);

/* routine for checking Lorentz violation */
int XLALSimLorentzInvarianceViolationTerm(COMPLEX16FrequencySeries **hptilde, COMPLEX16FrequencySeries **hctilde, REAL8 m1, REAL8 m2, REAL8 r, LALDict *LALparams);

/** Incoplete type for waveform generator */
struct tagLALSimInspiralGenerator;
typedef struct tagLALSimInspiralGenerator LALSimInspiralGenerator;


/* legacy generator templates */
extern const LALSimInspiralGenerator lalEOBNRv2HMGeneratorTemplate;
extern const LALSimInspiralGenerator lalEOBNRv2HM_ROMGeneratorTemplate;
extern const LALSimInspiralGenerator lalEOBNRv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalEOBNRv2_ROMGeneratorTemplate;
extern const LALSimInspiralGenerator lalEccentricFDGeneratorTemplate;
extern const LALSimInspiralGenerator lalEccentricTDGeneratorTemplate;
extern const LALSimInspiralGenerator lalHGimriGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomAGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomBGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomCGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomDGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomD_NRTidalGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomD_NRTidalv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomHMGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomNSBHGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomPGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomPv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomPv2_NRTidalGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomPv2_NRTidalv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomPv3HMGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomPv3GeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomTHMGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomTPHMGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomTPGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomTGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomXASGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomXHMGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomXPHMGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomXO4aGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomXPGeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomXAS_NRTidalv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalIMRPhenomXP_NRTidalv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalLackey_Tidal_2013_SEOBNRv2_ROMGeneratorTemplate;
extern const LALSimInspiralGenerator lalNRHybSur3dq8GeneratorTemplate;
extern const LALSimInspiralGenerator lalNRSur4d2sGeneratorTemplate;
extern const LALSimInspiralGenerator lalNRSur7dq2GeneratorTemplate;
extern const LALSimInspiralGenerator lalNRSur7dq4GeneratorTemplate;
extern const LALSimInspiralGenerator lalNR_hdf5GeneratorTemplate;
extern const LALSimInspiralGenerator lalPhenSpinTaylorRDGeneratorTemplate;
extern const LALSimInspiralGenerator lalPhenSpinTaylorGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv1GeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv1_ROM_DoubleSpinGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv1_ROM_EffectiveSpinGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv2TGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv2_ROM_DoubleSpinGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv2_ROM_DoubleSpin_HIGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv2_ROM_EffectiveSpinGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv2_optGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv3GeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv3_optGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv3_opt_rk4GeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv3_pertGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4HMGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4HM_ROMGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4PHMGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4PGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4TGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4T_surrogateGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4GeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4_ROMGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4_ROM_NRTidalGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4_ROM_NRTidalv2GeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4_ROM_NRTidalv2_NSBHGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4_optGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv4HM_PAGeneratorTemplate;
extern const LALSimInspiralGenerator lalpSEOBNRv4HM_PAGeneratorTemplate;
extern const LALSimInspiralGenerator lalSEOBNRv5_ROMGeneratorTemplate;
extern const LALSimInspiralGenerator lalSpinDominatedWfGeneratorTemplate;
extern const LALSimInspiralGenerator lalSpinTaylorF2GeneratorTemplate;
extern const LALSimInspiralGenerator lalSpinTaylorT1GeneratorTemplate;
extern const LALSimInspiralGenerator lalSpinTaylorT4FourierGeneratorTemplate;
extern const LALSimInspiralGenerator lalSpinTaylorT4GeneratorTemplate;
extern const LALSimInspiralGenerator lalSpinTaylorT5FourierGeneratorTemplate;
extern const LALSimInspiralGenerator lalSpinTaylorT5GeneratorTemplate;
extern const LALSimInspiralGenerator lalTEOBResumSGeneratorTemplate;
extern const LALSimInspiralGenerator lalTEOBResum_ROMGeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorEtGeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorF2EccGeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorF2NLTidesGeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorF2RedSpinTidalGeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorF2RedSpinGeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorF2GeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorR2F4GeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorT1GeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorT2GeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorT3GeneratorTemplate;
extern const LALSimInspiralGenerator lalTaylorT4GeneratorTemplate;

/* new generator templates */
extern const LALSimInspiralGenerator lalPythonGeneratorTemplate;

extern const LALSimInspiralGenerator *lalSimInspiralGeneratorTemplates[NumApproximants];

LALSimInspiralGenerator *XLALCreateSimInspiralGenerator(const LALSimInspiralGenerator *generator, LALDict *params);
void XLALDestroySimInspiralGenerator(LALSimInspiralGenerator *generator);

LALSimInspiralGenerator *XLALSimInspiralChooseGenerator(Approximant approx, LALDict *params);

int XLALSimInspiralGeneratorAddConditioningForApproximant(LALSimInspiralGenerator *generator, int approximant);
int XLALSimInspiralGeneratorAddStandardConditioning(LALSimInspiralGenerator *generator);

/* warning: returns a shallow pointer */
const char *XLALSimInspiralGeneratorName(LALSimInspiralGenerator *generator);

int XLALSimInspiralGenerateTDWaveform(
    REAL8TimeSeries **hplus,
    REAL8TimeSeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *generator
);

int XLALSimInspiralGenerateTDModes(
    SphHarmTimeSeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *generator
);

int XLALSimInspiralGenerateFDWaveform(
    COMPLEX16FrequencySeries **hplus,
    COMPLEX16FrequencySeries **hcross,
    LALDict *params,
    LALSimInspiralGenerator *generator
);

int XLALSimInspiralGenerateFDModes(
    SphHarmFrequencySeries **hlm,
    LALDict *params,
    LALSimInspiralGenerator *generator
);

void XLALSimInspiralParseDictionaryToChooseTDWaveform(
    REAL8 *m1,                             /**< [out] mass of companion 1 (kg) */
    REAL8 *m2,                             /**< [out] mass of companion 2 (kg) */
    REAL8 *S1x,                            /**< [out] x-component of the dimensionless spin of object 1 */
    REAL8 *S1y,                            /**< [out] y-component of the dimensionless spin of object 1 */
    REAL8 *S1z,                            /**< [out] z-component of the dimensionless spin of object 1 */
    REAL8 *S2x,                            /**< [out] x-component of the dimensionless spin of object 2 */
    REAL8 *S2y,                            /**< [out] y-component of the dimensionless spin of object 2 */
    REAL8 *S2z,                            /**< [out] z-component of the dimensionless spin of object 2 */
    REAL8 *distance,                       /**< [out] distance of source (m) */
    REAL8 *inclination,                    /**< [out] inclination of source (rad) */
    REAL8 *phiRef,                         /**< [out] reference orbital phase (rad) */
    REAL8 *longAscNodes,                   /**< [out] longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 *eccentricity,                   /**< [out] eccentrocity at reference epoch */
    REAL8 *meanPerAno,                     /**< [out] mean anomaly of periastron */
    REAL8 *deltaT,                         /**< [out] sampling interval (s) */
    REAL8 *f_min,                          /**< [out] starting GW frequency (Hz) */
    REAL8 *f_ref,                          /**< [out] reference frequency (Hz) */
    LALDict *params                        /**< Input lal dictionary with ChooseTDwaveform parameters */    
);

void XLALSimInspiralParseDictionaryToChooseTDModes(
    REAL8 *phiRef,
    REAL8 *deltaT,
    REAL8 *m1,
    REAL8 *m2,
    REAL8 *S1x,
    REAL8 *S1y,
    REAL8 *S1z,
    REAL8 *S2x,
    REAL8 *S2y,
    REAL8 *S2z,
    REAL8 *f_min,
    REAL8 *f_ref,
    REAL8 *distance,
    INT4 *lmax,
    LALDict *params
);

void XLALSimInspiralParseDictionaryToChooseFDWaveform(
    REAL8 *m1,                             /**< [out] mass of companion 1 (kg) */
    REAL8 *m2,                             /**< [out] mass of companion 2 (kg) */
    REAL8 *S1x,                            /**< [out] x-component of the dimensionless spin of object 1 */
    REAL8 *S1y,                            /**< [out] y-component of the dimensionless spin of object 1 */
    REAL8 *S1z,                            /**< [out] z-component of the dimensionless spin of object 1 */
    REAL8 *S2x,                            /**< [out] x-component of the dimensionless spin of object 2 */
    REAL8 *S2y,                            /**< [out] y-component of the dimensionless spin of object 2 */
    REAL8 *S2z,                            /**< [out] z-component of the dimensionless spin of object 2 */
    REAL8 *distance,                       /**< [out] distance of source (m) */
    REAL8 *inclination,                    /**< [out] inclination of source (rad) */
    REAL8 *phiRef,                         /**< [out] reference orbital phase (rad) */
    REAL8 *longAscNodes,                   /**< [out] longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
    REAL8 *eccentricity,                   /**< [out] eccentrocity at reference epoch */
    REAL8 *meanPerAno,                     /**< [out] mean anomaly of periastron */
    REAL8 *deltaF,                         /**< [out] frequency interval (Hz) */
    REAL8 *f_min,                          /**< [out] starting GW frequency (Hz) */
    REAL8 *f_max,                          /**< [out] ending GW frequency (Hz) */
    REAL8 *f_ref,                          /**< [out] reference frequency (Hz) */
    LALDict *params                        /**< Input lal dictionary with ChooseTDwaveform parameters **/    
);

void XLALSimInspiralParseDictionaryToChooseFDModes(
    REAL8 *m1,
    REAL8 *m2,
    REAL8 *S1x,
    REAL8 *S1y,
    REAL8 *S1z,
    REAL8 *S2x,
    REAL8 *S2y,
    REAL8 *S2z,
    REAL8 *deltaF,
    REAL8 *f_min,
    REAL8 *f_max,
    REAL8 *f_ref,
    REAL8 *phiRef,
    REAL8 *distance,
    REAL8 *inclination,
    LALDict *params
);

/*
 * some helper routines for XLALSimInspiralTD
 */
int XLALSimInspiralTDFromTD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 distance, REAL8 inclination, REAL8 phiRef, REAL8 longAscNodes, REAL8 eccentricity, REAL8 meanPerAno, REAL8 deltaT, REAL8 f_min, REAL8 f_ref, LALDict *LALparams, Approximant approximant);
int XLALSimInspiralTDFromFD(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, REAL8 m1, REAL8 m2, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, REAL8 distance, REAL8 inclination, REAL8 phiRef, REAL8 longAscNodes, REAL8 eccentricity, REAL8 meanPerAno, REAL8 deltaT, REAL8 f_min, REAL8 f_ref, LALDict *LALparams, Approximant approximant);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRAL_H */
