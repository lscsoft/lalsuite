/*
*  Copyright (C) 2007 Stas Babak, David Churches, Drew Keppel, Duncan Brown, Jolien Creighton, David McKechan, Patrick Brady, Peter Shawhan, Reinhard Prix, B.S. Sathyaprakash, Anand Sengupta, Craig Robinson , Sean Seader, Thomas Cokelaer, Riccardo Sturani,  Laszlo Vereb
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
 * \defgroup LALInspiral_h LALInspiral_h
 * \ingroup CBC_inspiral
 */

/**
\author Churches, D. K ,  B. S. Sathyaprakash,  T. Cokelaer.
\file
\ingroup LALInspiral_h

\brief %Header file for the template generation codes.

\heading{Synopsis}
\code
#include <lal/LALInspiral.h>
\endcode

This header file covers routines that are used in template generation.

*/

#ifndef _LALINSPIRAL_H
#define _LALINSPIRAL_H

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <lal/LALGSL.h>

# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>
# include <lal/SimulateCoherentGW.h>
# include <lal/GeneratePPNInspiral.h>
# include <lal/LIGOMetadataTables.h>
# include <lal/LALDatatypes.h>
# include <lal/LALComplex.h>


#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( LALINSPIRALH, "$Id$" );

# define oneby3    0.333333333333333333333333333
# define twoby3    0.666666666666666666666666667
# define fourby3   1.333333333333333333333333333
# define fiveby3   1.666666666666666666666666667
# define sevenby3  2.333333333333333333333333333
# define eightby3  2.666666666666666666666666667
# define tenby3    3.333333333333333333333333333
# define elevenby3 3.666666666666666666666666666
# define threeby8  0.375
# define fiveby8   0.625
# define threeby4  0.75
# define sevenby8  0.875
#define  ninty4by3etc  18.687902694437592603 /* (94/3 -41/31*pi*pi) */

/**\name Error Codes */ /*@{*/
#define LALINSPIRALH_ENULL           1
#define LALINSPIRALH_EMEM            2
#define LALINSPIRALH_EDIV0           3
#define LALINSPIRALH_ESIZE           4
#define LALINSPIRALH_ECHOICE         5
#define LALINSPIRALH_EORDER          6
#define LALINSPIRALH_EAPPROXIMANT    7
#define LALINSPIRALH_EPSI0           8
#define LALINSPIRALH_EPSI3           9
#define LALINSPIRALH_EALPHA         10
#define LALINSPIRALH_EFCUTOFF       11
#define LALINSPIRALH_ENOWAVEFORM    12
#define LALINSPIRALH_ESTOPPED       13
#define LALINSPIRALH_EROOTINIT      14
#define LALINSPIRALH_EFLOWER        15
#define LALINSPIRALH_EVECTOR        16
#define LALINSPIRALH_EFLOWERINJ     17
#define LALINSPIRALH_EORDERMISSING  18
#define LALINSPIRALH_EBPERR         19
#define LALINSPIRALH_ESWITCH        20
#define LALINSPIRALH_EMASSCHOICE    21

#define LALINSPIRALH_MSGENULL         "Arguments contained an unexpected null pointer"
#define LALINSPIRALH_MSGEMEM          "Memory allocation error"
#define LALINSPIRALH_MSGEDIV0         "Division by zero"
#define LALINSPIRALH_MSGESIZE         "Invalid input range"
#define LALINSPIRALH_MSGECHOICE       "Invalid choice for an input parameter"
#define LALINSPIRALH_MSGEORDER        "unknown order specified"
#define LALINSPIRALH_MSGEAPPROXIMANT  "Invalid model"
#define LALINSPIRALH_MSGEPSI0         "psi0 must be > 0"
#define LALINSPIRALH_MSGEPSI3         "psi3 must be < 0"
#define LALINSPIRALH_MSGEALPHA        "alpha must be defined positive"
#define LALINSPIRALH_MSGEFCUTOFF      "fcutoff must be defined and > 0"
#define LALINSPIRALH_MSGENOWAVEFORM   "No Waveform generated"
#define LALINSPIRALH_MSGESTOPPED      "Waveform generation stopped"
#define LALINSPIRALH_MSGEROOTINIT     "Can't find good bracket for BisectionFindRoot"
#define LALINSPIRALH_MSGEFLOWER       "fLower too low in comparison to flso"
#define LALINSPIRALH_MSGEVECTOR       "Attempting to write beyond the end of vector"
#define LALINSPIRALH_MSGEFLOWERINJ    "flower for the injection must be greater than zero"
#define LALINSPIRALH_MSGEORDERMISSING "The PN order requested is not implemented for this approximant"
#define LALINSPIRALH_MSGEBPERR        "Error in band passing signal"
#define LALINSPIRALH_MSGESWITCH       "Unknown case in switch"
#define LALINSPIRALH_MSGEMASSCHOICE   "Improper choice for massChoice"
/*@}*/


/** These are the input structures needed to solve for the mass
    ratio \f$\eta\f$ given the chirptimes \f$(\tau_0,\, \tau_2)\f$ or
    \f$(\tau_0, \, \tau_4).\f$

Here, \c t2\ \f$ = \tau_2,\f$ \c A2\ \f$ = A_2 ({\tau_0}/{A_0})^{3/5},\f$  and
\c B2\ \f$=B_2\f$,
where \f$A_0 = 5/[256 (\pi f_{s} )^{8/3}],\f$ \f$A_2 = 3715 / [64512 (\pi f_s)^2],\f$
\f$B_2 = 4620/3715.\f$

Similarly, \c t4\ \f$ = \tau_4,\f$ \c A4\ \f$ = A_4 ({\tau_0}/{A_0})^{1/5},\f$
\c B4\ \f$=B_4\f$ and \c C4\ \f$=C_4,\f$ where
where \f$A_0 = 5/[256 (\pi f_{s} )^{8/3}],\f$
\f$A_4 = 5 \times 3058673/ [128 \times 1016064  (\pi f_s)^{4/3}],\f$
\f$B_4 = 5429 \times 1016064 /(1008 \times 3058673),\f$ and \f$C_4 = 617 \times
1016064/(144 \times 3058673).\f$
*//* @{ */
typedef struct
tagEtaTau02In
{
   REAL8 t2;
   REAL8 A2;
   REAL8 B2;
} EtaTau02In;

typedef struct
tagEtaTau04In
{
   REAL8 t4;
   REAL8 A4;
   REAL8 B4;
   REAL8 C4;
} EtaTau04In;
/* @} */


/** Enum that tells which post-Newtonian order is being used.
<ul>
<li> \c LAL_PNORDER_NEWTONIAN: Newtonain order, flux and enrgy both to the lowest order.</li>
<li> \c LAL_PNORDER_HALF: same as before</li>
<li> \c LAL_PNORDER_ONE: Both energy and flux to order \f$O(v^2)\f$ beyond the Newtonian order.</li>
<li> \c LAL_PNORDER_ONE_POINT_FIVE: Energy to order \f$O(v^2)\f$ and flux to order \f$O(v^3)\f$</li>
<li> \c LAL_PNORDER_TWO: Both energy and flux to order \f$O(v^4)\f$</li>
<li> \c LAL_PNORDER_TWO_POINT_FIVE: Energy to order \f$O(v^4)\f$ and flux to order \f$O(v^5)\f$</li>
<li> \c LAL_PNORDER_THREE: Both energy and flux to order \f$O(v^6)\f$</li>
<li> \c LAL_PNORDER_THREE_POINT_FIVE: Energy to order \f$O(v^6)\f$ and flux to order \f$O(v^7)\f$</li>
<li> \c LAL_PNORDER_PSEUDIO_FOUR: Need to describe</li>
</ul>
In all cases, the gravitational wave phase (also frequency and time)
as an expansion of the gauge invariant parameter \f$v\f$ is given up to
the order specified by flux.  Note that there are certain undetermined
parameters at \c LAL_PNORDER_THREE and
\c LAL_PNORDER_THREE_POINT_FIVE. The waveform generation codes use
a specific value of those parameters while generating the wave.
*/
typedef enum {
  LAL_PNORDER_NEWTONIAN,
  LAL_PNORDER_HALF,
  LAL_PNORDER_ONE,
  LAL_PNORDER_ONE_POINT_FIVE,
  LAL_PNORDER_TWO,
  LAL_PNORDER_TWO_POINT_FIVE,
  LAL_PNORDER_THREE,
  LAL_PNORDER_THREE_POINT_FIVE,
  LAL_PNORDER_PSEUDO_FOUR,
  LAL_PNORDER_NUM_ORDER
 } LALPNOrder;





/** Enum that specifies the PN approximant to be used in computing the waveform.
<ul>
<li> \c TaylorT1: Time domain Taylor approximant in which
	the energy and flux are both kept as Taylor expansions
	and a first order ordinary differential equation is solved
	for the GW phase as a function of \f$t.\f$ Outputs a time-domain wave.</li>
<li> \c TaylorT2: Time domain Taylor approximant in which
	the phase evolution \f$\varphi(t)\f$ is obtained by iteratively
	solving post-Newtonian expansions \f$\varphi(v)\f$ and \f$t(v).\f$ Outputs a time-domain wave.</li>
<li> \c TaylorT3: Time domain Taylor approximant in which
phase is explicitly given as a function of time. Outputs a time-domain wave.</li>
<li> \c TaylorF1: The stationary phase approximation that
correctly represents, in the Fourier domain, the waveform given
by \c TaylorT1 approximant (see Ref. [\ref dis2000] for details). Outputs a frequency-domain wave.</li>
<li> \c TaylorF2: The standard stationary phase approximation. Outputs a frequency-domain wave.</li>
<li> \c PadeT1: Time-domain P-approximant. Outputs a time-domain wave.</li>
<li> \c PadeF1: Frequency-domain P-approximant (not yet implemented).</li>
<li> \c EOB: Effective one-body waveform  Outputs a time-domain wave.</li>
<li> \c BCV: Detection template family of Buonanno, Chen and
                    Vallisneri [\ref BCV03]. Outputs a frequency-domain wave.</li>
<li> \c BCVSpin: Detection template family of Buonanno, Chen and
                    Vallisneri including  spin effects[\ref BCV03b]. Outputs a frequency-domain wave.</li>
<li> \c SpinTaylorT3 Spinning case T3 models</li>
<li> \c SpinTaylor Spinning case PN models (should replace SpinTaylorT3 in the future)</li>
<li> \c SpinTaylorFrameless Spinning case PN models (replace SpinTaylor by removing the coordinate singularity)</li>
<li> \c PhenSpinTaylorRD Phenomenological waveforms, interpolating between a T4 spin-inspiral and the ringdown.</li>
<li> \c SpinQuadTaylor Spinning case PN models with quadrupole-monopole and self-spin interaction.</li>
<li> \c FindChirpSP The stationary phase templates implemented by FindChirpSPTemplate in the findchirp package (equivalent to TaylorF2 at twoPN order).</li>
<li> \c GeneratePPN The time domain templates generated by LALGeneratePPNInspiral() in the inject package (equivalent to TaylorT3 at twoPN order).</li>
<li> \c FrameFile The waveform contains arbitrary data read from a frame file.</li>
</ul>
*/
typedef enum {
   TaylorT1,
   TaylorT2,
   TaylorT3,
   TaylorF1,
   TaylorF2,
   PadeT1,
   PadeF1,
   EOB,
   BCV,
   BCVSpin,
   SpinTaylorT3,
   SpinTaylorFrameless,
   SpinTaylor,
   PhenSpinTaylorRD,
   PhenSpinTaylorRDF,
   SpinQuadTaylor,
   FindChirpSP,
   FindChirpPTF,
   GeneratePPN,
   BCVC,
   FrameFile,
   AmpCorPPN,
   NumRel,
   Eccentricity,
   EOBNR,
   IMRPhenomA,
   IMRPhenomB,
   IMRPhenomFA,
   IMRPhenomFB,
   TaylorEt,
   TaylorT4,
   TaylorN,
   NumApproximants
 } Approximant;


/** This structure is one of the members of the \c InspiralTemplate structure.

A user can specify the parameters of a binary using any of the
following combination of \e masses:
<ul>
<li> \c m1Andm2: component masses</li>
<li> \c totalMassAndEta: total mass and symmetric mass ratio</li>
<li> \c totalMassUAndEta: total mass and eta but uniform distribution in totalMass</li>
<li> \c totalMassAndMu: total mass and reduced mass</li>
<li> \c t01: unused; shouldn't be used.</li>
<li> \c t02: chirptimes \f$\tau_0\f$ and \f$\tau_2\f$</li>
<li> \c t03: chirptimes \f$\tau_0\f$ and \f$\tau_3\f$, and</li>
<li> \c t04: chirptimes \f$\tau_0\f$ and \f$\tau_4\f$</li>
<li> \c psi0Andpsi3: BCV parameters \f$\psi_0\f$ and \f$\psi_3\f$</li>
</ul>
The LALRandomInspiralSignal uses that structure as an input. Since the injected
waveform are not necessarely wanted to be random, we also provide the following
options
<ul>
<li> \c bhns: One of the mass is a Neutron star and the other a black
hole. (m1 \f$\in\f$ [minMass-3] and m2 \f$\in\f$ [3-maxMass]).</li>
<li> \c fixedMasses: The two masses are given by the input parameter structure.</li>
<li> \c fixedPsi: The two psi values are given by the input parameter structure.</li>
<li> \c fixedTau: The two tau values are given by the input parameter structure.</li>
</ul>
*/
typedef enum {
   m1Andm2,
   totalMassAndEta,
   totalMassUAndEta,
   totalMassAndMu,
   t01,
   t02,
   t03,
   t04,
   psi0Andpsi3,
   bhns,
   fixedMasses,
   fixedPsi,
   fixedTau,
   massesAndSpin,
   minmaxTotalMass,
   spinOnly
 } InputMasses;



/** Enumeration to specify which component will be used in the waveform
 * generation. Their combination also can be used by the bitwise or.
 **/
typedef enum {
	LAL_NOInter = 0,
	LAL_SOInter = 1, ///< Spin-orbit interaction
	LAL_SSInter = LAL_SOInter << 1, ///< Spin-spin interaction
	LAL_SSselfInter = LAL_SSInter << 1, ///< Spin-spin-self interaction
	LAL_QMInter = LAL_SSselfInter << 1, ///< quadrupole-monopole interaction
	LAL_AllInter = LAL_SOInter | LAL_SSInter | LAL_SSselfInter | LAL_QMInter	///< all interactions
} LALSpinInteraction;

/** The inspiral waveform parameter structure containing information about the waveform to be generated.

<ul>
<li> <tt> ieta:</tt> parameter that tells whether the symmetric mass ratio \f$\eta\f$
	  should be set to zero in the PN expansions of GW flux and binding energy.
	  If <tt>ieta=0</tt> \f$\eta\f$ will be set to zero, otherwise the appropriate
	  value of \f$\eta\f$ from the given parameters will be used.

  </li><li> <tt> level:</tt> Flag used in heirarical serached to indicate if this is a coarse or a fine template
  </li><li> <tt> *segmentIdVec:</tt> Vector of segment that have been filtered against this template needed for the LDAS implementation of the inspiral search.
  </li><li> <tt> number:</tt> Unique ID number for this template needed for the LDAS implementation of the inspiral search.
  </li><li> <tt> minMatch:</tt> The minimal match specified by the user when the bank that contains this template was created.
  </li><li> <tt> nStartPad:</tt> Number of leading elements in the signal generation to be set to zero (input). If template is requested, that value must be set to zero. In the injection routines related to inject package, that nStartPad is set to zero. However, for injection performed using the inspiral package, that value can be set to non zero.
  </li><li> <tt> nEndPad:</tt> Number of trailing bins to be set to zero, the
  resulting waveform will have at least this many bins zero at the end, probably
  more since we always deal with an integer power of 2 array (input).
  </li><li> <tt> mass1:</tt>  Mass of the primary in solar mass (input/output).
  </li><li> <tt> mass2:</tt>  Mass of the secondary in solar mass
  (\c mass1 need not be larger than \c mass2 (input/output).
  </li><li> <tt> spin1[3]:</tt> Spin vector of the primary (currently not in use)
  </li><li> <tt> spin2[3]:</tt> Spin vector of the secondary (currently not in use)
  </li><li> <tt> sourceTheta:</tt> Co-latitute in the direction to the source.
  </li><li> <tt> sourcePhi:</tt> Azimuth angle in the direction to the source.
  </li><li> <tt> orbitTheta0:</tt> Initial co-latitute of the orbit.
  </li><li> <tt> orbitPhi0:</tt> Initial azimuth angle of the orbit.
  </li><li> <tt> inclination:</tt> Inclination of the orbit  (currently not in use)
  </li><li> <tt> distance:</tt> Distance to the binary in seconds
  </li><li> <tt> psi0:</tt> BCV parameter \f$\psi_0.\f$
  </li><li> <tt> psi3:</tt> BCV parameter \f$\psi_3.\f$
  </li><li> <tt> alpha:</tt> BCV amplitude correction factor \f$\alpha f_\textrm{cut}^{2/3}\f$
  </li><li> <tt> eccentricity:</tt> initial eccentricity of the orbit  (currently not in use)
  </li><li> <tt> totalMass:</tt> total mass of the binary \f$m=m_1+m_2\f$ in solar mass (input/output).
  </li><li> <tt> eta:</tt> symmetric mass ratio \f$\eta=m_1m_2/m^2.\f$ (input/output).
  </li><li> <tt> chirpMass:</tt> chirp mass of the binary \f$=\eta^{3/5} m\f$ in solar mass (output).
  </li><li> <tt> t0:</tt> Newtonain chirp time in seconds (input/output).
  </li><li> <tt> t2:</tt> first post-Newtonian chirp time in seconds (input/output).
  </li><li> <tt> t3:</tt> 1.5 post-Newtonian chirp time in seconds (input/output).
  </li><li> <tt> t4:</tt> second post-Newtonian chirp time in seconds (output).
  </li><li> <tt> t5:</tt> 2.5 post-Newtonian chirp time in seconds (output).
  </li><li> <tt> t6:</tt> third post-Newtonian chirp time in seconds (output).
  </li><li> <tt> t7:</tt> 3.5 post-Newtonian chirp time in seconds (output).
  </li><li> <tt> tC:</tt> total chirp time seconds (output).
  </li><li> <tt> mu:</tt> reduced mass (in solar mass) (input/output)
  </li><li> <tt> fLower:</tt> lower frequency cutoff of the detector in Hz (input)
  </li><li> <tt> fCutoff:</tt> upper frequency cutoff in Hz to be used in generating the waveform.
  If the last stable orbit frequency is smaller than the upper cutoff it will be used
  in terminating the waveform instead of fCutoff (input).
  </li><li> <tt> tSampling:</tt> Sampling rate in Hz (input)
  </li><li> <tt> startPhase:</tt> starting phase of the waveform in radians (input)
  </li><li> <tt> startTime:</tt> starting time of the waveform (in sec); if different from zero, the
  waveform will start with an instantaneous frequency different from fLower and reach
  fLower at time (approximately) zero (input, not used in Stationary phase approximation)
  </li><li> <tt> signalAmplitude:</tt> dimensionless amplitude of the signal (input, currently unused.)
  </li><li> <tt> rInitial:</tt> initial radial separation of the two, in units of total mass
  bodies (used only in EOB waveforms) (output)
  </li><li> <tt> vInitial:</tt> initial velocity parameter, in units of the speed of light (output)
  </li><li> <tt> rFinal:</tt> final 'separation' between the bodies, in units of total mass (output)
  </li><li> <tt> vFinal:</tt> final velocity parameter, in units of the speed of light (output)
  </li><li> <tt> fFinal:</tt> final frequency reached, in units of Hz (output)
  </li><li> <tt> rLightRing:</tt> radial coordinate at the light ring, in units of total mass (output)
  </li><li> <tt> OmegaS:</tt> The 3PN (unknown) parameter; calculated to be equal to zero
  by Damour, Jaranowski and Schaffer (input).
  </li><li> <tt> Theta:</tt> The 3PN unknown flux parameter; likely to be around unity;
  most waveform generation routines take theta to be zero. Robustness of the EOB waveform
	  has been demonstrated for \f$-2 < \f$ \c Theta \f$< 2.\f$ (input)
  </li><li> <tt> massChoice:</tt> The pair of (mass) parameters given (see structure
		  defining this member for more details) (input).
  </li><li> <tt> order:</tt> Post-Newtonain order to be used in generating the wave (input).
  </li><li> <tt> approximant:</tt> Post-Newtonain approximant to be used in generating the wave (input).
  </li><li> <tt> tagInspiralTemplate *next:</tt> Linked list to the next coarse bank template
  (currently not filled by inspiral or bank codes)
  </li><li> <tt> tagInspiralTemplate *fine:</tt> Linked list to the next fine bank template
  (currently not filled by inspiral or bank codes)</li>
</ul>
*/
typedef struct
tagInspiralTemplate
{
/*  Parameters needed to generate Taylor/Pade waveforms */
  Approximant approximant;
  LALPNOrder order;
  LALPNOrder ampOrder;
  REAL8 mass1;
  REAL8 mass2;
  REAL8 fCutoff;
  REAL8 fLower;
  REAL8 tSampling;
  REAL8 distance;
  REAL8 signalAmplitude;
  REAL8 startPhase;
  REAL8 startTime;
  INT4  ieta;

/* Additional parameters for EOB waveforms */

  REAL8 Theta;
  REAL8 Zeta2;

/* Parameters for BCV1 template */

  REAL8 alpha;
  REAL8 psi0;
  REAL8 psi3;

/* Additional parameters for BCV2 template */

  REAL8 beta;
  REAL8 alpha1;
  REAL8 alpha2;
  REAL8 alpha3;
  REAL8 alpha4;
  REAL8 alpha5;
  REAL8 alpha6;

/* Parameters for spinning BH waveform */

  REAL8 inclination;
  REAL8 orbitTheta0;
  REAL8 orbitPhi0;
  REAL8 spin1[3];
  REAL8 spin2[3];
  REAL8 sourceTheta;
  REAL8 sourcePhi;
  REAL8 polarisationAngle;

/* Spin parameters for the PTF template */
  REAL8 chi; /* dimensionless spin of black hole (i.e. mass1) */
  REAL8 kappa; /* cosine of angle between spin of mass1 and orb ang mom */

/* Parameters which are currently might be used */

  REAL8 eccentricity;

/* Paramters which are computed using LALInspiralParameterCalc */

  REAL8 chirpMass;
  REAL8 eta;
  REAL8 totalMass;
  REAL8 fFinal;
  REAL8 t0;
  REAL8 t2;
  REAL8 t3;
  REAL8 t4;
  REAL8 t5;
  REAL8 t6;
  REAL8 t7;
  REAL8 tC;

/* Note that tc and fFinal are computed during waveform generation!!! */

  REAL4 minMatch;
  REAL8 mu;
  INT4  level;
  INT4  number;
  INT4  nStartPad;
  INT4  nEndPad;
  REAL8 OmegaS;
  REAL8 vFinal;
/*  REAL8 vInitial;
  REAL8 rFinal;
  REAL8 rInitial;
  REAL8 rLightRing;*/
  InputMasses massChoice;
  INT4Vector *segmentIdVec;
  LIGOTimeGPS end_time;
  EventIDColumn *event_id;
  CHAR ifo[LIGOMETA_IFO_MAX];

  /* Gamma[] is a vector that stores the upper triangular part of the metric in
   * the space of parameters. For time domain searches, Gamma[0,...,5] stores
   * the following information :
   *    Gamma[0] -> (tc,tc) metric component
   *    Gamma[1] -> (tc,t0) metric component
   *    Gamma[2] -> (tc,t3) metric component
   *    Gamma[3] -> (t0,t0) metric component
   *    Gamma[4] -> (t0,t3) metric component
   *    Gamma[5] -> (t3,t3) metric component
   * For spinBCV searches, (in 4 dimensions) Gamma[0,...,9] would be required.
   */
  REAL4  Gamma[10];
  REAL4  qmParameter[2];
  LALSpinInteraction	spinInteraction;

  InputAxis axisChoice;
  UINT4 fixedStep;
  UINT4 inspiralOnly;

  struct tagInspiralTemplate *next;
  struct tagInspiralTemplate *fine;
} InspiralTemplate;


/** This is a structure needed by the inner workings of the inspiral wave generation code.
 */
typedef struct
tagInspiralToffInput
{
   REAL8 tN;
   REAL8 t2;
   REAL8 t3;
   REAL8 t4;
   REAL8 t5;
   REAL8 t6;
   REAL8 t7;
   REAL8 tl6;
   REAL8 piM;
   REAL8 tc;
   REAL8 t;
 } InspiralToffInput;


/** This structure is needed to solve the differential equation
    giving the evolution of the orbital angular momentum and the
    spin angular momenta in the case of spinning black hole binaries.
<ul>
   <li>	\c v: parameter of 'integration': v=sqrt(M/r)
  </li><li> \c magS1: The constant spin magnitude of the primary.
  </li><li> \c magS2: The constant spin magnitude of the secondary.
  </li><li> <tt>NCap[3]:</tt> Source direction (unit vector) in detector coordinate system.
  </li><li> <tt>spin1[3]:</tt> Spin of the larger body.
  </li><li> \c M: Total mass of the binary (in seconds).
  </li><li> \c fourM1Plus: = \f$(4 m_1+3 m_2)/(2 m_1 M^3)\f$ (all masses expressed in seconds).
  </li><li> \c fourM2Plus: = \f$(4 m_2+3 m_1)/(2 m_2 M^3)\f$ (all masses expressed in seconds).
  </li><li> \c oneBy2Mcube: = \f$1/(2 M^3)\f$
  </li><li> \c threeBy2Mcube:  = \f$3/(2 M^3)\f$
  </li><li> \c thirtytwoBy5etc:=  \f$(32/5) \eta^2 M\f$</li>
</ul>
*/
typedef struct
tagInspiralACSTParams
{
   REAL8 v;
   REAL8 magS1;
   REAL8 magS2;
   REAL8 NCap[3];
   REAL8 spin1[3];
   REAL8 M;
   REAL8 fourM1Plus;
   REAL8 fourM2Plus;
   REAL8 oneBy2Mcube;
   REAL8 threeBy2Mcube;
   REAL8 thirtytwoBy5etc;
}  InspiralACSTParams;







/** This structure contains various post-Newtonian and P-approximant expansion
    coefficients; the meanings of the coefficients is indicated as comments
    before each list.

<ul>
<li> {Energy, flux, phase, time and frequency functions:} The following
	functions are generic function definitions that will be used in
	template generation. The function <tt>LALInspiralChooseModel,</tt>
	which is called by wave generation interface code, points these
	functions to the appropriate specific functions depending on the
	choices made by the user.</li>


<li> \c expnFunc: Structure to hold the pointers to the generic
	functions defined above.</li>


<li> <tt>TofVIn</tt> and <tt>TofVIntegrandIn:</tt> Structures needed to
	compute the time elapsed
	from/to the starting epoch of the waveform when the velocity
	parameter was \f$v_0,\f$ to/from the current epoch when velocity
	parameter is \f$v.\f$</li>


<li> <tt>InspiralPhaseIn</tt> and <tt>PhiofVIntegrandIn:</tt> Structures used
	to compute the phase of the signal from the `beginning', when the
	veolcity parameter is \f$v_0,\f$ to a time when the velocity parameter
	has evolved to a user input value \f$v.\f$</li>


<li> \c InspiralDerivativesIn: Structure used as an input to compute
	the derivatives needed in solving the phasing formula when the
	\c approximant is <tt>TaylorT1, TaylorP1</tt> or <tt>EOB.</tt></li>


<li> \c rk4GSLIntegrator: Structure containing steps and controls
for the GSL Runge-Kutta solver</li>

<li> \c rk4In: Structure used as an input to Runge-Kutta solver.</li>
</ul>
*/
typedef struct
tagexpnCoeffs {
   int ieta;
   /* coefficients in the Taylor expansion of new energy function*/
   REAL8 eTaN, eTa1, eTa2, eTa3;
   /* coefficients in the Pade expression of new energy function*/
   REAL8 ePaN, ePa1, ePa2, ePa3;
   /* coefficients in the Taylor expansion of usual energy function*/
   REAL8 ETaN, ETa1, ETa2, ETa3;
   /* coefficients in the Taylor expansion of the derivative of the
    usual energy function*/
   REAL8 dETaN, dETa1, dETa2, dETa3;

   /* Taylor expansion coefficients of energy flux*/
   REAL8 FTaN, FTa1, FTa2, FTa3, FTa4, FTa5, FTa6, FTa7, FTa8, FTl6, FTl8;
   /* Taylor expansion coefficients of factored flux*/
   REAL8 fTaN, fTa1, fTa2, fTa3, fTa4, fTa5, fTa6, fTa7, fTa8;
   /* Coefficients of the corresponding P-approximant*/
   REAL8 fPaN, fPa1, fPa2, fPa3, fPa4, fPa5, fPa6, fPa7, fPa8;

   /* Taylor expansion coefficents in t(v)*/
   REAL8 tvaN, tva2, tva3, tva4, tva5, tva6, tva7, tvl6;
   /* Taylor expansion coefficents in phi(v)*/
   REAL8 pvaN, pva2, pva3, pva4, pva5, pva6, pva7, pvl6;
   /* Taylor expansion coefficents in phi(t)*/
   REAL8 ptaN, pta2, pta3, pta4, pta5, pta6, pta7, ptl6;
   /* Taylor expansion coefficents in f(t)*/
   REAL8 ftaN, fta2, fta3, fta4, fta5, fta6, fta7, ftl6;
   /* Taylor expansion coefficents in psi(f) in the Fourier phase*/
   REAL8 pfaN, pfa2, pfa3, pfa4, pfa5, pfa6, pfa7, pfl5, pfl6;
   /* Taylor expansion for the spinning case */
   REAL8 ST[9], thetahat ;


   /* sampling rate and interval*/
   REAL8 samplingrate, samplinginterval;
   /* symmetric mass ratio, total mass, component masses*/
   REAL8 eta, totalmass, m1, m2;
   /* unknown 3PN parameters, euler constant*/
   REAL8 lambda, theta, EulerC, omegaS, zeta2;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase.*/
   REAL8 f0, fn, t0, tn, v0, vn, vf, vlso, flso, phiC;

   /* last stable orbit and pole defined by various Taylor and P-approximants*/
   REAL8 vlsoT0, vlsoT2, vlsoT4, vlsoT6;
   REAL8 vlsoP0, vlsoP2, vlsoP4, vlsoP6;
   REAL8 vpoleP4, vpoleP6;
}  expnCoeffs;

typedef REAL8 EnergyFunction(
   REAL8 v,
   expnCoeffs *ak);

typedef REAL8 FluxFunction(
   REAL8 v,
   expnCoeffs *ak);

typedef void (TestFunction)(
   REAL8Vector *vector1,
   REAL8Vector *vector2,
   void *params);

typedef void (InspiralPhasing2)(
   LALStatus *status,
   REAL8 *phase,
   REAL8 v,
   expnCoeffs *ak);

typedef void (InspiralPhasing3)(
   LALStatus *status,
   REAL8 *f,
   REAL8 td,
   expnCoeffs *ak);

typedef void (InspiralFrequency3)(
   LALStatus *status,
   REAL8 *f,
   REAL8 td,
   expnCoeffs *ak);

typedef void (InspiralTiming2) (
   LALStatus *status,
   REAL8 *toff,
   REAL8 f,
   void *params);

typedef struct
tagexpnFunc
{
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   InspiralTiming2 *timing2;
   InspiralPhasing2 *phasing2;
   InspiralPhasing3 *phasing3;
   InspiralFrequency3 *frequency3;
} expnFunc;


typedef struct
tagTofVIn
{
   REAL8 t;
   REAL8 v0;
   REAL8 t0;
   REAL8 vlso;
   REAL8 totalmass;
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
} TofVIn;

typedef struct
tagTofVIntegrandIn
{
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
} TofVIntegrandIn;


typedef struct
tagInspiralDerivativesIn
{
   REAL8 totalmass;
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
} InspiralDerivativesIn;







typedef struct
tagrk4In
{
   TestFunction *function;
   REAL8 x;
   REAL8Vector *y;
   REAL8Vector *dydx;
   REAL8Vector *yt;
   REAL8Vector *dym;
   REAL8Vector *dyt;
   REAL8 h;
   INT4 n;
} rk4In;

typedef struct
tagrk4GSLIntegrator
{
   const gsl_odeiv_step_type *type;
   gsl_odeiv_step *step;
   gsl_odeiv_control *control;
   gsl_odeiv_evolve *evolve;
   REAL8 *y;
   rk4In *input;
} rk4GSLIntegrator;

typedef struct
tagInspiralPhaseIn
{
   REAL8 v0;
   REAL8 phi0;
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
} InspiralPhaseIn;

typedef struct
tagPhiofVIntegrandIn
{
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
}  PhiofVIntegrandIn;

typedef struct
tagInspiralInit
{
  UINT4      nbins;
  expnCoeffs ak;
  expnFunc   func;

}  InspiralInit;

typedef enum
{
  INSPIRAL_TAPER_NONE,
  INSPIRAL_TAPER_START,
  INSPIRAL_TAPER_END,
  INSPIRAL_TAPER_STARTEND,
  INSPIRAL_TAPER_NUM_OPTS
}  InspiralApplyTaper;

/* Function prototypes */

/* --- HERE ARE SOME USEFUL PROTOTYPE FOR LENGTH, PARAMETER CALCULATION... --- */

void LALInspiralParameterCalc (
     LALStatus *status,
     InspiralTemplate *params);

void LALInspiralRestrictedAmplitude(
     LALStatus *status,
     InspiralTemplate  *params);





void LALInspiralWaveLength (
     LALStatus *status,
     UINT4 *n,
     InspiralTemplate params);






void LALInspiralChooseModel(
     LALStatus *status,
     expnFunc *func,
     expnCoeffs *ak,
     InspiralTemplate *params);





void LALInspiralSetup (
     LALStatus *status,
     expnCoeffs *ak,
     InspiralTemplate *params);




void
LALInspiralInit(
	LALStatus        *status,
	InspiralTemplate *params,
	InspiralInit     *paramsInit);





void LALInspiralWaveTaper(
     LALStatus    *status,
     REAL4Vector  *signalvec,
     UINT4       bookends
     );

int XLALInspiralWaveTaper(
                   REAL4Vector         *signalvec,
                   InspiralApplyTaper  bookends);

/* --- HERE ARE THE WAVEFORMS/MODELS PROTOTYPES --- */




void LALInspiralAmplitudeCorrectedWave(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralAmplitudeCorrectedWaveTemplates(
     LALStatus *status,
     REAL4Vector *filter1,
     REAL4Vector *filter2,
     InspiralTemplate *params);

void
LALInspiralAmplitudeCorrectedWaveForInjection(
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *params,
   PPNParamStruc  *ppnParams);





void LALInspiralWave(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralWaveTemplates(
     LALStatus *status,
     REAL4Vector *filter1,
     REAL4Vector *filter2,
     InspiralTemplate *params);

void
LALInspiralWaveForInjection(
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *params,
   PPNParamStruc  *ppnParams);





void LALInspiralWave1(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralWave1Templates(
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALInspiralWave1ForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams
     );

void LALInspiralEccentricity(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralEccentricityTemplates(
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALInspiralEccentricityForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams
     );






void LALInspiralWave2(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralWave2Templates (
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALInspiralWave2ForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams
     );





void LALInspiralWave3 (
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralWave3Templates (
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALInspiralWave3ForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams);




void LALInspiralStationaryPhaseApprox1 (
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);




void LALInspiralStationaryPhaseApprox2 (
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);





void LALEOBWaveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALEOBWaveformTemplates(
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALEOBWaveformForInjection(
     LALStatus *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams);





void LALBCVWaveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALTaylorEtWaveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALTaylorT4Waveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALTaylorT4WaveformForInjection(
     LALStatus        *status,
     CoherentGW       *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams);

void LALBCVSpinWaveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);





void LALTaylorNWaveform (
     LALStatus        *status,
     REAL4Vector      *signalvec,
     InspiralTemplate *params);





void LALInspiralSpinModulatedWave(
     LALStatus        *status,
     REAL4Vector      *signalvec,
     InspiralTemplate *in);


void LALInspiralSpinModulatedWaveForInjection(
     LALStatus *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams
     );





void
LALSTPNWaveformForInjection (
			    LALStatus        *status,
			    CoherentGW       *waveform,
			    InspiralTemplate *params,
			    PPNParamStruc  *ppnParams);

void
LALSTPNWaveformFramelessForInjection (
                             LALStatus        *status,
                             CoherentGW       *waveform,
                             InspiralTemplate *params,
                             PPNParamStruc    *ppnParams
                            );

void
LALSTPNWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec1,
                REAL4Vector      *signalvec2,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
		REAL4Vector      *shift,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                );
void
LALSTPNWaveformTemplates (
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   ) ;

void LALSTPNWaveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

/* Phen-Spin waveform functions*/





void LALSpinInspiralDerivatives(
                         REAL8Vector *values,
                         REAL8Vector *dvalues,
                         void *mparams );

void LALPSpinInspiralRD(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALPSpinInspiralRDTemplates (
    LALStatus        *status,
    REAL4Vector      *signalvec1,
    REAL4Vector      *signalvec2,
    InspiralTemplate *params
   );

void LALPSpinInspiralRDFreqDom (
				LALStatus        *status,
				REAL4Vector      *signalvec,
				InspiralTemplate *params);

void LALPSpinInspiralRDForInjection(
                           LALStatus        *status,
                           CoherentGW       *waveform,
                           InspiralTemplate *params,
                           PPNParamStruc  *ppnParams);

void LALPSpinInspiralRDEngine (
				LALStatus        *status,
				REAL8Vector      *signalvec1,
				REAL8Vector      *signalvec2,
				REAL8Vector      *h,
				REAL8Vector      *f,
				REAL8Vector      *phi,
				UINT4            *countback,
				InspiralTemplate *params,
				InspiralInit     *paramsInit
				);

/* Phenomenological waveform generation functions */

void LALBBHPhenWaveFreqDom ( LALStatus        *status,
			     REAL4Vector      *signalvec,
			     InspiralTemplate *params);

void LALBBHPhenWaveFreqDomTemplates( LALStatus        *status,
				     REAL4Vector      *signalvec1,
				     REAL4Vector      *signalvec2,
				     InspiralTemplate *params);

void LALBBHPhenWaveTimeDom ( LALStatus        *status,
			     REAL4Vector      *signalvec,
			     InspiralTemplate *insp_template);

void LALBBHPhenWaveTimeDomTemplates( LALStatus        *status,
				     REAL4Vector      *signalvec1,
				     REAL4Vector      *signalvec2,
				     InspiralTemplate *params);

void LALBBHPhenTimeDomEngine( LALStatus        *status,
			      REAL4Vector      *signalvec1,
			      REAL4Vector      *signalvec2,
			      REAL4Vector      *h,
			      REAL4Vector      *a,
			      REAL4Vector      *f,
			      REAL8Vector      *phiOut,
				  UINT4            *countback,
			      InspiralTemplate *params);

void LALBBHPhenWaveTimeDomForInjection (LALStatus        *status,
					CoherentGW       *waveform,
					InspiralTemplate *params,
					PPNParamStruc    *ppnParams);


/* --- OTHER PROTOTYPES --- */





void LALEtaTau02(
     LALStatus *status,
     REAL8 *x,
     REAL8 eta,
     void  *in);





void LALEtaTau04(
     LALStatus *status,
     REAL8 *x,
     REAL8 eta,
     void  *in);





void LALInspiralDerivatives (
     REAL8Vector *vec1,
     REAL8Vector *vec2,
     void *params);





void LALInspiralVelocity (
     LALStatus *status,
     REAL8 *v,
     TofVIn *params);





void LALInspiralPhasing1 (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     void *params);





void LALInspiralPhiofVIntegrand (
     LALStatus *status,
     REAL8 *,
     REAL8,
     void *);






void LALInspiralPhasing2_0PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);

void LALInspiralPhasing2_1PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);

void LALInspiralPhasing2_2PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);

void LALInspiralPhasing2_3PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);

void LALInspiralPhasing2_4PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);

void LALInspiralPhasing2_5PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);

void LALInspiralPhasing2_6PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);

void LALInspiralPhasing2_7PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     expnCoeffs *ak);






void LALInspiralPhasing3_0PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralPhasing3_1PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralPhasing3_2PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralPhasing3_3PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralPhasing3_4PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralPhasing3_5PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralPhasing3_6PN (
     LALStatus *,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralPhasing3_7PN (
     LALStatus *status,
     REAL8 *phase,
     REAL8 td,
     expnCoeffs *ak);





void LALInspiralTofV (
     LALStatus *,
     REAL8 *,
     REAL8,
     void *);






void LALInspiralTofVIntegrand (
     LALStatus *status,
     REAL8 *,
     REAL8,
     void *);





void LALInspiralTiming2_0PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);

void LALInspiralTiming2_1PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);

void LALInspiralTiming2_2PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);

void LALInspiralTiming2_3PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);

void LALInspiralTiming2_4PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);


void LALInspiralTiming2_5PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);

void LALInspiralTiming2_6PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);

void LALInspiralTiming2_7PN (
     LALStatus *,
     REAL8 *toff,
     REAL8 f,
     void *params);





void LALInspiralFrequency3_0PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralFrequency3_1PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralFrequency3_2PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralFrequency3_3PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralFrequency3_4PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralFrequency3_5PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralFrequency3_6PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

void LALInspiralFrequency3_7PN (
     LALStatus *status,
     REAL8 *frequency,
     REAL8 td,
     expnCoeffs *ak);

REAL4 LALInspiralHPlusPolarization(
     REAL8 phase,
     REAL8 v,
     InspiralTemplate *params );

REAL4 LALInspiralHCrossPolarization(
     REAL8 phase,
     REAL8 v,
     InspiralTemplate *params );





rk4GSLIntegrator * XLALRungeKutta4Init(
		   INT4 n,
                   rk4In *input);

void LALRungeKutta4(
     LALStatus *,
     REAL8Vector *,
     rk4GSLIntegrator *,
     void *);

void XLALRungeKutta4Free(
     rk4GSLIntegrator *integrator);

/* --- PARSING PROTOTYPE FOR INSPIRALTEMPLATE STRCUTURE --- */




void
LALInspiralITStructureParseParameters(
	LALStatus *status,
	UINT4 argc,
	CHAR **argv,
	InspiralTemplate *params);

void
LALInspiralITStructureSetDefault(
	LALStatus *status,
	InspiralTemplate *params);

void
LALInspiralITStructurePrint(
	LALStatus *status,
	InspiralTemplate  params);

void
LALInspiralITStructureHelp(void);

/* --- TEST PROTOTYPES --- */

INT4 XLALInspiralRingdownWave (
	REAL4Vector			*rdwave1,
	REAL4Vector			*rdwave2,
	InspiralTemplate		*params,
	REAL4VectorSequence		*inspwave1,
	REAL4VectorSequence		*inspwave2,
	COMPLEX8Vector			*modefreqs,
	UINT4				nmodes
	);

INT4 XLALGenerateWaveDerivatives (
	REAL4Vector		*dwave,
	REAL4Vector		*ddwave,
	REAL4Vector		*wave,
	InspiralTemplate	*params
	);

INT4 XLALGenerateQNMFreq(
	COMPLEX8Vector		*modefreqs,
	InspiralTemplate	*params,
	UINT4			l,
	UINT4			m,
	UINT4			nmodes
	);

INT4 XLALFinalMassSpin(
	REAL8			*finalMass,
	REAL8			*finalSpin,
	InspiralTemplate	*params
	);

INT4 XLALInspiralAttachRingdownWave (
        REAL4Vector 	 *Omega,
        REAL4Vector 	 *signalvec1,
        REAL4Vector  	 *signalvec2,
        InspiralTemplate *params);


INT4 XLALPSpinInspiralRingdownWave (
       REAL8Vector             *rdwave,
       InspiralTemplate        *params,
       REAL8Vector             *inspwave,
       COMPLEX8Vector          *modefreqs,
       UINT4                   nmodes
       );

INT4 XLALGenerateWaveDerivative (
	REAL8Vector		*dwave,
	REAL8Vector	        *wave,
	REAL8                    dt
       );

INT4 XLALPSpinGenerateQNMFreq (
       COMPLEX8Vector          *modefreqs,
       InspiralTemplate        *params,
       UINT4                   l,
       INT4                    m,
       UINT4                   nmodes,
       REAL8                   finalMass,
       REAL8                   finalSpin
       );

INT4 XLALPSpinFinalMassSpin(
       REAL8                   *finalMass,
       REAL8                   *finalSpin,
       InspiralTemplate        *params,
       REAL8                   energy,
       REAL8                   *LNhvec
       );

INT4 XLALPSpinInspiralAttachRingdownWave (
       REAL8Vector       *signalvec,
       InspiralTemplate  *params,
       UINT4             *attpos,
       UINT4             nmodes,
       UINT4             l,
       INT4              m,
       REAL8             finalMass,
       REAL8             finalSpin
);

int XLALInspiralGetApproximantString( CHAR        *output,
                                      UINT4       length,
                                      Approximant approx,
                                      LALPNOrder  order
                                    );

int XLALBandPassInspiralTemplate(
        REAL4Sequence  *sequence,
        REAL4          fLow,
        REAL4          fHigh,
        REAL4          fSampling
        );

int XLALInspiralGenerateIIRSet(
	REAL8Vector         *amp,
	REAL8Vector         *phase,
	double                  epsilon,
	double                  alpha,
	double                  beta,
	COMPLEX16Vector     **a1,
	COMPLEX16Vector     **b0,
	INT4Vector          **delay
	);

int XLALInspiralIIRSetResponse(
	COMPLEX16Vector     *a1,
	COMPLEX16Vector     *b0,
	INT4Vector          *delay,
	COMPLEX16Vector     *response
	);

int XLALInspiralGenerateIIRSetFourierTransform(
	int                 j,
	int                 jmax,
	COMPLEX16           a1,
	COMPLEX16           b0,
	int                 delay,
	COMPLEX16           *hfcos,
	COMPLEX16           *hfsin
	);

int XLALInspiralCalculateIIRSetInnerProduct(
	COMPLEX16Vector    *a1,
	COMPLEX16Vector    *b0,
	INT4Vector         *delay,
	REAL8Vector        *psd,
	double             *ip
	);

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRAL_H */
