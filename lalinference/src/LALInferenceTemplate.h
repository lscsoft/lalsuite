/*
 *
 *  LALInference:             Bayesian Followup        
 *  LALInferenceTemplate.h:   Template generation functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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
 *  \file LALInferenceTemplate.h
 *  \brief Main header file for LALInference signal template generating functions.
 *  \ingroup LALInference
 *
 *  All template functions have a parameter
 *    \param[in,out] IFOdata used for both specifiying the signal parameter values 
 *                           and returning the waveform template.
 *
 *  Signal parameter values are passed to the template generating functions via 
 *  the \c IFOdata->modelParams parameter.
 *
 *  Signal templates are either in time domain or frequency domain, and subsequent functions 
 *  (e.g. for likelihood computation) may then act accordingly by 
 *  checking the \c IFOdata->modelDomain parameter.
 * 
 *  The actual template waveforms then are stored in 
 *  the \c IFOdata->freqModelhPlus / \c IFOdata->freqModelhCross
 *  or \c IFOdata->timeModelhPlus / \c IFOdata->timeModelhCross slots.
 */


#ifndef LALInferenceTemplate_h
#define LALInferenceTemplate_h

#include <lal/LALInference.h>


/** De-bugging function writing a (frequency-domain) signal template to a CSV file.
 *  File contains real & imaginary parts of plus & cross components.
 *  Template amplitude is scaled to 1 Mpc luminosity distance.                        
 */
void LALInferenceDumptemplateFreqDomain(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                            LALInferenceTemplateFunction *template, const char *filename);


/** De-bugging function writing a (time-domain) signal template to a CSV file.
 *  File contains time series of plus & cross components.
 *  Template amplitude is scaled to 1 Mpc luminosity distance.                        
 */
void LALInferenceDumptemplateTimeDomain(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                            LALInferenceTemplateFunction *template, const char *filename);


/** Template function to generate LAL's "parametrized post-Newtonian" (PPN) inspiral waveform.
 *  Internally uses LAL's LALGeneratePPNInspiral() function. 
 *  Signal amplitude is scaled to 1 Mpc luminosity distance.
 */
void LALInferenceLALTemplateGeneratePPN(LALInferenceIFOData *IFOdata);


/** Generates a 2.0PN / 2.5PN stationary phase approximation inspiral template in frequency domain.
 * 
 *  Computations are done 
 *  following  Tanaka/Tagoshi (2000), Phys.Rev.D 62(8):082001
 *  or Christensen/Meyer (2001), Phys.Rev.D 64(2):022001.
 *  By supplying the optional \c IFOdata->modelParams "PNOrder"
 *  parameter, one may request a 2.0PN (instead of 2.5PN)
 *  template.
 *  Signal's amplitude corresponds to a luminosity distance
 *  of 1 Mpc; re-scaling will need to be taken care of e.g.
 *  in the calling likelihood function.
 */
void LALInferenceTemplateStatPhase(LALInferenceIFOData *IFOdata);


/** Returns a frequency-domain 'null' template
 *  (all zeroes, implying no signal present).
 */
void LALInferenceTemplateNullFreqdomain(LALInferenceIFOData *IFOdata);


/** Returns a time-domain 'null' template
 *  (all zeroes, implying no signal present).
 */
void LALInferenceTemplateNullTimedomain(LALInferenceIFOData *IFOdata);


/** Wrapper function to call LAL functions for waveform generation.
 *  Will always return frequency-domain templates (numerically FT'ed
 *  in case the LAL function returns time-domain).
 *
 * Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *   \c "chirpmass"        (REAL8, units of solar masses)
 * 
 *   \c "massratio"        (symmetric mass ratio:  0 < eta <= 0.25, REAL8)
 * 
 *   \c "phase"            (here: 'startPhase', not coalescence phase; REAL8, radians)
 * 
 *   \c "time"             (coalescence time, or equivalent/analog/similar; REAL8, GPS sec.)
 * 
 *   \c "inclination"      (inclination angle, REAL8, radians)
 * 
 *   \c "LAL_APPROXIMANT"  (INT4 value corresponding to `enum approximant' definition in `LALInspiral.h'.
 *                         Templates that (seem to) work by now are:
 *                         TaylorF2, TaylorT1, TaylorT2, TaylorT3, BCV, IMRPhenomA, EOB, EOBNR)
 * 
 *   \c "LAL_PNORDER"      (INT4 value corresponding to `enum LALPNOrder' definition in `LALInspiral.h'.)
 */
void LALInferenceTemplateLAL(LALInferenceIFOData *IFOdata);


/** Generate 3.5PN phase / 2.5PN amplitude time-domain inspiral templates.
 *
 * following
 * 
 *   Blanchet et al. (2001),   gr-qc/0104084
 * 
 *   Blanchet at al. (2002),   PRD 65(6):061501,   gr-qc/0105099
 * 
 *   Blanchet at al. (2005),   PRD 71(12):129902
 * 
 *   Arun et al. (2004),       CQG 21(15):3771
 * 
 *   Arun et al. (2004),       CQG 22(14):3115
 * 
 *   Blanchet et al. (2004),   PRL 93(9):091101
 *
 *  This is essentially the implementation that was also used in the paper by
 *  Roever/Meyer/Guidi/Vicere/Christensen (2007), CQG 24(19):S607.
 *
 *  Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *    \c "chirpmass"        (REAL8, units of solar masses)
 * 
 *    \c "massratio"        (symmetric mass ratio:  0 < eta <= 0.25, REAL8)
 * 
 *    \c "phase"            (coalescence phase; REAL8, radians)
 * 
 *    \c "time"             (coalescence time; REAL8, GPS seconds)
 * 
 *    \c "inclination"      (inclination angle, REAL8, radians)
 */
void LALInferenceTemplate3525TD(LALInferenceIFOData *IFOdata);


/** Sine-Gaussian (burst) template.
 * 
 *  The (plus-) waveform is given by:
 *    \f[ a \exp(-((t - \mu) / \sigma)^2) \times \sin(2 \pi f t - phi) \f]
 *
 *  Note that by setting f=0, phi=pi/2 you get a "plain" Gaussian template.
 *
 *  Signal is (by now?) linearly polarised, i.e., the cross-component remains zero.    
 *
 *  Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *    \c "time"       (the \f[ \mu \f[ parameter of the Gaussian part; REAL8, GPS seconds)
 * 
 *    \c "sigma"      (width, the \f[ \sigma \f[ parameter of the Gaussian part; REAL8, seconds)
 * 
 *    \c "frequency"  (frequency \f[ f \f[ of the sine part; REAL8, Hertz)
 * 
 *    \c "phase"      (phase \f[ \phi \f[ (at time \f[ \mu \f[); REAL8, radians)
 * 
 *    \c "amplitude"  (amplitude \f$ a \f$, REAL8)
 *
 */
void LALInferenceTemplateSineGaussian(LALInferenceIFOData *IFOdata);


/** Damped Sinusoid template.
 *
 *  The (plus-) waveform is an exponentially decaying sine wave:
 *    \f[ a \exp((t-time) / \tau) \times  sin(2 \pi f (t-time)) \f]
 *  where "time" is the time parameter denoting the instant at which the signal starts.
 *
 *  Signal is (by now?) linearly polarised, i.e., the cross-component remains zero.    
 *
 *  Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *    \c "time"       (the instant at which the signal starts; REAL8, GPS seconds)
 * 
 *    \c "tau"        (width parameter \f$ \tau \f$; REAL8, seconds)
 * 
 *    \c "frequency"  (frequency \f$ f \f$ of the sine part; REAL8, Hertz)
 * 
 *    \c "amplitude"  (amplitude \f$ a \f$, REAL8)
 */
void LALInferenceTemplateDampedSinusoid(LALInferenceIFOData *IFOdata);


/** Sinc function (burst) template.
 *
 *  The (plus-) waveform is a sinc function of given frequency:
 *    \f[ a sinc(2 \pi f (t-time)) = a \sin(2 \pi f (t-time)) / (2 \pi f (t-time)) \f]
 *  where "time" is the time parameter denoting the signal's central peak location.
 *
 *  Signal is (by now?) linearly polarised, i.e., the cross-component remains zero.    
 *
 *  Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *    \c "time"       (the instant at which the signal peaks; REAL8, GPS seconds)
 * 
 *    \c "frequency"  (frequency \f$ f \f$ of the sine part; REAL8, Hertz)
 * 
 *    \c "amplitude"  (amplitude \f$ a \f$, REAL8)
 */
void LALInferenceTemplateSinc(LALInferenceIFOData *IFOdata);


/** LALSTPN template
 *
 *  Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *   \c "m1"           (mass of object 1; REAL8, solar mass)
 * 
 *   \c "m2"           (mass of object 1; REAL8, solar mass)
 * 
 *   \c "inclination"  (inclination angle; REAL8, radians)
 * 
 *   \c "coa_phase"    (phase angle; REAL8, radians)
 * 
 *   \c "spin1x"       (x component of the spin of object 1; REAL8)
 * 
 *   \c "spin1y"       (y component of the spin of object 1; REAL8)
 * 
 *   \c "spin1z"       (z component of the spin of object 1; REAL8)
 * 
 *   \c "spin2x"       (x component of the spin of object 2; REAL8)
 * 
 *   \c "spin2y"       (y component of the spin of object 2; REAL8)
 * 
 *   \c "spin2z"       (z component of the spin of object 2; REAL8)
 * 
 *   \c "shift0"       (shift offset; REAL8, radians)
 * 
 *   \c "time"         (coalescence time, or equivalent/analog/similar; REAL8, GPS seconds)
 * 
 *   \c "PNorder"      (Phase PN order; REAL8)
 */
void LALInferenceTemplateLALSTPN(LALInferenceIFOData *IFOdata);


/** Trivial h(t) = A*sin(Omega*t) template
 *
 *  Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *    \c "A"       (dimensionless amplitude, REAL8)
 * 
 *    \c "Omega"   (frequency; REAL8, radians/sec)
 */
void LALInferenceTemplateASinOmegaT(LALInferenceIFOData *IFOdata);


/** LALGenerateInspiral wrapper.
 *
 *  Required ( \c IFOdata->modelParams ) parameters are:
 * 
 *    \c "m1"           (mass of object 1; REAL8, solar mass)
 * 
 *    \c "m2"           (mass of object 1; REAL8, solar mass)
 * 
 *    \c "inclination"  (inclination angle; REAL8, radians)
 * 
 *    \c "coa_phase"    (phase angle; REAL8, radians)
 * 
 *    \c "spin1x"       (x component of the spin of object 1; REAL8) (if SpinTaylor approx)
 * 
 *    \c "spin1y"       (y component of the spin of object 1; REAL8) (if SpinTaylor approx)
 * 
 *    \c "spin1z"       (z component of the spin of object 1; REAL8) (if SpinTaylor approx)
 * 
 *    \c "spin2x"       (x component of the spin of object 2; REAL8) (if SpinTaylor approx)
 * 
 *    \c "spin2y"       (y component of the spin of object 2; REAL8) (if SpinTaylor approx)
 * 
 *    \c "spin2z"       (z component of the spin of object 2; REAL8) (if SpinTaylor approx)
 * 
 *    \c "shift0"       (shift offset; REAL8, radians)
 * 
 *    \c "time"         (coalescence time, or equivalent/analog/similar; REAL8, GPS seconds)
 * 
 *    \c "PNorder"      (Phase PN order; REAL8)
 */
void LALInferenceTemplateLALGenerateInspiral(LALInferenceIFOData *IFOdata);


/** Template function for PhenSpinTaylorRingDown waveforms. 
 *
 *  (untested!)
 */
void LALInferenceTemplatePSTRD(LALInferenceIFOData *IFOdata);

#endif
