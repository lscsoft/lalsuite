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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */
#ifndef LALInferenceTemplate_h
#define LALInferenceTemplate_h

#include <lal/LALInference.h>

#ifdef SWIG // SWIG interface directives
SWIGLAL(
    FUNCTION_POINTER(
        LALInferenceTemplate3525TD,
        LALInferenceTemplateASinOmegaT,
        LALInferenceTemplateDampedSinusoid,
        LALInferenceTemplateLAL,
        LALInferenceTemplateLALGenerateInspiral,
        LALInferenceTemplateNullFreqdomain,
        LALInferenceTemplateNullTimedomain,
        LALInferenceTemplatePSTRD,
        LALInferenceTemplateSinc,
        LALInferenceTemplateSineGaussian,
        LALInferenceTemplateStatPhase,
        LALInferenceTemplateXLALSimInspiralChooseWaveform
    )
);
#endif

/**
 * \defgroup LALInferenceTemplate_h Header LALInferenceTemplate.h
 * \ingroup lalinference_general
 *
 * \brief Main header file for LALInference signal template generating functions.
 *
 * All template functions have a parameter
 * \param[in,out] IFOdata used for both specifiying the signal parameter values
 * and returning the waveform template.
 *
 * Signal parameter values are passed to the template generating functions via
 * the \c IFOdata->modelParams parameter.
 *
 * Signal templates are either in time domain or frequency domain, and subsequent functions
 * (e.g. for likelihood computation) may then act accordingly by
 * checking the \c IFOdata->modelDomain parameter.
 *
 * The actual template waveforms then are stored in
 * the \c IFOdata->freqModelhPlus / \c IFOdata->freqModelhCross
 * or \c IFOdata->timeModelhPlus / \c IFOdata->timeModelhCross slots.
 * The template's amplitude is (for physical models) scaled to
 * 1 Mpc luminosity distance.
 *
 */
/** @{ */

/**
 * De-bugging function writing a (frequency-domain) signal template to a CSV file.
 * File contains real & imaginary parts of plus & cross components.
 * Template amplitude is (usually) scaled to 1 Mpc luminosity distance.
 */
void LALInferenceDumptemplateFreqDomain(LALInferenceVariables *currentParams, LALInferenceModel *model,
                                        const char *filename);


/**
 * De-bugging function writing a (time-domain) signal template to a CSV file.
 * File contains time series of plus & cross components.
 * Template amplitude is (usually) scaled to 1 Mpc luminosity distance.
 */
void LALInferenceDumptemplateTimeDomain(LALInferenceVariables *currentParams, LALInferenceModel *model,
                                        const char *filename);


/**
 * Returns a frequency-domain 'null' template
 * (all zeroes, implying no signal present).
 */
void LALInferenceTemplateNullFreqdomain(LALInferenceModel *model);


/**
 * Returns a time-domain 'null' template
 * (all zeroes, implying no signal present).
 */
void LALInferenceTemplateNullTimedomain(LALInferenceModel *model);


/**
 * Sine-Gaussian (burst) template.
 *
 * The (plus-) waveform is given by:
 * \f[ s(t) = a \times \exp(-((t - \mu) / \sigma)^2) \times \sin(2 \pi f (t-\mu) - \phi) \f]
 *
 * Note that by setting f=0, phi=pi/2 you get a "plain" Gaussian template.
 *
 * Signal is (by now?) linearly polarised, i.e., the cross-component remains zero.
 *
 * Required ( \c IFOdata->modelParams ) parameters are:
 * - \c "time"       (REAL8, the \f$ \mu \f$ parameter of the Gaussian part, in GPS seconds)
 * - \c "sigma"      (REAL8, width, the \f$ \sigma \f$ parameter of the Gaussian part, seconds)
 * - \c "frequency"  (REAL8, frequency \f$ f \f$ of the sine part, Hertz)
 * - \c "phase"      (REAL8, phase \f$ \phi \f$ (at time \f$ \mu \f$), radians)
 * - \c "amplitude"  (REAL8, amplitude \f$ a \f$)
 */
void LALInferenceTemplateSineGaussian(LALInferenceModel *model);

void LALInferenceROQWrapperForXLALSimInspiralChooseFDWaveformSequence(LALInferenceModel *model);
/**
 * Damped Sinusoid template.
 *
 * The (plus-) waveform is an exponentially decaying sine wave:
 * \f[ s(t) = a \times \exp((t-time) / \tau) \times  sin(2 \pi f (t-time)) \f]
 * where "time" is the time parameter denoting the instant at which the signal starts.
 *
 * Signal is (by now?) linearly polarised, i.e., the cross-component remains zero.
 *
 * Required ( \c IFOdata->modelParams ) parameters are:
 * - \c "time"       (REAL8, the instant at which the signal starts, in GPS seconds)
 * - \c "tau"        (REAL8, width parameter \f$ \tau \f$, seconds)
 * - \c "frequency"  (REAL8, frequency \f$ f \f$ of the sine part, Hertz)
 * - \c "amplitude"  (REAL8, amplitude \f$ a \f$)
 */
void LALInferenceTemplateDampedSinusoid(LALInferenceModel *model);


/**
 * Sinc function (burst) template.
 *
 * The (plus-) waveform is a sinc function of given frequency:
 * \f[ s(t) = a \times sinc(2 \pi f (t-time)) = a \times \sin(2 \pi f (t-time)) / (2 \pi f (t-time)) \f]
 * where "time" is the time parameter denoting the signal's central peak location.
 *
 * Signal is (by now?) linearly polarised, i.e., the cross-component remains zero.
 *
 * Required ( \c IFOdata->modelParams ) parameters are:
 * - \c "time"       (REAL8, the instant at which the signal peaks, in GPS seconds)
 * - \c "frequency"  (REAL8, frequency \f$ f \f$ of the sine part, Hertz)
 * - \c "amplitude"  (REAL8, amplitude \f$ a \f$)
 */
void LALInferenceTemplateSinc(LALInferenceModel *model);


/**
 * Trivial h(t) = A*sin(Omega*t) template.
 *
 * Required ( \c IFOdata->modelParams ) parameters are:
 * - \c "A"       (REAL8, dimensionless amplitude)
 * - \c "Omega"   (REAL8, frequency, radians/sec)
 */
void LALInferenceTemplateASinOmegaT(LALInferenceModel *model);


/**
 * "XLALSimInspiralChooseWaveform{TD,FD}" wrapper.
 *
 * Required ( \c IFOdata->modelParams ) parameters are:
 * - \c "m1"           (REAL8, mass of object 1, solar masses)
 * - \c "m2"           (REAL8, mass of object 1, solar masses)
 * - \c "inclination"  (REAL8, inclination angle, radians)
 * - \c "coa_phase"    (REAL8, phase angle, radians)
 * - \c "spin1x"       (REAL8, x component of the spin of object 1)
 * - \c "spin1y"       (REAL8, y component of the spin of object 1)
 * - \c "spin1z"       (REAL8, z component of the spin of object 1)
 * - \c "spin2x"       (REAL8, x component of the spin of object 2)
 * - \c "spin2y"       (REAL8, y component of the spin of object 2)
 * - \c "spin2z"       (REAL8, z component of the spin of object 2)
 * - \c "shift0"       (REAL8, shift offset, radians)
 * - \c "time"         (REAL8, coalescence time, or equivalent/analog/similar, GPS seconds)
 * - \c "PNorder"      (REAL8, Phase PN order)
 * - \c "
 *
 * THIS IMPLEMENTATION IS NOT THREAD SAFE !!! (previous inclination value is stored as a static)
 *
 */
void LALInferenceTemplateXLALSimInspiralChooseWaveform(LALInferenceModel *model);

void LALInferenceTemplateXLALSimBurstChooseWaveform(LALInferenceModel *model);

void LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated(LALInferenceModel *model);

void LALInferenceTemplateXLALSimBurstSineGaussianF(LALInferenceModel *model);


/** @} */

#endif
