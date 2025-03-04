/*
 *
 *  LALInferenceLikelihood.c:   Likelihood functions for LALInference codes
 *  LALInferenceLikelihood.h:   header file
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
#ifndef LALInferenceLikelihood_h
#define LALInferenceLikelihood_h

#include <lal/LALInference.h>


#ifdef SWIG // SWIG interface directives
SWIGLAL(
	FUNCTION_POINTER(
			 LALInferenceUndecomposedFreqDomainLogLikelihood,
			 LALInferenceNoiseOnlyLogLikelihood,
			 LALInferenceZeroLogLikelihood,
			 LALInferenceFreqDomainLogLikelihood,
			 LALInferenceNullLogLikelihood,
			 LALInferenceFreqDomainStudentTLogLikelihood,
			 LALInferenceCorrelatedAnalyticLogLikelihood,
			 LALInferenceBimodalCorrelatedAnalyticLogLikelihood,
			 LALInferenceRosenbrockLogLikelihood,
			 LALInferenceMarginalisedPhaseLogLikelihood,
			 LALInferenceMarginalisedTimeLogLikelihood
			 )
);
#endif


extern const char *LALInferenceAnalyticNamesCBC[15];

extern const REAL8 LALInferenceAnalyticMeansCBC[15];

/* Scaling used for the CBC analytic likelihood parameters */
extern const REAL8 scaling[15];

/* Covariance matrix for use in analytic likelihoods */
extern const REAL8 CM[15][15];


/**
 * \defgroup LALInferenceLikelihood_h Header LALInferenceLikelihood.h
 * \ingroup lalinference_general
 *
 * \brief Header file for likelihood functions used by LALInference codes
 *
 * LALInferenceLikelihood contains all the necessary routines to compute the likelihood
 * from a template (computed with LALInferenceTemplate) and the data (initialised with LALInferenceReadData).
 *
 * Likelihood functions follow the basic naming convention: LALInference<type_of>LogLikelihood()
 *
 * Takes as input:
 * - a pointer to a LALInferenceVariable structure containing the parameters to compute the likelihood for,
 * - a pointer to a LALInferenceIFOData structure containing the linked list of interferometer data,
 * - a pointer to the LALInferenceModel model structure, containing the template function and model.
 *
 * Outputs as a REAL8 the natural logarithm value of the likelihood, as defined by:
 *
 * \f[
 * Likelihood(\vec{x}|\vec{\lambda},M)=\exp(-\tfrac{1}{2}<\vec{x}-\vec{h_M}(\vec{\lambda})|\vec{x}-\vec{h_M}(\vec{\lambda})>)
 * \f]
 *
 * where: \f$<x|y>=4Re\left ( \int \frac{\tilde{x}\,\tilde{y}^*}{S_f}\, df \right )\f$
 *
 * Note that the likelihood is reported unnormalised.
 *
 */
/** @{ */

/***********************************************************//**
 * (log-) likelihood function.
 * Returns the non-normalised logarithmic likelihood.
 *
 * Required (`currentParams') parameters are:
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)
 *   - "distance"        (REAL8, Mpc, >0)
 *   - "time"            (REAL8, GPS sec.)
 ***************************************************************/
REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

/**
 * For testing purposes (for instance sampling the prior),
 * likelihood that returns 0.0 = log(1) every
 * time.  Activated with the --zeroLogLike command flag.
 */
REAL8 LALInferenceZeroLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);


/**
 * Computes the <x|y> overlap in the Fourier domain.
 */
REAL8 LALInferenceComputeFrequencyDomainOverlap(LALInferenceIFOData * data,
        COMPLEX16Vector * freqData1, COMPLEX16Vector * freqData2);
/**
 * Computes the complex <x|y> overlap
 */
COMPLEX16 LALInferenceComputeFrequencyDomainComplexOverlap(LALInferenceIFOData * dataPtr,
                                                           COMPLEX16Vector * freqData1,
                                                           COMPLEX16Vector * freqData2);

/**
 * Identical to LALInferenceFreqDomainNullLogLikelihood, but returns the likelihood of a null template.
 * Used for normalising.
 */
REAL8 LALInferenceNullLogLikelihood(LALInferenceIFOData *data);

/***********************************************************//**
 * Student-t (log-) likelihood function
 * as described in Roever/Meyer/Christensen (2011):
 *   "Modelling coloured residual noise
 *   in gravitational-wave signal processing."
 *   Classical and Quantum Gravity, 28(1):015010.
 *   http://dx.doi.org/10.1088/0264-9381/28/1/015010
 *   http://arxiv.org/abs/0804.3853
 * Returns the non-normalised logarithmic likelihood.
 *
 * Required (`currentParams') parameters are:
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)
 *   - "distance"        (REAL8, Mpc, > 0)
 *   - "time"            (REAL8, GPS sec.)
 *
 * This function is essentially the same as the
 * "UndecomposedFreqDomainLogLikelihood()" function.
 * The additional parameter to be supplied is the (REAL8)
 * degrees-of-freedom parameter (nu) for each Ifo.
 * The additional "df" argument gives the corresponding
 * d.f. parameter for each element of the "*data" list.
 * The names of "df" must match the "->name" slot of
 * the elements of "data".
 *
 * (TODO: allow for d.f. parameter to vary with frequency,
 *        i.e., to be a set of vectors corresponding to
 *        frequencies)
 ***************************************************************/
REAL8 LALInferenceFreqDomainStudentTLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

/**
 * An analytic likeilhood that is a correlated Gaussian in 15
 * dimensions.
 */
REAL8 LALInferenceCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

/**
 * An analytic likeilhood that is two correlated Gaussians in 15
 * dimensions.
 */
REAL8 LALInferenceBimodalCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

/**
 * 15-D Rosenbrock log(L) function (see Eq (3) of
 * http://en.wikipedia.org/wiki/Rosenbrock_function .
 */
REAL8 LALInferenceRosenbrockLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

REAL8 LALInferenceMarginalisedPhaseLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

REAL8 LALInferenceMarginalisedTimePhaseLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

/** Compute delta-log-likelihood for given distance min, max and OptimalSNR and d_inner_h when evaluated at 1Mpc
  * cosmology: 0 = Euclidean distance prior , 1 = uniform in comoving volume
    margphi: 0 = use gaussian likelihood, 1 = phase-marginalised bessel likelihood */
double LALInferenceMarginalDistanceLogLikelihood(double dist_min, double dist_max, double OptimalSNR, double d_inner_h, int cosmology, int margphi);


/**
 * Returns the log-likelihood marginalised over the time dimension
 * from the prior min to the prior max.  See
 * https://dcc.ligo.org/LIGO-T1400460 for the details.
 */
REAL8 LALInferenceMarginalisedTimeLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);

/**
 * Initialisation function which reads runState->commaneLine and sets up the
 * likelihood function accordingly. Can choose between Gaussian, Student-t, marginalised
 * phase likelihoods
 */
void LALInferenceInitLikelihood(LALInferenceRunState *runState);

/** Get the intrinsic parameters from currentParams */
LALInferenceVariables LALInferenceGetInstrinsicParams(LALInferenceVariables *currentParams);

/** fast SineGaussian likelihood for LIB */
REAL8 LALInferenceFastSineGaussianLogLikelihood(LALInferenceVariables *currentParams,
                                                        LALInferenceIFOData *data,
                                                        LALInferenceModel *model);

/** Calculate the SNR across the network */
void LALInferenceNetworkSNR(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceModel *model);
/** @} */

#endif
