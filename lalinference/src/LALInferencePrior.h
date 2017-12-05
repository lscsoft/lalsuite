/*
 *
 *  LALInference:          LAL Inference library
 *  LALInferencePrior.h:   Collection of common priors
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
#ifndef LALInferencePrior_h
#define LALInferencePrior_h

#include <lal/LALInference.h>

#ifdef SWIG // SWIG interface directives
SWIGLAL(
    FUNCTION_POINTER(
        LALInferenceInspiralSkyLocPrior,
        LALInferenceInspiralPriorNormalised,
        LALInferenceAnalyticNullPrior
    )
);
#endif

/**
 * \defgroup LALInferencePrior_h Header LALInferencePrior.h
 * \ingroup lalinference_general
 * \brief Collection of commonly used Prior functions and utilities
 */
/*@{*/

#ifdef SWIG /* SWIG interface directives */
SWIGLAL(INPUT_SCALARS(REAL8*, min, max));
#endif

/**
 * Initialize the prior based on command line arguments.
*/
void LALInferenceInitCBCPrior(LALInferenceRunState *runState);

/**
 * Initialize the LIB prior based on command line arguments.
*/
void LALInferenceInitLIBPrior(LALInferenceRunState *runState);


/**
 * Return the log Prior for the glitch amplitude
*/
REAL8 logGlitchAmplitudeDensity(REAL8 A, REAL8 Q, REAL8 f);

/**
 * Return the logarithmic prior density of the variables specified, for the non-spinning/spinning inspiral signal case.
 */
REAL8 LALInferenceInspiralPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model);

/**
 * Convert the hypercube parameter to physical parameters, for the non-spinning/spinning inspiral signal case.
 */
UINT4 LALInferenceInspiralCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model, double *Cube, void *context);

/**
 * Apply cyclic and reflective boundaries to \c parameter to bring it
 * back within the allowed prior ranges that are specified in \c
 * priorArgs.  LALInferenceCyclicReflectiveBound() should not be
 * called after any multi-parameter update step in a jump proposal,
 * as this violates detailed balance.
 *
 * \param parameter [in] Pointer to an array of parameters
 * \param priorArgs [in] Pointer to an array of prior ranges
 */
void LALInferenceCyclicReflectiveBound(LALInferenceVariables *parameter, LALInferenceVariables *priorArgs);

/**
 * \brief Rotate initial phase if polarisation angle is cyclic around ranges
 *
 * If the polarisation angle parameter \f$\psi\f$ is cyclic about its upper and
 * lower ranges of \f$-\pi/4\f$ to \f$\pi/4\f$ then the transformation for
 * crossing a boundary requires the initial phase parameter \f$\phi_0\f$ to be
 * rotated through \f$\pi\f$ radians. The function assumes the value of
 * \f$\psi\f$ has been rescaled to be between 0 and \f$2\pi\f$ - this is a
 * requirement of the covariance matrix routine \c LALInferenceNScalcCVM
 * function.
 *
 * This is particularly relevant for pulsar analyses.
 *
 * \param parameter [in] Pointer to an array of parameters
 */
void LALInferenceRotateInitialPhase( LALInferenceVariables *parameter );

/**
 * Return the logarithmic prior density of the variables as specified for the sky localisation project
 * (see: https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/SkyLocComparison#priors ),
 * for the non-spinning/spinning inspiral signal case.
 */
REAL8 LALInferenceInspiralSkyLocPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model);

/**
 * Convert the hypercube parameter to physical parameters, for the prior density of the variables as
 * specified for the sky localisation project
 * (see: https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/SkyLocComparison#priors ),
 * for the non-spinning/spinning inspiral signal case.
 */
UINT4 LALInferenceInspiralSkyLocCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model, double *Cube, void *context);

/**
 * Function to add the minimum and maximum values for the uniform prior onto the \c priorArgs.
 */
void LALInferenceAddMinMaxPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 *min, REAL8 *max, LALInferenceVariableType type);

/**
 * Get the minimum and maximum values of the uniform prior from the \c priorArgs list, given a name.
 */
void LALInferenceGetMinMaxPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 *min, REAL8 *max);

/**
 * Function to remove the minimum and maximum values for the uniform prior onto the \c priorArgs.
 */
void LALInferenceRemoveMinMaxPrior(LALInferenceVariables *priorArgs, const char *name);


/**
 * Function to add the mu and sigma values for the Gaussian prior onto the \c priorArgs.
 */
void LALInferenceAddGaussianPrior(LALInferenceVariables *priorArgs,
                                  const char *name, REAL8 *mu, REAL8 *sigma,
                                  LALInferenceVariableType type);

/**
 * Get the mu and sigma values of the Gaussian prior from the \c priorArgs list, given a name.
 */
void LALInferenceGetGaussianPrior(LALInferenceVariables *priorArgs,
                                  const char *name, REAL8 *mu, REAL8 *sigma);


/**
 * Function to remove the mu and sigma values for the Gaussian prior onto the \c priorArgs.
 */
void LALInferenceRemoveGaussianPrior(LALInferenceVariables *priorArgs, const char *name);


/**
 * \brief Add a Gaussian Mixture Model prior
 *
 * Add a Gaussian Mixture Model prior defined by a number of multi-variate Gaussian modes, each with a
 * specified set of means, standard deviations, covariance matrices and weights (where weights are the
 * relative probabilities for each mode). The minumum and maximum allowed prior range for each
 * parameter should also be supplied, although if the array pointers are NULL these ranges will default
 * to +/-infinity.
 *
 * The \c name input should be a colon separated list of all the parameters in the multivariate GMM,
 * e.g. "H0:COSIOTA". The number of parameters in this list will be checked against the number of
 * means supplied for each mode, and the shape of the covariances for each mode to make sure that
 * they are consistent. If just one parameter is supplied (e.g. "H0") then this will just be a
 * one-dimensional GMM.
 *
 * Internally the function will convert the covariance matrices into correlation matrices and inverse
 * correlation matrices for use later (provided they are positive-definite). This will avoid
 * dynamic range/numerical precision issue with using covariances of parameters spanning a large range
 * of values. The standard deviations of each parameter will also be extracted from the covariance
 * matrices and stored, along with the determinants of the covariance matrices.
 */
void LALInferenceAddGMMPrior( LALInferenceVariables *priorArgs, const char *name,
                              REAL8Vector ***mus, gsl_matrix ***covs, REAL8Vector **weights,
                              REAL8Vector **minrange, REAL8Vector **maxrange );

/**
 * \brief Check for a Gaussian Mixture Model prior
 *
 * Check if the single parameter given by \c name has a Gaussian Mixture model prior. If the
 * parameter was within a multivariate GMM prior then it will be found.
 */
int LALInferenceCheckGMMPrior(LALInferenceVariables *priorArgs, const char *name);

/**
 * Remove a Gaussian Mixture Model prior
 */
void LALInferenceRemoveGMMPrior( LALInferenceVariables *priorArgs, const char *name );

/**
 * \brief Get the parameters defining a Gaussian Mixture Model prior
 *
 * For a single parameter given by \c name it will check if that parameter has a GMM prior (even if
 * it is within a multivariate GMM prior). Arrays of the following values for each GMM mode will be
 * returned: means of each parameter; standard deviations of each parameter; a correlation matrix;
 * and inverse correlation matrix; the weight (relative probability) of the mode; and, the determinant
 * of the covariance matrix. The minimum and maximum ranges for each parameter are returned.
 * The position (index) of the parameter \c name within a multivariate GMM will is returned. Finally,
 * the combined name of the prior (i.e. including all parameters) is returned.
 */
void LALInferenceGetGMMPrior( LALInferenceVariables *priorArgs, const char *name,
                              REAL8Vector ***mus, REAL8Vector ***sigmas, gsl_matrix ***cors, gsl_matrix ***invcors,
                              REAL8Vector **weights, REAL8Vector **minrange, REAL8Vector **maxrange,
                              REAL8Vector **dets, UINT4 *idx, CHAR **fullname );

/**
 * \brief Add a log-uniform prior
 *
 * Add a prior uniform in the log, i.e. PDF(x)~1/x
 * \f[p(h|h_{\rm min}, h_{\rm max}, I) = \frac{1/h}{\log{(h_{\rm max}/h_{\rm min})}},\f]
 * where \f$h_{\rm min}\f$ and \f$h_{\rm max}\f$ limit the domain of the PDF.
 * The function has no support outside this range.
 *
 * This function adds \c xmin  and \c xmax values for the Fermi-Dirac prior to the \c priorArgs.
 */
void LALInferenceAddLogUniformPrior(LALInferenceVariables *priorArgs,
                                    const char *name, REAL8 *xmin, REAL8 *xmax,
                                    LALInferenceVariableType type);

/**
 * Get the xmin and xmax values of the log-uniform prior from the \c priorArgs list, given a name.
 */
void LALInferenceGetLogUniformPrior(LALInferenceVariables *priorArgs,
                                    const char *name, REAL8 *xmin, REAL8 *xmax);

/**
 * Function to remove the min and max values for the log-uniform prior from the \c priorArgs.
 */
void LALInferenceRemoveLogUniformPrior(LALInferenceVariables *priorArgs, const char *name);

/**
 * \brief Add a Fermi-Dirac prior
 *
 * Add a prior defined by the Fermi-Dirac PDF
 * \f[p(h|\sigma, r, I) = \frac{1}{\sigma\log{\left(1+e^{r} \right)}}\left(e^{((h/\sigma) - r)} + 1\right)^{-1},\f]
 * where \f$r = \mu/\sigma\f$ to give a more familiar form of the function.
 *
 * This function adds \c sigma  and \c r values for the Fermi-Dirac prior onto the \c priorArgs.
 */
void LALInferenceAddFermiDiracPrior(LALInferenceVariables *priorArgs,
                                    const char *name, REAL8 *sigma, REAL8 *r,
                                    LALInferenceVariableType type);

/**
 * Get the r and sigma values of the Fermi-Dirac prior from the \c priorArgs list, given a name.
 */
void LALInferenceGetFermiDiracPrior(LALInferenceVariables *priorArgs,
                                    const char *name, REAL8 *sigma, REAL8 *r);


/**
 * Function to remove the r and sigma values for the Fermi-Dirac prior onto the \c priorArgs.
 */
void LALInferenceRemoveFermiDiracPrior(LALInferenceVariables *priorArgs, const char *name);

/** Check for types of standard prior */
/** Check for a uniform prior (with minimum and maximum) */
int LALInferenceCheckMinMaxPrior(LALInferenceVariables *priorArgs, const char *name);
/** Check for a Gaussian prior (with a mean and variance) */
int LALInferenceCheckGaussianPrior(LALInferenceVariables *priorArgs, const char *name);
/** Check for a log-uniform prior (with xmin and xmax parameters) */
int LALInferenceCheckLogUniformPrior(LALInferenceVariables *priorArgs, const char *name);
/** Check for a Fermi-Dirac prior (with a r and sigma parameter) */
int LALInferenceCheckFermiDiracPrior(LALInferenceVariables *priorArgs, const char *name);


/**
 * Function to add a correlation matrix and parameter index for a prior
 * defined as part of a multivariate Gaussian distribution onto the \c
 * priorArgs. The correlation coefficient matrix must be a gsl_matrix and the
 * index for the given parameter in the matrix must be supplied. The mean
 * and standard deviation the named parameter must also be supplied.
 */
void LALInferenceAddCorrelatedPrior( LALInferenceVariables *priorArgs,
                                     const char *name, gsl_matrix **cor,
                                     REAL8 *mu, REAL8 *sigma, UINT4 *idx );

/**
 * Get the correlation coefficient matrix and index for a parameter from the
 * \c priorArgs list.
 */
void LALInferenceGetCorrelatedPrior( LALInferenceVariables *priorArgs,
                                     const char *name, gsl_matrix **cor, gsl_matrix **invcor,
                                     REAL8 *mu, REAL8 *sigma, UINT4 *idx );

/**
 * Remove the correlation coefficient matrix and index for a parameter from the
 * \c priorArgs list.
 */
void LALInferenceRemoveCorrelatedPrior( LALInferenceVariables *priorArgs );

/**
 * Check for the existance of a correlation coefficient matrix and index for
 * a parameter from the \c priorArgs list.
 */
int LALInferenceCheckCorrelatedPrior( LALInferenceVariables *priorArgs,
                                      const char *name );

/** Draw variables from the prior ranges */
void LALInferenceDrawFromPrior( LALInferenceVariables *output,
                                LALInferenceVariables *priorArgs,
                                gsl_rng *rdm );

/** Draw an individual variable from its prior range */
void LALInferenceDrawNameFromPrior( LALInferenceVariables *output,
                                    LALInferenceVariables *priorArgs,
                                    char *name, LALInferenceVariableType type,
                                    gsl_rng *rdm );

/* Switch reads true if parameters lie within Malmquist prior */
UINT4 within_malmquist(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model);

/** Prior that is 1 everywhere in component mass space. */
REAL8 LALInferenceAnalyticNullPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model);

/** Analytic null prior converted from hypercube */
UINT4 LALInferenceAnalyticCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model, double *Cube, void *context);

/** Prior that is 1 everywhere. */
REAL8 LALInferenceNullPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model);

/**
 * Computes the numerical normalization of the mass prior \f$p(\mathcal{M}) \sim
 * \mathcal{M}^{-11/6}\f$ applying all cuts in the mass plane implied by the
 * various component, total, and chirp mass limits, and the mass
 * ratio limits.  Returns the integral of \f$\mathcal{M}^{-11/6}\f$ over the allowed
 * ranges in mass.
 */
REAL8 LALInferenceComputePriorMassNorm(const double MMin, const double MMax, const double MTotMax,
                    const double McMin, const double McMax,
                    const double massRatioMin, const double massRatioMax, const char *massRatioName);

/**
 * Prior that checks for minimum and maximum prior range specified in runState->priorArgs
 * and returns 0.0 if sample lies inside the boundaries, -DBL_MAX otherwise.
 * Can be used with MinMaxPrior functions.
 * Ignores variables which are not REAL8 or do not have min and max values set.
 */
REAL8 LALInferenceFlatBoundedPrior(LALInferenceRunState *runState, LALInferenceVariables *params);

/**
 * Utility CubeToPrior functions for psd-fit and both calibration models
 */
UINT4 LALInferenceCubeToPSDScaleParams(LALInferenceVariables *priorParams, LALInferenceVariables *params, INT4 *idx, double *Cube, void *context);
UINT4 LALInferenceCubeToConstantCalibrationPrior(LALInferenceRunState *runState, LALInferenceVariables *params, INT4 *idx, double *Cube, void *context);

/**
 * Prior that converts from a Cube parameter in [0,1] to the flat prior bounded by x1 and x2.
 */
REAL8 LALInferenceCubeToFlatPrior(double r, double x1, double x2);

/**
 * Prior that converts from a Cube parameter in [0,1] to the flat in log prior bounded by x1 and x2.
 */
REAL8 LALInferenceCubeToLogFlatPrior(double r, double x1, double x2);

/**
 * Prior that converts from a Cube parameter in [0,1] to the power prior bounded by x1 and x2 with power p.
 */
REAL8 LALInferenceCubeToPowerPrior(double p, double r, double x1, double x2);

/**
 * Prior that converts from a Cube parameter in [0,1] to the Gaussian prior with given mean and standard deviation.
 */
REAL8 LALInferenceCubeToGaussianPrior(double r, double mean, double sigma);

/**
 * Prior that converts from a Cube parameter in [0,1] to the sine prior with given
 * min (x1) and max (x2) values
 */
REAL8 LALInferenceCubeToSinPrior(double r, double x1, double x2);

/* Simple burst prior (only checks for dec and (log)hrss*/
REAL8 LALInferenceSineGaussianPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model);

/* return the log of the Fermi-Dirac prior */
REAL8 LALInferenceFermiDiracPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 value);

/**
 * \brief Calculate the log probability for the Gaussian Mixture Model prior
 */
REAL8 LALInferenceGMMPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 value);

/* Return the log Prior for a parameter that has a prior that is uniform in log space */
REAL8 LALInferenceLogUniformPrior( LALInferenceVariables *priorArgs, const char *name, REAL8 value );

/*@}*/

#endif
