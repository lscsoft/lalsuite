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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */
#ifndef LALInferenceLikelihood_h
#define LALInferenceLikelihood_h

#include <lal/LALInference.h>

/**
 * \defgroup LALInferenceLikelihood_h Header LALInferenceLikelihood.h
 * \ingroup pkg_LALInference
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
 * - a pointer to the LALInferenceTemplateFunction template function to be used.
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
/*@{*/

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
REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, 
                              LALInferenceTemplateFunction templt);

REAL8 LALInferenceNoiseOnlyLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceTemplateFunction templt);

/**
 * For testing purposes (for instance sampling the prior),
 * likelihood that returns 0.0 = log(1) every
 * time.  Activated with the --zeroLogLike command flag.
 */
REAL8 LALInferenceZeroLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceTemplateFunction templt);

/***********************************************************//**
 * (log-) likelihood function.                                 
 * Returns the non-normalised logarithmic likelihood.          
 * Slightly slower but cleaner than							   
 * UndecomposedFreqDomainLogLikelihood().          `		   
 *
 * Required (`currentParams') parameters are:                  
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        
 *   - "distance"        (REAL8, Mpc, >0)                      
 *   - "time"            (REAL8, GPS sec.)                     
 ***************************************************************/
REAL8 LALInferenceFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data,
                              LALInferenceTemplateFunction templt);

/***********************************************************//**
 * Chi-Square function.                                        
 * Returns the chi square of a template:                       
 * chisq= p * sum_i (dx_i)^2, with dx_i  =  <s,h>_i  - <s,h>/p
 * 
 * Required (`currentParams') parameters are:                  
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        
 *   - "distance"        (REAL8, Mpc, >0)                      
 *   - "time"            (REAL8, GPS sec.)                     
 ***************************************************************/
REAL8 LALInferenceChiSquareTest(LALInferenceVariables *currentParams, LALInferenceIFOData * data,
                              LALInferenceTemplateFunction templt);

/***********************************************************//**
 * Frequency-domain single-IFO response computation.           
 * Computes response for a given template.                    
 * Will re-compute template only if necessary                  
 * (i.e., if previous, as stored in data->freqModelhCross,     
 * was based on different parameters or template function).    
 * Carries out timeshifting for a given detector               
 * and projection onto this detector.                          
 * Result stored in freqResponse, assumed to be correctly      
 * initialized												   
 *
 * Required (`currentParams') parameters are:                  
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        
 *   - "distance"        (REAL8, Mpc, >0)                      
 *   - "time"            (REAL8, GPS sec.)                     
 ***************************************************************/				
void LALInferenceComputeFreqDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr,
                              LALInferenceTemplateFunction templt, COMPLEX16Vector *freqWaveform);

/**
 * Computes the <x|y> overlap in the Fourrier domain.
 */
REAL8 LALInferenceComputeFrequencyDomainOverlap(LALInferenceIFOData * data,
        COMPLEX16Vector * freqData1, COMPLEX16Vector * freqData2);

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
REAL8 LALInferenceFreqDomainStudentTLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data,
                                      LALInferenceTemplateFunction templt);

/**
 * An analytic likeilhood that is a correlated Gaussian in 15
 * dimensions.
 */
REAL8 LALInferenceCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams,
                                                  LALInferenceIFOData *data,
                                                  LALInferenceTemplateFunction templt);

/**
 * An analytic likeilhood that is two correlated Gaussians in 15
 * dimensions.
 */
REAL8 LALInferenceBimodalCorrelatedAnalyticLogLikelihood(LALInferenceVariables *currentParams,
                                                  LALInferenceIFOData *data,
                                                  LALInferenceTemplateFunction templt);

/**
 * 15-D Rosenbrock log(L) function (see Eq (3) of
 * http://en.wikipedia.org/wiki/Rosenbrock_function .
 */
REAL8 LALInferenceRosenbrockLogLikelihood(LALInferenceVariables *currentParams,
                                          LALInferenceIFOData *data,
                                          LALInferenceTemplateFunction templt);

REAL8 LALInferenceMarginalisedPhaseLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data,LALInferenceTemplateFunction templt);

/**
 * Returns the log-likelihood marginalised over the time dimension.
 * The boundaries of the time marginalization are the data segment;
 * any prior on time is ignored.
 */
REAL8 LALInferenceMarginalisedTimeLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, LALInferenceTemplateFunction templt);

/**
 * Initialisation function which reads runState->commaneLine and sets up the
 * likelihood function accordingly. Can choose between Gaussian, Student-t, marginalised
 * phase likelihoods
 */
void LALInferenceInitLikelihood(LALInferenceRunState *runState);

/** Get the intrinsic parameters from currentParams */
LALInferenceVariables LALInferenceGetInstrinsicParams(LALInferenceVariables *currentParams);

/** Switch for controlling the removal of lines */
INT4 LALInferenceLineSwitch(INT4 lineFlag, INT4 Nlines, INT4 *lines_array, INT4 *widths_array, INT4 i);

/*@}*/

#endif
