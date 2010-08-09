/* Nested Sampling algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */

#include "LALInference.h"

/* NestedSamplingAlgorithm implements the nested sampling algorithm,
 see e.g. Sivia "Data Analysis: A Bayesian Tutorial, 2nd edition */
void NestedSamplingAlgorithm(LALInferenceRunState *runState);

/* NestedSamplingOneIteration advances the state of the algorithm
 by one iteration of the nested sampling algorithm */
void NestedSamplingOneIteration(LALInferenceRunState *runState);
