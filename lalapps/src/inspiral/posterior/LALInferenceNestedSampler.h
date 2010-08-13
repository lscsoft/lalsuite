/* Nested Sampling algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */

#include "LALInference.h"

/* logadd(a,b) = log(exp(a) + exp(b)) using Stirling's approximation */
double logadd(double a,double b);

/* NestedSamplingAlgorithm implements the nested sampling algorithm,
 see e.g. Sivia "Data Analysis: A Bayesian Tutorial, 2nd edition */
void NestedSamplingAlgorithm(LALInferenceRunState *runState);

void LALInferenceProposalDifferentialEvolution(LALInferenceRunState *runState,
											   LALVariables *parameter);

void LALInferenceProposalNS(LALInferenceRunState *runState, LALVariables *parameter);
void LALInferenceProposalMultiStudentT(LALInferenceRunState *runState, LALVariables *parameter);
void LALInferenceCyclicReflectiveBound(LALVariables *parameter, LALVariables *priorArgs);
void calcCVM(gsl_matrix *cvm, LALVariables **Live, UINT4 Nlive);
double logadd(double a,double b);

void NestedSamplingOneStep(LALInferenceRunState *runState);
