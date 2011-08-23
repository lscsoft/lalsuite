/* PTMCMC algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */

#include <lal/LALInference.h>

#define COVMATRIXNAME "covarianceMatrix"
#define UNCORRSAMPNAME "uncorrelatedSample"
#define SIGMAVECTORNAME "sigmaJump"

void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState);
void PTMCMCOneStep(LALInferenceRunState *runState);
void PTMCMCAdaptationOneStep(LALInferenceRunState *runState);
//REAL8 PTUniformLALPrior(LALInferenceRunState *runState, LALInferenceVariables *params);
// Von Neumann rejection sampler for the prior !!
void VNRPriorOneStep(LALInferenceRunState *runState);
