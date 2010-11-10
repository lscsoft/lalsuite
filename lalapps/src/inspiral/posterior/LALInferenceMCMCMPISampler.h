/* PTMCMC algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */

#include "LALInference.h"

void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState);
void PTMCMCOneStep(LALInferenceRunState *runState);
void PTMCMCAdaptationOneStep(LALInferenceRunState *runState);
//REAL8 PTUniformLALPrior(LALInferenceRunState *runState, LALVariables *params);
void PTMCMCLALProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALProposaltemp(LALInferenceRunState *runState, LALVariables *proposedParams);

//REAL8 GaussianLikelihood(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template);
//REAL8 UnityLikelihood(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template);
//REAL8 PTUniformGaussianPrior(LALInferenceRunState *runState, LALVariables *params);
//void PTMCMCGaussianProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALAdaptationProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALAdaptationSingleProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
