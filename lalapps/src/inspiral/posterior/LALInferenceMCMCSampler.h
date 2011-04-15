/* PTMCMC algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */

#include "LALInference.h"

#define COVMATRIXNAME "covarianceMatrix"
#define UNCORRSAMPNAME "uncorrelatedSample"
#define SIGMAVECTORNAME "sigmaJump"

void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState);
void PTMCMCOneStep(LALInferenceRunState *runState);
void PTMCMCAdaptationOneStep(LALInferenceRunState *runState);
//REAL8 PTUniformLALPrior(LALInferenceRunState *runState, LALVariables *params);
void PTMCMCLALProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALBlockProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALSingleProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
//void PTMCMCLALProposaltemp(LALInferenceRunState *runState, LALVariables *proposedParams);

/* Makes correlated jumps if the square root (Cholesky decomp) of a
   covariance matrix has been specified; if not, falls back to PTMCMCLALBlockProposal. */
void PTMCMCLALBlockCorrelatedProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALSingleAdaptProposal(LALInferenceRunState *runState, LALVariables *proposedParams);

/* Flip inclination about observing plane (iota -> Pi - iota) */
void PTMCMCLALInferenceInclinationFlip(LALInferenceRunState *runState, LALVariables *proposedParams);

/* Rotate one or both spins about L. */
void PTMCMCLALInferenceRotateSpins(LALInferenceRunState *runState, LALVariables *proposedParams);

/* Increment orbital phase by Pi */
void PTMCMCLALInferenceOrbitalPhaseJump(LALInferenceRunState *runState, LALVariables *proposedParams);

/* Choose a random covariance matrix eigenvector to jump along. */
void PTMCMCLALInferenceCovarianceEigenvectorJump(LALInferenceRunState *runState, LALVariables *proposedParams);

/* Jump around by 0.01 radians in angle on the sky */
void PTMCMCLALInferenceSkyLocWanderJump(LALInferenceRunState *runState, LALVariables *proposedParams);

//REAL8 GaussianLikelihood(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template);
//REAL8 UnityLikelihood(LALVariables *currentParams, LALIFOData * data, LALTemplateFunction *template);
//REAL8 PTUniformGaussianPrior(LALInferenceRunState *runState, LALVariables *params);
//void PTMCMCGaussianProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALAdaptationProposal(LALInferenceRunState *runState, LALVariables *proposedParams);
void PTMCMCLALAdaptationSingleProposal(LALInferenceRunState *runState, LALVariables *proposedParams);

void PTMCMCLALInferenceRotateSky(LALInferenceRunState *state,LALVariables *parameter);
INT4 PTMCMCLALInferenceReflectDetPlane(LALInferenceRunState *state,LALVariables *parameter);

/* Jump in iota and distance so that approximate magnitude remains
   constant in one of the detectors. */
void PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump(LALInferenceRunState *state, LALVariables *proposedParams);

/* Differential evolution */
void PTMCMCLALInferenceDifferentialEvolutionFull(LALInferenceRunState *state, LALVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionNames(LALInferenceRunState *state, LALVariables *proposedParams, const char *names[]);
void PTMCMCLALInferenceDifferentialEvolutionMasses(LALInferenceRunState *state, LALVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionAmp(LALInferenceRunState *state, LALVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionSpins(LALInferenceRunState *state, LALVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionSky(LALInferenceRunState *state, LALVariables *proposedParams);
