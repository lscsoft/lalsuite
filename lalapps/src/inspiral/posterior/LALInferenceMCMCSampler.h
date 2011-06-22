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
void PTMCMCLALProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALBlockProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALSingleProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
//void PTMCMCLALProposaltemp(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Makes correlated jumps if the square root (Cholesky decomp) of a
   covariance matrix has been specified; if not, falls back to PTMCMCLALBlockProposal. */
void PTMCMCLALBlockCorrelatedProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALSingleAdaptProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Flip inclination about observing plane (iota -> Pi - iota) */
void PTMCMCLALInferenceInclinationFlip(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Rotate one or both spins about L. */
void PTMCMCLALInferenceRotateSpins(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Increment orbital phase by Pi */
void PTMCMCLALInferenceOrbitalPhaseJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Choose a random covariance matrix eigenvector to jump along. */
void PTMCMCLALInferenceCovarianceEigenvectorJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Jump around by 0.01 radians in angle on the sky */
void PTMCMCLALInferenceSkyLocWanderJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

//REAL8 GaussianLikelihood(LALInferenceVariables *currentParams, LALIFOData * data, LALTemplateFunction *template);
//REAL8 UnityLikelihood(LALInferenceVariables *currentParams, LALIFOData * data, LALTemplateFunction *template);
//REAL8 PTUniformGaussianPrior(LALInferenceRunState *runState, LALInferenceVariables *params);
//void PTMCMCGaussianProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALAdaptationProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALAdaptationSingleProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

void PTMCMCLALInferenceRotateSky(LALInferenceRunState *state,LALInferenceVariables *parameter);
INT4 PTMCMCLALInferenceReflectDetPlane(LALInferenceRunState *state,LALInferenceVariables *parameter);

/* Jump in iota and distance so that approximate magnitude remains
   constant in one of the detectors. */
void PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump(LALInferenceRunState *state, LALInferenceVariables *proposedParams);

/* Differential evolution */
void PTMCMCLALInferenceDifferentialEvolutionFull(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionNames(LALInferenceRunState *state, LALInferenceVariables *proposedParams, const char *names[]);
void PTMCMCLALInferenceDifferentialEvolutionMasses(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionAmp(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionSpins(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionSky(LALInferenceRunState *state, LALInferenceVariables *proposedParams);

/*draws a value from the prior, uniformly in individual parameters used for jumps.*/
void PTMCMCLALInferenceDrawFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
/*draws a value from the prior, using Von Neumann rejection sampling.*/
void PTMCMCLALInferenceDrawUniformlyFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
// Von Neumann rejection sampler for the prior !!
void VNRPriorOneStep(LALInferenceRunState *runState);
