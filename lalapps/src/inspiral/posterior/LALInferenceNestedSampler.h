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

void calcCVM(gsl_matrix **cvm, LALVariables **Live, UINT4 Nlive);
double logadd(double a,double b);

void NestedSamplingOneStep(LALInferenceRunState *runState);

REAL8 mean(REAL8 *array,int N);
REAL8 sample_logt(int Nlive,gsl_rng *RNG);
REAL8 ang_dist(REAL8 a1, REAL8 a2);
REAL8 ang_var(LALVariables **list,const char *pname, int N);

/* Functions for proposal distributions */

void LALInferenceProposalNS(LALInferenceRunState *runState, LALVariables *parameter);
void LALInferenceProposalMultiStudentT(LALInferenceRunState *runState, LALVariables *parameter);
//Declared in LALInferencePrior.h instead:
//void LALInferenceCyclicReflectiveBound(LALVariables *parameter, LALVariables *priorArgs);
INT4 LALInferenceReflectDetPlane(
							 LALInferenceRunState *state,
							 LALVariables *parameter
								 );
void LALInferenceRotateSky(
						   LALInferenceRunState *state,
						   LALVariables *parameter
						   );
void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3]);
void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude);
void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude);



/* Returns a REAL4vector drawn from the multivariate normal distribution
 with covariance matrix matrix, mean zero, of dimension dim. randParam is
 an already initialised RandomParams structure */

void XLALMultiNormalDeviates(REAL4Vector *vector,
							 gsl_matrix *matrix,
							 UINT4 dim,
							 RandomParams *randParam);

/* Returns a REAL4vector drawn from the multivariate student-t distribution
 with n degrees of freedom, with covariance matrix matrix, mean zero,
 of dimension dim. randParam is an already initialised RandomParams structure */

void
XLALMultiStudentDeviates(REAL4Vector  *vector,
						 gsl_matrix   *matrix,
						 UINT4         dim,
						 UINT4         n,
						 RandomParams *randParam);


/* Check that the gsl_matrix is positive definite. dim = number of dimensions */
UINT4 LALInferenceCheckPositiveDefinite(gsl_matrix *matrix, UINT4 dim);
