#include "LALInference.h"

REAL8 LALInferenceInspiralPrior(LALInferenceRunState *runState, LALVariables *variables);
void LALInferenceCyclicReflectiveBound(LALVariables *parameter, LALVariables *priorArgs);
REAL8 LALInferenceInspiralPriorNormalised(LALInferenceRunState *runState, LALVariables *params);
double innerIntegrand(double M2, void *viData);
double outerIntegrand(double M1, void *voData);
double computePriorMassNorm(const double MMin, const double MMax, const double MTotMax, const double McMin, const double McMax, const double etaMin, const double etaMax);
