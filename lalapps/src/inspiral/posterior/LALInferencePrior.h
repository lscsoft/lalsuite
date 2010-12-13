#include "LALInference.h"

REAL8 LALInferenceInspiralPrior(LALInferenceRunState *runState, LALVariables *variables);
void LALInferenceCyclicReflectiveBound(LALVariables *parameter, LALVariables *priorArgs);
REAL8 LALInferenceInspiralPriorNormalised(LALInferenceRunState *runState, LALVariables *params);
