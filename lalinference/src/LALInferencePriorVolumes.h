/**
 * \file LALInferencePriorVolumes.h
 * \brief Header file for initialisation functions used by LALInference codes
 *
 */
#ifndef LALInferencePriorVolumes_h
#define LALInferencePriorVolumes_h
#include <lal/LALInference.h>

double LALInferenceMassPriorVolume(LALInferenceRunState *state);
double LALInferenceMassDistancePriorVolume(LALInferenceRunState *state);

#endif
