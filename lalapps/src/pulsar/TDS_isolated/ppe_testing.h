/* ******************************************************************************
  Matt Pitkin - 2011

  ppe_testing.h

  Header file for ppe_testing.c

*******************************************************************************/

/*
  Author:
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin
 *
 * \brief Header file for functions used in testing the parameter
 * estimation code for known pulsar searches using the nested sampling
 * algorithm.
 */

#ifndef _PPE_TESTING_H
#define _PPE_TESTING_H

#include "pulsar_parameter_estimation_nested.h"

#ifdef __cplusplus
extern "C" {
#endif

/** The usage format for the test case of performing the analysis on a
 * one-dimensional grid. */
#define USAGEGRID \
"Usage: %s [options]\n\n"\
" --grid              perform the posterior evalution on a 1D grid over the\n\
                     parameter given by --gridpar\n"\
" --gridpar           The parameter over which to perform the 1D posterior\n\
                     evaluation\n"\
" --gridmin           The lower end of the range over which to evaluate the\n\
                     paramter given by --gridpar\n"\
" --gridmax           The upper end of the range over which to evaluate the\n\
                     paramter given by --gridpar\n"\
" --gridsteps         The number of points in the grid\n"\
"\n"

/* testing functions */
void gridOutput( LALInferenceRunState *runState );

REAL8 test_gaussian_log_likelihood( LALInferenceVariables *vars,
                                    LALInferenceIFOData *data,
                                    LALInferenceTemplateFunction get_model );

void outputPriorSamples( LALInferenceRunState *runState );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_TESTING_H */
