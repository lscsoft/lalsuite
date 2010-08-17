/*******************************************************************************
  Matt Pitkin - 08/08/07

  pulsar_parameter_estimation.h

  Header file for pulsar_parameter_estimation.c

*******************************************************************************/

/*
  Author:
  $Id$
*/

#ifndef _PULSAR_PARAMETER_ESTIMATION_H
#define _PULSAR_PARAMETER_ESTIMATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>

/* LAL headers */
#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/Random.h>
#include <lal/LALString.h>
#include <lal/SFTutils.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/MatrixUtils.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/LALRCSID.h>
#include <lal/ComputeFstat.h>
#include <lal/TimeSeries.h>

#include <lalapps.h>

#include "LALInference.h"
#include "LALInferenceNestedSampler.h"


#ifdef __cplusplus
extern "C" {
#endif

/** define macros */
#define ROUND(x) (floor(x+0.5))

/* macro to perform logarithmic addition log(exp(x) + exp(y)) */
#define PLUS(x,y) ( x>y ? x+log(1.+exp(y-x)) : y+log(1.+exp(x-y)) )

/* macro to calculate the area of a trapezium for on logs */
/* logarea = log(0.5) + log(width) + log(exp(logHeight1) + exp(logHeight2)) */
#define LOG_TRAPEZIUM(h1, h2, w) ( -0.693147180559945 + log(w) + PLUS(h1, h2) )

#define MAXLENGTH 1000000
#define MAXPARAMS 35

#define SIXTH 0.166666666666666666666666666666666666666667L
#define TWENTYFOURTH 0.04166666666666666666666666666666666666666667L

/** define functions */
	
/* function to return the (REAL8) log factorial of an integer */
REAL8 log_factorial(INT4 num);

UINT4Vector *get_chunk_lengths( LALIFOData *data, INT4 chunkMax );

REAL8Vector *get_phase_model( BinaryPulsarParams params, LALIFOData *data );

void add_initial_variables( LALVariables *ini, BinaryPulsarParams pars );
	
REAL8Vector * sum_data( LALIFOData *data );

void get_amplitude_model( BinaryPulsarParams pars, LALIFOData *data );

REAL8 priorFunction(LALInferenceRunState *runState, LALVariables *params);

LALIFOData *readPulsarData(int argc, char *argv[]);

void initialiseAlgorithm(LALInferenceRunState *runState);

void setupLookupTables(LALInferenceRunState *runState, LALSource *source);

void initialiseProposal(LALInferenceRunState *runState);

void setupFromParFile(LALInferenceRunState *runState);

void setupLivePointsArray(LALInferenceRunState *runState);

REAL8 pulsar_log_likelihood( LALVariables *vars, LALIFOData *data,
  LALTemplateFunction *get_pulsar_model );

void get_pulsar_model( LALIFOData *data );

void response_lookup_table(REAL8 t0, LALDetAndSource detAndSource,
  INT4 timeSteps, INT4 psiSteps, gsl_matrix *LUfplus, gsl_matrix *LUfcross);
	
#ifdef __cplusplus
}
#endif

#endif /* _PULSAR_PARAMETER_ESTIMATION_H */
