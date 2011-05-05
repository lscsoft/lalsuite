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

#include <lal/LALInference.h>
#include <lal/LALInferenceNestedSampler.h>


#ifdef __cplusplus
extern "C" {
#endif

/** define macros */
#define ROUND(x) (floor(x+0.5))

#define MAXLENGTH 1000000

#define SIXTH 0.166666666666666666666666666666666666666667L
#define TWENTYFOURTH 0.04166666666666666666666666666666666666666667L

/* default values */
#define CHUNKMIN 5
#define CHUNKMAX 30
#define PSIBINS 50
#define TIMEBINS 1440

/** define functions */

/* initialisation functions */
void initialiseAlgorithm( LALInferenceRunState *runState );

void readPulsarData( LALInferenceRunState *runState );

void setupFromParFile( LALInferenceRunState *runState );

void setupLookupTables(LALInferenceRunState *runState, LALSource *source);

UINT4 add_initial_variables( LALInferenceVariables *ini, LALInferenceVariables *scaleFac,
                             LALInferenceVariables *priorArgs, 
                             BinaryPulsarParams pars ); 
  
UINT4 add_variable_scale_prior( LALInferenceVariables *var, LALInferenceVariables *scale, 
                                LALInferenceVariables *prior, const char *name, 
                                REAL8 value, REAL8 sigma);

void initialiseProposal( LALInferenceRunState *runState );

void setupLivePointsArray( LALInferenceRunState *runState );

/* likelihood and prior */
REAL8 pulsar_log_likelihood( LALInferenceVariables *vars, LALInferenceIFOData *data,
                             LALInferenceTemplateFunction *get_pulsar_model );
                             
REAL8 priorFunction( LALInferenceRunState *runState, LALInferenceVariables *params );

/* model functions */
void get_pulsar_model( LALInferenceIFOData *data );

REAL8Vector *get_phase_model( BinaryPulsarParams params, LALInferenceIFOData *data );

void get_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data );
  
REAL8 noise_only_model( LALInferenceIFOData *data );
  
/* helper functions */
UINT4Vector *get_chunk_lengths( LALInferenceIFOData *data, INT4 chunkMax );

REAL8Vector * sum_data( LALInferenceIFOData *data );

void response_lookup_table( REAL8 t0, LALDetAndSource detAndSource,
                            INT4 timeSteps, INT4 psiSteps, gsl_matrix *LUfplus,
                            gsl_matrix *LUfcross );
                            
REAL8 log_factorial(UINT4 num);

void rescaleOutput( LALInferenceRunState *runState );

#ifdef __cplusplus
}
#endif

#endif /* _PULSAR_PARAMETER_ESTIMATION_H */
