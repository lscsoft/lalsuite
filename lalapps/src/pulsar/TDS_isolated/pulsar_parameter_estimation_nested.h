/*******************************************************************************
  Matt Pitkin, John Veitch - 2011

  pulsar_parameter_estimation_nested.h

  Header file for pulsar_parameter_estimation_nested.c

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
#include <sys/time.h>

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
#include <lal/LALNoiseModels.h>
#include <lal/Units.h>

#include <lalapps.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>

#ifdef __cplusplus
extern "C" {
#endif

/** define macros */
#define ROUND(x) (floor(x+0.5))

/* macro to perform logarithmic addition log(exp(x) + exp(y)) */
#define LOGPLUS(x,y) ( x>y ? x+log(1.+exp(y-x)) : y+log(1.+exp(x-y)) )
#define SQUARE(x) ( (x) * (x) );

#define MAXLENGTH 1000000

/* default values */
#define CHUNKMIN 5
#define CHUNKMAX 0
#define PSIBINS 500
#define TIMEBINS 2880

/* INCREASE THESE VALUES IF ADDING ADDITIONAL PARAMETERS */
/* number of amplitude parameters e.g. h0, phi0, psi, ciota */
#define NUMAMPPARS 4
CHAR amppars[NUMAMPPARS][VARNAME_MAX] = { "h0", "phi0", "psi", "cosiota" };

/* number of frequency parameters e.g. f0 */
#define NUMFREQPARS 7
CHAR freqpars[NUMFREQPARS][VARNAME_MAX] = { "f0", "f1", "f2", "f3", "f4", "f5",
                                            "pepoch" };

/* number of sky position parameters e.g. ra, dec */ 
#define NUMSKYPARS 5
CHAR skypars[NUMSKYPARS][VARNAME_MAX] = { "ra", "pmra", "dec", "pmdec",
                                          "posepoch" };

/* number of binary parameters e.g. e, x */       
#define NUMBINPARS 33
CHAR binpars[NUMBINPARS][VARNAME_MAX] = { "Pb", "e", "eps1", "eps2", "T0",
                                          "Tasc", "x", "w0", "Pb2", "e2", "T02",
                                          "x2", "w02", "Pb3", "e3", "T03", "x3",
                                          "w03", "xpbdot", "eps1dot", "eps2dot",
                                          "wdot", "gamma", "Pbdot", "xdot",
                                          "edot", "s", "dr", "dth", "a0", "b0",
                                          "M", "m2" };

/** define functions */

/* initialisation functions */
void initialiseAlgorithm( LALInferenceRunState *runState );

void readPulsarData( LALInferenceRunState *runState );

void readDoublePulsarData( LALInferenceRunState *runState );

void setSignalModelType( LALInferenceRunState *runState );

void setupFromParFile( LALInferenceRunState *runState );

void setupLookupTables(LALInferenceRunState *runState, LALSource *source);

void add_initial_variables( LALInferenceVariables *ini, 
                            LALInferenceVariables *scaleFac,
                            LALInferenceVariables *priorArgs, 
                            BinaryPulsarParams pars ); 
  
void add_variable_scale_prior( LALInferenceVariables *var, 
                               LALInferenceVariables *scale, 
                               LALInferenceVariables *prior, const char *name, 
                               REAL8 value, REAL8 sigma );

void initialiseProposal( LALInferenceRunState *runState );

void setupLivePointsArray( LALInferenceRunState *runState );

/* likelihood and prior */
REAL8 pulsar_log_likelihood( LALInferenceVariables *vars, 
                             LALInferenceIFOData *data,
                             LALInferenceTemplateFunction *get_pulsar_model );

REAL8 pulsar_double_log_likelihood( LALInferenceVariables *vars, 
                             LALInferenceIFOData *data,
                             LALInferenceTemplateFunction *get_pulsar_model );
                             
REAL8 priorFunction( LALInferenceRunState *runState, 
                     LALInferenceVariables *params );

/* model functions */
void get_pulsar_model( LALInferenceIFOData *data );

REAL8 rescale_parameter( LALInferenceIFOData *data, const CHAR *parname );

void get_triaxial_pulsar_model( BinaryPulsarParams params, 
                                LALInferenceIFOData *data );

void get_pinsf_pulsar_model( BinaryPulsarParams params, 
                                LALInferenceIFOData *data );

REAL8Vector *get_phase_model( BinaryPulsarParams params, 
                              LALInferenceIFOData *data );

REAL8Vector *get_ssb_delay( BinaryPulsarParams pars, 
                            LIGOTimeGPSVector *datatimes,
                            EphemerisData *ephem,
                            LALDetector *detector,
                            REAL8 interptime );
                            
REAL8Vector *get_bsb_delay( BinaryPulsarParams pars,
                            LIGOTimeGPSVector *datatimes,
                            REAL8Vector *dts );                
                              
void get_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data );

void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data );
  
REAL8 noise_only_model( LALInferenceIFOData *data );

/* software injection functions */
void injectSignal( LALInferenceRunState *runState );

/* helper functions */
UINT4Vector *get_chunk_lengths( LALInferenceIFOData *data, INT4 chunkMax );

UINT4Vector *chop_n_merge( LALInferenceIFOData *data, INT4 chunkMin, 
                           INT4 chunkMax );

COMPLEX16Vector *subtract_running_median( COMPLEX16Vector *data );

UINT4Vector *chop_data( COMPLEX16Vector *data, INT4 chunkMin );

UINT4 find_change_point( COMPLEX16Vector *data, REAL8 *logodds, INT4 chunkMin );

void rechop_data( UINT4Vector *segs, INT4 chunkMax, INT4 chunkMin );

void merge_data( COMPLEX16Vector *data, UINT4Vector *segs );

REAL8Vector * sum_data( LALInferenceIFOData *data );

void response_lookup_table( REAL8 t0, LALDetAndSource detAndSource,
                            INT4 timeSteps, INT4 psiSteps, gsl_matrix *LUfplus,
                            gsl_matrix *LUfcross );

void rescaleOutput( LALInferenceRunState *runState );

INT4 count_csv( CHAR *csvline );

INT4 recognised_parameter( CHAR *parname );

REAL8 calculate_time_domain_snr( LALInferenceIFOData *data );

void get_loudest_snr( LALInferenceRunState *runState );

/* testing functions */
void gridOutput( LALInferenceRunState *runState );

REAL8 test_gaussian_log_likelihood( LALInferenceVariables *vars,
                                    LALInferenceIFOData *data,
                                    LALInferenceTemplateFunction *get_model );

#ifdef __cplusplus
}
#endif

#endif /* _PULSAR_PARAMETER_ESTIMATION_H */
