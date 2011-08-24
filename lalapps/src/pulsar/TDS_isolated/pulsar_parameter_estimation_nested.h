/*******************************************************************************
  Matt Pitkin, John Veitch - 2011

  pulsar_parameter_estimation_nested.h

  Header file for pulsar_parameter_estimation_nested.c

*******************************************************************************/

/*
  Author:
  $Id$
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Header file for the parameter estimation code for known pulsar
 * searches using the nested sampling algorithm.
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

/** Macro to round a value to the nearest integer. */
#define ROUND(x) (floor(x+0.5))

/** Macro to perform addition of two values within logarithm space \f$ \log{(e^x
+ e^y)}\f$. */
#define LOGPLUS(x,y) ( x>y ? x+log(1.+exp(y-x)) : y+log(1.+exp(x-y)) )

/** Macro that gives the integer number of times that \c x goes in to \c y. */
#define FACTOR(x,y) ((INT4)floor(x/y))

/** Macro to square a value. */
#define SQUARE(x) ( (x) * (x) );

/** The maximum allowable length of the input data stream. Note: this may be
 * removed in the future with memory allocated dynamically. */
#define MAXLENGTH 1000000

/* default values */
/** Default value of the minimum length into which the data can be split. */
#define CHUNKMIN 5
/** Default value of the maximum length into which the data can be split. */
#define CHUNKMAX 0
/** Default number of bins in polarisation angle \f$ \psi \f$ for the time vs.
 * \f$ \psi \f$ antenna pattern lookup table. */
#define PSIBINS 500
/** Default number of bins in time (over one sidereal day) for the time vs.
 * \f$ \psi \f$ antenna pattern lookup table. */
#define TIMEBINS 2880

/** The total number of 'amplitude' parameters that can define a signal e.g.
 * gravitational wave amplitude from a triaxial star \f$ h_0 \f$, initial phase
 * of the signal \f$ \phi_0 \f$, polarisation angle \f$ psi \f$, and cosine of
 * the inclination angle \f$ \cos{\iota} \f$. 
 * 
 * Note: These should be increased if additional model parameters are added.
 */ 
#define NUMAMPPARS 7
/** A list of the amplitude parameters. The names given here are those that are
 * recognised within the code. */
CHAR amppars[NUMAMPPARS][VARNAME_MAX] = { "h0", "phi0", "psi", "cosiota", "h1",
                                          "lambda", "theta" };

/** The total number of frequency parameters that can defined a signal e.g.
 * the signal frequency and its time derivatives, and the frequency (period)
 * epoch. */
#define NUMFREQPARS 7
/** A list of the frequency parameters. The names given here are those that are
 * recognised within the code. */
CHAR freqpars[NUMFREQPARS][VARNAME_MAX] = { "f0", "f1", "f2", "f3", "f4", "f5",
                                            "pepoch" };

/** The total number of sky position parameters that can define a signal e.g.
 * right ascension, declination, proper motion and the positional epoch. */ 
#define NUMSKYPARS 5
/** A list of the sky position parameters. The names given here are those that
 * are recognised within the code. */
CHAR skypars[NUMSKYPARS][VARNAME_MAX] = { "ra", "pmra", "dec", "pmdec",
                                          "posepoch" };
     
/** The total number of binary system parameters that can define a signal e.g.
 * binary period, orbital eccentricity, projected semi-major axis, time of
 * periastron and angle of periastron. */
#define NUMBINPARS 33
/** A list of the binary system parameters. The names given here are those that
 * are recognised within the code. */
CHAR binpars[NUMBINPARS][VARNAME_MAX] = { "Pb", "e", "eps1", "eps2", "T0",
                                          "Tasc", "x", "w0", "Pb2", "e2", "T02",
                                          "x2", "w02", "Pb3", "e3", "T03", "x3",
                                          "w03", "xpbdot", "eps1dot", "eps2dot",
                                          "wdot", "gamma", "Pbdot", "xdot",
                                          "edot", "s", "dr", "dth", "a0", "b0",
                                          "M", "m2" };


/* initialisation functions */
void initialiseAlgorithm( LALInferenceRunState *runState );

void readPulsarData( LALInferenceRunState *runState );

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
                             
REAL8 priorFunction( LALInferenceRunState *runState, 
                     LALInferenceVariables *params );

/* model functions */
void get_pulsar_model( LALInferenceIFOData *data );

REAL8 rescale_parameter( LALInferenceIFOData *data, const CHAR *parname );

void pulsar_model( BinaryPulsarParams params, 
                                LALInferenceIFOData *data );

REAL8Vector *get_phase_model( BinaryPulsarParams params, 
                              LALInferenceIFOData *data,
                              REAL8 freqFactor );

REAL8Vector *get_ssb_delay( BinaryPulsarParams pars, 
                            LIGOTimeGPSVector *datatimes,
                            EphemerisData *ephem,
                            LALDetector *detector,
                            REAL8 interptime );
                            
REAL8Vector *get_bsb_delay( BinaryPulsarParams pars,
                            LIGOTimeGPSVector *datatimes,
                            REAL8Vector *dts );                
                              
void get_triaxial_amplitude_model( BinaryPulsarParams pars, 
                                   LALInferenceIFOData *data );

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

void sumData( LALInferenceRunState *runState );

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
