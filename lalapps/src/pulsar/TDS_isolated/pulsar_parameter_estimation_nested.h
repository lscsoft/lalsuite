/* ******************************************************************************
  Matt Pitkin, John Veitch - 2011

  pulsar_parameter_estimation_nested.h

  Header file for pulsar_parameter_estimation_nested.c

*******************************************************************************/

/*
  Author:
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
#include <dirent.h>
#include <complex.h>
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
#include <lal/ComputeFstat.h>
#include <lal/TimeSeries.h>
#include <lal/LALNoiseModels.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/StringVector.h>
#include <lal/XLALGSL.h>

#include <lalapps.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceProposal.h>

#include <lal/LALSimNoise.h>

#ifdef HAVE_LIBLALXML
#include <lal/LALInferenceXML.h>
#endif

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/** Macro to round a value to the nearest integer. */
#define ROUND(x) (floor(x+0.5))

/** Macro to perform addition of two values within logarithm space \f$ \log{(e^x + e^y)}\f$. */
#define LOGPLUS(x,y) ( x>y ? x+log(1.+exp(y-x)) : y+log(1.+exp(x-y)) )

/** Macro that gives the integer number of times that \c x goes in to \c y. */
#define FACTOR(x,y) ((INT4)floor(x/y))

/** Macro to square a value. */
#define SQUARE(x) ( (x) * (x) );

/**
 * The maximum allowable length of the input data stream. Note: this may be
 * removed in the future with memory allocated dynamically.
 */
#define MAXLENGTH 1000000

/* default values */
/** Default value of the minimum length into which the data can be split. */
#define CHUNKMIN 5
/** Default value of the maximum length into which the data can be split. */
#define CHUNKMAX 0
/**
 * Default number of bins in polarisation angle \f$ \psi \f$ for the time vs.
 * \f$ \psi \f$ antenna pattern lookup table.
 */
#define PSIBINS 1000
/**
 * Default number of bins in time (over one sidereal day) for the time vs.
 * \f$ \psi \f$ antenna pattern lookup table.
 */
#define TIMEBINS 2880

/**
 * The total number of 'amplitude' parameters that can define a signal e.g. gravitational wave amplitude from a
 * triaxial star emitting from the \f$l=m=2\f$ mode we have \f$ h_0 \f$, initial phase of the signal \f$ \phi_0 \f$,
 * polarisation angle \f$ psi \f$, and cosine of the inclination angle \f$ \cos{\iota} \f$. Or, more generally for
 * emission from \f$l=2\f$ and \f$m=1,2\f$ instead of \f$ h_0 \f$ and \f$ \phi_0 \f$ there can be complex amplitude and
 * phase parameters \f$C_{22}\f$, \f$C_{21}\f$, \f$\phi_{22}\f$ and \f$\phi_{21}\f$.
 *
 * Note: These should be increased if additional model parameters are added.
 */
#define NUMAMPPARS 8

/**
 * The total number of frequency parameters that can defined a signal e.g.
 * the signal frequency and its time derivatives, and the frequency (period)
 * epoch.
 */
#define NUMFREQPARS 8

/**
 * The total number of sky position parameters that can define a signal e.g.
 * right ascension, declination, proper motion and the positional epoch.
 */
#define NUMSKYPARS 5

/**
 * The total number of binary system parameters that can define a signal e.g.
 * binary period, orbital eccentricity, projected semi-major axis, time of
 * periastron and angle of periastron.
 */
#define NUMBINPARS 33

/**
 * A list of the amplitude parameters. The names given here are those that are
 * recognised within the code.
 */
static const CHAR amppars[NUMAMPPARS][VARNAME_MAX] = { "h0", "phi0", "psi",
"cosiota", "C22", "C21", "phi22", "phi21" };

/**
 * A list of the frequency parameters. The names given here are those that are
 * recognised within the code.
 */
static const CHAR freqpars[NUMFREQPARS][VARNAME_MAX] = { "f0", "f1", "f2", "f3",
"f4", "f5", "pepoch", "cgw" };

/**
 * A list of the sky position parameters. The names given here are those that
 * are recognised within the code.
 */
static const CHAR skypars[NUMSKYPARS][VARNAME_MAX] = { "ra", "pmra", "dec",
"pmdec", "posepoch" };

/**
 * A list of the binary system parameters. The names given here are those that
 * are recognised within the code.
 */
static const CHAR binpars[NUMBINPARS][VARNAME_MAX] = { "Pb", "e", "eps1",
"eps2", "T0", "Tasc", "x", "w0", "Pb2", "e2", "T02", "x2", "w02", "Pb3", "e3",
"T03", "x3", "w03", "xpbdot", "eps1dot", "eps2dot", "wdot", "gamma", "Pbdot",
"xdot", "edot", "s", "dr", "dth", "a0", "b0", "M", "m2" };

/**
 * A flag to specify if phase parameters are being searched over and
 * therefore the pulsar model requires phase evolution to be re-calculated (0 =
 * no, 1 = yes).
 */
extern UINT4 varyphase;

/**
 * A flag to specify if the sky position will be searched over, and therefore
 * whether the solar system barycentring needs recalculating (0 = no, 1 = yes).
 */
extern UINT4 varyskypos;

/**
 * A flag to specify if the binary system parameters will be searched over,
 * and therefore whether the binary system barycentring needs recalculating (0 =
 * no, 1 = yes)
 */
extern UINT4 varybinary;

extern LALStringVector *corlist;

/* initialisation functions */
void initialiseAlgorithm( LALInferenceRunState *runState );

void readPulsarData( LALInferenceRunState *runState );

void setupFromParFile( LALInferenceRunState *runState );

void setupLookupTables(LALInferenceRunState *runState, LALSource *source);

void add_initial_variables( LALInferenceVariables *ini, LALInferenceVariables *scaleFac, BinaryPulsarParams pars );

void add_variable_scale( LALInferenceVariables *var, LALInferenceVariables *scale, const char *name, REAL8 value );

void initialisePrior( LALInferenceRunState *runState );

void initialiseProposal( LALInferenceRunState *runState );

void add_correlation_matrix( LALInferenceVariables *ini,
                             LALInferenceVariables *priors, REAL8Array *corMat,
                             LALStringVector *parMat );

/* software injection functions */
void injectSignal( LALInferenceRunState *runState );

/* helper functions */
UINT4Vector *get_chunk_lengths( LALInferenceIFOData *data, INT4 chunkMax );

UINT4Vector *chop_n_merge( LALInferenceIFOData *data, INT4 chunkMin, INT4 chunkMax );

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

TimeCorrectionType XLALAutoSetEphemerisFiles( CHAR *efile, CHAR *sfile,
                                              CHAR *tfile,
                                              BinaryPulsarParams pulsar,
                                              INT4 gpsstart, INT4 gpsend );

void samples_prior( LALInferenceRunState *runState );

#ifdef __cplusplus
}
#endif

#endif /* _PULSAR_PARAMETER_ESTIMATION_H */
