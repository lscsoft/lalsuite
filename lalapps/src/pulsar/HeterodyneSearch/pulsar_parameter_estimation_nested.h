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
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Header file for the parameter estimation code for known pulsar
 * searches using the nested sampling algorithm.
 */

#ifndef _PULSAR_PARAMETER_ESTIMATION_NESTED_H
#define _PULSAR_PARAMETER_ESTIMATION_NESTED_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <time.h>
#include <sys/time.h>

/* LAL headers */
#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
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
#include <lal/FileIO.h>

#include <lalapps.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceGenerateROQ.h>

#include <lal/LALSimNoise.h>

#ifdef HAVE_LIBLALXML
#include <lal/LALInferenceXML.h>
#endif

/* check whether openmp is enabled and if so include omp.h */
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>

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
#define SQUARE(x) ( (x) * (x) )

/**
 * The maximum allowable length of the input data stream. Note: this may be
 * removed in the future with memory allocated dynamically.
 */
#define MAXLENGTH 1000000

/**
 * The maximum line length (in characters) of a heterodyned data file.
 */
#define PPEN_MAXLINELENGTH 1024

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
#define NUMAMPPARS 24

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
#define NUMBINPARS 34

/** The maximum number of different detectors allowable. */
#define MAXDETS 6

/** The usage format for the code.  */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --detectors         all IFOs with data to be analysed e.g. H1,H2\n\
                     (delimited by commas) (if generating fake data these\n\
                     should not be set)\n"\
" --par-file          pulsar parameter (.par) file (full path)\n"\
" --cor-file          pulsar TEMPO-fit parameter correlation matrix\n"\
" --input-files       full paths and file names for the data for each\n\
                     detector and model harmonic in the list (must be in the\n\
                     same order) delimited by commas. These files can be gzipped.\n\
                     If not set you can generate fake data (see --fake-data below)\n"\
" --sample-interval   (REAL8) the time interval bewteen samples (default to 60 s)\n"\
" --outfile           name of output data file [required]\n"\
" --non-fixed-only    output only the non-fixed (i.e. variable) parameters\n\
                     specified in the prior and .par files.\n"\
" --gzip              gzip the output text file\n"\
" --outXML            name of output XML file [not required]\n"\
" --chunk-min         (INT4) minimum stationary length of data to be used in\n\
                     the likelihood e.g. 5 mins\n"\
" --chunk-max         (INT4) maximum stationary length of data to be used in\n\
                     the likelihood e.g. 30 mins\n"\
" --time-bins         (INT4) no. of time bins in the time-psi lookup table\n"\
" --prior-file        file containing the parameters to search over and\n\
                     their upper and lower ranges\n"\
" --ephem-earth       Earth ephemeris file\n"\
" --ephem-sun         Sun ephemeris file\n"\
" --ephem-timecorr    Einstein delay time correction ephemeris file\n"\
" --harmonics         (CHAR) the signal model frequency harmonics that you want\n\
                     to use (delimited by commas). Currently this can be either\n\
                     the 'triaxial' model for which you use 2 (the default\n\
                     value) or 1,2 for a model with emission at both the rotation\n\
                     frequency and twice the rotation frequency.\n"\
" --biaxial           Set this if the waveform model parameters spcify a biaxial star\n"\
" --gaussian-like     Set this if a Gaussian likelihood is to be used. If the input\n\
                     file contains a column specifying the noise standard deviation of\n\
                     the data then that will be used in the Gaussian likelihood function,\n\
                     otherwise the noise variance will be calculated from the data.\n"\
" --nonGR             Set to allow non-GR polarisation modes and/or a variable\n\
                     speed of gravitational waves.\n"\
" --randomise         Set this to randomise the data (through permutations of the\n\
                     time stamps) for use in Monte-Carlo studies. NOTE: this will not\n\
                     work if using the code to create injections.\n"\
"\n"\
" Nested sampling parameters:\n"\
" --Nlive             (INT4) no. of live points for nested sampling\n"\
" --Nmcmc             (INT4) no. of for MCMC used to find new live points\n\
                     (if not specified an adaptive number of points is used)\n"\
" --Nmcmcinitial      (INT4) no. of MCMC points to use in the initial resampling of\n\
                     the prior (default is to use MAXMCMC).\n"\
" --Nruns             (INT4) no. of parallel runs\n"\
" --tolerance         (REAL8) tolerance of nested sampling integrator\n"\
" --randomseed        seed for random number generator\n"\
"\n"\
" MCMC proposal parameters:\n"\
" --covariance        (UINT4) relative weight of using covariance matrix\n\
                     of the live points as the proposal (DEFAULT = 0,\n\
                     e.g. 0%%)\n"\
" --temperature       (REAL8) temperature for covariance proposal\n\
                     distribution (DEFAULT = 0.1)\n"\
" --kDTree            (UINT4) relative weight of using a k-D tree of the live\n\
                     points to use as a proposal (DEFAULT = 0, e.g. 0%%)\n"\
" --kDNCell           (INT4) maximum number of samples in a k-D tree cell\n"\
" --kDUpdateFactor    (REAL8) how often the k-D tree gets updated as a\n\
                     factor of the number of live points\n"\
" --diffev            (UINT4) relative weight of using differential evolution\n\
                     of the live points as the proposal (DEFAULT = 0, e.g.\n\
                     0%%)\n"\
" --freqBinJump       (UINT4) relative weight of using jumps to adjacent\n\
                     frequency bins as a proposal (DEFAULT = 0, e.g. this is\n\
                     not required unless searching over frequency)\n"\
" --ensembleStretch   (UINT4) relative weight of the ensemble stretch\n\
                     proposal (DEFAULT = 15, e.g. 37.5%%)\n"\
" --ensembleWalk      (UINT4) relative weight of the ensemble walk\n\
                     proposal (DEFAULT = 15, e.g. 37.5%%)\n"\
" --uniformprop       (UINT4) relative weight of a proposal that draws\n\
                     points uniformly from any flat prior ranges.\n\
                     (DEFAULT = 10, e.g. 25%%)\n"\
"\n"\
" Reduced order quadrature parameters:\n"\
" --roq               Set this to use reduced order quadrature to compute the\n\
                     likelihood.\n"\
" --ntraining         (UNIT4) The number of training models used to generate an\n\
                     orthonormal basis of waveform models.\n"\
" --roq-tolerance     (REAL8) The tolerance used during the basis generation\n\
                     (DEFAULT = 1e-11)\n"\
" --test-basis        If this is set then the reduced basis set will be tested\n\
                     against another set of waveforms to check they really are\n\
                     within the required tolerance.\n"\
" --output-weights    (CHAR) If this is set then the weights will be output to\n\
                     the (binary) file that is named and the programme will\n\
                     exit. These could be read in later instead of being\n\
                     regenerated. This allows the ROQ to be generated on a\n\
                     machine with a large amount of RAM, whilst the full\n\
                     parameter estimation can run on a machine with less RAM.\n"\
" --input-weights     (CHAR) A binary file containing all the weights in a\n\
                     defined format. If this is present then the RQO will\n\
                     not be recalculated.\n"\
"\n"\
" Signal injection parameters:\n"\
" --inject-file       a pulsar parameter (par) file containing the parameters\n\
                     of a signal to be injected. If this is given a signal\n\
                     will be injected\n"\
" --inject-output     a filename to which the injected signal will be\n\
                     output if specified\n"\
" --fake-data         a list of IFO's for which fake data will be generated\n\
                     e.g. H1,L1 (delimited by commas). Unless the --fake-psd\n\
                     flag is set the power spectral density for the data will\n\
                     be generated from the noise models in LALNoiseModels.\n\
                     For Advanced detectors (e.g Advanced LIGO) prefix the\n\
                     name with an A (e.g. AH1 for an Advanced LIGO detector\n\
                     at Hanford). The noise will be white across the data\n\
                     band width.\n"\
" --fake-psd          if you want to generate fake data with specific power\n\
                     spectral densities for each detector giving in\n\
                     --fake-data then they should be specified here delimited\n\
                     by commas (e.g. for --fake-data H1,L1 then you could use\n\
                     --fake-psd 1e-48,1.5e-48) where values are single-sided\n\
                     PSDs in Hz^-1.\n"\
" --fake-starts       the start times (in GPS seconds) of the fake data for\n\
                     each detector separated by commas (e.g.\n\
                     910000000,910021000). If not specified these will all\n\
                     default to 900000000\n"\
" --fake-lengths      the length of each fake data set (in seconds) for each\n\
                     detector separated by commas. If not specified these\n\
                     will all default to 86400 (i.e. 1 day)\n"\
" --fake-dt           the data sample rate (in seconds) for the fake data for\n\
                     each detector. If not specified this will default to\n\
                     60s.\n"\
" --scale-snr         give a (multi-detector) SNR value to which you want to\n\
                     scale the injection. This is 1 by default.\n"\
"\n"\
" Flags for using a Nested sampling file as a prior:\n"\
" --sample-files     a list of (comma separated) file containing the nested\n\
                    samples from a previous run of the code (these should\n\
                    contain samples in ascending likelihood order and be\n\
                    accompanied by a file containg a list of the parameter\n\
                    names for each column with the suffix _params.txt). If\n\
                    this is set this WILL be used as the only prior, but the\n\
                    prior ranges set in the --prior-file and --par-file are\n\
                    still needed (and should be consistent with the variable\n\
                    parameters in the nested sample file).\n"\
" --sample-nlives    a list (comma separated) of the number of live point\n\
                    that where used when creating each nested sample file.\n"\
" --prior-cell       The number of samples to use in a k-d tree cell for the\n\
                    prior (the default will be 8).\n"\
" --Npost            The (approximate) number of posterior samples to be\n\
                    generated from each nested sample file (default = 1000)\n"\
"\n"\
" Legacy code flags:\n"\
" --oldChunks        Set if using fixed chunk sizes for dividing the data as\n\
                    in the old code, rather than the calculating chunks\n\
                    using the change point method.\n"\
" --jones-model      Set if using both 1 and 2 multiples of the frequency and\n\
                    requiring the use of the original signal model parameters\n\
                    from Jones, MNRAS, 402 (2010).\n"\
"\n"\
" Benchmarking:\n"\
" --time-it          Set if wanting to time the various parts of the code.\n\
                    A output file with the \"outfile\" filename appended with\n\
                    \"_timings\" will contain the timings.\n"\
" --sampleprior      (UINT4) Set this to be a number of samples generated from\n\
                    the prior. The nested sampling will not be performed.\n"\
"\n"

/**
 * A list of the amplitude parameters. The names given here are those that are
 * recognised within the code.
 */
static const CHAR amppars[NUMAMPPARS][VARNAME_MAX] = { "H0", "PHI0", "PSI",
"COSIOTA", "C22", "C21", "PHI22", "PHI21", "HSCALARB", "HSCALARL", "HVECTORX",
"HVECTORY", "PSIVECTOR", "PHI0VECTOR", "PSISCALAR", "PHI0SCALAR", "PSITENSOR",
"PHI0TENSOR", "I21", "I31", "LAMBDA", "COSTHETA", "IOTA", "THETA" };

/**
 * A list of the frequency parameters. The names given here are those that are
 * recognised within the code.
 */
static const CHAR freqpars[NUMFREQPARS][VARNAME_MAX] = { "F0", "F1", "F2", "F3",
"F4", "F5", "PEPOCH", "CGW" };

/**
 * A list of the sky position parameters. The names given here are those that
 * are recognised within the code.
 */
static const CHAR skypars[NUMSKYPARS][VARNAME_MAX] = { "RA", "PMRA", "DEC",
"PMDEC", "POSEPOCH" };

/**
 * A list of the binary system parameters. The names given here are those that
 * are recognised within the code.
 */
static const CHAR binpars[NUMBINPARS][VARNAME_MAX] = { "PB", "ECC", "EPS1",
"EPS2", "T0", "TASC", "A1", "OM", "PB_2", "ECC_2", "T0_2", "A1_2", "OM_2", "PB_3", "ECC_3",
"T0_3", "A1_3", "OM_3", "XPBDOT", "EPS1DOT", "EPS2DOT", "OMDOT", "GAMMA", "PBDOT",
"XDOT", "EDOT", "SINI", "DR", "DTHETA", "A0", "B0", "MTOT", "M2", "FB" };

extern LALStringVector *corlist;

extern UINT4 verbose_output;

#ifdef __cplusplus
}
#endif

#endif /* _PULSAR_PARAMETER_ESTIMATION__NESTED_H */
