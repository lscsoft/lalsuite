/*
*  Copyright (C) 2012 Matthew Pitkin, Colin Gill, John Veitch
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Parameter estimation code for known pulsar searches using the nested sampling algorithm.

\heading{Description}
This code is used to perform parameter estimation and evidence calculation in targeted/semi-targeted searches for
gravitational waves from known pulsars. It may also be used to follow-up on signal candidates from semi-coherent all-sky
searches for unknown sources.

It uses the Bayesian technique of 'Nested Sampling' to sample over a defined prior parameter space (unknown signal
parameters such as the gravitational wave amplitude). These samples can then be used to create posterior probability
density functions of the unknown parameters. The samples are also used to calculate the Bayesian evidence, also known as
the marginal likelihood, for a given signal model. This can be compared with other models, in particular the model that
the data is described by Gaussian noise alone.

As input the code requires time domain data that has been heterodyned using the known (or close to) phase evolution of
the pulsar. The time domain input should consist of a three column text file containing the GPS time stamp of the data
point, the real part of the heterodyned data and the imaginary part of the heterodyned data, e.g.
900000000.000000  1.867532e-24  -7.675329e-25
900000060.000000  2.783651e-24  3.654386e-25
...

Most commonly such data will have a sample rate of 1/60 Hz, giving a bandwidth of the same amount, but the code can
accept any rate, or downsample data (by averaging) by a given factor.

The code also requires that you specify which parameters are to be searched over, and the prior ranges over these. Any
of the signal parameters can be searched over, including frequency, sky position and binary system parameters, although
the bandwidth of the data and search efficiency need to be taken into account.

The 'Nested Sampling' algorithm (developed by [\ref Skilling2006]) used is that defined in LALinferenceNestedSampler
(see [\ref VeitchVecchio2010]). It is essentially an efficient way to perform the integral
\f[
Z = \int^{\mathbf{\theta}} p(d|\mathbf{\theta}) p(\mathbf{\theta}) \mathrm{d}\mathbf{\theta},
\f]
where \f$ \mathbf{\theta} \f$ is a vector of parameters, \f$ p(d|\mathbf{\theta}) \f$ is the likelihood of the data
given the parameters, and \f$ p(\mathbf{\theta}) \f$ is the prior on the parameters. It does this by changing the
multi-dimensional integral over N parameters into a one-dimensional integral
\f[
Z = \int^X L(X) \mathrm{d}X \approx \sum_i L(X_i) \Delta{}X_i,
\f]
where \f$ L(X) \f$ is the likelihood, and \f$ X \f$ is the prior mass. The algorithm will draw a number (\f$ N \f$) of
samples (live points) from the parameter priors, calculate the likelihood for each point and find the lowest likelihood
value. The lowest likelihood value will be added to the summation in the above equation, with \f$ \log{\Delta{}X_i}
\approx 1/N \f$ coming from the fact that the prior would be normalised to unity and therefore each point should occupy
an equal fraction and at each iteration the prior volume will decrease geometrically (for \f$\log{\Delta{}X_0} = 0\f$).
A new point is then drawn from the prior with the criterion that it has a higher likelihood than the previous lowest
point and substitutes that point. To draw the new point a Markov Chain Monte Carlo (MCMC) procedure is used - there are
two methods used to sample points within this: i) drawing from a proposal distributions based on the covariance matrix
if the current live points (although to keep things computationally efficient this no updated at every iteration), ii)
picking a point via differential evolution (two random live points are selected and a new point half way between the two
is created). The probability of using either method is currently set at 80\% and 20\% respectively. The procedure is
continued until a stopping criterion is reached, which in this case is that the remaining prior volume is less than the
\c tolerance value set (see below). The implementation of this can be seen in [\ref VeitchVecchio2010].

\heading{Usage}
The usage format is given below and can also be found by running the code with
\code
lalapps_pulsar_parameter_estimation_nested --help
\endcode

An example of running the code on to search over the four unknown parameters \f$ h_0 \f$, \f$ \phi_0 \f$, \f$ \psi \f$
and \f$ \cos{\iota} \f$, for pulsar J0534-2200, given heterodyned time domain data from the H1 detector in the file
\c finehet_J0534-2200_H1, is:
\code
lalapps_pulsar_parameter_estimation_nested --detectors H1 --par-file
J0534-2200.par --input-files finehet_J0534-2200_H1 --outfile ns_J0534-2200
--prior-file prior_J0534-2200.txt --ephem-earth
lscsoft/share/lalpulsar/earth05-09.dat --ephem-sun
lscsoft/share/lalpulsar/sun05-09.dat --model-type triaxial --Nlive 1000 --Nmcmc
100 --Nruns 1 --tolerance 0.25
\endcode
The \c par-file is a TEMPO-style file containing the parameters of the pulsar used to perform the heterodyne (the
frequency parameters are the rotation frequency and therefore not necessarily the gravitational wave frequency) e.g.
\code
RA      12:54:31.87523895
DEC     -54:12:43.6572033
PMRA    1.7
PMDEC   2.8
POSEPOCH 54320.8531
F0  123.7896438753
F1  4.592e-15
PEPOCH 54324.8753
\endcode
The \c prior-file is a text file containing a list of the parameters to be searched over, the prior type ("uniform" or
"gaussian" and their given lower/mean and upper/standard deviation ranges e.g.
\code
h0 uniform 0 1e-21
phi0 uniform 0 6.283185307179586
cosiota uniform -1 1
psi uniform -0.785398163397448 0.785398163397448
\endcode
Note that if searching over frequency parameters the ranges specified in the \c prior-file should be given in terms of
the pulsars rotation frequency and not necessarily the gravitational wave frequency e.g. for a triaxial star emitting
gravitational waves at 100 Hz (which will be at twice the rotation frequency) if you wanted to search over 99.999 to
100.001 Hz then you should used
\code
f0 49.9995 50.0005
\endcode

An example of running the code as above, but this time on fake data created using the Advanced LIGO design noise curves
and with a signal injected into the data is:
\code
lalapps_pulsar_parameter_estimation_nested --fake-data AH1 --inject-file
fake.par --par-file fake.par --outfile ns_fake --prior-file prior_fake.txt
--ephem-earth lscsoft/share/lalpulsar/earth05-09.dat --ephem-sun
lscsoft/share/lalpulsar/sun05-09.dat --model-type triaxial --Nlive 1000 --Nmcmc
100 --Nruns 1 --tolerance 0.25
\endcode
In this case the \c inject-file parameter file must contain the values of \c h0, \c phi0, \c psi and \c cosiota,
otherwise these will be set to zero by default. The parameter files given for \c inject-file and \c par-file do not
have to be the same - the injection can be offset from the 'heterodyned' values, which will be reflected in the data. If
an \c inject-output file is also specified then the fake data containing the signal, and a fake signal-only data set,
will be output.
 */

#include "pulsar_parameter_estimation_nested.h"
#include "ppe_models.h"
#include "ppe_likelihood.h"
#include "ppe_testing.h"

/** The maximum number of different detectors allowable. */
#define MAXDETS 6

/* global variables */
/** An array to contain the log of factorials up to a certain number. */
REAL8 *logfactorial = NULL;

UINT4 verbose_output = 0;
UINT4 varyphase = 0;
UINT4 varyskypos = 0;
UINT4 varybinary = 0;

LALStringVector *corlist = NULL;

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
                     same order) delimited by commas. If not set you can\n\
                     generate fake data (see --fake-data below)\n"\
" --sample-interval   (REAL8) the time interval bewteen samples (default to 60 s)\n"\
" --downsample-factor (INT4) factor by which to downsample the input data\n\
                     (default is for no downsampling and this is NOT\n\
                     applied to fake data)\n"\
" --outfile           name of output data file [required]\n"\
" --outXML            name of output XML file [not required]\n"\
" --chunk-min         (INT4) minimum stationary length of data to be used in\n\
                     the likelihood e.g. 5 mins\n"\
" --chunk-max         (INT4) maximum stationary length of data to be used in\n\
                     the likelihood e.g. 30 mins\n"\
" --psi-bins          (INT4) no. of psi bins in the time-psi lookup table\n"\
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
"\n"\
" Nested sampling parameters:\n"\
" --Nlive             (INT4) no. of live points for nested sampling\n"\
" --Nmcmc             (INT4) no. of for MCMC used to find new live points\n\
                     (if not specified an adaptive number of points is used)\n"\
" --Nruns             (INT4) no. of parallel runs\n"\
" --tolerance         (REAL8) tolerance of nested sampling integrator\n"\
" --randomseed        seed for random number generator\n"\
"\n"\
" MCMC proposal parameters:\n"\
" --covariance        (REAL8) relative weigth of using covariance matrix\n\
                     of the live points as the proposal (DEFAULT = 14,\n\
                     e.g. 70%%)\n"\
" --temperature       (REAL8) temperature for covariance proposal\n\
                     distribution (DEFAULT = 0.1)\n"\
" --kDTree            (REAL8) relative weigth of using a k-D tree of the live\n\
                     points to use as a proposal (DEFAULT = 3, e.g. 15%%)\n"\
" --kDNCell           (INT4) maximum number of samples in a k-D tree cell\n"\
" --kDUpdateFactor    (REAL8) how often the k-D tree gets updated as a\n\
                     factor of the number of live points\n"\
" --diffev            (REAL8) relative weight of using differential evolution\n\
                     of the live points as the proposal (DEFAULT = 3, e.g.\n\
                     15%%)\n"\
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
" --Npost            The (approxiamate) number of posterior samples to be\n\
                    generated from each nested sample file (default = 1000)\n"\
"\n"\
" Phase parameter search speed-up factors:\n"\
" --mismatch         Maximum allowed phase mismatch between consecutive\n\
                    models (if phase mismatch is small do not update phase\n\
                    correction)\n"\
" --mm-factor        (INT4) Downsampling factor for the phase models used\n\
                    to calculate the mismatch\n"\
 "\n"\
" Legacy code flags:\n"\
" --oldChunks        Set if using fixed chunk sizes for dividing the data as\n\
                    in the old code, rather than the calculating chunks\n\
                    using the change point method.\n"\
"\n"

INT4 main( INT4 argc, CHAR *argv[] ){
  ProcessParamsTable *param_table;
  LALInferenceRunState runState;
  REAL8 logZnoise = 0.;

  /* set error handler to abort in main function */
  XLALSetErrorHandler( XLALExitErrorHandler );

  /* Get ProcParamsTable from input arguments */
  param_table = LALInferenceParseCommandLine( argc, argv );
  runState.commandLine = param_table;

  /* Initialise the algorithm structures from the command line arguments */
  /* Include setting up random number generator etc */
  initialiseAlgorithm( &runState );

  /* read in data */
  readPulsarData( &runState );

  /* set algorithm to use Nested Sampling */
  runState.algorithm = &LALInferenceNestedSamplingAlgorithm;
  runState.evolve = &LALInferenceNestedSamplingOneStep;

  /* set likelihood function */
  runState.likelihood = &pulsar_log_likelihood;

  /* set prior function */
  runState.prior = &priorFunction;

  /* set signal model/template */
  runState.templt = get_pulsar_model;

  /* set output style (change this when the code if fixed for using XML) */
  runState.logsample = LALInferenceLogSampleToFile;

  /* Generate the lookup tables and read parameters from par file */
  setupFromParFile( &runState );

  /* add injections if requested */
  injectSignal( &runState );

  /* create sum square of the data to speed up the likelihood calculation */
  sumData( &runState );

  /* Initialise the prior distribution given the command line arguments */
  initialisePrior( &runState );

  gridOutput( &runState );

  /* get noise likelihood and add as variable to runState */
  logZnoise = noise_only_model( &runState );
  LALInferenceAddVariable( runState.algorithmParams, "logZnoise", &logZnoise, LALINFERENCE_REAL8_t,
                           LALINFERENCE_PARAM_FIXED );

  /* Create live points array and fill initial parameters */
  LALInferenceSetupLivePointsArray( &runState );

  /* output the live points sampled from the prior */
  outputPriorSamples( &runState );

  /* Initialise the MCMC proposal distribution */
  initialiseProposal( &runState );

  /* Call the nested sampling algorithm */
  runState.algorithm( &runState );

  /* get SNR of highest likelihood point */
  get_loudest_snr( &runState );

  /* re-read in output samples and rescale appropriately */
  rescaleOutput( &runState );

  return 0;
}

/******************************************************************************/
/*                      INITIALISATION FUNCTIONS                              */
/******************************************************************************/

/** \brief Initialises the nested sampling algorithm control
 *
 * Memory is allocated for the parameters, priors and proposals. The nested sampling control parameters are set: the
 * number of live points \c Nlive, the number of points for each MCMC \c Nmcmc, the number of independent runs within
 * the algorithm \c Nruns, and the stopping criterion \c tolerance.
 *
 * The random number generator is initialise (the GSL Mersenne Twister algorithm \c gsl_rng_mt19937) using either a user
 * defined seed \c randomseed, the system defined \c /dev/random file, or the system clock time.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void initialiseAlgorithm( LALInferenceRunState *runState )
{
  ProcessParamsTable *ppt = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;
  REAL8 tmp;
  INT4 tmpi;
  INT4 randomseed;
  UINT4 verbose = 0;

  FILE *devrandom = NULL;
  struct timeval tv;

  /* print out help message */
  ppt = LALInferenceGetProcParamVal( commandLine, "--help" );
  if(ppt){
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  /* Initialise parameters structure */
  runState->algorithmParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
  runState->priorArgs = XLALCalloc( 1, sizeof(LALInferenceVariables) );
  runState->proposalArgs = XLALCalloc( 1, sizeof(LALInferenceVariables) );

  ppt = LALInferenceGetProcParamVal( commandLine, "--verbose" );
  if( ppt ) {
    LALInferenceAddVariable( runState->algorithmParams, "verbose", &verbose , LALINFERENCE_UINT4_t,
                             LALINFERENCE_PARAM_FIXED );
    verbose_output = 1;
  }

  /* Number of live points */
  ppt = LALInferenceGetProcParamVal( commandLine, "--Nlive" );
  if( ppt ){
    tmpi = atoi( LALInferenceGetProcParamVal(commandLine, "--Nlive")->value );
    LALInferenceAddVariable( runState->algorithmParams,"Nlive", &tmpi, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
  }
  else{
   XLALPrintError("Error... Number of live point must be specified.\n");
   XLAL_ERROR_VOID(XLAL_EIO);
  }

  /* Number of points in MCMC chain */
  ppt = LALInferenceGetProcParamVal( commandLine, "--Nmcmc" );
  if( ppt ){
    tmpi = atoi( LALInferenceGetProcParamVal(commandLine, "--Nmcmc")->value );
    LALInferenceAddVariable( runState->algorithmParams, "Nmcmc", &tmpi, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
  }

  /* set sloppiness! */
  ppt = LALInferenceGetProcParamVal(commandLine,"--sloppyfraction");
  if( ppt ) { tmp = atof(ppt->value); }
  else { tmp = 0.0; }
  LALInferenceAddVariable( runState->algorithmParams, "sloppyfraction", &tmp, LALINFERENCE_REAL8_t,
                           LALINFERENCE_PARAM_OUTPUT );

  /* Optionally specify number of parallel runs */
  ppt = LALInferenceGetProcParamVal( commandLine, "--Nruns" );
  if( ppt ) {
    tmpi = atoi( ppt->value );
    LALInferenceAddVariable( runState->algorithmParams, "Nruns", &tmpi, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
  }

  /* Tolerance of the Nested sampling integrator */
  ppt = LALInferenceGetProcParamVal( commandLine, "--tolerance" );
  if( ppt ){
    tmp = strtod( ppt->value, (char **)NULL );
    LALInferenceAddVariable( runState->algorithmParams, "tolerance", &tmp, LALINFERENCE_REAL8_t,
                             LALINFERENCE_PARAM_FIXED );
  }

  /* Set up the random number generator */
  gsl_rng_env_setup();
  runState->GSLrandom = gsl_rng_alloc( gsl_rng_mt19937 );

  /* (try to) get random seed from command line: */
  ppt = LALInferenceGetProcParamVal( commandLine, "--randomseed" );
  if ( ppt != NULL ) { randomseed = atoi( ppt->value ); }
  else { /* otherwise generate "random" random seed: */
    if ( (devrandom = fopen("/dev/random","r")) == NULL ) {
      gettimeofday( &tv, 0 );
      randomseed = tv.tv_sec + tv.tv_usec;
    }
    else {
      if( fread(&randomseed, sizeof(randomseed), 1, devrandom) != 1 ){
        fprintf(stderr, "Error... could not read random seed\n");
        exit(3);
      }
      fclose( devrandom );
    }
  }

  gsl_rng_set( runState->GSLrandom, randomseed );

  return;
}

/** \brief Reads in the input gravitational wave data files, or creates fake data sets.
 *
 * The function will check whether data files are being input of fake data is to be generated. If using real data the \c
 * detectors command line input must list the names of the detectors from which each of the data sets comes, with names
 * separated by commas - allowed detector names are H1, H2, L1, V1, G1 or T1, for LIGO Hanford Observatory 1, LIGO
 * Hanford Observatory 2, LIGO Livingston Observatory, Virgo, GEO600/HF and TAMA300 respectively. If requiring fake data
 * to be generated then the \c fake-data command line argument must list the detectors for which fake data is required
 * (again separated by commas) - these can be the same names as above, although if an 'A' is appended to the LIGO of
 * Virgo detector name then the Advanced detector is assumed (for use if generating data from designed sensitivity
 * curves).
 *
 * If generating fake data the noise floor can either be created using models of the detector design sensitivities (if
 * just \c fake-data is set), or with noise levels defined using the \c fake-psd command line argument. If using \c
 * fake-psd it should list the signle-sided power spectral density required for each detector listed in \c fake-data
 * (again separated by commas). If using the design sensitivity models the \c par-file argument will be used to find the
 * noise at the correct frequency, which is here assumed to be twice the rotation frequency. The start time (in GPS
 * seconds), data length (in seconds) and sample interval (in seconds) can be specified for each fake data set by
 * supplying comma separated lists to the \c fake-starts, \c fake-lengths and \c fake-dt command line arguments. By
 * default these values are GPS 900000000 (13th July 2008 at 15:59:46 UTC), 86400 seconds (1 solar day) and 60 seconds
 * (1/60 Hz sample rate) respectively. The fake data is drawn from a normal distribution using \c XLALNormalDeviates and
 * scaled with the appropriate PSD.
 *
 * The number of data files required to be read in, or number of fake data sets generated depends on the pulsar model
 * type, which is specified by the number of frequency harmonics given by the \c harmonics command line argument. This
 * should be a list of comma separated values giving the frequency of the signal harmonics to be included. E.g.
 * <ol>
 * <li>For a model with \f$l=m=2\f$ (i.e. a triaxial star with a signal defined in e.g. [\ref DupuisWoan2005]), which
 * purely emits at twice the rotation frequency, the \c harmonics input would just be \c 2. This requires that for each
 * specified detector there is <b>one</b> input file containing data heterodyned at twice the rotation frequency of the
 * pulsar.</li>
 * <li>For a model including the two harmonics \f$l=2\f$, \f$m=1,2\f$ (see e.g. [\ref Jones2010]), which produces
 * emission at both the rotation frequency <i>and</i> twice the rotation frequency, the \c harmonics input would be \c
 * 1,2. This requires that for each specified detector there are two input files containing data heterodyned and the
 * rotation frequency <i>and</i> twice the rotation frequency (these must be in the same order as the harmonics
 * are listed).</li>
 * </ol>
 * The default model, if none is specified, is the triaxial source model with emission at just twice the rotation
 * frequency. At the moment only these two models above can be used, although this could be expanded in the future.
 *
 * If creating fake data for a specific model then the number of PSDs, start time, lengths and sample intervals
 * specified must be equivalent to the number of input files that would have been required e.g. if using the pinned
 * superfluid model and requiring data for H1 and L1 then four fake PSDs would be required (the first pair at the
 * pulsars rotation frequency and twice that in H1, and the seconds pair at the pulsars rotation frequency and twice
 * that in L1). These most be specified in the same order as the detectors.
 *
 * Any new models added can require and arbitrary number of inputs for a given detector, however the heterodyned
 * frequencies of each input must be hardcoded into the function.
 *
 * If using real data the files must be specified in the \c input-files command line argument - these should be comma
 * separated for multiple files and be in the same order at the associated detector from which they came given by the
 * \c detectors command. The data can be downsampled by a factor given by the \c downsample-factor command line
 * argument, but by default no down-sampling is applied. The downsampling is performed by averaging points.
 *
 * The function also check that valid Earth and Sun ephemeris files (from the lalpulsar suite) are set with the \c
 * ephem-earth and \c ephem-sun arguments, and that a valid output file for the nested samples is set via the \c outfile
 * argument.
 *
 * The \c log_factorial array is also filled in with values of the log of the factorial of all numbers up to the maximum
 * length of the data. This is used in likelihood calculations.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void readPulsarData( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL, *ppt2 = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;

  CHAR *detectors = NULL;
  CHAR *inputfile = NULL;

  CHAR *filestr = NULL;

  CHAR *efile = NULL, *sfile = NULL, *tfile = NULL;
  TimeCorrectionType ttype;

  CHAR *tempdets = NULL;
  CHAR *tempdet = NULL;

  REAL8 *fpsds = NULL;
  CHAR *fakestarts = NULL, *fakelengths = NULL, *fakedt = NULL;
  REAL8 *fstarts = NULL, *flengths = NULL, *fdt = NULL;

  CHAR dets[MAXDETS][256];
  INT4 numDets = 0, i = 0, numPsds = 0, numLengths = 0, numStarts = 0;
  INT4 numDt = 0, count = 0;
  UINT4 maxlen = 0;

  LALInferenceIFOData *ifodata = NULL;
  LALInferenceIFOData *prev = NULL;

  UINT4 seed = 0; /* seed for data generation */
  RandomParams *randomParams = NULL;

  CHAR *harmonics = NULL;
  REAL8Vector *modelFreqFactors = NULL;
  INT4 ml = 1, downs = 1;

  CHAR *parFile = NULL;
  BinaryPulsarParams pulsar;

  runState->data = NULL;

  /* check pulsar model required by getting the frequency harmonics */
  ppt = LALInferenceGetProcParamVal( commandLine, "--harmonics" );
  if ( ppt ) { harmonics = XLALStringDuplicate( ppt->value ); }
  else { harmonics = XLALStringDuplicate( "2" ); } /* default to using harmonic at twice the rotation rate */

  CHAR *tmpharms = NULL, *tmpharm = NULL, harmval[256];

  ml = count_csv( harmonics );
  modelFreqFactors = XLALCreateREAL8Vector( ml );

  tmpharms = XLALStringDuplicate( harmonics );
  for( i = 0; i < ml; i++ ){
    tmpharm = strsep( &tmpharms, "," );
    XLALStringCopy( harmval, tmpharm, strlen(tmpharm)+1 );

    modelFreqFactors->data[i] = atof(harmval);
  }

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--par-file" );
  if( ppt == NULL ) { fprintf(stderr,"Must specify --par-file!\n"); exit(1); }
  parFile = XLALStringDuplicate( ppt->value );

  /* get the pulsar parameters to give a value of f */
  XLALReadTEMPOParFile( &pulsar, parFile );

  /* get the detectors - must */
  ppt = LALInferenceGetProcParamVal( commandLine, "--detectors" );
  ppt2 = LALInferenceGetProcParamVal( commandLine, "--fake-data" );
  if( ppt && !ppt2 ){
    detectors = XLALStringDuplicate( ppt->value );

    /* count the number of detectors from command line argument of comma separated vales and set their names */
    tempdets = XLALStringDuplicate( detectors );

    if ( (numDets = count_csv( detectors )) > MAXDETS ){
      fprintf(stderr, "Error... too many detectors specified. Increase MAXDETS to be greater than %d if necessary.\n",
MAXDETS);
      exit(0);
    }

    for( i = 0; i < numDets; i++ ){
      tempdet = strsep( &tempdets, "," );
      XLALStringCopy( dets[i], tempdet, strlen(tempdet)+1 );
    }

  }
  /*Get psd values for generating fake data.*/
  /*=========================================================================*/
  /*if using fake data and no detectors are specified*/
  else if( ppt2 && !ppt ){
    detectors = XLALStringDuplicate( ppt2->value );

    fpsds = XLALCalloc( MAXDETS * ml, sizeof(REAL8) );

    if ( (numDets = count_csv( detectors )) > MAXDETS ){
      fprintf(stderr, "Error... too many detectors specified. Increase MAXDETS to be greater than %d if necessary.\n",
MAXDETS);
      exit(0);
    }
    /*------------------------------------------------------------------------*/
    /* get noise psds if specified */
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-psd");
    if( ppt ){
      CHAR *psds = NULL, *tmppsds = NULL, *tmppsd = NULL, psdval[256];

      psds = XLALStringDuplicate( ppt->value );

      tmppsds = XLALStringDuplicate( psds );
      tempdets = XLALStringDuplicate( detectors );

      /* count the number of PSDs (comma seperated values) to compare to number of detectors */
      if( (numPsds = count_csv( psds )) != ml*numDets ){
        fprintf(stderr, "Error... for %d harmonics the number of PSDs for fake data must be %d times the number of \
detectors specified (no. dets = %d)\n", ml, ml, numDets);
        exit(0);
      }

      for( i = 0; i < ml*numDets; i++ ){
        CHAR *tmpstr = NULL;

        tmppsd = strsep( &tmppsds, "," );
        XLALStringCopy( psdval, tmppsd, strlen(tmppsd)+1 );
        fpsds[i] = atof(psdval);

        /* set detector */
        if ( i%ml == 0 ){
          tempdet = strsep( &tempdets, "," );

          if( (tmpstr = strstr(tempdet, "A")) != NULL ){ /* have advanced */
            XLALStringCopy( dets[FACTOR(i,ml)], tmpstr+1, strlen(tmpstr)+1 );
          }
          else { XLALStringCopy( dets[FACTOR(i,ml)], tempdet, strlen(tempdet)+1 ); }
        }
      }
    }
    /*------------------------------------------------------------------------*/
    else{ /* get PSDs from model functions and set detectors */
      REAL8 pfreq = 0.;

      /* putting in pulsar frequency at f here */
      pfreq = pulsar.f0;

      tempdets = XLALStringDuplicate( detectors );

      for( i = 0; i < ml*numDets; i++ ){
        CHAR *tmpstr = NULL;
        REAL8 psdvalf = 0.;

        numPsds++;

        if( i%ml == 0 ) { tempdet = strsep( &tempdets, "," ); }

        if( (tmpstr = strstr(tempdet, "A")) != NULL ){ /* have Advanced */
          XLALStringCopy( dets[FACTOR(i,ml)], tmpstr+1, strlen(tmpstr)+1 );

          if( !strcmp(dets[FACTOR(i,ml)], "H1") || !strcmp(dets[FACTOR(i,ml)], "L1") ||
              !strcmp(dets[FACTOR(i,ml)], "H2") ){ /* ALIGO */
            psdvalf = XLALSimNoisePSDaLIGOZeroDetHighPower( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "V1") ){ /* AVirgo */
            psdvalf = XLALSimNoisePSDAdvVirgo( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else{
            fprintf(stderr, "Error... trying to use Advanced detector that is not available!\n");
            exit(0);
          }
        }
        else{ /* initial detector */
          XLALStringCopy( dets[FACTOR(i,ml)], tempdet, strlen(tempdet)+1 );

          if( !strcmp(dets[FACTOR(i,ml)], "H1") || !strcmp(dets[FACTOR(i,ml)], "L1") ||
              !strcmp(dets[FACTOR(i,ml)], "H2") ){ /* Initial LIGO */
            psdvalf = XLALSimNoisePSDiLIGOSRD( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );

            /* divide H2 psds by 2 */
            if( !strcmp(dets[FACTOR(i,ml)], "H2") ) { psdvalf /= 2.; }
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "V1") ){ /* Initial Virgo */
            psdvalf = XLALSimNoisePSDVirgo( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "G1") ){ /* GEO 600 */
            psdvalf = XLALSimNoisePSDGEO( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "T1") ){ /* TAMA300 */
            psdvalf = XLALSimNoisePSDTAMA( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else{
            fprintf(stderr, "Error... trying to use detector that is not available!\n");
            exit(0);
          }
        }

        fpsds[i] = psdvalf;
      }
    }
    /*generate the fake data timestamps.*/
    /*====================================================================*/

    fstarts = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-starts");
    if( ppt ){
      CHAR *tmpstarts = NULL, *tmpstart = NULL, startval[256];
      fakestarts = XLALStringDuplicate( ppt->value );

      if( (numStarts = count_csv( fakestarts )) != numDets*ml ){
        fprintf(stderr, "Error... for %d harmonics the number of start times for fake data must be %d times the number \
of detectors specified (no. dets = %d)\n", ml, ml, numDets);
        exit(0);
      }

      tmpstarts = XLALStringDuplicate( fakestarts );
      for( i = 0; i < ml*numDets; i++ ){
        tmpstart = strsep( &tmpstarts, "," );
        XLALStringCopy( startval, tmpstart, strlen(tmpstart)+1 );

        fstarts[i] = atof(startval);
      }
    }
    else{ /* set default GPS 900000000 - 13th July 2008 at 15:59:46 */
      for(i = 0; i < ml*numDets; i++ ){ fstarts[i] = 900000000.; }
    }

    flengths = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal( commandLine, "--fake-lengths" );
    if( ppt ){
      CHAR *tmplengths = NULL, *tmplength = NULL, lengthval[256];
      fakelengths = XLALStringDuplicate( ppt->value );

      if( (numLengths = count_csv( fakelengths )) != numDets*ml ){
        fprintf(stderr, "Error... for %d harmonics the number of data lengths for fake data must be %d times the \
number of detectors specified (no. dets = %d)\n", ml, ml, numDets);
        exit(0);
      }

      tmplengths = XLALStringDuplicate( fakelengths );
      for( i = 0; i < ml*numDets; i++ ){
        tmplength = strsep( &tmplengths, "," );
        XLALStringCopy( lengthval, tmplength, strlen(tmplength)+1 );
        flengths[i] = atof(lengthval);
      }
    }
    else{ /* set default (86400 seconds or 1 day) */
      for(i = 0; i < ml*numDets; i++ ) { flengths[i] = 86400.; }
    }

    fdt = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal( commandLine, "--fake-dt" );
    if( ppt ){
      CHAR *tmpdts = NULL, *tmpdt = NULL, dtval[256];
      fakedt = XLALStringDuplicate( ppt->value );

      if( (numDt = count_csv( fakedt )) != ml*numDets ){
        fprintf(stderr, "Error... for %d harmonics the number of sample time steps for fake data must be %d times the \
number of detectors specified (no. dets =\%d)\n", ml, ml, numDets);
        exit(0);
      }

      tmpdts = XLALStringDuplicate( fakedt );

      for( i = 0; i < ml*numDets; i++ ){
        tmpdt = strsep( &tmpdts, "," );
        XLALStringCopy( dtval, tmpdt, strlen(tmpdt)+1 );
        fdt[i] = atof(dtval);
      }
    }
    else{ /* set default (60 sesonds) */
      for(i = 0; i < ml*numDets; i++) { fdt[i] = 60.; }
    }

  }
  /*psds set and timestamps set.*/
  /*====================================================================*/
  else{
    fprintf(stderr, "Error... --detectors OR --fake-data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  ppt = LALInferenceGetProcParamVal( commandLine,"--input-files" );
  if( ppt ) { inputfile = XLALStringDuplicate( ppt->value ); }

  if ( inputfile == NULL && !ppt2 ){
    fprintf(stderr, "Error... an input file or fake data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  /* get the output directory */
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( !ppt ){
    fprintf(stderr, "Error... --outfile needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  /* count the number of input files (by counting commas) and check it's equal to twice the number of detectors */
  if ( !ppt2 ){ /* if using real data */
    count = count_csv( inputfile );

    if ( count != ml*numDets ){
      fprintf(stderr, "Error... for %d harmonics the number of input files given must be %d times the number of \
detectors specified (no. dets =\%d)\n", ml, ml, numDets);
      exit(0);
    }
  }

  /* set random number generator in case when that fake data is used */
  ppt = LALInferenceGetProcParamVal( commandLine, "--randomseed" );
  if ( ppt != NULL ) { seed = atoi( ppt->value ); }
  else { seed = 0; } /* will be set from system clock */

  /* check if we want to down-sample the input data by a given factor */
  ppt = LALInferenceGetProcParamVal( commandLine, "--downsample-factor" );
  if ( ppt != NULL ) { downs = atoi( ppt->value ); }
  else { downs = 1; } /* no downsampling */

  /* reset filestr if using real data (i.e. not fake) */
  if ( !ppt2 ) { filestr = XLALStringDuplicate( inputfile ); }

  /* read in data, needs to read in two sets of data for each ifo for pinsf model */
  for( i = 0, prev=NULL ; i < ml*numDets ; i++, prev=ifodata ){
    CHAR *datafile = NULL;
    REAL8 times = 0;
    LIGOTimeGPS gpstime;
    REAL8 dataValsRe = 0., dataValsIm = 0.;
    REAL8Vector *temptimes = NULL;
    INT4 j = 0, k = 0, datalength = 0;
    ProcessParamsTable *ppte = NULL, *ppts = NULL, *pptt = NULL;

    FILE *fp = NULL;

    count = 0;

    /* initialise random number generator */
    /* Moved into det loop so same random seed can be used with */
    /* different detector combos and still get same noise realisation */
    randomParams = XLALCreateRandomParams( seed+i );

    ifodata = XLALCalloc( 1, sizeof(LALInferenceIFOData) );
    ifodata->modelParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
    ifodata->modelDomain = LAL_SIM_DOMAIN_TIME;
    ifodata->next = NULL;
    ifodata->dataParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );

    /* add the pulsar model */
    /* LALInferenceAddVariable( ifodata->dataParams, "modeltype", &modeltype, LALINFERENCE_string_t,
                             LALINFERENCE_PARAM_FIXED ); */

    /* add data sample interval */
    /*LALInferenceAddVariable( ifodata->dataParams, "dt", &fdt[i],
                             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );*/

    /* add frequency factors variable */
    LALInferenceAddVariable( ifodata->dataParams, "freqfactors", &modelFreqFactors, LALINFERENCE_REAL8Vector_t,
                             LALINFERENCE_PARAM_FIXED );

    if( i == 0 ) { runState->data = ifodata; }
    if( i > 0 ) { prev->next = ifodata; }

    /* set detector */
    ifodata->detector = XLALGetSiteInfo( dets[FACTOR(i,ml)] );
    snprintf(ifodata->name, sizeof(char)*DETNAMELEN, "%s", dets[FACTOR(i,ml)]);

    /* set dummy initial time */
    gpstime.gpsSeconds = 0;
    gpstime.gpsNanoSeconds = 0;

    /* allocate time domain data - will be dynamically allocated as data read*/
    ifodata->compTimeData = NULL;
    ifodata->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 1 );

    /* allocate time domain model */
    ifodata->compModelData = NULL;
    ifodata->compModelData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0., 1., &lalSecondUnit, 1 );

    /*============================ GET DATA ==================================*/
    /* get i'th filename from the comma separated list */
    if ( !ppt2 ){ /* if using real data read in from the file */
      UINT4 counter = 0;

      datafile = strsep(&filestr, ",");

      /* open data file */
      if( (fp = fopen(datafile, "r")) == NULL ){
        fprintf(stderr, "Error... can't open data file %s!\n", datafile);
        exit(0);
      }

      j=0;

      /* read in data */
      temptimes = XLALCreateREAL8Vector( 1 );

      /* variables to aid downsampling of data */
      REAL8 tmpre = 0., tmpim = 0., timetmp = 0., dtcur = 0., dtprev = 0.;
      REAL8 tnow = 0., tprev = 0.;

      /* read in data - ignore lines starting with # or % */
      long offset;
      while(!feof(fp)){
        offset = ftell(fp);
        int testchar;

        if( (testchar = fgetc(fp)) == EOF ){ break; }
        else{ fseek(fp, offset, SEEK_SET); } /* return to start of line */

        /* test for comment characters */
        if( ( testchar == '%' ) || ( testchar == '#' ) ){
          /* skip to end of line */
          /* some lines for testing the comment check */
          // char *line = NULL;
          // size_t len = 0;
          // ssize_t testread;

          // if( (testread = getline(&line, &len, fp)) == -1 ){ break; }
          // fprintf(stderr, "%s", line);

          if ( fscanf(fp, "%*[^\n]") == EOF ){ break; }
        }
        else{ /* read in data */
          int rc = fscanf(fp, "%lf%lf%lf", &times, &dataValsRe, &dataValsIm);
          if( rc == EOF ){ break; }
          else if( rc != 3 ){ continue; } /* ignore the line */
        }

        j++;

        tnow = times;

        /* downsample by averaging if required */
        if ( j%downs != 0 ){
          if ( j > 1 ) { dtcur = tnow - tprev; }

          tmpre += dataValsRe;
          tmpim += dataValsIm;
          timetmp += times;
          tprev = tnow;

          /* if timesteps between points aren't equal then reset and move on i.e. we only want to use contiguous data
           * segments with equal time spacings */
          if ( j > 2 && dtcur != dtprev ){
            timetmp = times;
            tmpre = dataValsRe;
            tmpim = dataValsIm;
            dtcur = dtprev = 0.;
            j = 1;
          }

          dtprev = dtcur;

          continue;
        }
        else{
          /* if downsampling is not occuring just set individual values */
          if ( !tmpre && !tmpim && !timetmp ){
            tmpre = dataValsRe;
            tmpim = dataValsIm;
            timetmp = times;
          }
          else{ /* if downsampling get averages */
            if( j > 1 ) dtcur = tnow - tprev;

            /* check for contiguous segments */
            if( j > 2 && dtcur != dtprev ){
              timetmp = times;
              tmpre = dataValsRe;
              tmpim = dataValsIm;
              dtcur = dtprev = 0.;
              j = 1;
              continue;
            }

            /* add on final point and average */
            tmpre = (tmpre + dataValsRe) / (REAL8)downs;
            tmpim = (tmpim + dataValsIm ) / (REAL8)downs;
            timetmp = (timetmp + times) / (REAL8)downs;
          }
        }

        counter++;

        /* dynamically allocate more memory */
        ifodata->compTimeData = XLALResizeCOMPLEX16TimeSeries( ifodata->compTimeData, 0, counter );

        ifodata->compModelData = XLALResizeCOMPLEX16TimeSeries( ifodata->compModelData, 0, counter );

        temptimes = XLALResizeREAL8Vector( temptimes, counter );

        temptimes->data[counter-1] = timetmp;
        ifodata->compTimeData->data->data[counter-1] = tmpre + I*tmpim;

        tmpre = tmpim = timetmp = 0.;
        dtcur = dtprev = 0.;
        j = 0;
      }

      fclose(fp);

      datalength = counter;

      /* allocate data time stamps */
      ifodata->dataTimes = NULL;
      ifodata->dataTimes = XLALCreateTimestampVector( datalength );

      /* fill in time stamps as LIGO Time GPS Vector */
      REAL8 sampledt = INFINITY; /* sample interval */
      for ( k = 0; k < datalength; k++ ) {
        XLALGPSSetREAL8( &ifodata->dataTimes->data[k], temptimes->data[k] );

        if ( k > 0 ){
          /* get sample interval from the minimum time difference in the data */
          if ( temptimes->data[k] - temptimes->data[k-1] < sampledt ) {
            sampledt = temptimes->data[k] - temptimes->data[k-1];
          }
        }
      }

      ifodata->compTimeData->epoch = ifodata->dataTimes->data[0];
      ifodata->compModelData->epoch = ifodata->dataTimes->data[0];

      /* add data sample interval */
      ppt = LALInferenceGetProcParamVal( commandLine, "--sample-interval" );
      if( ppt ){
        sampledt = atof( ppt->value );
      }
      LALInferenceAddVariable( ifodata->dataParams, "dt", &sampledt, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

      XLALDestroyREAL8Vector( temptimes );
    }
    else{ /* set up fake data */
      datalength = flengths[i] / fdt[i];

      /* temporary real and imaginary data vectors */
      REAL4Vector *realdata = NULL;
      REAL4Vector *imagdata = NULL;

      REAL8 psdscale = 0.;

      /* allocate data time stamps */
      ifodata->dataTimes = NULL;
      ifodata->dataTimes = XLALCreateTimestampVector( (UINT4)datalength );

      /* add data sample interval */
      LALInferenceAddVariable( ifodata->dataParams, "dt", &fdt[0], LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

      /* resize the data and model times series */
      ifodata->compTimeData = XLALResizeCOMPLEX16TimeSeries( ifodata->compTimeData, 0, datalength );

      ifodata->compModelData = XLALResizeCOMPLEX16TimeSeries( ifodata->compModelData, 0, datalength );

      /* create data drawn from normal distribution with zero mean and unit variance */
      realdata = XLALCreateREAL4Vector( datalength );
      imagdata = XLALCreateREAL4Vector( datalength );

      XLALNormalDeviates( realdata, randomParams );
      XLALNormalDeviates( imagdata, randomParams );

      /* converts single sided psd into double sided psd, and then into a time domain noise standard deviation */
      psdscale = sqrt( ( fpsds[i] / 2.) / ( 2. * fdt[i] ) );

      /* create time stamps and scale data with the PSD */
      for( k = 0; k < datalength; k++ ){
        /* set time stamp */
        XLALGPSSetREAL8( &ifodata->dataTimes->data[k], fstarts[i] + fdt[i] * (REAL8)k );

        ifodata->compTimeData->data->data[k] = (REAL8)realdata->data[k] * psdscale +
          I * (REAL8)imagdata->data[k] * psdscale;
      }

      ifodata->compTimeData->epoch = ifodata->dataTimes->data[0];
      ifodata->compModelData->epoch = ifodata->dataTimes->data[0];

      XLALDestroyREAL4Vector( realdata );
      XLALDestroyREAL4Vector( imagdata );
    }

    /* set ephemeris data */
    ifodata->ephem = XLALMalloc( sizeof(EphemerisData) );

    /* get ephemeris files */
    ppte = LALInferenceGetProcParamVal( commandLine, "--ephem-earth" );
    ppts = LALInferenceGetProcParamVal( commandLine, "--ephem-sun" );
    pptt = LALInferenceGetProcParamVal( commandLine, "--ephem-timecorr" );
    if( ppte && ppts ){
      efile = XLALStringDuplicate( ppte->value );
      sfile = XLALStringDuplicate( ppts->value );

      if ( pptt ){
        tfile = XLALStringDuplicate( pptt->value );

        if ( pulsar.units != NULL ){
          if( !strcmp(pulsar.units, "TDB") ) { ttype = TIMECORRECTION_TDB; }
          else { ttype = TIMECORRECTION_TCB; } /* default to TCB otherwise */
        }
        else { ttype = TIMECORRECTION_TCB; }
      }
      else{
        tfile = NULL;
        ttype = TIMECORRECTION_ORIGINAL;
      }
    }
    else{ /* try getting files automatically */
      CHAR efiletmp[1024], sfiletmp[1024], tfiletmp[1024];

      if( !( ttype = XLALAutoSetEphemerisFiles( efiletmp, sfiletmp, tfiletmp, pulsar,
            ifodata->dataTimes->data[0].gpsSeconds, ifodata->dataTimes->data[datalength-1].gpsSeconds ) ) ){
        fprintf(stderr, "Error... not been able to set ephemeris files!\n");
        exit(3);
      }

      efile = XLALStringDuplicate(efiletmp);
      sfile = XLALStringDuplicate(sfiletmp);
      tfile = XLALStringDuplicate(tfiletmp);
    }

    /* check ephemeris files exist and if not output an error message */
    if( fopen(sfile, "r") == NULL || fopen(efile, "r") == NULL ){
      fprintf(stderr, "Error... ephemeris files not, or incorrectly, defined!\n");
      exit(3);
    }

    /* set up ephemeris information */
    XLAL_CHECK_VOID( (ifodata->ephem = XLALInitBarycenter( efile, sfile ) ) != NULL, XLAL_EFUNC );
    if( tfile ){ XLAL_CHECK_VOID( (ifodata->tdat = XLALInitTimeCorrections( tfile ) ) != NULL, XLAL_EFUNC ); }
    else { ifodata->tdat = NULL; }
    ifodata->ttype = ttype;

    XLALDestroyRandomParams( randomParams );

    /* get maximum data length */
    if ( ifodata->compTimeData->data->length > maxlen ) { maxlen = ifodata->compTimeData->data->length; }
  }

  /* set global variable logfactorial */
  logfactorial = XLALCalloc( maxlen+1, sizeof(REAL8) );
  for ( i = 2; i < (INT4)(maxlen+1); i++ ) { logfactorial[i] = logfactorial[i-1] + log((REAL8)i); }

  /* free memory */
  XLALFree( fdt );
  XLALFree( flengths );
  XLALFree( fstarts );
  XLALFree( fpsds );
}


/** \brief Reads in the parameters of the pulsar being searched for
 *
 * This function reads in a pulsars parameters from the specified TEMPO-style .par file given by \c par-file using \c
 * XLALReadTEMPOParFile. This file must be specified and should contain at least the pulsars frequency, right
 * ascension and declination (any value not included will be zero by default). The file should contain the parameters
 * with which the detector data was heterodyned, as these are used to produce a signal phase template based on this
 * assumption.
 *
 * A example .par file may look like
 * \code
 * RA 12:31:56.17643
 * DEC 43:21:35.2531
 * F0 100.78634 1 0.00005
 * F1 2.34e-15
 * PEPOCH 54323.785634
 * \endcode
 * which shows several parameters mostly defined by the parameter name and a parameter value. However, the \c F0 value
 * contains 4 items. If a parameter has a \c 1 as the third entry then it means that this was a parameter that was fit
 * by TEMPO with the entry after the \c 1 being the 1 standard deviation error on that parameter. For parameters where
 * an error is present the code will attempt to search over that parameter using a Gaussian prior defined by the
 * 1\f$\sigma\f$ error value. Other parameters will be set as fixed by default. These can be overridden by the prior
 * file values described in \c initialisePrior().
 *
 * Based on the defined sky position defined in the par file a lookup table of the detector antenna response over time
 * and polarisation will be set by \c setupLookupTables().
 *
 * The function \c add_initial_variables() is used to pass the parameter values from the .par file to the algorithm.
 *
 * Using the parameters from the .par file the phase template, including the solar system and binary system barycentring
 * time delays will be setup. These define the phase template used to perform the initial heterodyne, which is used as
 * the reference in cases when phase parameters (other than the initial phase) are being searched over.
 *
 * Values used for scaling the parameters (to avoid dynamic range issues) are initialised although will be set as
 * default values.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 *
 * \sa setupLookupTables
 * \sa add_initial_variables
 * \sa get_phase_model
 * \sa add_correlation_matrix
 */
void setupFromParFile( LALInferenceRunState *runState )
/* Read the PAR file of pulsar parameters and setup the code using them */
/* Generates lookup tables also */
{
  LALSource psr;
  BinaryPulsarParams pulsar;
  REAL8Vector *phase_vector = NULL;
  LALInferenceIFOData *data = runState->data;
  LALInferenceVariables *scaletemp;
  ProcessParamsTable *ppt = NULL;
  UINT4 mmfactor = 0;
  REAL8 mm = 0;

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--par-file" );
  if( ppt == NULL ) { fprintf(stderr,"Must specify --par-file!\n"); exit(1); }
  CHAR *parFile = ppt->value;

  /* check if we needed a downsampled time stamp series */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--mm-factor" );
  if( ppt != NULL ) mmfactor = atoi( ppt->value );

  /* get mismatch value if required */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--mismatch" );
  if( ppt != NULL ) mm = atof( ppt->value );

  /* get the pulsar parameters */
  XLALReadTEMPOParFile( &pulsar, parFile );
  psr.equatorialCoords.longitude = pulsar.ra;
  psr.equatorialCoords.latitude = pulsar.dec;
  psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* Setup lookup tables for amplitudes */
  setupLookupTables( runState, &psr );

  runState->currentParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );

  scaletemp = XLALCalloc( 1, sizeof(LALInferenceVariables) );

  /* Add initial (unchanging) variables for the model, initial (unity) scale factors, from the par file */
  add_initial_variables( runState->currentParams, scaletemp, pulsar );

  /* Setup initial phase, and barycentring delays */
  while( data ){
    REAL8Vector *freqFactors = NULL;
    UINT4 j = 0;

    freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams,
                                                            "freqfactors" );

    for( j = 0; j < freqFactors->length; j++ ){
      UINT4 i = 0;
      LALInferenceVariableItem *scaleitem = scaletemp->head;
      REAL8Vector *dts = NULL, *bdts = NULL;

      dts = get_ssb_delay( pulsar, data->dataTimes, data->ephem, data->tdat, data->ttype, data->detector, 0. );

      LALInferenceAddVariable( data->dataParams, "ssb_delays", &dts, LALINFERENCE_REAL8Vector_t,
                               LALINFERENCE_PARAM_FIXED );

      bdts = get_bsb_delay( pulsar, data->dataTimes, dts, data->ephem );

      LALInferenceAddVariable( data->dataParams, "bsb_delays", &bdts, LALINFERENCE_REAL8Vector_t,
                               LALINFERENCE_PARAM_FIXED );

      phase_vector = get_phase_model( pulsar, data, freqFactors->data[j], 0 );

      data->timeData = NULL;
      data->timeData = XLALCreateREAL8TimeSeries( "", &data->dataTimes->data[0], 0., 1., &lalSecondUnit,
                                                  phase_vector->length );

      for ( i=0; i<phase_vector->length; i++ ) { data->timeData->data->data[i] = phase_vector->data[i]; }

      /* add the scale factors from scaletemp into the data->dataParams structure */
      for( ; scaleitem; scaleitem = scaleitem->next ){
        LALInferenceAddVariable( data->dataParams, scaleitem->name, scaleitem->value, scaleitem->type,
                                 scaleitem->vary );
      }

      /* get down sampled time stamps if required and set mismatch */
      if ( mmfactor != 0 && mm != 0. ){
        LIGOTimeGPSVector *downst = XLALCreateTimestampVector( floor(phase_vector->length/mmfactor) );
        UINT4 k = 0;

        /* array to contain down-sampled phases */
        REAL8Vector *dsphase = XLALCreateREAL8Vector( floor(phase_vector->length/mmfactor) );

        if ( downst->length < 2 ){
          XLALPrintError("Error, downsampled time stamp factor to high!\n");
          XLAL_ERROR_VOID(XLAL_EFAILED);
        }

        for( k = 1; k < downst->length+1; k++ ){
          downst->data[k-1] = data->dataTimes->data[(k-1)*mmfactor];
          dsphase->data[k-1] = 0.;
        }

        LALInferenceAddVariable( data->dataParams, "downsampled_times", &downst, LALINFERENCE_void_ptr_t,
                                 LALINFERENCE_PARAM_FIXED );

        LALInferenceAddVariable( data->dataParams, "ds_phase", &dsphase, LALINFERENCE_REAL8Vector_t,
                                 LALINFERENCE_PARAM_FIXED );

        LALInferenceAddVariable( data->dataParams, "mismatch", &mm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
      }

      data = data->next;
    }
  }

  return;
}


/** \brief Sets the time vs polarisation angle antenna response lookup table
 *
 * This function sets up a lookup table in time vs. polarisation angle \f$\psi\f$ for each detector from which data
 * exists (either real or fake). The time ranges over one sidereal day and the polarisation angle range from
 * \f$-\pi/4\f$ to \f$\pi/4\f$. The number of bins for the grid over these two parameters can be specified on the
 * command line via \c time-bins and \c psi-bins respectively, but if these are not given then default values are used.
 * The data times as a fraction of a sidereal day from the start time will also be calculated.
 *
 * The function will by default also call \c chop_n_merge() for each data set, which will split the data into chunks
 * during which it can be considered Gaussian and stationary. The command line arguments \c chunk-min and \c chunk-max
 * can be used to specify hardwire minimum and maximum lengths of chunks that are allowable. By default the maximum
 * chunk length is 0, which corresponds to no maximum value being set. If the \c --oldChunks flag is set then data will
 * be split as in the older version of the parameter estimation code, in which the chunk length is fixed, except for the
 * possibility of shorter segments at the end of science segments.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 * \param source [in] A pointer to a LALSource variable containing the source location
 *
 * \sa chop_n_merge
 *
 */
void setupLookupTables( LALInferenceRunState *runState, LALSource *source ){
  /* Set up lookup tables */
  /* Using psi bins, time bins */
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;

  INT4 chunkMin, chunkMax;

  /* Get chunk min and chunk max */
  ppt = LALInferenceGetProcParamVal( commandLine, "--chunk-min" );
  if( ppt ) { chunkMin = atoi( ppt->value ); }
  else { chunkMin = CHUNKMIN; } /* default minimum chunk length */

  ppt = LALInferenceGetProcParamVal( commandLine, "--chunk-max" );
  if( ppt ) { chunkMax = atoi( ppt->value ); }
  else { chunkMax = CHUNKMAX; } /* default maximum chunk length */

  LALInferenceIFOData *data = runState->data;

  gsl_matrix *LUfplus = NULL;
  gsl_matrix *LUfcross = NULL;

  REAL8 t0;
  LALDetAndSource detAndSource;

  ppt = LALInferenceGetProcParamVal( commandLine, "--psi-bins" );
  INT4 psiBins;
  if( ppt ) { psiBins = atoi( ppt->value ); }
  else { psiBins = PSIBINS; } /* default psi bins */

  ppt = LALInferenceGetProcParamVal( commandLine, "--time-bins" );
  INT4 timeBins;
  if( ppt ) { timeBins = atoi( ppt->value ); }
  else { timeBins = TIMEBINS; } /* default time bins */

  while(data){
    UINT4Vector *chunkLength = NULL;
    REAL8Vector *sidDayFrac = NULL;
    UINT4 i = 0;

    LALInferenceAddVariable( data->dataParams, "psiSteps", &psiBins, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( data->dataParams, "timeSteps", &timeBins, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );

    t0 = XLALGPSGetREAL8( &data->dataTimes->data[0] );

    sidDayFrac = XLALCreateREAL8Vector( data->dataTimes->length );

    /* set the time in sidereal days since the first data point (mod 1 sidereal
       day) */
    for( i = 0; i < data->dataTimes->length; i++ ){
      sidDayFrac->data[i] = fmod( XLALGPSGetREAL8(&data->dataTimes->data[i]) - t0, LAL_DAYSID_SI );
    }

    LALInferenceAddVariable( data->dataParams, "siderealDay", &sidDayFrac, LALINFERENCE_REAL8Vector_t,
                             LALINFERENCE_PARAM_FIXED );

    detAndSource.pDetector = data->detector;
    detAndSource.pSource = source;

    LUfplus = gsl_matrix_alloc( psiBins, timeBins );

    LUfcross = gsl_matrix_alloc( psiBins, timeBins );

    response_lookup_table( t0, detAndSource, timeBins, psiBins, LUfplus, LUfcross );

    LALInferenceAddVariable( data->dataParams, "LU_Fplus", &LUfplus, LALINFERENCE_gslMatrix_t,
                             LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( data->dataParams, "LU_Fcross", &LUfcross, LALINFERENCE_gslMatrix_t,
                             LALINFERENCE_PARAM_FIXED );

    LALInferenceAddVariable( data->dataParams, "chunkMin", &chunkMin, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( data->dataParams, "chunkMax", &chunkMax, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );

    ppt = LALInferenceGetProcParamVal( commandLine, "--oldChunks" );
    if ( ppt ){ /* use old style quasi-fixed data chunk lengths */
      /* if a chunk max wasn't set use 30 mins by default */
      if ( !LALInferenceGetProcParamVal( commandLine, "--chunk-max" ) ){
        chunkMax = 30;
        LALInferenceSetVariable( data->dataParams, "chunkMax", &chunkMax );
      }

      chunkLength = get_chunk_lengths( data, chunkMax );
    }
    /* use new change points analysis to get chunks */
    else { chunkLength = chop_n_merge( data, chunkMin, chunkMax ); }

    LALInferenceAddVariable( data->dataParams, "chunkLength", &chunkLength, LALINFERENCE_UINT4Vector_t,
                             LALINFERENCE_PARAM_FIXED );

    data = data->next;
  }

  return;
}


/** \brief Set up all the allowed variables for a known pulsar search
 *
 * This functions sets up all possible variables that are possible in a known pulsar search. Parameter values read in
 * from a .par file and passed in via the \c pars variable will be set. Scale factors will be initialised for all
 * variables (so that they exist) although they will be set to 1.
 *
 * \param ini [in] A pointer to a \c LALInferenceVariables type that will be filled in with pulsar parameters
 * \param scaleFac [in] A pointer to a \c LALInferenceVariables type that will be initialised to hold scale factors for
 * each corresponding pulsar parameter
 * \param pars [in] A \c BinaryPulsarParams type containing pulsar parameters read in from a TEMPO-style .par file
 *
 * \sa add_variable_scale_prior
 */
void add_initial_variables( LALInferenceVariables *ini,  LALInferenceVariables *scaleFac, BinaryPulsarParams pars ){
  /* include a scale factor of 1 scaling values and if the parameter file contains an uncertainty then set the prior to
   * be Gaussian with the uncertainty as the standard deviation */

  /* amplitude model parameters for l=m=2 harmonic emission */
  add_variable_scale( ini, scaleFac, "h0", pars.h0 );
  add_variable_scale( ini, scaleFac, "phi0", pars.phi0 );
  add_variable_scale( ini, scaleFac, "cosiota", pars.cosiota );
  add_variable_scale( ini, scaleFac, "psi", pars.psi );

  /* amplitude model parameters for l=2, m=1 and 2 harmonic emission */
  /* add_variable_scale( ini, scaleFac, "I21", pars.I21 );
  add_variable_scale( ini, scaleFac, "I31", pars.I31 );
  add_variable_scale( ini, scaleFac, "lambda", pars.lambda );
  add_variable_scale( ini, scaleFac, "costheta", pars.costheta ); */

  /* amplitude model parameters in phase and amplitude form */
  add_variable_scale( ini, scaleFac, "C22", pars.C22 );
  add_variable_scale( ini, scaleFac, "C21", pars.C21 );
  add_variable_scale( ini, scaleFac, "phi22", pars.phi22 );
  add_variable_scale( ini, scaleFac, "phi21", pars.phi21 );

  /***** phase model parameters ******/
  /* frequency */
  add_variable_scale( ini, scaleFac, "f0", pars.f0 );
  add_variable_scale( ini, scaleFac, "f1", pars.f1 );
  add_variable_scale( ini, scaleFac, "f2", pars.f2 );
  add_variable_scale( ini, scaleFac, "f3", pars.f3 );
  add_variable_scale( ini, scaleFac, "f4", pars.f4 );
  add_variable_scale( ini, scaleFac, "f5", pars.f5 );
  add_variable_scale( ini, scaleFac, "pepoch", pars.pepoch );

  /* sky position */
  add_variable_scale( ini, scaleFac, "ra", pars.ra );
  add_variable_scale( ini, scaleFac, "pmra", pars.pmra );
  add_variable_scale( ini, scaleFac, "dec", pars.dec );
  add_variable_scale( ini, scaleFac, "pmdec", pars.pmdec );
  add_variable_scale( ini, scaleFac, "posepoch", pars.posepoch );

  /* only add binary system parameters if required */
  if ( pars.model ){
    LALInferenceAddVariable( ini, "model", &pars.model, LALINFERENCE_string_t, LALINFERENCE_PARAM_FIXED );

    add_variable_scale( ini, scaleFac, "Pb", pars.Pb );
    add_variable_scale( ini, scaleFac, "e", pars.e );
    add_variable_scale( ini, scaleFac, "eps1", pars.eps1 );
    add_variable_scale( ini, scaleFac, "eps2", pars.eps2 );
    add_variable_scale( ini, scaleFac, "T0", pars.T0 );
    add_variable_scale( ini, scaleFac, "Tasc", pars.Tasc );
    add_variable_scale( ini, scaleFac, "x", pars.x );
    add_variable_scale( ini, scaleFac, "w0", pars.w0 );

    add_variable_scale( ini, scaleFac, "Pb2", pars.Pb2 );
    add_variable_scale( ini, scaleFac, "e2", pars.e2 );
    add_variable_scale( ini, scaleFac, "T02", pars.T02 );
    add_variable_scale( ini, scaleFac, "x2", pars.x2 );
    add_variable_scale( ini, scaleFac, "w02", pars.w02 );

    add_variable_scale( ini, scaleFac, "Pb3", pars.Pb3 );
    add_variable_scale( ini, scaleFac, "e3", pars.e3 );
    add_variable_scale( ini, scaleFac, "T03", pars.T03 );
    add_variable_scale( ini, scaleFac, "x3", pars.x3 );
    add_variable_scale( ini, scaleFac, "w03", pars.w03 );

    add_variable_scale( ini, scaleFac, "xpbdot", pars.xpbdot );
    add_variable_scale( ini, scaleFac, "eps1dot", pars.eps1dot );
    add_variable_scale( ini, scaleFac, "eps2dot", pars.eps2dot );
    add_variable_scale( ini, scaleFac, "wdot", pars.wdot );
    add_variable_scale( ini, scaleFac, "gamma", pars.gamma );
    add_variable_scale( ini, scaleFac, "Pbdot", pars.Pbdot );
    add_variable_scale( ini, scaleFac, "xdot", pars.xdot );
    add_variable_scale( ini, scaleFac, "edot", pars.edot );

    add_variable_scale( ini, scaleFac, "s", pars.s );
    add_variable_scale( ini, scaleFac, "dr", pars.dr );
    add_variable_scale( ini, scaleFac, "dth", pars.dth );
    add_variable_scale( ini, scaleFac, "a0", pars.a0 );
    add_variable_scale( ini, scaleFac, "b0", pars.b0 );
    add_variable_scale( ini, scaleFac, "M", pars.M );
    add_variable_scale( ini, scaleFac, "m2", pars.m2 );
  }
}


/** \brief Adds variables, scale factors and priors
 *
 * This function adds a variable with a name and a value. For all parameters a scale factor and scale minimum range will
 * be set. These are just initialised to 1 and 0 respectively and will be set in \c initialisePrior for any parameters
 * that require them.
 *
 * \param var [in] Pointer to \c LALInferenceVariables type to contain parameter information
 * \param scale [in] Pointer to \c LALInferenceVariables type to contain parameter scaling information
 * \param name [in] string containing the parameter name
 * \param value [in] the value of the parameter
 */
void add_variable_scale( LALInferenceVariables *var, LALInferenceVariables *scale, const CHAR *name, REAL8 value ){
  REAL8 scaleVal = 1., scaleMin = 0.;
  CHAR scaleName[VARNAME_MAX] = "", scaleMinName[VARNAME_MAX] = "";

  /* add the variable */
  LALInferenceAddVariable( var, name, &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  /* add the initial scale factor of 1 */
  sprintf( scaleName, "%s_scale", name );
  LALInferenceAddVariable( scale, scaleName, &scaleVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  /* add initial scale offset of zero */
  sprintf( scaleMinName, "%s_scale_min", name );
  LALInferenceAddVariable( scale, scaleMinName, &scaleMin, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
}


/** \brief Sets up the parameters to be searched over and their prior ranges
 *
 * This function sets up any parameters that you require the code to search over and specifies the prior range and type
 * for each. This information is contained in a prior file specified by the command line argument \c prior-file. This
 * file should contain four columns: the first has the name of a parameter to be searched over; the second has the prior
 * type (e.g. "uniform" for a prior that is flat over the given range or "gaussian"); the third has the lower limit, or
 * mean, of the prior, for "uniform" and "gaussian" priors respectively; and the fourth has the upper limit, or standard
 * deviation, for "uniform" and "gaussian" priors respectively. E.g.
 * \code
 h0 uniform 0 1e-21
 phi0 uniform 0 6.6.283185307179586
 cosiota uniform -1 1
 psi uniform -0.785398163397448 0.785398163397448
 \endcode
 *
 * Any parameter specified in the file will have its vary type set to \c LALINFERENCE_PARAM_LINEAR, (except
 * \f$\phi_0\f$, which if it is defined to have a prior covering \f$2\pi\f$ wraps around at the edges of its range and
 * has a \c LALINFERENCE_PARAM_CIRCULAR type). Parameters, and their priors, with linear variable type are scaled such
 * that parameter \f$x\f$ with priors in the range \f$[a, b]\f$ will become \f$(x - a) / (b - a)\f$. As such the new
 * prior ranges will cover from 0 to 1. The parameter scale factor is set the to value of \f$(b - a)\f$ and the
 * minimum scale range is set to \f$a\f$  - this allows the true parameter value to be reconstructed. For parameters
 * with Gaussian priors the scale factor is applied differently, so as to give a Gaussian with zero mean and unit
 * variance.
 *
 * If a parameter correlation matrix is given by the \c cor-file command then this is used to construct a multi-variate
 * Gaussian prior for the given parameters (it is assumed that this file is created using TEMPO and the parameters it
 * contains are the same as those for which a standard deviation is defined in the par file). This overrules the
 * Gaussian priors that will have been set for these parameters. Due to the scalings applied to the parameters this
 * correlation coefficient matrix does not have to be converted into the true covariance matrix for use when calculating
 * the prior. Note that these multi-variate Gaussian priors will not overrule values given in the proposal file.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void initialisePrior( LALInferenceRunState *runState )
{
  CHAR *propfile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  FILE *fp=NULL;

  CHAR tempPar[VARNAME_MAX] = "", tempPrior[VARNAME_MAX] = "";
  REAL8 low, high;

  LALInferenceIFOData *data = runState->data;

  /* parameters in correlation matrix */
  LALStringVector *corParams = NULL;
  REAL8Array *corMat = NULL;

  INT4 phidef = 0; /* check if phi0 is in prior file */

  ppt = LALInferenceGetProcParamVal( commandLine, "--prior-file" );
  if( ppt ) { propfile = XLALStringDuplicate( LALInferenceGetProcParamVal(commandLine,"--prior-file")->value ); }
  else{
    fprintf(stderr, "Error... --prior-file is required.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }

  /* open file */
  if( (fp = fopen(propfile, "r")) == NULL ){
    fprintf(stderr, "Error... Could not open prior file %s.\n", propfile);
    exit(3);
  }

  while(fscanf(fp, "%s %s %lf %lf", tempPar, tempPrior, &low, &high) != EOF){
    REAL8 tempVar;
    LALInferenceVariableType type;
    INT4 isthere = 0, i = 0;

    REAL8 scale = 0., scaleMin = 0.;
    LALInferenceVariableType scaleType;
    CHAR tempParScale[VARNAME_MAX] = "";
    CHAR tempParScaleMin[VARNAME_MAX] = "";
    CHAR tempParPrior[VARNAME_MAX] = "";

    LALInferenceIFOData *datatemp = data;

    LALInferenceParamVaryType varyType;

    phidef = 0.;

    if( high < low ){
      fprintf(stderr, "Error... In %s the %s parameters ranges are wrongly set.\n", propfile, tempPar);
      exit(3);
    }

    sprintf(tempParScale, "%s_scale", tempPar);
    sprintf(tempParScaleMin, "%s_scale_min", tempPar);
    sprintf(tempParPrior, "%s_gaussian_mean", tempPar);

    tempVar = *(REAL8*)LALInferenceGetVariable( runState->currentParams, tempPar );
    type = LALInferenceGetVariableType( runState->currentParams, tempPar );

    /* remove variable value */
    LALInferenceRemoveVariable( runState->currentParams, tempPar );

    if ( !strcmp(tempPrior, "uniform") ){
      scale = high - low; /* the prior range */
      scaleMin = low;     /* the lower limit of the prior range */
    }
    else if( !strcmp(tempPrior, "gaussian") ){
      scale = high;   /* the standard deviation of the Gaussian prior */
      scaleMin = low; /* the mean of the Gaussian prior */
    }
    else{
      fprintf(stderr, "Error... prior type '%s' not recognised\n", tempPrior);
      exit(3);
    }

    /* if phi0 is defined and covers the 2pi range set phidef so re-parameterisation can take place  */
    if( ( !strcmp(tempPar, "phi0") || !strcmp(tempPar, "phi22") || !strcmp(tempPar, "phi21") )
      && !strcmp(tempPrior, "uniform")  ){
      if ( scale/LAL_TWOPI > 0.99 && scale/LAL_TWOPI < 1.01 ) { phidef = 1; }
    }

    /* if psi is covering the range -pi/4 to pi/4, i.e. as used in the triaxial
     * model, set psidef, so re-parameterisation can take place */
    //if( !strcmp(tempPar, "psi") )
    //  if ( scale/LAL_PI_2 > 0.99 && scale/LAL_PI_2 < 1.01 ) psidef = 1;

    /* set the scale factor to be the width of the prior */
    while( datatemp ){
      scaleType = LALInferenceGetVariableType( datatemp->dataParams, tempParScale );
      LALInferenceRemoveVariable( datatemp->dataParams, tempParScale );
      LALInferenceRemoveVariable( datatemp->dataParams, tempParScaleMin );

      LALInferenceAddVariable( datatemp->dataParams, tempParScale, &scale, scaleType, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( datatemp->dataParams, tempParScaleMin, &scaleMin, scaleType, LALINFERENCE_PARAM_FIXED );

      datatemp = datatemp->next;
    }

    /* scale variable and priors */
    tempVar = (tempVar - scaleMin) / scale;
    low = 0.;
    high = (high - scaleMin) / scale;

    /* re-add variable */
    if( !strcmp(tempPar, "phi0") && phidef ) { varyType = LALINFERENCE_PARAM_CIRCULAR; }
    else { varyType = LALINFERENCE_PARAM_LINEAR; }

    LALInferenceAddVariable( runState->currentParams, tempPar, &tempVar, type, varyType );

    /* Add the prior variables */
    if ( !strcmp(tempPrior, "uniform") ){
      LALInferenceAddMinMaxPrior( runState->priorArgs, tempPar, &low, &high, type );
    }
    else if( !strcmp(tempPrior, "gaussian") ){
      LALInferenceAddGaussianPrior( runState->priorArgs, tempPar, &low, &high, LALINFERENCE_REAL8_t );
    }

    /* if there is a phase parameter defined in the proposal then set varyphase to 1 */
    for ( i = 0; i < NUMAMPPARS; i++ ){
      if ( !strcmp(tempPar, amppars[i]) ){
        isthere = 1;
        break;
      }
    }
    if ( !isthere ) { varyphase = 1; }

    /* check if there are sky position parameters that will be searched over */
    for ( i = 0; i < NUMSKYPARS; i++ ){
      if ( !strcmp(tempPar, skypars[i]) ){
        varyskypos = 1;
        break;
      }
    }

    /* check if there are any binary parameters that will be searched over */
    for ( i = 0; i < NUMBINPARS; i++ ){
      if ( !strcmp(tempPar, binpars[i]) ){
        varybinary = 1;
        break;
      }
    }
  }

  /* If source physical parameters (rather than amplitude/phase parameters) have been defined in the prior file then
   * convert to amplitude/phase parameters. Only do this for the l=m=2 model of emission at 2f for a triaxial star.
   * For other modes the phase/amplitude parameters most be defined in the prior file. */
  REAL8Vector *freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "freqfactors" );

  if ( freqFactors->length == 1 && freqFactors->data[0] == 2. ){
    if ( LALInferenceGetVariableVaryType( runState->currentParams, "h0" ) == LALINFERENCE_PARAM_LINEAR &&
         LALInferenceGetVariableVaryType( runState->currentParams, "C22" ) == LALINFERENCE_PARAM_FIXED ){
      REAL8 h0 = *(REAL8*)LALInferenceGetVariable( runState->currentParams, "h0" );
      REAL8 h0scale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "h0_scale" );
      REAL8 h0min = *(REAL8*)LALInferenceGetVariable( data->dataParams, "h0_scale_min" );
      REAL8 C22 = h0 / 2.;
      REAL8 C22scale = h0scale / 2.;
      REAL8 C22scalemin = h0min / 2.;
      low = 0., high = 1.;

      LALInferenceRemoveVariable( runState->currentParams, "C22" );
      LALInferenceAddVariable( runState->currentParams, "C22", &C22, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR );

      LALInferenceIFOData *datatemp = data;
      while( datatemp ){
        LALInferenceRemoveVariable( datatemp->dataParams, "C22_scale" );
        LALInferenceRemoveVariable( datatemp->dataParams, "C22_scale_min" );

        LALInferenceAddVariable( datatemp->dataParams, "C22_scale", &C22scale, LALINFERENCE_REAL8_t,
                                 LALINFERENCE_PARAM_FIXED );
        LALInferenceAddVariable( datatemp->dataParams, "C22_scale_min", &C22scalemin, LALINFERENCE_REAL8_t,
                                 LALINFERENCE_PARAM_FIXED );

        datatemp = datatemp->next;
      }

      if ( LALInferenceCheckGaussianPrior( runState->priorArgs, "h0" ) ){
        LALInferenceRemoveGaussianPrior( runState->priorArgs, "h0" );
        LALInferenceAddGaussianPrior( runState->priorArgs, "C22", &low, &high, LALINFERENCE_REAL8_t );
      }
      else{
        LALInferenceRemoveMinMaxPrior( runState->priorArgs, "h0" );
        LALInferenceAddMinMaxPrior( runState->priorArgs, "C22", &low, &high, LALINFERENCE_REAL8_t );
      }
    }

    if ( LALInferenceGetVariableVaryType( runState->currentParams, "phi0" ) != LALINFERENCE_PARAM_FIXED &&
         LALInferenceGetVariableVaryType( runState->currentParams, "phi22" ) == LALINFERENCE_PARAM_FIXED ){
      REAL8 phi0 = *(REAL8*)LALInferenceGetVariable( runState->currentParams, "phi0" );
      REAL8 phi0scale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "phi0_scale" );
      REAL8 phi0min = *(REAL8*)LALInferenceGetVariable( data->dataParams, "phi0_scale_min" );
      REAL8 phi22;
      REAL8 phi22scale = phi0scale;
      REAL8 phi22scalemin = phi0min;
      low = 0., high = 1.;

      phi22 = phi0*phi0scale + phi0min - LAL_PI;
      phi22 = phi22 - LAL_TWOPI*floor(phi22/LAL_TWOPI); /* modulo 2*pi */
      phi22 /= phi22scale;

      LALInferenceParamVaryType vary = LALInferenceGetVariableVaryType( runState->currentParams, "phi0" );

      LALInferenceRemoveVariable( runState->currentParams, "phi22" );
      LALInferenceAddVariable( runState->currentParams, "phi22", &phi22, LALINFERENCE_REAL8_t, vary );

      LALInferenceIFOData *datatemp = data;
      while( datatemp ){
        LALInferenceRemoveVariable( datatemp->dataParams, "phi22_scale" );
        LALInferenceRemoveVariable( datatemp->dataParams, "phi22_scale_min" );

        LALInferenceAddVariable( datatemp->dataParams, "phi22_scale", &phi22scale, LALINFERENCE_REAL8_t,
                                 LALINFERENCE_PARAM_FIXED );
        LALInferenceAddVariable( datatemp->dataParams, "phi22_scale_min", &phi22scalemin, LALINFERENCE_REAL8_t,
                                 LALINFERENCE_PARAM_FIXED );

        datatemp = datatemp->next;
      }

      if ( LALInferenceCheckGaussianPrior( runState->priorArgs, "phi0" ) ){
        LALInferenceRemoveGaussianPrior( runState->priorArgs, "phi0" );
        LALInferenceAddGaussianPrior( runState->priorArgs, "phi22", &low, &high, LALINFERENCE_REAL8_t );
      }
      else{
        LALInferenceRemoveMinMaxPrior( runState->priorArgs, "phi0" );
        LALInferenceAddMinMaxPrior( runState->priorArgs, "phi22", &low, &high, LALINFERENCE_REAL8_t );
      }
    }

    /* remove superfluous h0 parameter */
    if ( LALInferenceGetVariableVaryType( runState->currentParams, "h0" ) != LALINFERENCE_PARAM_FIXED &&
         LALInferenceGetVariableVaryType( runState->currentParams, "C22" ) != LALINFERENCE_PARAM_FIXED ){
      LALInferenceRemoveVariable( runState->currentParams, "h0" );

      LALInferenceIFOData *datatemp = data;
      while( datatemp ){
        LALInferenceRemoveVariable( datatemp->dataParams, "h0_scale" );
        LALInferenceRemoveVariable( datatemp->dataParams, "h0_scale_min" );
        datatemp = datatemp->next;
      }

      if ( LALInferenceCheckGaussianPrior( runState->priorArgs, "h0" ) ){
        LALInferenceRemoveGaussianPrior( runState->priorArgs, "h0" );
      }
      else { LALInferenceRemoveMinMaxPrior( runState->priorArgs, "h0" ); }
    }

    /* remove superfluous phi0 parameter */
    if ( LALInferenceGetVariableVaryType( runState->currentParams, "phi0" ) != LALINFERENCE_PARAM_FIXED &&
         LALInferenceGetVariableVaryType( runState->currentParams, "phi22" ) != LALINFERENCE_PARAM_FIXED ){
      LALInferenceRemoveVariable( runState->currentParams, "phi0" );

      LALInferenceIFOData *datatemp = data;
      while( datatemp ){
        LALInferenceRemoveVariable( datatemp->dataParams, "phi0_scale" );
        LALInferenceRemoveVariable( datatemp->dataParams, "phi0_scale_min" );
        datatemp = datatemp->next;
      }

      if ( LALInferenceCheckGaussianPrior( runState->priorArgs, "phi0" ) ){
        LALInferenceRemoveGaussianPrior( runState->priorArgs, "phi0" );
      }
      else { LALInferenceRemoveMinMaxPrior( runState->priorArgs, "phi0" ); }
    }
  }

  /* if phi0 and psi have been given in the prop-file and defined at the limits of their range (for a triaxial star, so
   * -pi/4 <= psi <= pi/4 and 0 <= phi0 <= 2pi) then remove them and add the phi0' and psi' coordinates */
  //if( phidef && psidef && !strcmp( "triaxial", *(CHAR **)LALInferenceGetVariable( runState->data->dataParams,
  //                                       "modeltype" ) ) ){
  //  LALInferenceIFOData *datatemp = data;
  //
  //  REAL8 phi0 = *(REAL8*)LALInferenceGetVariable( runState->currentParams, "phi0" );
  //  REAL8 phi0scale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "phi0_scale" );
  //  REAL8 phi0min = *(REAL8*)LALInferenceGetVariable( data->dataParams, "phi0_scale_min" );
  //  REAL8 psi = *(REAL8*)LALInferenceGetVariable( runState->currentParams,  "psi" );
  //  REAL8 psiscale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "psi_scale" );
  //  REAL8 psimin = *(REAL8*)LALInferenceGetVariable( data->dataParams, "psi_scale_min" );
  //  REAL8 theta = atan2(1.,2.);
  //  REAL8 primemin = -LAL_PI_2*cos(theta);
  //  REAL8 primescale = 2.*fabs(primemin);
  //  REAL8 phi0prime = 0., psiprime = 0.;
  //
  //  phi0 = phi0*phi0scale + phi0min;
  //  psi = psi*psiscale + psimin;
  //
  //  /* convert to phi0' and psi' */
  //  phi0_psi_transform( phi0, psi, &phi0prime, &psi );
  //
  //  /* scale phi0' and psi' */
  //  phi0prime = (phi0prime - primemin)/primescale;
  //  psiprime = (psiprime - primemin)/primescale;
  //
  //  /* remove phi0 and psi */
  //  LALInferenceRemoveVariable( runState->currentParams, "phi0" );
  //  LALInferenceRemoveVariable( runState->currentParams, "psi" );
  //
  //  /* add new variables */
  //  LALInferenceAddVariable( runState->currentParams, "phi0prime", &phi0prime, LALINFERENCE_REAL8_t,
  //                           LALINFERENCE_PARAM_CIRCULAR );
  //
  //  LALInferenceAddVariable( runState->currentParams, "psiprime", &psiprime, LALINFERENCE_REAL8_t,
  //                           LALINFERENCE_PARAM_CIRCULAR );
  //
  //  /* remove old scale factors and add new ones */
  //  while( datatemp ){
  //    LALInferenceRemoveVariable( datatemp->dataParams, "phi0_scale" );
  //    LALInferenceRemoveVariable( datatemp->dataParams, "phi0_scale_min" );
  //
  //    LALInferenceRemoveVariable( datatemp->dataParams, "psi_scale" );
  //    LALInferenceRemoveVariable( datatemp->dataParams, "psi_scale_min" );
  //
  //    LALInferenceAddVariable( datatemp->dataParams, "phi0prime_scale", &primescale, LALINFERENCE_REAL8_t,
  //                             LALINFERENCE_PARAM_FIXED );
  //    LALInferenceAddVariable( datatemp->dataParams, "phi0prime_scale_min", &primemin, LALINFERENCE_REAL8_t,
  //                             LALINFERENCE_PARAM_FIXED );
  //    LALInferenceAddVariable( datatemp->dataParams, "psiprime_scale", &primescale, LALINFERENCE_REAL8_t,
  //                             LALINFERENCE_PARAM_FIXED );
  //    LALInferenceAddVariable( datatemp->dataParams, "psiprime_scale_min", &primemin, LALINFERENCE_REAL8_t,
  //                             LALINFERENCE_PARAM_FIXED );
  //
  //    datatemp = datatemp->next;
  //  }
  //
  //  /* change prior */
  //  low = 0.;
  //  high = 1.;
  //  LALInferenceRemoveMinMaxPrior( runState->priorArgs, "phi0" );
  //  LALInferenceAddMinMaxPrior( runState->priorArgs, "phi0prime", &low, &high, LALINFERENCE_REAL8_t );
  //
  //  LALInferenceRemoveMinMaxPrior( runState->priorArgs, "psi" );
  //  LALInferenceAddMinMaxPrior( runState->priorArgs, "psiprime", &low, &high, LALINFERENCE_REAL8_t );
  //}

  /* now check for a parameter correlation coefficient matrix file */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--cor-file" );
  if( ppt ){
    CHAR *corFile = XLALStringDuplicate( ppt->value );
    UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
    dims->data[0] = 1;
    dims->data[1] = 1;

    corMat = XLALCreateREAL8Array( dims );

    corParams = XLALReadTEMPOCorFile( corMat, corFile );

    /* if the correlation matrix is given then add it as the prior for values with Gaussian errors specified in the par
     * file */
    add_correlation_matrix( runState->currentParams, runState->priorArgs, corMat, corParams );

    XLALDestroyUINT4Vector( dims );
  }

  /* check if using a previous nested sampling file as a prior */
  samples_prior( runState );

  return;
}


/** \brief Initialise the MCMC proposal distribution for sampling new points
 *
 * There are various proposal distributions that can be used to sample new live points via an MCMC. A combination of
 * different ones can be used to help efficiency for awkward posterior distributions. Here the proposals that can be
 * used are:
 *   \c covariance Drawing from a multi-variate Gaussian described by the covariance matrix of the current live points,
 * with the spread of the distribution controlled by the \c temperature. One parameter is evolved during a single draw.
 *   \c diffev Drawing a new point by differential evolution of two randomly chosen live points. All parameters are
 * evolved during a single draw.
 *   \c kDTree Drawing points from a distributions created from a k-D tree of the current live points, with
 * probabilities of each leaf being inversely their volume. All parameters are evolved during a single draw.
 *
 * Note: also add ability to jump between frequency modes.
 *
 * This function sets up the relative weights with which each of above distributions is used.
 *
 * \param runState [in] A pointer to the run state
*/
void initialiseProposal( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL;
  UINT4 covfrac = 0, defrac = 0, kdfrac = 0;
  REAL8 temperature = 0.;
  const CHAR *defaultPropName = NULL;
  defaultPropName = XLALStringDuplicate( "none" );

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--covariance" );
  if( ppt ) { covfrac = atoi( ppt->value ); }
  else { covfrac = 14; } /* default value */

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--diffev" );
  if( ppt ) { defrac = atoi( ppt->value ); }
  else { defrac = 3; } /* default value */

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--kDTree" );
  if( ppt ) { kdfrac = atoi( ppt->value ); }
  else { kdfrac = 3; } /* default value */

  if( !covfrac && !defrac && !kdfrac ){
    XLALPrintError("All proposal weights are zero!\n");
    XLAL_ERROR_VOID(XLAL_EFAILED);
  }

  runState->proposalStats = NULL;
  if(!runState->proposalStats) runState->proposalStats = XLALCalloc(1,sizeof(LALInferenceVariables));

  /* add proposals */
  if( covfrac ){
    LALInferenceAddProposalToCycle( runState, covarianceEigenvectorJumpName, &LALInferenceCovarianceEigenvectorJump,
                                    covfrac );
  }

  if( defrac ){
    LALInferenceAddProposalToCycle( runState, differentialEvolutionFullName, &LALInferenceDifferentialEvolutionFull,
                                    defrac );
  }

  if( kdfrac ){
    /* set the maximum number of points in a kd-tree cell if given */
    ppt = LALInferenceGetProcParamVal( runState->commandLine, "--kDNCell" );
    if( ppt ){
      INT4 kdncells = atoi( ppt->value );

      LALInferenceAddVariable( runState->proposalArgs, "KDNCell", &kdncells, LALINFERENCE_INT4_t,
                               LALINFERENCE_PARAM_FIXED );
    }

    LALInferenceAddProposalToCycle( runState, KDNeighborhoodProposalName, &LALInferenceKDNeighborhoodProposal, kdfrac );

    LALInferenceSetupkDTreeNSLivePoints( runState );
  }

  LALInferenceRandomizeProposalCycle( runState );
  /* set temperature */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--temperature" );
  if( ppt ) { temperature = atof( ppt->value ); }
  else { temperature = 0.1; }

  LALInferenceAddVariable( runState->proposalArgs, "temperature", &temperature, LALINFERENCE_REAL8_t,
                           LALINFERENCE_PARAM_FIXED );

  /* add default proposal name */
  LALInferenceAddVariable( runState->proposalArgs, LALInferenceCurrentProposalName, &defaultPropName,
                           LALINFERENCE_string_t, LALINFERENCE_PARAM_OUTPUT );

  /* set proposal */
  runState->proposal = LALInferenceDefaultProposal;
}


/** \brief Adds a correlation matrix for a multi-variate Gaussian prior
 *
 * If a TEMPO-style parameter correlation coefficient file has been given, then this function will use it to set the
 * prior distribution for the given parameters. It is assumed that the equivalent par file contained standard
 * deviations for all parameters given in the correlation matrix file, but if  the correlation matrix contains more
 * parameters they will be ignored.
 */
void add_correlation_matrix( LALInferenceVariables *ini, LALInferenceVariables *priors, REAL8Array *corMat,
                             LALStringVector *parMat ){
  UINT4 i = 0, j = 0, k = 0;
  LALStringVector *newPars = NULL;
  gsl_matrix *corMatg = NULL;
  UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
  UINT4 corsize = corMat->dimLength->data[0];
  UINT4 corshrink = corsize;

  /* loop through parameters and find ones that have Gaussian priors set - these should match with parameters in the
   * correlation coefficient matrix */
  for ( i = 0; i < parMat->length; i++ ){
    UINT4 incor = 0;
    LALInferenceVariableItem *checkPrior = ini->head;

    for( ; checkPrior ; checkPrior = checkPrior->next ){
      if( LALInferenceCheckGaussianPrior(priors, checkPrior->name) ){
        /* ignore parameter name case */
        if( !strcasecmp(parMat->data[i], checkPrior->name) ){
          incor = 1;

          /* add parameter to new parameter string vector */
          newPars = XLALAppendString2Vector( newPars, parMat->data[i] );
          break;
        }
      }
    }

    /* if parameter in the corMat did not match one with a Gaussian defined prior, then remove it from the matrix */
    if ( incor == 0 ){
      /* remove the ith row and column from corMat, and the ith name from parMat */
      /* shift rows up */
      for ( j = i+1; j < corsize; j++ )
        for ( k = 0; k < corsize; k++ )
          corMat->data[(j-1)*corsize + k] = corMat->data[j*corsize + k];

      /* shift columns left */
      for ( k = i+1; k < corsize; k++ )
        for ( j = 0; j < corsize; j++ )
          corMat->data[j*corsize + k-1] = corMat->data[j*corsize + k];

      /* resize array */
      corshrink--;
    }
  }

  XLALDestroyUINT4Vector( dims );

  /* return new parameter string vector as old one */
  XLALDestroyStringVector( parMat );
  parMat = newPars;

  /* copy the corMat into a gsl_matrix */
  corMatg = gsl_matrix_alloc( parMat->length, parMat->length );
  for ( i = 0; i < parMat->length; i++ )
    for ( j = 0; j < parMat->length; j++ )
      gsl_matrix_set( corMatg, i, j, corMat->data[i*corsize + j] );

  /* re-loop over parameters removing Gaussian priors on those in the parMat and replacing with a correlation matrix */
  for ( i = 0; i < parMat->length; i++ ){
    LALInferenceVariableItem *checkPrior = ini->head;

    /* allocate global variable giving the list of the correlation matrix parameters */
    corlist = XLALAppendString2Vector( corlist, parMat->data[i] );

    for( ; checkPrior ; checkPrior = checkPrior->next ){
      if( LALInferenceCheckGaussianPrior(priors, checkPrior->name) ){
        if( !strcasecmp(parMat->data[i], checkPrior->name) ){
          /* remove the Gaussian prior */
          LALInferenceRemoveGaussianPrior( priors, checkPrior->name );

          /* replace it with the correlation matrix as a gsl_matrix */
          LALInferenceAddCorrelatedPrior( priors, checkPrior->name, &corMatg, &i );

          break;
        }
      }
    }
  }
}


/*------------------- END INITIALISATION FUNCTIONS ---------------------------*/


/******************************************************************************/
/*                       SOFTWARE INJECTION FUNCTIONS                         */
/******************************************************************************/

/** \brief Inject a simulated signal into the data
 *
 * This function will create an simulated signal (of the required model) to inject into the data from multiple
 * detectors. The parameters of the signal to be injected must be specified in a TEMPO-stype .par file given with the
 * \c inject-file command line argument. The parameters do not have to be the same as those in the .par file controlling
 * the analysis (although should ideally contain a signal within the bandwidth of the data).
 *
 * If a signal of a specific signal-to-noise ratio is required then the \c scale-snr command line argument can be used
 * to give the multi-detector SNR to which the signal needs to be scaled.
 *
 * The injected signal can be output if \c inject-output is set. Two files will be output: one containing the signal
 * only, and one containing the signal plus noise. These will both be in the format of a standard data input file. The
 * files will have names given by the \c inject-output value, with a prefix of the detector name, and a suffix of of \c
_* signal_only, respectively.
 *
 * \param runState [in] the program information structure
 *
 * \sa calculate_time_domain_snr
 */
void injectSignal( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;

  CHAR *injectfile = NULL, *snrfile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;

  BinaryPulsarParams injpars;

  FILE *fpsnr = NULL; /* output file for SNRs */
  INT4 ndats = 0, j = 1, k = 0;

  REAL8Vector *freqFactors = NULL;
  REAL8 snrmulti = 0.;
  REAL8 snrscale = 0;

  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( !ppt ){
    fprintf(stderr, "Error... no output file specified!\n");
    exit(0);
  }

  snrfile = XLALStringDuplicate( ppt->value );
  snrfile = XLALStringAppend( snrfile, "_SNR" );

  if( (fpsnr = fopen(snrfile, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open output SNR file!\n");
    exit(0);
  }

  ppt = LALInferenceGetProcParamVal( commandLine, "--inject-file" );
  if( ppt ){
    injectfile = XLALStringDuplicate( ppt->value );

    /* check that the file exists */
    if ( fopen(injectfile, "r") == NULL ){
      fprintf(stderr, "Error... Injection specified, but the injection parameter file %s is wrong.\n", injectfile);
      exit(3);
    }

    /* read in injection parameter file */
    XLALReadTEMPOParFile( &injpars, injectfile );

    /* make sure that we have parameters in terms of amplitude and phase parameters */
    invert_source_params( &injpars );
  }
  else{
    fclose( fpsnr );
    return;
  }

  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "freqfactors" );

  /* get the SNR scale factor if required */
  ppt = LALInferenceGetProcParamVal( commandLine, "--scale-snr" );
  if ( ppt ) { snrscale = atof( ppt->value ); }

  /* create signal to inject */
  while( data ){
    /* for injection always attempt to include the signal phase model even if the search is not going to be over
     * phase */
    UINT4 varyphasetmp = varyphase;
    varyphase = 1;

    pulsar_model( injpars, data );

    /* reset varyphase to its original value */
    varyphase = varyphasetmp;

    data = data->next;

    /* If model uses more than one data stream need to advance data on to next, so this loop only runs once if there
     * is only 1 det */
    for ( k = 1; k < (INT4)freqFactors->length; k++ ){ data = data->next; }
  }

  /* reset data to head */
  data = runState->data;

  fprintf(fpsnr, "# Injected SNR\n");

  /* calculate SNRs */
  while ( data ){
    REAL8 snrval = calculate_time_domain_snr( data );

    snrmulti += SQUARE(snrval);

    //fprintf(fpsnr, "freq_factor: %lf, non-scaled snr: %le\n", freqFactors->data[ndats%(INT4)freqFactors->length],
    //        snrval);

    /* if not scaling print out individual detector/datastream SNRs */
    if( snrscale == 0. ){
      fprintf(fpsnr, "%s\t%.3lf\t%le\n", data->name, freqFactors->data[ndats%(INT4)freqFactors->length], snrval);
    }

    data = data->next;

    ndats++;
  }

  /* get overall multi-detector SNR */
  snrmulti = sqrt(snrmulti);

  /* only need to print out multi-detector snr if the were multiple detectors or data streams */
  if ( snrscale == 0. ){
    if ( ndats > 1 ) { fprintf(fpsnr, "Coherent\t%le\n", snrmulti); }
  }
  else{
    /* rescale the signal and calculate the SNRs */
    data = runState->data;
    snrscale /= snrmulti;

    injpars.C22 *= snrscale;
    injpars.C21 *= snrscale;

    /* recreate the signal with scaled amplitude */
    while( data ){
      UINT4 varyphasetmp = varyphase;
      varyphase = 1;

      pulsar_model( injpars, data );

      /* reset varyphase to its original value */
      varyphase = varyphasetmp;

      data = data->next;

      for ( k = 1; k < (INT4)freqFactors->length; k++ ) { data = data->next; }
    }

    data = runState->data;

    /* get new snrs */
    snrmulti = 0;
    ndats = 0;

    while( data ){
      REAL8 snrval = 0.;

      /* recalculate the SNR */
      snrval = calculate_time_domain_snr( data );

      snrmulti += SQUARE(snrval);

      fprintf(fpsnr, "%s\t%.3lf\t%le\n", data->name, freqFactors->data[ndats%(INT4)freqFactors->length], snrval);

      data = data->next;

      ndats++;
    }

    snrmulti = sqrt( snrmulti );
    //fprintf(stderr, "scaled multi data snr: %le\n", snrmulti);

    if( ndats > 1 ){
      fprintf(fpsnr, "Coherent\t%le\n", snrmulti);
    }
  }

  fclose( fpsnr );

  /* reset data to head */
  data = runState->data;

  /* add signal to data */
  while( data ){
    FILE *fp = NULL, *fpso = NULL;
    ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal( commandLine, "--inject-output" );
    INT4 i = 0, length = data->dataTimes->length;

    /* check whether to output the data */
    if ( ppt2 ){
      /* add the site prefix to the start of the output name */
      CHAR *outfile = NULL;
      CHAR *signalonly = NULL; /* file containing only signal and no noise */
      CHAR suffix[5];
      INT4 sf;

      outfile = XLALStringDuplicate( ppt2->value );

      /* append detector name to file */
      outfile = XLALStringAppend( outfile, "_" );
      outfile = XLALStringAppend( outfile, data->detector->frDetector.prefix );

      /* append the harmonic frequency of the signal in the file */
      sf = sprintf(suffix, "_%.1lf", freqFactors->data[(INT4)fmod(j,freqFactors->length)]);
      outfile = XLALStringAppend( outfile, suffix );

      if ( (fp = fopen(outfile, "w")) == NULL || !sf ){
        fprintf(stderr, "Non-fatal error... unable to open file %s to output injection\n", outfile);
      }

      signalonly = XLALStringDuplicate( outfile );
      signalonly = XLALStringAppend( signalonly, "_signal_only" );

      if ( (fpso = fopen(signalonly, "w")) == NULL ){
        fprintf(stderr, "Non-fatal error... unable to open file %s to output injection\n", signalonly);
      }
    }

    /* add the signal to the data */
    for ( i = 0; i < length; i++ ){
      data->compTimeData->data->data[i] += data->compModelData->data->data[i];

      /* write out injection to file */
      if( fp != NULL && fpso != NULL ){
        /* print out data - time stamp, real and imaginary parts of data (injected signal + noise) */
        fprintf(fp, "%.5lf\t%le\t%le\n", XLALGPSGetREAL8( &data->dataTimes->data[i] ),
                creal(data->compTimeData->data->data[i]), cimag(data->compTimeData->data->data[i]) );

        /* print signal only data - time stamp, real and imaginary parts of signal */
        fprintf(fpso, "%.5lf\t%le\t%le\n", XLALGPSGetREAL8( &data->dataTimes->data[i] ),
                creal(data->compModelData->data->data[i]), cimag(data->compModelData->data->data[i]) );
      }
    }

    if ( fp != NULL ) { fclose( fp ); }
    if ( fpso != NULL ) { fclose( fpso ); }

    data = data->next;
    j++;
  }
}

/*-------------------- END OF SOFTWARE INJECTION FUNCTIONS -------------------*/


/******************************************************************************/
/*                            HELPER FUNCTIONS                                */
/******************************************************************************/

/** \brief Split the data into segments
 *
 * This function is deprecated to \c chop_n_merge, but gives the functionality of the old code.
 *
 * It cuts the data into as many contiguous segments of data as possible of length \c chunkMax. Where contiguous is
 * defined as containing consecutive point within 180 seconds of each other. The length of segments that do not fit into
 * a \c chunkMax length are also included.
 *
 * \param data [in] a data structure
 * \param chunkMax [in] the maximum length of a data chunk/segment
 *
 * \return A vector of chunk/segment lengths
 */
UINT4Vector *get_chunk_lengths( LALInferenceIFOData *data, INT4 chunkMax ){
  INT4 i = 0, j = 0, count = 0;
  INT4 length;

  REAL8 t1, t2;

  UINT4Vector *chunkLengths = NULL;

  length = data->dataTimes->length;

  chunkLengths = XLALCreateUINT4Vector( length );

  REAL8 dt = *(REAL8*)LALInferenceGetVariable( data->dataParams, "dt" );

  /* create vector of data segment length */
  while( 1 ){
    count++; /* counter */

    /* break clause */
    if( i > length - 2 ){
      /* set final value of chunkLength */
      chunkLengths->data[j] = count;
      j++;
      break;
    }

    i++;

    t1 = XLALGPSGetREAL8( &data->dataTimes->data[i-1 ]);
    t2 = XLALGPSGetREAL8( &data->dataTimes->data[i] );

    /* if consecutive points are within two sample times of each other count as in the same chunk */
    if( t2 - t1 > 2.*dt || count == chunkMax ){
      chunkLengths->data[j] = count;
      count = 0; /* reset counter */

      j++;
    }
  }

  chunkLengths = XLALResizeUINT4Vector( chunkLengths, j );

  return chunkLengths;
}


/* function to use change point analysis to chop up and remerge the data to find stationary chunks (i.e. lengths of data
 * which look like they have the same statistics e.g. the same standard deviation) */
/** \brief Chops and remerges data into stationary segments
 *
 * This function finds segments of data that appear to be stationary (have the same standard deviation).
 *
 * The function first attempts to chop up the data into as many stationary segments as possible. The splitting may not
 * be optimal, so it then tries remerging consecutive segments to see if the merged segments show more evidence of
 * stationarity. <b>[NOTE: Remerging is currently turned off and will make very little difference to the algorithm]</b>.
 * It then, if necessary, chops the segments again to make sure there are none greater than the required \c chunkMax.
 * The default \c chunkMax is 0, so this rechopping will not normally happen.
 *
 * This is all performed on data that has had a running median subtracted, to try and removed any underlying trends in
 * the data (e.g. those caused by a strong signal), which might affect the calculations (which assume the data is
 * Gaussian with zero mean).
 *
 * If the \c verbose flag is set then a list of the segments will be output to a file called \c data_segment_list.txt,
 * with a prefix of the detector name.
 *
 * \param data [in] A data structure
 * \param chunkMin [in] The minimum length of a segment
 * \param chunkMax [in] The maximum length of a segment
 *
 * \return A vector of segment/chunk lengths
 *
 * \sa subtract_running_median
 * \sa chop_data
 * \sa merge_data
 * \sa rechop_data
 */
UINT4Vector *chop_n_merge( LALInferenceIFOData *data, INT4 chunkMin, INT4 chunkMax ){
  UINT4 j = 0;

  UINT4Vector *chunkLengths = NULL;
  UINT4Vector *chunkIndex = NULL;

  COMPLEX16Vector *meddata = NULL;

  /* subtract a running median value from the data to remove any underlying trends (e.g. caused by a string signal) that
   * might affect the chunk calculations (which can assume the data is Gaussian with zero mean). */
  meddata = subtract_running_median( data->compTimeData->data );

  chunkIndex = chop_data( meddata, chunkMin );

  /* DON'T BOTHER WITH THE MERGING AS IT WILL MAKE VERY LITTLE DIFFERENCE */
  /* merge_data( meddata, chunkIndex ); */

  /* if a maximum chunk length is defined then rechop up the data, to segment any chunks longer than this value */
  if ( chunkMax > chunkMin ) { rechop_data( chunkIndex, chunkMax, chunkMin ); }

  chunkLengths = XLALCreateUINT4Vector( chunkIndex->length );

  /* go through segments and turn into vector of chunk lengths */
  for ( j = 0; j < chunkIndex->length; j++ ){
    if ( j == 0 ) { chunkLengths->data[j] = chunkIndex->data[j]; }
    else { chunkLengths->data[j] = chunkIndex->data[j] - chunkIndex->data[j-1]; }
  }

  /* if verbose print out the segment end indices to a file */
  if ( verbose_output ){
    FILE *fpsegs = NULL;

    CHAR *outfile = NULL;

    /* set detector name as prefix */
    outfile = XLALStringDuplicate( data->detector->frDetector.prefix );

    outfile = XLALStringAppend( outfile, "data_segment_list.txt" );

    if ( (fpsegs = fopen(outfile, "w")) == NULL ){
      fprintf(stderr, "Non-fatal error open file to output segment list.\n");
      return chunkLengths;
    }

    for ( j = 0; j < chunkIndex->length; j++ ) { fprintf(fpsegs, "%u\n", chunkIndex->data[j]); }

    /* add space at the end so that you can separate lists from different detector data streams */
    fprintf(fpsegs, "\n");

    fclose( fpsegs );
  }

  return chunkLengths;
}


/** \brief Subtract the running median from complex data
 *
 * This function uses \c gsl_stats_median_from_sorted_data to subtract a running median, calculated from the 30
 * consecutive point around a set point, from the data. At the start of the data running median is calculated from
 * 30-15+(i-1) points, and at the end it is calculated from 15+(N-i) points, where i is the point index and N is the
 * total number of data points.
 *
 * \param data [in] A complex data vector
 *
 * \return A complex vector containing data with the running median removed
 */
COMPLEX16Vector *subtract_running_median( COMPLEX16Vector *data ){
  COMPLEX16Vector *submed = NULL;
  UINT4 length = data->length, i = 0, j = 0, n = 0;
  UINT4 RANGE = 30; /* perform running median with 30 data points */
  UINT4 N = (UINT4)floor(RANGE/2);
  INT4 sidx = 0;

  submed = XLALCreateCOMPLEX16Vector( length );

  for ( i = 1; i < length+1; i++ ){
    REAL8 *dre = NULL;
    REAL8 *dim = NULL;

    /* get median of data within RANGE */
    if ( i < N ){
      n = N+i;
      sidx = 0;
    }
    else if ( i > length - N ){
      n = length - i + N;
      sidx = (i-N)-1;
    }
    else{
      n = RANGE;
      sidx = i-N;
    }

    dre = XLALCalloc( n, sizeof(REAL8) );
    dim = XLALCalloc( n, sizeof(REAL8) );

    for ( j = 0; j < n; j++ ){
      dre[j] = creal(data->data[j+sidx]);
      dim[j] = cimag(data->data[j+sidx]);
    }

    /* sort data */
    gsl_sort( dre, 1, n );
    gsl_sort( dim, 1, n );

    /* get median and subtract from data*/
    submed->data[i-1] = ( creal(data->data[i-1]) - gsl_stats_median_from_sorted_data( dre, 1, n ) )
      + I * ( cimag(data->data[i-1]) - gsl_stats_median_from_sorted_data( dim, 1, n ) );

    XLALFree( dre );
    XLALFree( dim );
  }

  return submed;
}


/** \brief Chops the data into stationary segments based on Bayesian change point analysis
 *
 * This function splits data into two (and recursively runs on those two segments) if it is found that the odds ratio
 * for them being from two independent Gaussian distributions is greater than a certain threshold.
 *
 * The threshold is for the natural logarithm of the odds ratio is empirically set to be:
 * \f[
 * T = 0.57\ln{N} + 2.71,
 * \f]
 * where \f$N\f$ is the length of the data set. This comes from a fit to the threshold value required to give a 1%
 * chance of splitting actual Gaussian data (drawn from one distribution) for data of various lengths. The first two
 * terms come from a fit to odds ratios for a Monte Carlo of Gaussian noise (with real and imaginary components) of
 * various lengths, and the final term comes from an offset to give the 1% false alarm rate.
 *
 * \param data [in] A complex data vector
 * \param chunkMin [in] The minimum allowed segment length
 *
 * \return A vector of segment lengths
 *
 * \sa find_change_point
 */
UINT4Vector *chop_data( COMPLEX16Vector *data, INT4 chunkMin ){
  UINT4Vector *chunkIndex = NULL;

  UINT4 length = data->length;

  REAL8 logodds = 0;
  UINT4 changepoint = 0;

  REAL8 threshold = 0.; /* may need tuning or setting globally */

  chunkIndex = XLALCreateUINT4Vector( 1 );

  changepoint = find_change_point( data, &logodds, chunkMin );

  threshold = 0.57*log(length) + 2.71;

  if ( logodds > threshold ){
    UINT4Vector *cp1 = NULL;
    UINT4Vector *cp2 = NULL;

    COMPLEX16Vector *data1 = XLALCreateCOMPLEX16Vector( changepoint );
    COMPLEX16Vector *data2 = XLALCreateCOMPLEX16Vector( length - changepoint );

    UINT4 i = 0, l = 0;

    /* fill in data */
    for (i = 0; i < changepoint; i++) { data1->data[i] = data->data[i]; }

    for (i = 0; i < length - changepoint; i++) { data2->data[i] = data->data[i+changepoint]; }

    cp1 = chop_data( data1, chunkMin );
    cp2 = chop_data( data2, chunkMin );

    l = cp1->length + cp2->length;

    chunkIndex = XLALResizeUINT4Vector( chunkIndex, l );

    /* combine new chunks */
    for (i = 0; i < cp1->length; i++) { chunkIndex->data[i] = cp1->data[i]; }

    for (i = 0; i < cp2->length; i++) { chunkIndex->data[i+cp1->length] = cp2->data[i] + changepoint; }

    /* free memory */
    XLALDestroyCOMPLEX16Vector( data1 );
    XLALDestroyCOMPLEX16Vector( data2 );
  }
  else{ chunkIndex->data[0] = length; }

  return chunkIndex;
}


/** \brief Find a change point in complex data
 *
 * This function is based in the Bayesian Blocks algorithm of [\ref Scargle1998] that finds "change points" in data -
 * points at which the statistics of the data change. It is based on calculating evidence, or odds, ratios. The
 * function first computes the marginal likelihood (or evidence) that the whole of the data is described by a single
 * Gaussian (with mean of zero). This comes from taking a Gaussian likelihood function and analytically marginalising
 * over the standard deviation (using a prior on the standard deviation of \f$1/\sigma\f$), giving (see [\ref
 * DupuisWoan2005]) a Students-t distribution (see
 * <a href="https://wiki.ligo.org/foswiki/pub/CW/PulsarParameterEstimationNestedSampling/studentst.pdf">here</a>).
 * Following this the data is split into two segments (with lengths greater than, or equal to the minimum chunk length)
 * for all possible combinations, and the joint evidence for each of the two segments consisting of independent
 * Gaussian (basically multiplying the above equation calculated for each segment separately) is calculated and the
 * split point recorded. However, the value required for comparing to that for the whole data set, to give the odds
 * ratio, is the evidence that having any split is better than having no split, so the individual split data evidences
 * need to be added incoherently to give the total evidence for a split. The index at which the evidence for a single
 * split is maximum (i.e. the most favoured split point) that which is returned.
 *
 * \param data [in] a complex data vector
 * \param logodds [in] a pointer to return the natural logarithm of the odds ratio/Bayes factor
 * \param minlength [in] the minimum chunk length
 *
 * \return The position of the change point
 */
UINT4 find_change_point( COMPLEX16Vector *data, REAL8 *logodds, INT4 minlength ){
  UINT4 changepoint = 0, i = 0;
  UINT4 length = data->length, lsum = 0;

  REAL8 datasum = 0.;

  REAL8 logsingle = 0., logtot = -LAL_REAL8_MAX;
  REAL8 logdouble = 0., logdouble_min = -LAL_REAL8_MAX;
  REAL8 logratio = 0.;

  REAL8Vector *sumforward = NULL, *sumback = NULL;

  /* check that data is at least twice the minimum length, if not return an odds ratio of zero (log odds = -inf [or
   * close to that!]) */
  if ( length < (UINT4)(2*minlength) ){
    logratio = -LAL_REAL8_MAX;
    memcpy(logodds, &logratio, sizeof(REAL8));
    return 0;
  }

  /* calculate the sum of the data squared */
  for (i = 0; i < length; i++) { datasum += SQUARE( cabs(data->data[i]) ); }

  /* calculate the evidence that the data consists of a Gaussian data with a
     single standard deviation */
  logsingle = -2 + logfactorial[length-1] - (REAL8)length * log( datasum );

  /* to speed up process calculate data sums first */
  lsum = length - 2*minlength + 1;
  sumforward = XLALCreateREAL8Vector( lsum );
  sumback = XLALCreateREAL8Vector( lsum );

  sumforward->data[0] = 0.;
  sumback->data[0] = 0.;

  for ( i = 0; i < length - minlength; i++ ){
    if ( i < (UINT4)minlength ){
      sumforward->data[0] += SQUARE( cabs(data->data[i]) );
      sumback->data[0] += SQUARE( cabs(data->data[length-(i+1)]) );
    }
    else{
      sumforward->data[i+1-minlength] = sumforward->data[i-minlength] + SQUARE( cabs(data->data[i]) );
      sumback->data[i+1-minlength] = sumback->data[i-minlength] + SQUARE( cabs(data->data[length-(i+1)]) );
    }
  }

  /* go through each possible change point and calculate the evidence for the data consisting of two independent
   * Gaussian's either side of the change point. Also calculate the total evidence for any change point.
   * Don't allow single points, so start at the second data point. */
  for (i = 0; i < lsum; i++){
    UINT4 ln1 = i+minlength, ln2 = (length-i-minlength);

    REAL8 log_1 = 0., log_2 = 0.;

    /* get log evidences for the individual segments */
    log_1 = -2 + logfactorial[ln1-1] - (REAL8)ln1 * log( sumforward->data[i] );
    log_2 = -2 + logfactorial[ln2-1] - (REAL8)ln2 * log( sumback->data[lsum-i-1] );

    /* get evidence for the two segments */
    logdouble = log_1 + log_2;

    /* add to total evidence for a change point */
    logtot = LOGPLUS(logtot, logdouble);

    /* find maximum value of logdouble and record that as the change point */
    if ( logdouble > logdouble_min ){
      changepoint = ln1;
      logdouble_min = logdouble;
    }
  }

  /* get the log odds ratio of segmented versus non-segmented model */
  logratio = logtot - logsingle;
  memcpy(logodds, &logratio, sizeof(REAL8));

  XLALDestroyREAL8Vector( sumforward );
  XLALDestroyREAL8Vector( sumback );

  return changepoint;
}


/** \brief Chop up the data into chunks smaller the the maximum allowed length
 *
 * This function chops any chunks that are greater than \c chunkMax into chunks smaller than, or equal to \c chunkMax,
 * and greater than \c chunkMin. On some occasions this might result in a segment smaller than \c chunkMin, but these
 * are ignored in the likelihood calculation anyway.
 *
 * \param chunkIndex [in] a vector of segment split positions
 * \param chunkMax [in] the maximum allowed segment/chunk length
 * \param chunkMin [in] the minimum allowed segment/chunk length
 */
void rechop_data( UINT4Vector *chunkIndex, INT4 chunkMax, INT4 chunkMin ){
  INT4 i = 0, j = 0, count = 0;
  INT4 length = chunkIndex->length;
  INT4 endIndex = (INT4)chunkIndex->data[length-1];
  UINT4 startindex = 0, chunklength = 0;

  UINT4Vector *newindex = NULL;
  newindex = XLALCreateUINT4Vector( ceil((REAL8)endIndex / (REAL8)chunkMax ) );

  /* chop any chunks that are greater than chunkMax into chunks smaller than, or equal to chunkMax, and greater than
   * chunkMin */
  for ( i = 0; i < length; i++ ){
    if ( i == 0 ) { startindex = 0; }
    else { startindex = chunkIndex->data[i-1]+1; }

    chunklength = chunkIndex->data[i] - startindex;

    if ( chunklength > (UINT4)chunkMax ){
      INT4 remain = chunklength % chunkMax;

      /* cut segment into as many chunkMin chunks as possible */
      for ( j = 0; j < floor(chunklength / chunkMax); j++ ){
        newindex->data[count] = startindex + (j+1)*chunkMax;
        count++;
      }

      /* last chunk values */
      if ( remain != 0 ){
        /* set final value */
        newindex->data[count] = startindex + j*chunkMax + remain;

        if ( remain < chunkMin ){
          /* split the last two cells into one that is chunkMin long and one that is (chunkMax+remainder)-chunkMin long
           * - this may leave a cell shorter than chunkMin, but we'll have to live with that! */
          INT4 n1 = (chunkMax + remain) - chunkMin;

          /* reset second to last value two values */
          newindex->data[count-1] = newindex->data[count] - chunkMin;

          if ( n1 < chunkMin && verbose_output ){
            fprintf(stderr, "Non-fatal error... segment no. %d is %d long, which is less than chunkMin = %d.\n",
                    count, n1, chunkMin);
          }
        }

        count++;
      }
    }
    else{
      newindex->data[count] = chunkIndex->data[i];
      count++;
    }
  }

  chunkIndex = XLALResizeUINT4Vector( chunkIndex, count );

  for ( i = 0; i < count; i++ ) { chunkIndex->data[i] = newindex->data[i]; }

  XLALDestroyUINT4Vector( newindex );
}


/** \brief Merge adjacent segments
 *
 * This function will attempt to remerge adjacent segments if statistically favourable (as calculated by the odds
 * ratio). For each pair of adjacent segments the joint likelihood of them being from two independent distributions is
 * compared to the likelihood that combined they are from one distribution. If the likelihood is highest for the
 * combined segments they are merged.
 *
 * \param data [in] A complex data vector
 * \param segs [in] A vector of split segment indexes
 */
void merge_data( COMPLEX16Vector *data, UINT4Vector *segs ){
  UINT4 j = 0;
  REAL8 threshold = 0.; /* may need to be passed to function in the future, or defined globally */

  /* loop until stopping criterion is reached */
  while( 1 ){
    UINT4 ncells = segs->length;

    UINT4 mergepoint = 0;
    REAL8 logodds = 0., minl = -LAL_REAL8_MAX;

    for (j = 1; j < ncells; j++){
      REAL8 summerged = 0., sum1 = 0., sum2 = 0.;
      UINT4 i = 0, n1 = 0, n2 = 0, nm = 0;
      UINT4 cellstarts1 = 0, cellends1 = 0, cellstarts2 = 0, cellends2 = 0;
      REAL8 log_merged = 0., log_individual = 0.;

      /* get the evidence for merged adjacent cells */
      if( j == 1 ) { cellstarts1 = 0; }
      else { cellstarts1 = segs->data[j-2]; }

      cellends1 = segs->data[j-1];

      cellstarts2 = segs->data[j-1];
      cellends2 = segs->data[j];

      n1 = cellends1 - cellstarts1;
      n2 = cellends2 - cellstarts2;
      nm = cellends2 - cellstarts1;

      for( i = cellstarts1; i < cellends1; i++ ) { sum1 += SQUARE( cabs(data->data[i]) ); }

      for( i = cellstarts2; i < cellends2; i++ ) { sum2 += SQUARE( cabs(data->data[i]) ); }

      summerged = sum1 + sum2;

      /* calculated evidences */
      log_merged = -2 + logfactorial[nm-1] - (REAL8)nm * log( summerged );

      log_individual = -2 + logfactorial[n1-1] - (REAL8)n1 * log( sum1 );
      log_individual += -2 + logfactorial[n2-1] - (REAL8)n2 * log( sum2 );

      logodds = log_merged - log_individual;

      if ( logodds > minl ){
        mergepoint = j - 1;
        minl = logodds;
      }
    }

    /* set break criterion */
    if ( minl < threshold ) { break; }
    else{ /* merge cells */
      /* remove the cell end value between the two being merged and shift */
      for( UINT4 i=0; i < ncells-(mergepoint+1); i++ ){
        segs->data[mergepoint+i] = segs->data[mergepoint+i+1];
      }

      segs = XLALResizeUINT4Vector( segs, ncells - 1 );
    }
  }
}


/** \brief Calculates the sum of the square of the data
 *
 * This function calculates the sum of the square of the data:
 * \f[
 * s = \sum_i^N \Re{d_i}^2 + \Im{d_i}^2,
 * \f]
 * for each stationary segment given in the \c chunkLength vector. These value are used in the likelihood calculation in
 * \c pulsar_log_likelihood and are precomputed here to speed that calculation up. The vector of value is output
 * in a \c sumData parameter of the \c data structure.
 *
 * \param runState [in] The analysis information structure
 */
void sumData( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;

  while( data ){
    REAL8Vector *sumdat = NULL;

    INT4 chunkLength = 0, length = 0, i = 0, j = 0, count = 0;
    COMPLEX16 B;

    UINT4Vector *chunkLengths;

    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams, "chunkLength" );

    length = data->dataTimes->length + 1 - chunkLengths->data[chunkLengths->length - 1];

    sumdat = XLALCreateREAL8Vector( chunkLengths->length );

    for( i = 0 ; i < length ; i+= chunkLength ){
      chunkLength = chunkLengths->data[count];

      sumdat->data[count] = 0.;

      for( j = i ; j < i + chunkLength ; j++){
        B = data->compTimeData->data->data[j];

        /* sum up the data */
        sumdat->data[count] += (creal(B)*creal(B) + cimag(B)*cimag(B));
      }

      count++;
    }

    LALInferenceAddVariable( data->dataParams, "sumData", &sumdat, LALINFERENCE_REAL8Vector_t,
                             LALINFERENCE_PARAM_FIXED );

    data = data->next;
  }

  return;
}


/** \brief Creates a lookup table of the detector antenna pattern
 *
 * This function creates a 2D lookup table of the 'plus' and 'cross' antenna patterns for a given detector orientation
 * and source sky position. The lookup table spans one sidereal day in time (this being the period over which the
 * antenna pattern changes) and goes between \f$\pm\pi/2\f$ radians in \f$\psi\f$ (this is the full range of \f$\psi\f$
 * that should be required for and model).
 *
 * Note: the may want to be converted into an XLAL function and moved into LAL at some point.
 *
 * \param t0 [in] initial GPS time of the data
 * \param detNSource [in] structure containing the detector and source orientations and locations
 * \param timeSteps [in] the number of grid bins to use in time
 * \param psiSteps [in] the number of grid bins to use in polarisation angle \f$\psi\f$
 * \param LUfplus [in] a matrix into which the 'plus' antenna pattern lookup table will be output
 * \param LUfcross [in] a matrix into which the 'cross' antenna pattern lookup table will be output
 */
void response_lookup_table( REAL8 t0, LALDetAndSource detNSource, INT4 timeSteps, INT4 psiSteps, gsl_matrix *LUfplus,
                            gsl_matrix *LUfcross ){
  LIGOTimeGPS gps;
  REAL8 T = 0;

  REAL8 fplus = 0., fcross = 0.;
  REAL8 psteps = (REAL8)psiSteps;
  REAL8 tsteps = (REAL8)timeSteps;

  INT4 i = 0, j = 0;

  for( i = 0 ; i < psiSteps ; i++ ){
    detNSource.pSource->orientation = -(LAL_PI_2) + (REAL8)i*(LAL_PI) / ( psteps - 1. );

    for( j = 0 ; j < timeSteps ; j++ ){
      T = t0 + (REAL8)j*LAL_DAYSID_SI / tsteps;

      XLALGPSSetREAL8(&gps, T);

      XLALComputeDetAMResponse( &fplus, &fcross, detNSource.pDetector->response,
                                detNSource.pSource->equatorialCoords.longitude,
                                detNSource.pSource->equatorialCoords.latitude,
                                detNSource.pSource->orientation, XLALGreenwichMeanSiderealTime( &gps ) );

      gsl_matrix_set( LUfplus, i, j, fplus );
      gsl_matrix_set( LUfcross, i, j, fcross );
    }
  }
}


/** \brief Rescale the values output by the Nested Sampling algorithm
 *
 * This function reads in the file of samples output from the Nested Sampling algorithm (in the file specified by \c
 * outfile) and scales them back to their true values. It removes the string variable "model" from the output and shifts
 * the logPrior and logLikelihood values to the end of the parameter list.
 *
 * Note: The output may soon be in an XML format, so this function will need to be amended.
 *
 * \param runState [in] The analysis information structure
 */
void rescaleOutput( LALInferenceRunState *runState ){
  /* Open original output output file */
  CHAR *outfile, outfiletmp[256] = "";
  CHAR outfilepars[256] = "", outfileparstmp[256] = "";
  FILE *fp = NULL, *fptemp = NULL, *fppars = NULL, *fpparstmp = NULL;

  LALStringVector *paramsStr = NULL;

  ProcessParamsTable *ppt1 = LALInferenceGetProcParamVal( runState->commandLine, "--outfile" );

  if( ppt1 ){
    outfile = ppt1->value;

    /* set temporary file for re-writing out samples */
    sprintf(outfiletmp, "%s_tmp", outfile);

    /* open output file */
    if( (fp = fopen(outfile, "r")) == NULL ){
      XLALPrintError("Error... cannot open output file %s.\n", outfile);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    /* open temporary output file for reading */
    if( (fptemp = fopen(outfiletmp, "w")) == NULL ){
      XLALPrintError("Error... cannot open temporary output file %s.\n", outfile);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    /* open file for printing out list of parameter names - this should already
      exist */
    sprintf(outfilepars, "%s_params.txt", outfile);
    if( (fppars = fopen(outfilepars, "r")) == NULL ){
      XLALPrintError("Error... cannot open parameter name output file %s.\n", outfilepars);
      XLAL_ERROR_VOID(XLAL_EIO);
    }
    /* read in the parameter names and remove the "model" value */
    sprintf(outfileparstmp, "%s_params.txt_tmp", outfile);
    if( (fpparstmp = fopen(outfileparstmp, "w")) == NULL ){
      XLALPrintError("Error... cannot open parameter name output file %s.\n", outfileparstmp);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    CHAR v[128] = "";
    while( fscanf(fppars, "%s", v) != EOF ){
      paramsStr = XLALAppendString2Vector( paramsStr, v );

      /* re-output everything but the "model" value to a temporary file */
      if( strcmp(v, "model") ) { fprintf(fpparstmp, "%s\t", v); }
    }

    fclose(fppars);
    fclose(fpparstmp);

    /* move the temporary file name to the standard outfile_param name */
    rename( outfileparstmp, outfilepars );

    while ( 1 ){
      UINT4 i = 0;

      /* scan through line, get value and reprint out scaled value to temporary file */
      for( i = 0; i < paramsStr->length; i++ ){
        CHAR scalename[VARNAME_MAX] = "";
        CHAR scaleminname[VARNAME_MAX] = "";
        REAL8 scalefac = 1., scalemin = 0.;
        CHAR value[128];

        if( fscanf(fp, "%s", value) == EOF ) break;

        sprintf(scalename, "%s_scale", paramsStr->data[i]);
        sprintf(scaleminname, "%s_scale_min", paramsStr->data[i]);

        if ( LALInferenceCheckVariable( runState->data->dataParams, scalename ) &&
          LALInferenceCheckVariable( runState->data->dataParams, scaleminname ) ){
          scalefac = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, scalename );
          scalemin = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, scaleminname );

          fprintf(fptemp, "%.12le\t", atof(value)*scalefac + scalemin);
        }
        else if( strcmp(paramsStr->data[i], "model") ){
          fprintf(fptemp, "%.12le\t", atof(value));
        }
      }

      if( feof(fp) ) break;

      /* print out the last two items to be the logPrior and logLikelihood */
      fprintf(fptemp, "\n");
    }

    fclose(fp);
    fclose(fptemp);

    XLALDestroyStringVector( paramsStr );

    /* move the temporary file name to the standard outfile name */
    rename( outfiletmp, outfile );
  }
/* if we have XML enabled */
#ifdef HAVE_LIBLALXML
  ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal( runState->commandLine, "--outXML" );
  LALInferenceVariables *output_array = NULL;
  UINT4 N_output_array = 0, i = 0;
  CHAR *outVOTable = NULL;

  if ( !ppt2 && !ppt1 ){
    XLALPrintError("Must specify either --outfile or --outXML\n");
    XLAL_ERROR_VOID( XLAL_EIO );
  }

  /* rescale parameters held in array and recreate XML output - we don't need
     to remove any variables. */
  if( ppt2 ){
    outVOTable = ppt2->value;

    if( LALInferenceCheckVariable(runState->algorithmParams,"outputarray")
      && LALInferenceCheckVariable(runState->algorithmParams,"N_outputarray")){
      output_array = *(LALInferenceVariables **)LALInferenceGetVariable( runState->algorithmParams, "outputarray" );
      N_output_array = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams, "N_outputarray" );
    }

    /* loop through output array and rescale values accordingly */
    for( i = 0; i < N_output_array; i++ ){
      LALInferenceVariableItem *scaleitem = NULL;

      scaleitem = output_array[i].head;

      /* loop through tmparr parameters and scale if necessary */
      for( ; scaleitem; scaleitem = scaleitem->next ){
        CHAR scalename[VARNAME_MAX] = "";
        CHAR scaleminname[VARNAME_MAX] = "";
        REAL8 scalefac = 1., scalemin = 0., value = 0;

        sprintf(scalename, "%s_scale", scaleitem->name);
        sprintf(scaleminname, "%s_scale_min", scaleitem->name);

        /* check if scale values are present */
        if ( LALInferenceCheckVariable( runState->data->dataParams, scalename ) &&
          LALInferenceCheckVariable( runState->data->dataParams, scaleminname ) ){
          scalefac = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, scalename );
          scalemin = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, scaleminname );

          /* get the value and scale it */
          value = *(REAL8 *)LALInferenceGetVariable( &output_array[i], scaleitem->name );
          value = value*scalefac + scalemin;

          /* reset the value */
          LALInferenceSetVariable( &output_array[i], scaleitem->name, &value );

          /* change type to be REAL8 */
          scaleitem->type = LALINFERENCE_REAL8_t;
        }
      }
    }

    if( output_array && outVOTable && N_output_array > 0 ){
      xmlNodePtr votable = XLALInferenceVariablesArray2VOTTable( output_array, N_output_array, "Nested Samples");
      xmlNewProp( votable, CAST_CONST_XMLCHAR("utype"), CAST_CONST_XMLCHAR("lalinference:results:nestedsamples") );

      xmlNodePtr stateResource = XLALInferenceStateVariables2VOTResource(runState, "Run State Configuration");

      xmlNodePtr nestResource = XLALCreateVOTResourceNode("lalinference:results", "Nested sampling run", votable);

      if( stateResource ) { xmlAddChild( nestResource, stateResource ); }

      CHAR *xmlString = XLALCreateVOTStringFromTree( nestResource );

      /* Write to disk */
      if ( (fp = fopen(outVOTable, "w")) == NULL ){
        XLALPrintError("Error... can't open output XML file\n");
        XLAL_ERROR_VOID(XLAL_EIO);
      }

      fprintf(fp, "%s", xmlString);
      fclose(fp);
    }
  }

  if( output_array ) { XLALFree( output_array ); }
#else
  if ( !ppt1 ){
    XLALPrintError("Error... --outfile not defined!\n");
    XLAL_ERROR_VOID( XLAL_EIO );
  }
#endif

  return;
}


/** \brief Counts the number of comma separated values in a string
 *
 * This function counts the number of comma separated values in a given input string.
 *
 * \param csvline [in] Any string
 *
 * \return The number of comma separated value in the input string
 */
INT4 count_csv( CHAR *csvline ){
  CHAR *inputstr = NULL;
  INT4 count = 0;

  inputstr = XLALStringDuplicate( csvline );

  /* count number of commas */
  while(1){
    if( strsep(&inputstr, ",") == NULL ){
      XLALPrintError("Error... problem counting number of commas!\n");
      XLAL_ERROR( XLAL_EFUNC );
    }

    if ( inputstr == NULL ) { break; }

    count++;
  }

  return count+1;
}


/** \brief Checks if a given parameter is recognised
 *
 * This function checks whether a given parameter is one of the defined amplitude (\c amppars), frequency (\c freqpars),
 * sky location (\c skypars) or binary system (\c binpars) parameters given in the header file.
 *
 * \param parname [in] The name of a parameter
 *
 * \return true (1) if the parameter is recognised and false (0) if not
 */
INT4 recognised_parameter( CHAR *parname ){
  INT4 i = 0;

  for( i = 0; i < NUMAMPPARS; i++ ) { if (!strcmp(parname, amppars[i])) { return 1; } }

  for( i = 0; i < NUMFREQPARS; i++ ) { if (!strcmp(parname, freqpars[i])) { return 1; } }

  for( i = 0; i < NUMSKYPARS; i++ ) { if (!strcmp(parname, skypars[i])) { return 1; } }

  for( i = 0; i < NUMBINPARS; i++ ) { if (!strcmp(parname, binpars[i])) { return 1; } }

  return 0;
}


/** \brief Calculates the optimal matched filter signal-to-noise ratio for a given signal
 *
 * This function calculates the optimal matched filter signal-to-noise ratio (SNR) of a given signal model for a set of
 * detector data via:
 * \f[
 * \rho = \sqrt{\sum_{i=1}^N \frac{d_i^2}{\sigma^2}},
 * \f]
 * where \f$\{d\}\f$ is a time series of data, and \f$\sigma^2\f$ is its variance. As the data and model used here are
 * complex the real and imaginary SNRs are added in quadrature to give the total SNR.
 *
 * The data variance \f$\sigma^2\f$ is calculated on data that has had the running median subtracted in order to remove
 * any underlying trends (e.g. caused by a string signal). The variance is assumed constant over segments given in the
 * \c chunkLength vector and the SNR from each segment is added in quadrature.
 *
 * \param data [in] A data pointer containing detector data and the signal model
 *
 * \return The optimal matched filter signal-to-noise ratio
 */
REAL8 calculate_time_domain_snr( LALInferenceIFOData *data ){
  REAL8 snrval = 0., chunkLength;

  INT4 i = 0, j = 0, length = 0, cl = 0;

  COMPLEX16Vector *meddata = NULL;
  INT4 chunkMin = 0, count = 0;

  /* subtract a running median value from the data to remove any underlying
     trends (e.g. caused by a strong signal) */
  meddata = subtract_running_median( data->compTimeData->data );

  UINT4Vector *chunkLengths = NULL;
  chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams, "chunkLength" );
  chunkMin = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMin" );

  length = data->compTimeData->data->length;

  for ( i = 0; i < length; i+=chunkLength ){
    REAL8 snrcRe = 0., snrcIm = 0., variRe = 0., variIm = 0.;

    chunkLength = (REAL8)chunkLengths->data[count];

    /* skip section of data if its length is less than the minimum allowed
      chunk length */
    if( chunkLength < chunkMin ){
      count++;
      continue;
    }

    cl = i + (INT4)chunkLength;

    for( j = i ; j < cl ; j++ ){
      variRe += SQUARE( creal(meddata->data[j]) );
      variIm += SQUARE( cimag(meddata->data[j]) );

      /* calculate optimal signal power */
      snrcRe += SQUARE( creal(data->compModelData->data->data[j]) );
      snrcIm += SQUARE( cimag(data->compModelData->data->data[j]) );
    }

    variRe /= (chunkLength - 1.);
    variIm /= (chunkLength - 1.);

    /* add SNRs for each chunk in quadrature */
    snrval += ( snrcRe/variRe ) + ( snrcIm/variIm );

    count++;
  }

  snrval = sqrt( snrval );

  return snrval;
}


/** \brief Get the signal-to-noise ratio of the maximum likelihood signal
 *
 * The function uses the signal with the highest likelihood (which will be the final point in the live points array) and
 * calculates the optimal signal-to-noise ratio (SNR) for it. This is output to a file based on the \c outfile value,
 * but with \c _SNR appended to it. For multiple detector, and/or models with multiple data sets, the individual
 * detector/data set SNR values will be output, with the final value being the multi-detector SNR. If a fake signal has
 * been injected into the data this file will already contain the optimal SNR of the true signal.
 *
 * \param runState [in] The analysis information structure
 *
 * \sa calculate_time_domain_snr
 */
void get_loudest_snr( LALInferenceRunState *runState ){
  INT4 ndats = 0;
  INT4 Nlive = *(INT4 *)LALInferenceGetVariable( runState->algorithmParams, "Nlive" );
  REAL8 snrmulti = 0.;
  REAL8Vector *freqFactors = NULL;

  CHAR *snrfile = NULL;
  FILE *fpsnr = NULL;

  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;

  LALInferenceVariables *loudestParams = NULL;
  LALInferenceIFOData *data;

  loudestParams = calloc( 1, sizeof(LALInferenceVariables) );

  /* max likelihood point should have been sorted to be the final value */
  LALInferenceCopyVariables( runState->livePoints[Nlive-1], loudestParams );

  /* make sure that the signal model in runState->data is that of the loudest signal */
  REAL8 logLnew = runState->likelihood( loudestParams, runState->data, runState->templt );

  if ( logLnew != *(REAL8 *)LALInferenceGetVariable(
    runState->livePoints[Nlive-1], "logL" ) ){
    fprintf(stderr, "Error... maximum log likelihood problem!\n");
    exit(0);
  }

  LALInferenceClearVariables( loudestParams );

  /* setup output file */
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( !ppt ){
    fprintf(stderr, "Error... no output file specified!\n");
    exit(0);
  }

  snrfile = XLALStringDuplicate( ppt->value );
  snrfile = XLALStringAppend( snrfile, "_SNR" );

  /* append to previous injection SNR file if it exists */
  if( (fpsnr = fopen(snrfile, "a")) == NULL ){
    fprintf(stderr, "Error... cannot open output SNR file!\n");
    exit(0);
  }

  /* get SNR of loudest point and print out to file */
  data = runState->data;

  //FILE *fp = NULL;
  //fp = fopen("max_like_signal.txt", "w");

  fprintf(fpsnr, "# Recovered SNR\n");

  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "freqfactors" );

  while( data ){
    REAL8 snrval = calculate_time_domain_snr( data );

    //UINT4 length = data->compModelData->data->length;

    /* print out maxlikelihood template */
    //for ( UINT4 j=0; j < length; j++ ){
    //  fprintf(fp, "%lf\t%le\t%le\n",
    //          XLALGPSGetREAL8( &data->dataTimes->data[j] ),
    //          data->compModelData->data->data[j].re,
    //          data->compModelData->data->data[j].im );
    //}

    ndats++;

    snrmulti += SQUARE( snrval );

    /* print out SNR value */
    fprintf(fpsnr, "%s\t%.3lf\t%le\n", data->name, freqFactors->data[ndats%(INT4)freqFactors->length], snrval);

    data = data->next;
  }

  //fclose(fp);

  /* print out multi-detector/multi-datastream SNR value */
  if ( ndats > 1 ) { fprintf(fpsnr, "Coherent\t%le\n", sqrt(snrmulti)); }

  fclose( fpsnr );
}


/** \brief Automatically set the solar system ephemeris file based on environment variables and data time span
 *
 * This function will attempt to find Earth and Sun ephemeris files based on LAL environment variables (as set up by
 * <code> lalpulsar-user-env.(c)sh </code>) and a given start and end GPS time (presumably taken from the data that is
 * to be analysed). It requires \c LALPULSAR is installed and the \c LALPULSAR_PREFIX variable is set, which should mean
 * that ephemeris files are installed in the directory \c ${LALPULSAR_PREFIX}/share/lalpulsar/.
 *
 *
 * \param efile [in] a string that will return the Earth ephemeris file
 * \param sfile [in] a string that will return the Sun ephemeris file
 * \param tfile [in] a string that will return the time correction file
 * \param pulsar [in] the pulsar parameters read from a .par file
 * \param gpsstart [in] the GPS time of the start of the data
 * \param gpsend [in] the GPS time of the end of the data
 *
 * \return The TimeCorrectionType e.g. TDB or TCB
 */
TimeCorrectionType XLALAutoSetEphemerisFiles( CHAR *efile, CHAR *sfile, CHAR *tfile, BinaryPulsarParams  pulsar,
                                              INT4 gpsstart, INT4 gpsend ){
  /* set the times that the ephemeris files span */
  INT4 ephemstart = 630720013; /* GPS time of Jan 1, 2000, 00:00:00 UTC */
  INT4 ephemend = 1261872015; /* GPS time of Jan 1, 2020, 00:00:00 UTC */
  CHAR *eftmp = NULL, *sftmp = NULL, *tftmp = NULL;
  TimeCorrectionType ttype = TIMECORRECTION_NONE;

  CHAR *lalpath = NULL, *lalpulsarpath = NULL;

  if( gpsstart < ephemstart || gpsend < ephemstart || gpsstart > ephemend || gpsend > ephemend ){
    XLALPrintError("Start and end times are outside the ephemeris file ranges!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* first check that the path to the Ephemeris files is available in the
     environment variables */
  if((lalpath = getenv("LALPULSAR_PREFIX")) == NULL){
    XLALPrintError("LALPULSAR_PREFIX environment variable not set. Cannot automatically generate ephemeris files!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  lalpulsarpath = XLALStringDuplicate( lalpath );

  if ( (lalpulsarpath = XLALStringAppend(lalpulsarpath, "/share/lalpulsar/")) == NULL ) { XLAL_ERROR(XLAL_EFUNC); }

  eftmp = XLALStringDuplicate(lalpulsarpath);
  sftmp = XLALStringDuplicate(lalpulsarpath);
  tftmp = XLALStringDuplicate(lalpulsarpath);

  eftmp = XLALStringAppend(eftmp, "earth00-19-");
  sftmp = XLALStringAppend(sftmp, "sun00-19-");

  if( pulsar.ephem == NULL ){
    /* default to use DE405 */
    eftmp = XLALStringAppend(eftmp, "DE405");
    sftmp = XLALStringAppend(sftmp, "DE405");
  }
  else{
    if( !strcmp(pulsar.ephem, "DE405") || !strcmp(pulsar.ephem, "DE200") ||
        !strcmp(pulsar.ephem, "DE414") ){
      eftmp = XLALStringAppend(eftmp, pulsar.ephem);
      sftmp = XLALStringAppend(sftmp, pulsar.ephem);
    }
    else{
      XLALPrintError("Unknown ephemeris %s in par file\n", pulsar.ephem);
      XLAL_ERROR(XLAL_EFUNC);
    }
  }

  /* add .dat extension */
  eftmp = XLALStringAppend(eftmp, ".dat.gz");
  sftmp = XLALStringAppend(sftmp, ".dat.gz");

  if ( eftmp == NULL || sftmp == NULL ) { XLAL_ERROR(XLAL_EFUNC); }

  XLALStringCopy( efile, eftmp, strlen(eftmp)+1 );
  XLALStringCopy( sfile, sftmp, strlen(sftmp)+1 );

  /* double check that the files exist */
  if( fopen(sfile, "r") == NULL || fopen(efile, "r") == NULL ){
    XLALPrintError("Error... ephemeris files not, or incorrectly, defined!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  if( pulsar.units == NULL ){
    /* default to using TCB units */
    tftmp = XLALStringAppend(tftmp, "te405_2000-2019.dat.gz");
    ttype = TIMECORRECTION_TCB;
  }
  else{
    if ( !strcmp( pulsar.units, "TDB" ) ){
      tftmp = XLALStringAppend(tftmp, "tdb_2000-2019.dat.gz");
      ttype = TIMECORRECTION_TDB;
    }
    else{
      XLALPrintError("Error... unknown units %s in par file!\n", pulsar.units);
      XLAL_ERROR(XLAL_EFUNC);
    }
  }

  if ( tftmp == NULL ) { XLAL_ERROR(XLAL_EFUNC); }

  XLALStringCopy( tfile, tftmp, strlen(tftmp)+1 );

  if( fopen(tfile, "r") == NULL ){
    XLALPrintError("Error... time ephemeris files not, or incorrectly, defined!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  return ttype;
}


/*----------------------- END OF HELPER FUNCTIONS ----------------------------*/

/** \brief Read in an ascii text file of nested samples, convert to posterior samples and create k-d tree
 *
 * This function reads in an ascii file of nested samples, converted them into posterior samples and them add them to a
 * k-d tree. The file name containing the samples must have been given as the command line argument \c sample-file and
 * there must be an accompanying file with the names of each column with the same file name with _params.txt appended.
 *
 * It is assumed that the samples are in ascending log likelihood order. It is also assumed that variable values in the
 * file (and are not likelihood-like values) are consistent with those given that have prior ranges defined in the prior
 * file/par file (as these ranges will be used as bounds in a k-d tree created from this data).
 *
 * As it is assumed that the points read in are from a previous nested sampling run the number of live points used for
 * that run are also required to be given with the \c sample-nlive argument. This will be used during the conversion to
 * posterior samples.
 *
 * If given the k-d tree cell size for using the posterior as a prior can be set with the \c prior-cell argument, if not
 * set this defaults to 32.
 *
 * In the future this will be altered so as to also read in an XML file of samples.
 *
 * NOTE: add the ability to read in multiple files and combine the posterior samples
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void samples_prior( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL;

  UINT4 Ncell = 8; /* default prior cell size */

  UINT4 i = 0, k = 0, nsamps = 0, nnlive = 0, n = 0;
  UINT4Vector *nlive = NULL, *Nsamps = NULL;
  CHAR *nlivevals = NULL, *templives = NULL, *templive = NULL;

  CHAR *sampfile = NULL;
  CHAR *tempsamps = NULL, *tempsamp = NULL;

  LALStringVector *sampfilenames = NULL;

  LALInferenceVariables ***params = NULL;

  FILE *fp = NULL;

  const CHAR *fn = __func__;

  /* get names of nested sample file columns */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--sample-files" );
  if ( ppt != NULL ){
    sampfile = XLALStringDuplicate( ppt->value );

    /* count the number of sample files from the comma separated vales and set their names */
    tempsamps = XLALStringDuplicate( sampfile );

    nsamps = count_csv( tempsamps );

    for( i = 0; i < nsamps; i++ ){
      tempsamp = strsep( &tempsamps, "," );
      sampfilenames = XLALAppendString2Vector( sampfilenames, tempsamp );
    }
  }
  else return; /* no file so we don't use this function */

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--sample-nlives" );
  if ( ppt != NULL ){
    nlivevals = XLALStringDuplicate( ppt->value );

    templives = XLALStringDuplicate( nlivevals );

    nnlive = count_csv( templives );

    if( nnlive != nsamps ){
      XLALPrintError("%s: Number of live points not equal to number of posterior files!\n", fn, sampfile);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    for( i = 0; i < nnlive; i++ ){
      templive = strsep( &templives, "," );
      nlive = XLALResizeUINT4Vector( nlive, i+1 );
      nlive->data[i] = atoi( templive );
    }

    LALInferenceAddVariable( runState->algorithmParams, "numberlive", &nlive, LALINFERENCE_UINT4Vector_t,
                             LALINFERENCE_PARAM_FIXED );
  }
  else{
    fprintf(stderr, "Must set the number of live points used in the input nested samples file.\n\n");
    fprintf(stderr, USAGE, runState->commandLine->program);
    exit(0);
  }

  /* allocate memory for nested samples */
  params = XLALCalloc( nsamps, sizeof(LALInferenceVariables**) );

  /* loop over files, convert to posterior samples and combine them */
  for ( n = 0; n < nsamps; n++ ){
    CHAR *namefile = NULL, name[256];
    LALStringVector *paramNames = NULL;

    i = 0;

    /* initialise array as NULL */
    params[n] = NULL;

    namefile = XLALStringDuplicate( sampfilenames->data[n] );
    namefile = XLALStringAppend( namefile, "_params.txt" );

    /* check file exists */
    if ( fopen(namefile, "r") == NULL || fopen(sampfilenames->data[n], "r") == NULL ){
      XLALPrintError("%s: Cannot access either %s or %s!\n", fn, namefile, sampfilenames->data[n]);
      XLAL_ERROR_VOID(XLAL_EIO);
    }

    /* read in parameter names */
    fp = fopen( namefile, "r" );
    while( fscanf(fp, "%s", name) != EOF ){ paramNames = XLALAppendString2Vector( paramNames, name ); }

    fclose(fp);

    /* read in parameter values */
    fp = fopen( sampfilenames->data[n], "r" );
    while( !feof(fp) ){
      REAL8 ps[paramNames->length];
      UINT4 j = 0;

      for( j=0; j<paramNames->length; j++ ) { if( fscanf(fp, "%lf", &ps[j]) == EOF ) { break; } }

      if( feof(fp) ) { break; }

      /* dynamically allocate memory */
      params[n] = XLALRealloc( params[n], (i+1)*sizeof(LALInferenceVariables*) );
      params[n][i] = NULL;
      params[n][i] = XLALCalloc( 1, sizeof(LALInferenceVariables) );

      /* add variables */
      for( j=0; j<paramNames->length; j++ ){
        /* use vary type of this analyses parameters i.e. those set by the prior
           and par file, otherwise set the parameter to fixed */
        LALInferenceParamVaryType vary;

        if ( LALInferenceCheckVariable( runState->currentParams, paramNames->data[j] ) ){
          vary = LALInferenceGetVariableVaryType( runState->currentParams, paramNames->data[j] );
        }
        else { vary = LALINFERENCE_PARAM_FIXED; }

        LALInferenceAddVariable( params[n][i], paramNames->data[j], &ps[j], LALINFERENCE_REAL8_t, vary );
      }

      i++;
    }

    /* check that non-fixed, or output parameters actually do vary, otherwise
      complain */
    LALInferenceVariableItem *item1 = params[n][0]->head;

    while ( item1 ){
      UINT4 allsame = 0;

      for ( k=1; k<i; k++ ){
        LALInferenceVariableItem *item2 = LALInferenceGetItem( params[n][k], item1->name );

        if( item1->vary != LALINFERENCE_PARAM_FIXED && item1->vary != LALINFERENCE_PARAM_OUTPUT ){
          if ( *(REAL8*)item1->value != *(REAL8*)item2->value ) { allsame++; }
        }
      }

      if( ( item1->vary != LALINFERENCE_PARAM_FIXED && item1->vary != LALINFERENCE_PARAM_OUTPUT ) && allsame == 0 ){
        XLALPrintError("%s: Apparently variable parameter %s does not vary!\n", fn, item1->name );
        XLAL_ERROR_VOID(XLAL_EFUNC);
      }

      item1 = item1->next;
    }

    Nsamps = XLALResizeUINT4Vector( Nsamps, n+1 );
    Nsamps->data[n] = i;
  }

  LALInferenceAddVariable( runState->algorithmParams, "nestedsamples", &params, LALINFERENCE_void_ptr_t,
                           LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( runState->algorithmParams, "Nsamps", &Nsamps, LALINFERENCE_UINT4Vector_t,
                           LALINFERENCE_PARAM_FIXED );

  /* get cell size */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--prior-cell" );
  if ( ppt != NULL ) { Ncell = atoi( ppt->value ); }

  LALInferenceAddVariable( runState->priorArgs, "kDTreePriorNcell", &Ncell, LALINFERENCE_UINT4_t,
                           LALINFERENCE_PARAM_FIXED );

  /* convert samples to posterior */
  ns_to_posterior( runState );

  /* create k-d tree of the samples for use as a prior */
  create_kdtree_prior( runState );
}
