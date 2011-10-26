/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Parameter estimation code for known pulsar searches using the nested
 * sampling algorithm.

\heading{Description}
This code is used to perform parameter estimation and evidence calculation in
targeted/semi-targeted searches for gravitational waves from known pulsars. It
may also be used to follow-up on signal candidates from semi-coherent all-sky
searches for unknown sources.

It uses the Bayesian technique of 'Nested Sampling' to sample over a defined
prior parameter space (unknown signal parameters such as the gravitational wave
amplitude). These samples can then be used to create posterior probability
density functions of the unknown parameters. The samples are also used to
calculate the Bayesian evidence, also known as the marginal likelihood, for a
given signal model. This can be compared with other models, in particular the
model that the data is described by Gaussian noise alone.

As input the code requires time domain data that has been heterodyned using the
known (or close to) phase evolution of the pulsar. The time domain input should
consist of a three column text file containing the GPS time stamp of the data
point, the real part of the heterodyned data and the imaginary part of the
heterodyned data, e.g.
900000000.000000  1.867532e-24  -7.675329e-25
900000060.000000  2.783651e-24  3.654386e-25
...

Most commonly such data will have a sample rate of 1/60 Hz, giving a bandwidth
of the same amount.

The code also requires that you specify which parameters are to be searched
over, and the prior ranges over these. Any of the signal parameters can be
searched over, including frequency, sky position and binary system parameters,
although the bandwidth of the data and search efficiency need to be taken into
account.

The 'Nested Sampling' algorithm (developed by [\ref Skilling2006]) used is that
defined in LALinferenceNestedSampler (see [\ref VeitchVecchio2010]). It is
essentially an efficient way to perform the integral
\f[
Z = \int^{\mathbf{\theta}} p(d|\mathbf{\theta}) p(\mathbf{\theta}) {\rm
d}\mathbf{\theta},
\f]
where \f$ \mathbf{\theta} \f$ is a vector of parameters, \f$
p(d|\mathbf{\theta}) \f$ is the likelihood of the data given the parameters,
and \f$ p(\mathbf{\theta}) \f$ is the prior on the parameters. It does this by
changing the multi-dimensional integral over N parameters into a
one-dimensional integral
\f[
Z = \int^X L(X) {\rm d}X \approx \sum_i L(X_i) \Delta{}X_i,
\f]
where \f$ L(X) \f$ is the likelihood, and \f$ X \f$ is the prior mass. The
algorithm will draw a number (\f$ N \f$) of samples (live points) from the
parameter priors, calculate the likelihood for each point and find the lowest
likelihood value. The lowest likelihood value will be added to the summation in
the above equation, with \f$ \log{\Delta{}X_i} \approx 1/N \f$ coming from the
fact that the prior would be normalised to unity and therefore each point should
occupy an equal fraction and at each iteration the prior volume will decrease
geometrically (for \f$ \log{\Delta{}X_0} = 0 \f$). A new point is then drawn
from the prior with the criterion that it has a higher likelihood than the
previous lowest point and substitutes that point. To draw the new point a Markov
Chain Monte Carlo (MCMC) procedure is used - there are two methods used to
sample points within this: i) drawing from a proposal distributions based on the
covariance matrix if the current live points (although to keep things
computationally efficient this no updated at every iteration), ii) picking a
point via differential evolution (two random live points are selected and a new
point half way between the two is created). The probability of using either
method is currently set at 80\% and 20\% respectively. The procedure is
continued until a stopping criterion is reached, which in this case is that the
remaining prior volume is less than the \c tolerance value set (see below). The
implementation of this can be seen in [\ref VeitchVecchio2010].

\heading{Usage}
The usage format is given below and can also be found by running the code with
\code
lalapps_pulsar_parameter_estimation_nested --help
\endcode

An example of running the code on to search over the four unknown parameters
\f$ h_0 \f$, \f$ \phi_0 \f$, \f$ \psi \f$ and \f$ \cos{\iota} \f$, for pulsar
J0534-2200, given heterodyned time domain data from the H1 detector in the file
\c finehet_J0534-2200_H1, is:
\code
lalapps_pulsar_parameter_estimation_nested --detectors H1 --par-file
J0534-2200.par --input-files finehet_J0534-2200_H1 --outfile ns_J0534-2200
--prop-file prop_J0534-2200.txt --ephem-earth
lscsoft/share/lalpulsar/earth05-09.dat --ephem-sun
lscsoft/share/lalpulsar/sun05-09.dat --model-type triaxial --Nlive 1000 --Nmcmc
100 --Nruns 1 --tolerance 0.25
\endcode
The \c par-file is a TEMPO-style file containing the parameters of the
pulsar used to perform the heterodyne (the frequency parameters are the
rotation frequency and therefore not necessarily the gravitational wave
frequency) e.g.
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
The \c prop-file is a text file containing a list of the parameters to be
searched over and their given upper and lower prior ranges e.g.
\code
h0 0 1e-21
phi0 0 6.283185307179586
cosiota -1 1
psi -0.785398163397448 0.785398163397448
\endcode
Note that if searching over frequency parameters the ranges specified in the 
\c prop-file should be given in terms of the pulsars rotation frequency and not
necessarily the gravitational wave frequency e.g. for a triaxial star emitting
gravitational waves at 100 Hz (which will be at twice the rotation frequency)
if you wanted to search over 99.999 to 100.001 Hz then you should used
\code
f0 49.9995 50.0005
\endcode

An example of running the code as above, but this time on fake data created
using the Advanced LIGO design noise curves and with a signal injected into the
data is:
\code
lalapps_pulsar_parameter_estimation_nested --fake-data AH1 --inject-file
fake.par --par-file fake.par --outfile ns_fake --prop-file prop_fake.txt
--ephem-earth lscsoft/share/lalpulsar/earth05-09.dat --ephem-sun
lscsoft/share/lalpulsar/sun05-09.dat --model-type triaxial --Nlive 1000 --Nmcmc
100 --Nruns 1 --tolerance 0.25
\endcode
In this case the \c inject-file parameter file must contain the values of \c
h0, \c phi0, \c psi and \c cosiota, otherwise these will be set to zero by
default. The parameter files given for \c inject-file and \c par-file do not
have to be the same - the injection can be offset from the 'heterodyned' values,
which will be reflected in the data. If an \c inject-output file is also
specified then the fake data containing the signal, and a fake signal-only data
set, will be output.
 */

#include "pulsar_parameter_estimation_nested.h"

/** The inverse of the factorials of the numbers 0 to 6. */
static const REAL8 inv_fact[7] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0),
(1.0/120.0), (1.0/720.0) };

/** The maximum number of different detectors allowable. */
#define MAXDETS 6

RCSID("$Id$");

/* global variables */
/** The flag to specify verbose output. */
INT4 verbose = 0;

/** An array to contain the log of factorials up to a certain number. */
REAL8 *logfactorial = NULL;

/** A flag to specify if phase parameters are being searched over and
 * therefore the pulsar model requires phase evolution to be re-calculated (0 =
 * no, 1 = yes). */ 
UINT4 varyphase = 0; 

/** A flag to specify if the sky position will be searched over, and therefore
 * whether the solar system barycentring needs recalculating (0 = no, 1 = yes).
*/
UINT4 varyskypos = 0; 

/** A flag to specify if the binary system parameters will be searched over,
 * and therefore whether the binary system barycentring needs recalculating (0 =
 * no, 1 = yes) */
UINT4 varybinary = 0; 

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/** The usage format for the code.  */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --detectors         all IFOs with data to be analysed e.g. H1,H2\n\
                     (delimited by commas) (if generating fake data these\n\
                     should not be set)\n"\
" --par-file          pulsar parameter (.par) file (full path) \n"\
" --input-files       full paths and file names for the data for each\n\
                     detector in the list (must be in the same order)\n\
                     delimited by commas. If not set you can generate fake\n\
                     data (see --fake-data below)\n"\
" --outfile           name of output data file\n"\
" --chunk-min         (INT4) minimum stationary length of data to be used in\n\
                     the likelihood e.g. 5 mins\n"\
" --chunk-max         (INT4) maximum stationary length of data to be used in\n\
                     the likelihood e.g. 30 mins\n"\
" --psi-bins          (INT4) no. of psi bins in the time-psi lookup table\n"\
" --time-bins         (INT4) no. of time bins in the time-psi lookup table\n"\
" --prop-file         file containing the parameters to search over and their\n\
                     upper and lower ranges\n"\
" --ephem-earth       Earth ephemeris file\n"\
" --ephem-sun         Sun ephemeris file\n"\
" --model-type        (CHAR) the signal model that you want to use. Currently\n\
                     this can only be 'triaxial', which is the default\n\
                     value.\n"\
"\n"\
" Nested sampling parameters:\n"\
" --Nlive             (INT4) no. of live points for nested sampling\n"\
" --Nmcmc             (INT4) no. of for MCMC used to find new live points\n"\
" --Nruns             (INT4) no. of parallel runs\n"\
" --tolerance         (REAL8) tolerance of nested sampling integrator\n"\
" --randomseed        seed for random number generator\n"\
"\n"\
" Signal injection parameters:\n"\
" --inject-file       a pulsar parameter (par) file containing the parameters\n\
                     of a signal to be injected. If this is given a signal\n\
                     will be injected\n"\
" --inject-output     a filename to with the injected signal will be\n\
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
"\n"

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


INT4 main( INT4 argc, CHAR *argv[] ){
  ProcessParamsTable *param_table;
  LALInferenceRunState runState;
  REAL8 logZnoise = 0.;
  
  REAL8Vector *logLikelihoods = NULL;
  
  /* set error handler to abort in main function */
  lalDebugLevel = 7;
  XLALSetErrorHandler(XLALAbortErrorHandler);
  
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
  runState.template = get_pulsar_model;
 
  /* Generate the lookup tables and read parameters from par file */
  setupFromParFile( &runState );
  
  /* add injections if requested */
  injectSignal( &runState );
  
  /* create sum square of the data to speed up the likelihood calculation */
  sumData( &runState );
  
  /* Initialise the prior distribution given the command line arguments */
  /* Initialise the proposal distribution given the command line arguments */
  initialiseProposal( &runState );
  
  gridOutput( &runState );
  
  runState.proposal = LALInferenceProposalPulsarNS;
  
  /* get noise likelihood and add as variable to runState */
  logZnoise = noise_only_model( runState.data );
  LALInferenceAddVariable( runState.algorithmParams, "logZnoise", &logZnoise, 
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
               
  /* Create live points array and fill initial parameters */
  LALInferenceSetupLivePointsArray( &runState );
  
  logLikelihoods = *(REAL8Vector **)
    LALInferenceGetVariable( runState.algorithmParams, "logLikelihoods" );
  
  /* Call the nested sampling algorithm */
  runState.algorithm( &runState );
  printf("Done nested sampling algorithm\n");

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
 * Memory is allocated for the parameters, priors and proposals. The nested
 * sampling control parameters are set: the number of live points \c Nlive, the
 * number of points for each MCMC \c Nmcmc, the number of independent runs
 * within the algorithm \c Nruns, and the stopping criterion \c tolerance.
 *
 * The random number generator is initialise (the GSL Mersenne
 * Twister algorithm \c gsl_rng_mt19937) using either a user defined seed \c
 * randomseed, the system defined \c /dev/random file, or the system clock
 * time. 
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 */ 
void initialiseAlgorithm( LALInferenceRunState *runState )
/* Populates the structures for the algorithm control in runState, given the
 commandLine arguments. Includes setting up a random number generator.*/
{
  ProcessParamsTable *ppt = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;
  REAL8 tmp;
  INT4 tmpi;
  INT4 randomseed;
        
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
    verbose = 1;
    LALInferenceAddVariable( runState->algorithmParams, "verbose", &verbose , 
                             LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  }

  /* Number of live points */
  tmpi = atoi( LALInferenceGetProcParamVal(commandLine, "--Nlive")->value );
  LALInferenceAddVariable( runState->algorithmParams,"Nlive", &tmpi, LALINFERENCE_INT4_t, 
                           LALINFERENCE_PARAM_FIXED );
        
  /* Number of points in MCMC chain */
  tmpi = atoi( LALInferenceGetProcParamVal(commandLine, "--Nmcmc")->value );
  LALInferenceAddVariable( runState->algorithmParams, "Nmcmc", &tmpi, LALINFERENCE_INT4_t, 
                           LALINFERENCE_PARAM_FIXED );

  /* Optionally specify number of parallel runs */
  ppt = LALInferenceGetProcParamVal( commandLine, "--Nruns" );
  if(ppt) {
    tmpi = atoi( ppt->value );
    LALInferenceAddVariable( runState->algorithmParams, "Nruns", &tmpi, LALINFERENCE_INT4_t,
                             LALINFERENCE_PARAM_FIXED );
  }
        
  /* Tolerance of the Nested sampling integrator */
  ppt = LALInferenceGetProcParamVal( commandLine, "--tolerance" );
  if( ppt ){
    tmp = strtod( ppt->value, (char **)NULL );
    LALInferenceAddVariable( runState->algorithmParams, "tolerance", &tmp, 
                             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  }
        
  /* Set up the random number generator */
  gsl_rng_env_setup();
  runState->GSLrandom = gsl_rng_alloc( gsl_rng_mt19937 );
        
  /* (try to) get random seed from command line: */
  ppt = LALInferenceGetProcParamVal( commandLine, "--randomseed" );
  if ( ppt != NULL )
    randomseed = atoi( ppt->value );
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

/** \brief Reads in the input gravitational wave data files, or creates fake
 * data sets.
 * 
 * The function will check whether data files are being input of fake data is
 * to be generated. If using real data the \c detectors command line input must
 * list the names of the detectors from which each of the data sets comes, with
 * names separated by commas - allowed detector names are H1, H2, L1, V1, G1
 * or T1, for LIGO Hanford Observatory 1, LIGO Hanford Observatory 2, LIGO
 * Livingston Observatory, Virgo, GEO600/HF and TAMA300 respectively. If
 * requiring fake data to be generated then the \c fake-data command line
 * argument must list the detectors for which fake data is required (again
 * separated by commas) - these can be the same names as above, although if an
 * 'A' is appended to the LIGO of Virgo detector name then the Advanced
 * detector is assumed (for use if generating data from designed sensitivity
 * curves).
 * 
 * If generating fake data the noise floor can either be created using models
 * of the detector design sensitivities (if just \c fake-data is set), or with
 * noise levels defined using the \c fake-psd command line argument. If using \c
 * fake-psd it should list the signle-sided power spectral density required for
 * each detector listed in \c fake-data (again separated by commas). If using
 * the design sensitivity models the \c par-file argument will be used to find
 * the noise at the correct frequency, which is here assumed to be twice the
 * rotation frequency. The start time (in GPS seconds), data length (in seconds)
 * and sample interval (in seconds) can be specified for each fake data set by
 * supplying comma separated lists to the \c fake-starts, \c fake-lengths and \c
 * fake-dt command line arguments. By default these values are GPS 900000000 
 * (13th July 2008 at 15:59:46 UTC), 86400 seconds (1 solar day) and 60 seconds
 * (1/60 Hz sample rate) respectively. The fake data is drawn from a normal
 * distribution using \c XLALNormalDeviates and scaled with the appropriate PSD.
 * 
 * The number of data files required to be read in, or number of fake data sets
 * generated depends on the pulsar model type, specified by the \c model-type 
 * command line argument. At the moment two models are available (although this
 * can be expanded):
 * <ol>
 * <li>the triaxial model - this assumes that for each specified detector there
 * is <b>one</b> input file containing data heterodyned at twice the rotation
 * frequency of the pulsar (see e.g. that defined in [\ref DupuisWoan2005].</li>
 * <li>the pinned superfluid model - this assumes for each specified detector
 * there are two input files containing data heterodyned and the rotation 
 * frequency <i>and</i> twice the rotation frequency (see 
 * [\ref Jones2010]).</li>
 * </ol> 
 * The default model, if none is specified, is the triaxial model.
 * 
 * If creating fake data for a specific model then the number of PSDs, start
 * time, lengths and sample intervals specified must be equivalent to the number
 * of input files that would have been required e.g. if using the pinned 
 * superfluid model and requiring data for H1 and L1 then four fake PSDs would
 * be required (the first pair at the pulsars rotation frequency and twice that
 * in H1, and the seconds pair at the pulsars rotation frequency and twice that
 * in L1). These most be specified in the same order as the detectors.
 * 
 * Any new models added can require and arbitrary number of inputs for a given
 * detector, however the heterodyned frequencies of each input must be 
 * hardcoded into the function.
 * 
 * If using real data the files must be specified in the \c input-files command
 * line argument - these should be comma separated for mulitple files and be in
 * the same order at the associated detector from which they came given by the
 * \c detectors command.
 * 
 * The function also check that valid Earth and Sun ephemeris files (from the
 * lalpulsar suite) are set with the \c ephem-earth and \c ephem-sun arguments,
 * and that a valid output file for the nested samples is set via the \c outfile
 * argument.
 * 
 * The \c log_factorial array is also filled in with values of the log of the
 * factorial of all numbers up to the maximum length of the data. This is used
 * in likelihood calculations.
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void readPulsarData( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL, *ppt2 = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;
  
  CHAR *detectors = NULL;
  CHAR *outfile = NULL;
  CHAR *inputfile = NULL;
  
  CHAR *filestr = NULL;
  
  CHAR efile[1024];
  CHAR sfile[1024];
  
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
  
  CHAR *modeltype = NULL;
  REAL8Vector *modelFreqFactors = NULL;
  INT4 ml = 1;
  
  runState->data = NULL;
  
  /* check pulsar model required */
  ppt = LALInferenceGetProcParamVal( commandLine, "--model-type" );
  if ( ppt ){
    modeltype = XLALStringDuplicate( ppt->value );
  }
  else{
    /* default model type is triaxial */
    modeltype = XLALStringDuplicate( "triaxial" );
  }
  
  /* set parameters to read in data for different models - if adding a new
     model different parameters may be required. Different models may use sets
     of data at multiple frequencies. */
  modelFreqFactors = XLALCreateREAL8Vector( ml );
  if( !strcmp( modeltype, "triaxial" ) ){
    /* this model requires 1 set of data for each detector at a frequency
       of twice the pulsars rotation rate */
    modelFreqFactors->data[0] = 2.;
  }
  else if( !strcmp( modeltype, "pinsf" ) ){
    /* this model requires 2 sets of data for each detector at a frequency
       of once and twice the pulsars rotation rate */
    ml = 2;
    modelFreqFactors = XLALResizeREAL8Vector( modelFreqFactors, ml );
    modelFreqFactors->data[0] = 1.;
    modelFreqFactors->data[1] = 2.;
  }
  /* ADD EXTRA MODELS HERE */
  else{
    fprintf(stderr, "Error... model type %s is unknown!\n", modeltype);
    exit(0);
  }
  
  /* allocate memory for data sample intervals (in seconds) */
  
  
  /* get the detectors - must */
  ppt = LALInferenceGetProcParamVal( commandLine, "--detectors" );
  ppt2 = LALInferenceGetProcParamVal( commandLine, "--fake-data" );
  if( ppt && !ppt2 ){
   detectors = XLALStringDuplicate( ppt->value );
      
    /* count the number of detectors from command line argument of comma
       separated vales and set their names */
    tempdets = XLALStringDuplicate( detectors );
    
    if ( (numDets = count_csv( detectors )) > MAXDETS ){
      fprintf(stderr, "Error... too many detectors specified. Increase MAXDETS\
 to be greater than %d if necessary.\n", MAXDETS);
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
      fprintf(stderr, "Error... too many detectors specified. Increase MAXDETS\
 to be greater than %d if necessary.\n", MAXDETS);
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
      
      /* count the number of PSDs (comma seperated values) to compare to number
         of detectors */
      if( (numPsds = count_csv( psds )) != ml*numDets ){
        fprintf(stderr, "Error... for model type \"%s\" number of PSDs for fake\
 data must be %d times the number of detectors specified (no. dets = %d)\n",
                modeltype, ml, numDets);
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
        
          if( (tmpstr = strstr(tempdet, "A")) != NULL ) /* have advanced */
            XLALStringCopy( dets[FACTOR(i,ml)],
              tmpstr+1, strlen(tmpstr)+1 );
          else
            XLALStringCopy( dets[FACTOR(i,ml)],
              tempdet, strlen(tempdet)+1 );
        }
      }
    }
    /*------------------------------------------------------------------------*/
    else{ /* get PSDs from model functions and set detectors */
      CHAR *parFile = NULL;
      BinaryPulsarParams pulsar;
      REAL8 pfreq = 0.;
      
      ppt = LALInferenceGetProcParamVal( runState->commandLine, "--par-file" );
      if( ppt == NULL ) { 
        fprintf(stderr,"Must specify --par-file!\n"); exit(1);
      }
      parFile = XLALStringDuplicate( ppt->value );
        
      /* get the pulsar parameters to give a value of f */
      XLALReadTEMPOParFile( &pulsar, parFile );
      
      /* putting in pulsar frequency at f here */
      pfreq = pulsar.f0;
      
      tempdets = XLALStringDuplicate( detectors );
      
      for( i = 0; i < ml*numDets; i++ ){
        LALStatus status;
	
        CHAR *tmpstr = NULL;
        REAL8 psdvalf = 0.;
        
        numPsds++;
       
        if( i%ml == 0 ) tempdet = strsep( &tempdets, "," );
        
        if( (tmpstr = strstr(tempdet, "A")) != NULL ){ /* have Advanced */
          XLALStringCopy( dets[FACTOR(i,ml)], tmpstr+1, strlen(tmpstr)+1 );
          
          if( !strcmp(dets[FACTOR(i,ml)], "H1") || 
              !strcmp(dets[FACTOR(i,ml)], "L1") ||  
              !strcmp(dets[FACTOR(i,ml)], "H2") ){ /* ALIGO */
            LALAdvLIGOPsd( &status, &psdvalf, 
                           pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
            psdvalf *= 1.e-49; /* scale factor in LALAdvLIGOPsd.c */
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "V1") ){ /* AVirgo */
            LALEGOPsd( &status, &psdvalf,
                       pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else{
            fprintf(stderr, "Error... trying to use Advanced detector that is\
 not available!\n");
            exit(0);
          }
        }
        else{ /* initial detector */
          XLALStringCopy( dets[FACTOR(i,ml)], tempdet, strlen(tempdet)+1 );
          
          if( !strcmp(dets[FACTOR(i,ml)], "H1") || 
              !strcmp(dets[FACTOR(i,ml)], "L1") ||  
              !strcmp(dets[FACTOR(i,ml)], "H2") ){ /* Initial LIGO */
            psdvalf = 
              XLALLIGOIPsd( pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          
            /* divide H2 psds by 2 */
            if( !strcmp(dets[FACTOR(i,ml)], "H2") ){
              psdvalf /= 2.;
            }
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "V1") ){ /* Initial Virgo */
            LALVIRGOPsd( &status, &psdvalf,
                         pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "G1") ){ /* GEO 600 */
            LALGEOPsd( &status, &psdvalf,
                       pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
            psdvalf *= 1.e-46; /* scale factor in LALGEOPsd.c */
          }
          else if( !strcmp(dets[FACTOR(i,ml)], "T1") ){ /* TAMA300 */
            LALTAMAPsd( &status, &psdvalf,
                        pfreq*modelFreqFactors->data[(INT4)(fmod(i,ml))] );
            psdvalf *= 75. * 1.e-46; /* scale factor in LALTAMAPsd.c */
          }
          else{
            fprintf(stderr, "Error... trying to use detector that is\
 not available!\n");
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
        fprintf(stderr, "Error... for model type \"%s\" number of start times\
 for fake data must be %d times the number of detectors specified (no. dets =\
 %d)\n", modeltype, ml, numDets);
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
      for(i = 0; i < ml*numDets; i++ ){ 
        fstarts[i] = 900000000.;
      }
    }
      
    flengths = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-lengths");
    if( ppt ){
      CHAR *tmplengths = NULL, *tmplength = NULL, lengthval[256];
      fakelengths = XLALStringDuplicate( ppt->value );
      
      if( (numLengths = count_csv( fakelengths )) != numDets*ml ){
        fprintf(stderr, "Error... for model type \"%s\" number of data lengths\
 for fake data must be %d times the number of detectors specified (no. dets =\
 %d)\n", modeltype, ml, numDets);
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
      for(i = 0; i < ml*numDets; i++ ) flengths[i] = 86400.;
    }
    
    fdt = XLALCalloc( MAXDETS*ml, sizeof(REAL8) );
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-dt");
    if( ppt ){
      CHAR *tmpdts = NULL, *tmpdt = NULL, dtval[256];
      fakedt = XLALStringDuplicate( ppt->value );
      
      if( (numDt = count_csv( fakedt )) != ml*numDets ){
        fprintf(stderr, "Error... for model type \"%s\" number of sample time\
 steps for fake data must be %d times the number of detectors specified (no.\
 dets =\%d)\n", modeltype, ml, numDets);
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
      for(i = 0; i < ml*numDets; i++) fdt[i] = 60.;
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
  if( ppt ) inputfile = XLALStringDuplicate( ppt->value );
  
  if ( inputfile == NULL && !ppt2 ){
    fprintf(stderr, "Error... an input file or fake data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* get the output directory */
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( ppt ) outfile = XLALStringDuplicate( ppt->value );
  else{
    fprintf(stderr, "Error... --outfile needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
 
  /* count the number of input files (by counting commas) and check it's equal
     to twice the number of detectors */
  if ( !ppt2 ){ /* if using real data */
    count = count_csv( inputfile );
    
    if ( count != ml*numDets ){
      fprintf(stderr, "Error... for model type \"%s\" number of input files\
given must be %d times the number of detectors specified (no. dets =\%d)\n",
              modeltype, ml, numDets);
      exit(0);
    }
  }
  
  /* set random number generator in case when that fake data is used */
  ppt = LALInferenceGetProcParamVal( commandLine, "--randomseed" );
  if ( ppt != NULL ) seed = atoi( ppt->value );
  else seed = 0; /* will be set from system clock */
      
  /* reset filestr if using real data (i.e. not fake) */
  if ( !ppt2 ) filestr = XLALStringDuplicate( inputfile );
 
  /* read in data, needs to read in two sets of data for each ifo! */
  for( i = 0, prev=NULL ; i < ml*numDets ; i++, prev=ifodata ){
    CHAR *datafile = NULL;
    REAL8 times = 0;
    LIGOTimeGPS gpstime;
    COMPLEX16 dataVals;
    REAL8Vector *temptimes = NULL;
    INT4 j = 0, k = 0, datalength = 0;
    ProcessParamsTable *ppte = NULL, *ppts = NULL;
    
    FILE *fp = NULL;
    
    count = 0;
    
    /* initialise random number generator */
    /* Moved into det loop so same random seed can be used with */
    /* Different detector combos and still get same noise realisation */
    randomParams = XLALCreateRandomParams( seed+i );
  
    ifodata = XLALCalloc( 1, sizeof(LALInferenceIFOData) );
    ifodata->modelParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
    ifodata->modelDomain = LALINFERENCE_DOMAIN_TIME;
    ifodata->next = NULL;
    ifodata->dataParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
    
    /* add the pulsar model */
    LALInferenceAddVariable( ifodata->dataParams, "modeltype", &modeltype,
                             LALINFERENCE_string_t, LALINFERENCE_PARAM_FIXED );
    
    /* add data sample interval */
    /*LALInferenceAddVariable( ifodata->dataParams, "dt", &fdt[i],
                             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );*/
    
    /* add frequency factors variable */
    LALInferenceAddVariable( ifodata->dataParams, "freqfactors",
      &modelFreqFactors, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    
    if( i == 0 ) runState->data = ifodata;
    if( i > 0 ) prev->next = ifodata;

    /* set detector */
    ifodata->detector = XLALGetSiteInfo( dets[FACTOR(i,ml)] );
    
    /* set dummy initial time */
    gpstime.gpsSeconds = 0;
    gpstime.gpsNanoSeconds = 0;
    
    /* allocate time domain data - will be dynamically allocated as data read*/
    ifodata->compTimeData = NULL;
    ifodata->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                            1., &lalSecondUnit,
                                                            1 );
    
    /* allocate time domain model */
    ifodata->compModelData = NULL;
    ifodata->compModelData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                            1., &lalSecondUnit,
                                                            1 );
    
    /*============================ GET DATA ==================================*/
    /* get i'th filename from the comma separated list */
    if ( !ppt2 ){ /* if using real data read in from the file */
      datafile = strsep(&filestr, ",");
   
      /* open data file */
      if( (fp = fopen(datafile, "r")) == NULL ){
        fprintf(stderr, "Error... can't open data file %s!\n", datafile);
        exit(0);
      }

      j=0;

      /* read in data */
      temptimes = XLALCreateREAL8Vector( 1 );

      /* read in data */
      while(fscanf(fp, "%lf%lf%lf", &times, &dataVals.re, &dataVals.im) != EOF){
        /* dynamically allocate more memory */
        ifodata->compTimeData =
          XLALResizeCOMPLEX16TimeSeries( ifodata->compTimeData, 0, j+1 );

        ifodata->compModelData = 
          XLALResizeCOMPLEX16TimeSeries( ifodata->compModelData, 0, j+1 );
          
        temptimes = XLALResizeREAL8Vector( temptimes, j+1 );

        temptimes->data[j] = times;
        ifodata->compTimeData->data->data[j].re = dataVals.re;
        ifodata->compTimeData->data->data[j].im = dataVals.im;
      
        j++;
      }
    
      fclose(fp);
    
      datalength = j;
      
      /* allocate data time stamps */
      ifodata->dataTimes = NULL;
      ifodata->dataTimes = XLALCreateTimestampVector( datalength );
    
      /* fill in time stamps as LIGO Time GPS Vector */
      for ( k = 0; k<j; k++ )
        XLALGPSSetREAL8( &ifodata->dataTimes->data[k], temptimes->data[k] );
      
      ifodata->compTimeData->epoch = ifodata->dataTimes->data[0];
      ifodata->compModelData->epoch = ifodata->dataTimes->data[0];
      
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
      
      /* resize the data and model times series */
      ifodata->compTimeData =
        XLALResizeCOMPLEX16TimeSeries( ifodata->compTimeData, 0, datalength );
        
      ifodata->compModelData = 
        XLALResizeCOMPLEX16TimeSeries( ifodata->compModelData, 0, datalength );
        
      /* create data drawn from normal distribution with zero mean and unit
         variance */
      realdata = XLALCreateREAL4Vector( datalength );
      imagdata = XLALCreateREAL4Vector( datalength );
      
      XLALNormalDeviates( realdata, randomParams );
      XLALNormalDeviates( imagdata, randomParams );
      
      /* converts single sided psd into double sided psd, and then into a time
         domain noise standard deviation */
      psdscale = sqrt( ( fpsds[i] / 2.) / ( 2. * fdt[i] ) );
      
      /* create time stamps and scale data with the PSD */
      for( k = 0; k < datalength; k++ ){
        /* set time stamp */
        XLALGPSSetREAL8( &ifodata->dataTimes->data[k], 
                         fstarts[i] + fdt[i] * (REAL8)k );
        
        ifodata->compTimeData->data->data[k].re = (REAL8)realdata->data[k] *
          psdscale;
        ifodata->compTimeData->data->data[k].im = (REAL8)imagdata->data[k] *
          psdscale;
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
    if( ppte && ppts ){
      XLALStringCopy( efile, ppte->value, sizeof(efile) );
      XLALStringCopy( sfile, ppts->value, sizeof(sfile) );
    }
    else{ /* try getting files automatically */
      if( XLALAutoSetEphemerisFiles( efile, sfile,
        ifodata->dataTimes->data[0].gpsSeconds,
        ifodata->dataTimes->data[datalength-1].gpsSeconds ) ){
        fprintf(stderr, "Error... not been able to set ephemeris files!\n");
        exit(3);
      }
    }

    fprintf(stderr, "%s %s\n", efile, sfile);

    /* check ephemeris files exist and if not output an error message */
    if( access(sfile, F_OK) != 0 || access(efile, F_OK) != 0 ){
      fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
      exit(3);
    }
    
    /* set up ephemeris information */
    ifodata->ephem = XLALInitBarycenter( efile, sfile );
    XLALDestroyRandomParams( randomParams );
    
    /* get maximum data length */
    if ( ifodata->compTimeData->data->length > maxlen )
      maxlen = ifodata->compTimeData->data->length;
  }
  
  /* set global variable logfactorial */
  logfactorial = XLALCalloc( maxlen+1, sizeof(REAL8) );
  for ( i = 2; i < (INT4)(maxlen+1); i++ )
    logfactorial[i] = logfactorial[i-1] + log((REAL8)i);
  
  /* free memory */
  XLALFree( fdt );
  XLALFree( flengths );
  XLALFree( fstarts );
  XLALFree( fpsds );
}


/** \brief Reads in the parameters of the pulsar being searched for
 * 
 * This function reads in a pulsars parameters from the specified TEMPO-style 
 * .par file given by \c par-file using \c XLALReadTEMPOParFile. This file must
 * be specified and should contain at least the pulsars frequency, right
 * ascension and declination (any value not included will be zero by default).
 * The file should contain the parameters with which the detector data was
 * heterodyned, as these are used to produce a signal phase template based on
 * this assumption.
 * 
 * A example .par file may look like
 * \code
 * RA 12:31:56.17643
 * DEC 43:21:35.2531
 * F0 100.78634 1 0.00005
 * F1 2.34e-15
 * PEPOCH 54323.785634
 * \endcode
 * which shows several parameters mostly defined by the parameter name and a
 * parameter value. However, the \c F0 value contains 4 items. If a parameter
 * has a \c 1 as the third entry then it means that this was a parameter that
 * was fit by TEMPO with the entry after the \c 1 being the 1 standard
 * deviation error on that parameter. For parameters where an error is present
 * the code will attempt to search over that parameter using a Gaussian prior
 * defined by the 1\f$\sigma\f$ error value. Other parameters will be set as
 * fixed by default. These can be overridden by the proposal file values
 * described in \c initialiseProposal().
 * 
 * Based on the defined sky position defined in the par file a lookup table of
 * the detector antenna response over time and polarisation will be set by \c 
 * setupLookupTables().
 * 
 * The function \c add_initial_variables() is used to pass the parameter values
 * from the .par file to the algorithm.
 * 
 * Using the parameters from the .par file the phase template, including the
 * solar system and binary system barycentring time delays will be setup. These
 * define the phase template used to perform the initial heterodyne, which is
 * used as the reference in cases when phase parameters (other than the initial
 * phase) are being searched over.
 * 
 * Values used for scaling the parameters (to avoid dynamic range issues) are
 * initialised although will be set as default values. 
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 * 
 * \sa setupLookupTables
 * \sa add_initial_variables
 * \sa get_phase_model
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
  
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--par-file" );
  if( ppt == NULL ) { fprintf(stderr,"Must specify --par-file!\n"); exit(1); }
  char *parFile = ppt->value;
        
  /* get the pulsar parameters */
  XLALReadTEMPOParFile( &pulsar, parFile );
  psr.equatorialCoords.longitude = pulsar.ra;
  psr.equatorialCoords.latitude = pulsar.dec;
  psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  
  /* Setup lookup tables for amplitudes */
  setupLookupTables( runState, &psr );
  
  runState->currentParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
  
  scaletemp = XLALCalloc( 1, sizeof(LALInferenceVariables) );
  
  /* if no binary model set the value to "None" */
  if ( pulsar.model == NULL )
    pulsar.model = XLALStringDuplicate("None");
  
  /* Add initial (unchanging) variables for the model, initial (unity) scale
     factors, and any Gaussian priors defined from the par file */
  add_initial_variables( runState->currentParams, scaletemp,
                         runState->priorArgs, pulsar );
                                     
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
 
    
      dts = get_ssb_delay( pulsar, data->dataTimes, data->ephem, data->detector,
                           0. );
   
      bdts = get_bsb_delay( pulsar, data->dataTimes, dts );
    
      LALInferenceAddVariable( data->dataParams, "ssb_delays", &dts,
                              LALINFERENCE_REAL8Vector_t, 
                              LALINFERENCE_PARAM_FIXED );
    
      LALInferenceAddVariable( data->dataParams, "bsb_delays", &bdts,
                              LALINFERENCE_REAL8Vector_t, 
                              LALINFERENCE_PARAM_FIXED );

      phase_vector = get_phase_model( pulsar, data, freqFactors->data[j] );
      
      data->timeData = NULL;
      data->timeData = XLALCreateREAL8TimeSeries( "",
                                                  &data->dataTimes->data[0], 
                                                  0., 1., &lalSecondUnit,
                                                  phase_vector->length );
      
      for ( i=0; i<phase_vector->length; i++ )
        data->timeData->data->data[i] = phase_vector->data[i];
      
      /* add the scale factors from scaletemp into the data->dataParams
         structure */
      for( ; scaleitem; scaleitem = scaleitem->next ){
        LALInferenceAddVariable( data->dataParams, scaleitem->name, 
                                 scaleitem->value, scaleitem->type,
                                 scaleitem->vary );
      }
    
      data = data->next;
    }
  }
                                     
  return;
}


/** \brief Sets the time vs polarisation angle antenna reponse lookup table
 * 
 * This function sets up a lookup table in time vs. polarisation angle
 * \f$\psi\f$ for each detector from which data exists (either real or fake). 
 * The time ranges over one sidereal day and the polarisation angle range from
 * \f$-\pi/4\f$ to \f$\pi/4\f$. The number of bins for the grid over these two
 * parameters can be specified on the command line via \c time-bins and \c 
 * psi-bins respectively, but if these are not given then default values are
 * used. 
 * 
 * The function will also calls \c chop_n_merge() for each data set, which will 
 * split the data into chunks during which it can be considered Gaussian and
 * stationary. The command line arguments \c chunk-min and \c chunk-max can be
 * used to specify hardwire minimum and maximum lengths of chunks that are
 * allowable. By default the maximum chunk length is 0, which corresponds to no
 * maximum value being set.
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 * \param source [in] A pointer to a LALSource variable containing the source
 * location
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
  if( ppt ) chunkMin = atoi( ppt->value );
  else chunkMin = CHUNKMIN; /* default minimum chunk length */
    
  ppt = LALInferenceGetProcParamVal( commandLine, "--chunk-max" );
  if( ppt ) chunkMax = atoi( ppt->value );
  else chunkMax = CHUNKMAX; /* default maximum chunk length */
   
  LALInferenceIFOData *data = runState->data;

  gsl_matrix *LUfplus = NULL;
  gsl_matrix *LUfcross = NULL;

  REAL8 t0;
  LALDetAndSource detAndSource;
        
  ppt = LALInferenceGetProcParamVal( commandLine, "--psi-bins" );
  INT4 psiBins;
  if( ppt ) psiBins = atoi( ppt->value );
  else psiBins = PSIBINS; /* default psi bins */
                
  ppt = LALInferenceGetProcParamVal( commandLine, "--time-bins" );
  INT4 timeBins;
  if( ppt ) timeBins = atoi( ppt->value );
  else timeBins = TIMEBINS; /* default time bins */
    
  while(data){
    UINT4Vector *chunkLength = NULL;
    
    LALInferenceAddVariable( data->dataParams, "psiSteps", &psiBins, 
                             LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( data->dataParams, "timeSteps", &timeBins, 
                             LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    
    t0 = XLALGPSGetREAL8( &data->dataTimes->data[0] );
    detAndSource.pDetector = data->detector;
    detAndSource.pSource = source;
                
    LUfplus = gsl_matrix_alloc( psiBins, timeBins );

    LUfcross = gsl_matrix_alloc( psiBins, timeBins );
    
    response_lookup_table( t0, detAndSource, timeBins, psiBins, LUfplus, 
                           LUfcross );

    LALInferenceAddVariable( data->dataParams, "LU_Fplus", &LUfplus, 
                             LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( data->dataParams, "LU_Fcross", &LUfcross,
                             LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED );

    LALInferenceAddVariable( data->dataParams, "chunkMin", &chunkMin, 
                             LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( data->dataParams, "chunkMax", &chunkMax, 
                             LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );
    
    /* get chunk lengths of data */
    /* chunkLength = get_chunk_lengths( data, chunkMax ); */
    
    chunkLength = chop_n_merge( data, chunkMin, chunkMax );
    
    LALInferenceAddVariable( data->dataParams, "chunkLength", &chunkLength, 
                             LALINFERENCE_UINT4Vector_t, LALINFERENCE_PARAM_FIXED );

    data = data->next;
  }
        
  return;
}


/** \brief Set up all the allowed variables for a known pulsar search
 * 
 * This functions sets up all possible variables that are possible in a known
 * pulsar search. Parameter values read in from a .par file and passed in via
 * the \c pars variable will be set. For parameters for which the .par file
 * contained an uncertainty the \c priorArgs value will be set to give a
 * Gaussian prior. Scale factors will be initialised for all variables (so that
 * they exist) although they will be set to 1. This is all performed by the \c
 * add_variable_scale_prior() function.
 * 
 * \param ini [in] A pointer to a \c LALInferenceVariables type that will
 * be filled in with pulsar parameters
 * \param scaleFac [in] A pointer to a \c LALInferenceVariables type that will
 * be initialised to hold scale factors for each corresponding pulsar parameter
 * \param priorArgs [in] A pointer to a \c LALInferenceVariables type that will 
 * hold the prior type and ranges for each parameter to be searched over
 * \param pars [in] A \c BinaryPulsarParams type containing pulsar parameters
 * read in from a TEMPO-style .par file
 * 
 * \sa add_variable_scale_prior
 */
void add_initial_variables( LALInferenceVariables *ini, 
                            LALInferenceVariables *scaleFac,
                            LALInferenceVariables *priorArgs, 
                            BinaryPulsarParams pars ){ 
  /* include a scale factor of 1 scaling values and if the parameter file
     contains an uncertainty then set the prior to be Gaussian with the
     uncertainty as the standard deviation */
  
  /* amplitude model parameters */
  add_variable_scale_prior( ini, scaleFac, priorArgs, "h0", pars.h0,
                            pars.h0Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "phi0", pars.phi0,
                            pars.phi0Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "cosiota",
                            pars.cosiota, pars.cosiotaErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "psi", pars.psi,  
                            pars.psiErr );
  /* amplitude model parameters for pinned superfluid model*/
  add_variable_scale_prior( ini, scaleFac, priorArgs, "h1", pars.h1,  
                            pars.h1Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "lambda", pars.lambda,  
                            pars.lambdaErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "theta", pars.theta,  
                            pars.thetaErr );
    
  /* phase model parameters */
  
  /* frequency */
  add_variable_scale_prior( ini, scaleFac, priorArgs, "f0", pars.f0, 
                            pars.f0Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "f1", pars.f1, 
                            pars.f1Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "f2", pars.f2, 
                            pars.f2Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "f3", pars.f3, 
                            pars.f3Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "f4", pars.f4, 
                            pars.f4Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "f5", pars.f5, 
                            pars.f5Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "pepoch", pars.pepoch,
                            pars.pepochErr );
  
  /* sky position */
  add_variable_scale_prior( ini, scaleFac, priorArgs, "ra", pars.ra, 
                            pars.raErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "pmra", pars.pmra, 
                            pars.pmraErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "dec", pars.dec,
                            pars.decErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "pmdec", pars.pmdec,
                            pars.pmdecErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "posepoch", pars.posepoch,
                            pars.posepochErr );
  
  /* binary system parameters */
  LALInferenceAddVariable( ini, "model", &pars.model, LALINFERENCE_string_t, 
                           LALINFERENCE_PARAM_FIXED );
  
  add_variable_scale_prior( ini, scaleFac, priorArgs, "Pb", pars.Pb, 
                            pars.PbErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "e", pars.e, pars.eErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "eps1", pars.eps1,
                            pars.eps1Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "eps2", pars.eps2,
                            pars.eps2Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "T0", pars.T0, 
                            pars.T0Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "Tasc", pars.Tasc,
                            pars.TascErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "x", pars.x, pars.xErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "w0", pars.w0, 
                            pars.w0Err );

  add_variable_scale_prior( ini, scaleFac, priorArgs, "Pb2", pars.Pb2,
                            pars.Pb2Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "e2", pars.e2, 
                            pars.e2Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "T02", pars.T02,
                            pars.T02Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "x2", pars.x2, 
                            pars.x2Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "w02", pars.w02,
                            pars.w02Err );

  add_variable_scale_prior( ini, scaleFac, priorArgs, "Pb3", pars.Pb3,
                            pars.Pb3Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "e3", pars.e3, 
                            pars.e3Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "T03", pars.T03,
                            pars.T03Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "x3", pars.x3, 
                            pars.x3Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "w03", pars.w03,
                            pars.w03Err );

  add_variable_scale_prior( ini, scaleFac, priorArgs, "xpbdot", pars.xpbdot,
                            pars.xpbdotErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "eps1dot", pars.eps1dot,
                            pars.eps1dotErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "eps2dot", pars.eps2dot,
                            pars.eps2dotErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "wdot", pars.wdot,
                            pars.wdotErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "gamma", pars.gamma,
                            pars.gammaErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "Pbdot", pars.Pbdot,
                            pars.PbdotErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "xdot", pars.xdot,
                            pars.xdotErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "edot", pars.edot,
                            pars.edotErr );
 
  add_variable_scale_prior( ini, scaleFac, priorArgs, "s", pars.s, pars.sErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "dr", pars.dr, 
                            pars.drErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "dth", pars.dth,
                            pars.dthErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "a0", pars.a0, 
                            pars.a0Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "b0", pars.b0, 
                            pars.b0Err );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "M", pars.M, pars.MErr );
  add_variable_scale_prior( ini, scaleFac, priorArgs, "m2", pars.m2, 
                            pars.m2Err );
}


/** \brief Adds variables, scale factors and priors
 * 
 * This function adds a variable with a name and a value. If a non-zero standard
 * deviation value is given then a Gaussian prior will be set for that parameter
 * and it will be set as a parameter to be searched over by the algorithm
 * (using the \c LALINFERENCE_PARAM_LINEAR vary type) - the mean of the
 * Gaussian prior will be the given parameter value. If the standard deviation
 * is zero the parameter will be set as fixed (using the \c
 * LALINFERENCE_PARAM_FIXED vary type).
 * 
 * If any non-amplitude parameters have non-zero standard deviations (i.e. will
 * be searched over) then the \c varyphase variable will be set to 1. If any sky
 * position parameters have non-zero standard deviations then the \c varyskypos
 * variable will be set to 1. If any binary system parameters have non-zero
 * standard deviations then the \c varybinary variable will be set to 1. 
 * 
 * For all parameters a scale factor and scale minimum range will be set. These 
 * are just initialised to 1 and 0 respectively and will be set in \c
 * initialiseProposal for any parameters that require them.
 * 
 * \param var [in] Pointer to \c LALInferenceVariables type to contain
 * parameter information
 * \param scale [in] Pointer to \c LALInferenceVariables type to contain
 * parameter scaling information
 * \param prior [in] Pointer to \c LALInferenceVariables type to contain
 * prior range and type information
 * \param name [in] string containing the parameter name
 * \param value [in] the value of the parameter
 * \param sigma [in] the standard deviation of the parameter
 */
void add_variable_scale_prior( LALInferenceVariables *var, 
                               LALInferenceVariables *scale, 
                               LALInferenceVariables *prior, const CHAR *name, 
                               REAL8 value, REAL8 sigma ){
  REAL8 scaleVal = 1., scaleMin = 0.;
  LALInferenceParamVaryType vary;
  CHAR scaleName[VARNAME_MAX] = "", scaleMinName[VARNAME_MAX] = "";
  INT4 i = 0;
  
  /* if the sigma is non-zero then set a Gaussian prior */
  if ( sigma != 0. ){
    vary = LALINFERENCE_PARAM_LINEAR;
    INT4 isthere = 0;
    
    /* set the prior to a Gaussian prior with mean value and sigma */
    LALInferenceAddGaussianPrior( prior, name, (void *)&value, (void *)&sigma, 
                                  LALINFERENCE_REAL8_t );
    
    /* if the parameter is not one of the amplitude parameters then set
       global variable varyphase to 1, so that the phase evolution will be
       calculated */
    for ( i = 0; i < NUMAMPPARS; i++ ){
      if ( !strcmp(name, amppars[i]) ){
        isthere = 1;
        break;
      }
    }
    if ( !isthere ) varyphase = 1;
   
    /* check if there are sky position parameters that will be searched over */
    for ( i = 0; i < NUMSKYPARS; i++ ){
      if ( !strcmp(name, skypars[i]) ){
        varyskypos = 1;
        break;
      }
    }
    
    /* check if there are any binary parameters that will be searched over */
    for ( i = 0; i < NUMBINPARS; i++ ){
      if ( !strcmp(name, binpars[i]) ){
        varybinary = 1;
        break;
      }
    }
  }
  else vary = LALINFERENCE_PARAM_FIXED;
  
  /* add the variable */
  LALInferenceAddVariable( var, name, &value, LALINFERENCE_REAL8_t, vary );
  
  /* add the initial scale factor of 1 */
  sprintf( scaleName, "%s_scale", name );
  LALInferenceAddVariable( scale, scaleName, &scaleVal, LALINFERENCE_REAL8_t, 
                           LALINFERENCE_PARAM_FIXED );
  
  /* add initial scale offset of zero */
  sprintf( scaleMinName, "%s_scale_min", name );
  LALInferenceAddVariable( scale, scaleMinName, &scaleMin,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
}


/** \brief Sets up the parameters to be searched over and their prior ranges
 * 
 * This function sets up any parameters that you require the code searches over
 * (it will overwrite anything set via the .par file if required) and specifies
 * the prior range for each. This information is contained in a proposal file 
 * specified by the command line argument \c prop-file. This file should
 * contain three columns: the first has the name of a parameter to be searched
 * over; the second has the lower limit of that parameters prior range; and the
 * third has the upper limit of the prior range. E.g.
 * \code
 * h0 0 1e-21
 * phi0 0 6.6.283185307179586
 * cosiota -1 1
 * psi -0.785398163397448 0.785398163397448
 * \endcode
 * 
 * Any parameter specified in the file will have its vary type set to \c
 * LALINFERENCE_PARAM_LINEAR, (except \f$\phi_0\f$, which if it is defined
 * to have a prior covering \f$2\pi\f$ wraps around at the edges of its range
 * and has a \c LALINFERENCE_PARAM_CIRCULAR type). Parameters, and their priors,
 * which linear variable type are scaled such that parameter x with priors
 * in the range [a, b] will become (x - a) / (b - a). As such the new prior
 * ranges will cover from 0 to 1. The parameter scale factor 
 * is set the to value of (b - a) and the minimum scale range is
 * set to a - this allows the true parameter value to be reconstructed.
 * 
 * For parameters with Gaussian priors set from the .par file the scale factor
 * is applied differently, so as to give a Gaussian with zero mean and unit
 * variance.
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void initialiseProposal( LALInferenceRunState *runState )
{
  CHAR *propfile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  FILE *fp=NULL;
  
  CHAR tempPar[VARNAME_MAX] = "";
  REAL8 low, high;
  
  LALInferenceIFOData *data = runState->data;
  
  ppt = LALInferenceGetProcParamVal( commandLine, "--prop-file" );
  if( ppt ){
  propfile =
    XLALStringDuplicate( LALInferenceGetProcParamVal(commandLine,"--prop-file")->value );
  }
  else{
    fprintf(stderr, "Error... --prop-file is required.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* open file */
  if( (fp = fopen(propfile, "r")) == NULL ){
    fprintf(stderr, "Error... Could not open proposal file %s.\n", propfile);
    exit(3);
  }
  
  while(fscanf(fp, "%s %lf %lf", tempPar, &low, &high) != EOF){
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
    
    if( high < low ){
      fprintf(stderr, "Error... In %s the %s parameters ranges are wrongly \
set.\n", propfile, tempPar);
      exit(3);
    }
    
    sprintf(tempParScale, "%s_scale", tempPar);
    sprintf(tempParScaleMin, "%s_scale_min", tempPar);
    sprintf(tempParPrior, "%s_gaussian_mean", tempPar);

    tempVar = *(REAL8*)LALInferenceGetVariable( runState->currentParams, 
                                                tempPar );
    type = LALInferenceGetVariableType( runState->currentParams, tempPar );
    
    /* remove variable value */
    LALInferenceRemoveVariable( runState->currentParams, tempPar );
    
    /* if a Gaussian prior has already been defined (i.e. from the par file)
       remove this and overwrite with values from the propfile */
    if ( LALInferenceCheckVariable(runState->priorArgs, tempParPrior) )
      LALInferenceRemoveGaussianPrior( runState->priorArgs, tempPar );
    
    scale = high - low;
    scaleMin = low;
    
    /* don't rescale phi0, so that it varies over 0-2*pi - required by the 
       LALInferenceAngularVariance function when calculating the covariance 
       (unless the prop-file a phi0 range that is not 2*pi) */
    if( !strcmp(tempPar, "phi0") ){
      if ( scale/LAL_TWOPI > 0.99 && scale/LAL_TWOPI < 1.01 ){ 
        scale = 1.;
        scaleMin = 0.;
      }
    }
    
    /* if psi is covering the range -pi/4 to pi/4 scale it, so that it cover
       the 0 to 2pi range of a circular parameter */
    if( !strcmp(tempPar, "psi") ){
      if ( scale/LAL_PI_2 > 0.99 && scale/LAL_PI_2 < 1.01 ){ 
        scale = 0.25;
        scaleMin = -LAL_PI/4.;
      }
    }
    
    /* set the scale factor to be the width of the prior */
    while( datatemp ){
      scaleType = LALInferenceGetVariableType( datatemp->dataParams, 
                                               tempParScale );
      LALInferenceRemoveVariable( datatemp->dataParams, tempParScale );
      LALInferenceRemoveVariable( datatemp->dataParams, tempParScaleMin );
      
      LALInferenceAddVariable( datatemp->dataParams, tempParScale, &scale, 
                               scaleType, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( datatemp->dataParams, tempParScaleMin,
                               &scaleMin, scaleType, 
                               LALINFERENCE_PARAM_FIXED );
                               
      datatemp = datatemp->next;
    }
    
    /* scale variable and priors */
    tempVar = (tempVar - scaleMin) / scale;
    low = 0.;
    high = (high - scaleMin) / scale;
    
    /* re-add variable */    
    if( !strcmp(tempPar, "phi0") && scale == 1. ) 
      varyType = LALINFERENCE_PARAM_CIRCULAR;
    else if ( !strcmp(tempPar, "phi0") && scale == 0.25 ) 
      varyType = LALINFERENCE_PARAM_CIRCULAR;
    else varyType = LALINFERENCE_PARAM_LINEAR;
    
    LALInferenceAddVariable( runState->currentParams, tempPar, &tempVar, type,
                             varyType );
    
    /* Add the prior variables */
    LALInferenceAddMinMaxPrior( runState->priorArgs, tempPar, &low, 
                                &high, type );
    
    /* if there is a phase parameter defined in the proposal then set varyphase
       to 1 */
    for ( i = 0; i < NUMAMPPARS; i++ ){
      if ( !strcmp(tempPar, amppars[i]) ){
        isthere = 1;
        break;
      }
    }
    if ( !isthere ) varyphase = 1;
      
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
  
  /* check for any parameters with Gaussian priors and rescale to mean value */
  LALInferenceVariableItem *checkPrior = runState->currentParams->head;
  for( ; checkPrior ; checkPrior = checkPrior->next ){
    REAL8 scale = 0., scaleMin = 0;
    REAL8 tempVar;
    
    REAL8 mu, sigma;
    
    LALInferenceIFOData *datatemp = data;
    CHAR tempParPrior[VARNAME_MAX] = "";

    sprintf(tempParPrior,"%s_gaussian_mean",checkPrior->name);

    if ( LALInferenceCheckVariable(runState->priorArgs, tempParPrior) ){
      tempVar = *(REAL8 *)checkPrior->value;
      
      /* get the mean and standard deviation of the Gaussian prior */
      LALInferenceGetGaussianPrior( runState->priorArgs, checkPrior->name, 
                                    (void *)&mu, (void *)&sigma );
      
      /* set the scale factor to be the sigma value */
      scale = sigma;
      scaleMin = mu;
      tempVar = (tempVar - scaleMin) / scale;
      
      /* scale the parameter value and reset it */
      memcpy( checkPrior->value, &tempVar, 
              LALInferenceTypeSize[checkPrior->type] );
      
      mu -= scaleMin;
      sigma /= scale;
      
      /* remove the Gaussian prior values and reset as scaled values */
      LALInferenceRemoveGaussianPrior( runState->priorArgs, checkPrior->name );
      LALInferenceAddGaussianPrior( runState->priorArgs, checkPrior->name, 
                        (void *)&mu, (void *)&sigma, checkPrior->type );
        
      /* set scale factor in data structure */
      while( datatemp ){
        LALInferenceVariableType scaleType;
        CHAR tempParScale[VARNAME_MAX] = "";
        CHAR tempParScaleMin[VARNAME_MAX] = "";
        
        sprintf(tempParScale,"%s_scale",checkPrior->name);
        sprintf(tempParScaleMin,"%s_scale_min",checkPrior->name);
        
        scaleType = LALInferenceGetVariableType( datatemp->dataParams, 
                                                 tempParScale );
        LALInferenceRemoveVariable( datatemp->dataParams, tempParScale );
        LALInferenceRemoveVariable( datatemp->dataParams, tempParScaleMin );
        
        LALInferenceAddVariable( datatemp->dataParams, tempParScale, &scale, 
                                 scaleType, LALINFERENCE_PARAM_FIXED );
        LALInferenceAddVariable( datatemp->dataParams, tempParScaleMin,
                                 &scaleMin, scaleType, 
                                 LALINFERENCE_PARAM_FIXED );
      
        datatemp = datatemp->next;
      }
    }
  }
    
  return;
}


/*------------------- END INITIALISATION FUNCTIONS ---------------------------*/


/******************************************************************************/
/*                     LIKELIHOOD AND PRIOR FUNCTIONS                         */
/******************************************************************************/
/** \brief The log likelihood function
 * 
 * This function calculates natural logarithm of the likelihood of a signal
 * model (specified by a given set of parameters) given the data from a set of 
 * detectors.
 * 
 * The likelihood is the joint likelihood of chunks of data over which the noise
 * is assumed stationary and Gaussian. For each chunk a Gaussian likelihood for
 * the noise and data has been marginalised over the unknown noise standard 
 * deviation using a Jeffreys prior on the standard deviation. Given the
 * data consisting of independent real and imaginary parts this gives a
 * Students-t distribution for each chunk (of length \f$m\f$) with \f$m/2\f$
 * degrees of freedom:
 * \f[
 * p(\mathbf{\theta}|\mathbf{B}) = \prod_{j=1}^M \frac{(m_j-1)!}{2\pi^{m_j}} 
 * \left( \sum_{k=k_0}^{k_0+(m_j-1)} |B_k - y(\mathbf{\theta})_k|^2 
 * \right)^{-m_j},
 * \f]
 * where \f$\mathbf{B}\f$ is a vector of the complex data, 
 * \f$y(\mathbf{\theta})\f$ is the model for a set of parameters
 * \f$\mathbf{\theta}\f$, \f$M\f$ is the total number of independent data chunks
 * with lengths \f$m_j\f$ and \f$k_0 = \sum_{i=1}^j 1 + m_{i-1}\f$ (with
 * \f$m_0 = 0\f$) is the index of the first data point in each chunk. The
 * product of this for each detector will give the full joint likelihood. In
 * the calculation here the unnecessary proportionality factors are left out
 * (this would effect the actual value of the marginal likelihood/evidence, but
 * since we are only interested in evidence ratios/Bayes factors these
 * factors would cancel out anyway. See [\ref DupuisWoan2005] for a more
 * detailed description.
 *
 * In this function data in chunks smaller than a certain minimum length 
 * \c chunkMin are ignored.
 * 
 * \param vars [in] The parameter values
 * \param data [in] The detector data and initial signal phase template
 * \param get_model [in] The signal template/model function
 * 
 * \return The natural logarithm of the likelihood function
 */
REAL8 pulsar_log_likelihood( LALInferenceVariables *vars, 
                             LALInferenceIFOData *data,
                             LALInferenceTemplateFunction *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  UINT4 i = 0;
  CHAR *modeltype = NULL;/*need to check model type in this function*/
  
  modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );
  
  LALInferenceIFOData *datatemp1 = data, *datatemp2 = data;
  
  /* copy model parameters to data parameters */
  while( datatemp1 ){
    LALInferenceCopyVariables( vars, datatemp1->modelParams );
    datatemp1 = datatemp1->next;
  }
  
  /* get pulsar model */
  while( datatemp2 ){
    /*fprintf(bugtest,"getting model in log like func\n");*/
    get_model( datatemp2 );
    datatemp2 = datatemp2->next;
    /* If modeltype is pinsf need to advance data on to next, so this loop only
     runs once if there is only 1 det*/
    if ( !strcmp( modeltype, "pinsf" ) ) datatemp2 = datatemp2->next;
  }
  
  while ( data ){
    UINT4 j = 0, count = 0, cl = 0;
    UINT4 length = 0, chunkMin, chunkMax;
    REAL8 chunkLength = 0.;
    REAL8 logliketmp = 0.;

    REAL8 sumModel = 0., sumDataModel = 0.;
    REAL8 chiSquare = 0.;
    COMPLEX16 B, M;
  
    REAL8Vector *sumDat = NULL;
    UINT4Vector *chunkLengths = NULL;
    
    /*fprintf(bugtest,"calc log like for one set of data\n");*/
    
    sumDat = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, 
                                                       "sumData" );
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams,
                                                             "chunkLength" );
    chunkMin = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMin" );
    chunkMax = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMax" );
  
    length = data->compTimeData->data->length;
  
    for( i = 0 ; i < length ; i += chunkLength ){
      chunkLength = (REAL8)chunkLengths->data[count];
    
      /* skip section of data if its length is less than the minimum allowed
        chunk length */
      if( chunkLength < chunkMin ){
        count++;
        continue;
      }

      sumModel = 0.;
      sumDataModel = 0.;

      cl = i + (INT4)chunkLength;
    
      for( j = i ; j < cl ; j++ ){
        B.re = data->compTimeData->data->data[j].re;
        B.im = data->compTimeData->data->data[j].im;

        M.re = data->compModelData->data->data[j].re;
        M.im = data->compModelData->data->data[j].im;
      
        /* sum over the model */
        sumModel += M.re*M.re + M.im*M.im;
        
        /* sum over that data and model */
        sumDataModel += B.re*M.re + B.im*M.im;
        /*fprintf(bugtest,"B.re= %e, B.im= %e, M.re: %e, M.im: %e\n",B.re, B.im, M.re, M.im);*/
      }
 
      chiSquare = sumDat->data[count];
      chiSquare -= 2.*sumDataModel;
      chiSquare += sumModel;
      
      logliketmp -= chunkLength*log(chiSquare);
      
      count++;
    }
    loglike += logliketmp;
    data = data->next;
  }
  return loglike;
}


/** \brief The prior function
 *
 * This function calculates the natural logarithm of the prior for a set of
 * parameters. If the prior on a particular parameter is uniform over a given
 * range then \f$p(\theta) = 1/(\theta_{\rm max} - \theta_{\rm min})\f$. If the 
 * prior is Gaussian then the probability of that value given the mean and 
 * standard deviation of the Gaussian is calculated.
 * 
 * \param runState [in] A pointer to the LALInferenceRunState
 * \param params [in] The set of parameter values
 * 
 * \return The natural logarithm of the prior value for a set of parameters 
 */
REAL8 priorFunction( LALInferenceRunState *runState, LALInferenceVariables *params ){
  LALInferenceIFOData *data = runState->data;
  (void)runState;
  LALInferenceVariableItem *item = params->head;
  REAL8 min, max, mu, sigma, prior = 0, value = 0.;
  
  for(; item; item = item->next ){
    /* get scale factor */
    CHAR scalePar[VARNAME_MAX] = "";
    CHAR scaleMinPar[VARNAME_MAX] = "";
    CHAR priorPar[VARNAME_MAX] = "";
    REAL8 scale = 0., scaleMin = 0.;
    
    if( item->vary == LALINFERENCE_PARAM_FIXED || 
      item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }
    
    sprintf(scalePar, "%s_scale", item->name);
    scale = *(REAL8 *)LALInferenceGetVariable( data->dataParams, scalePar );
    
    sprintf(scaleMinPar, "%s_scale_min", item->name);
    scaleMin = *(REAL8 *)LALInferenceGetVariable( data->dataParams, 
                                                  scaleMinPar );
    
    if( item->vary == LALINFERENCE_PARAM_LINEAR || item->vary == LALINFERENCE_PARAM_CIRCULAR ){
      sprintf(priorPar, "%s_gaussian_mean", item->name);
      /* Check for a gaussian */
      if ( LALInferenceCheckVariable(runState->priorArgs, priorPar) ){
        LALInferenceGetGaussianPrior( runState->priorArgs, item->name, 
                                      (void *)&mu, (void *)&sigma );
      
       value = (*(REAL8 *)item->value) * scale + scaleMin;
       mu += scaleMin;
       sigma *= scale;
       prior -= log(sqrt(2.*LAL_PI)*sigma);
       prior -= (value - mu)*(value - mu) / (2.*sigma*sigma);
      }
      /* Otherwise use a flat prior */
      else{
	LALInferenceGetMinMaxPrior( runState->priorArgs, item->name, 
                                    &min, &max );
      
        if( (*(REAL8 *) item->value) < min || (*(REAL8 *)item->value) > max ){
          return -DBL_MAX;
        }
        else prior -= log( (max - min) * scale );
      }
    }
  }
  return prior; 
}

/*--------------- END OF LIKELIHOOD AND PRIOR FUNCTIONS ----------------------*/

/******************************************************************************/
/*                            MODEL FUNCTIONS                                 */
/******************************************************************************/

/** \brief Defines the pulsar model/template to use
 * 
 * This function is the wrapper for functions defining the pulsar model 
 * template to be used in the analysis. It also uses \c rescale_parameter to
 * scale any parameters back to there true values for use in the model and 
 * places them into a \c BinaryPulsarParams structure.
 * 
 * Note: Any additional models should be added into this function.
 * 
 * \param data [in] The data structure hold data and current parameter info
 * 
 * \sa rescale_parameter
 * \sa pulsar_model
 */
void get_pulsar_model( LALInferenceIFOData *data ){
  BinaryPulsarParams pars;
  
  /* set model parameters (including rescaling) */
  pars.h0 = rescale_parameter( data, "h0" );
  pars.cosiota = rescale_parameter( data, "cosiota" );
  pars.psi = rescale_parameter( data, "psi" );
  pars.phi0 = rescale_parameter( data, "phi0" );
  
  /*pinned superfluid parameters*/
  pars.h1 = rescale_parameter( data, "h1" );
  pars.lambda = rescale_parameter( data, "lambda" );
  pars.theta = rescale_parameter( data, "theta" );
 
  /* set the potentially variable parameters */
  pars.pepoch = rescale_parameter( data, "pepoch" );
  pars.posepoch = rescale_parameter( data, "posepoch" );
 
  pars.ra = rescale_parameter( data, "ra" );
  pars.pmra = rescale_parameter( data, "pmra" );
  pars.dec = rescale_parameter( data, "dec" );
  pars.pmdec = rescale_parameter( data, "pmdec" );
 
  pars.f0 = rescale_parameter( data, "f0" );
  pars.f1 = rescale_parameter( data, "f1" );
  pars.f2 = rescale_parameter( data, "f2" );
  pars.f3 = rescale_parameter( data, "f3" );
  pars.f4 = rescale_parameter( data, "f4" );
  pars.f5 = rescale_parameter( data, "f5" );
  
  /* binary system model - NOT pulsar model */
  pars.model = *(CHAR**)LALInferenceGetVariable( data->modelParams, "model" );

  /* binary parameters */
  if( pars.model != NULL ){
    pars.e = rescale_parameter( data, "e" );
    pars.w0 = rescale_parameter( data, "w0" );
    pars.Pb = rescale_parameter( data, "Pb" );
    pars.x = rescale_parameter( data, "x" );
    pars.T0 = rescale_parameter( data, "T0" );
    
    pars.e2 = rescale_parameter( data, "e2" );
    pars.w02 = rescale_parameter( data, "w02" );
    pars.Pb2 = rescale_parameter( data, "Pb2" );
    pars.x2 = rescale_parameter( data, "x2" );
    pars.T02 = rescale_parameter( data, "T02" );
   
    pars.e3 = rescale_parameter( data, "e3" );
    pars.w03 = rescale_parameter( data, "w03" );
    pars.Pb3 = rescale_parameter( data, "Pb3" );
    pars.x3 = rescale_parameter( data, "x3" );
    pars.T03 = rescale_parameter( data, "T03" );
    
    pars.xpbdot = rescale_parameter( data, "xpbdot" );
    pars.eps1 = rescale_parameter( data, "eps1" );
    pars.eps2 = rescale_parameter( data, "eps2" );
    pars.eps1dot = rescale_parameter( data, "eps1dot" );
    pars.eps2dot = rescale_parameter( data, "eps2dot" );
    pars.Tasc = rescale_parameter( data, "Tasc" );
   
    pars.wdot = rescale_parameter( data, "wdot" );
    pars.gamma = rescale_parameter( data, "gamma" );
    pars.Pbdot = rescale_parameter( data, "Pbdot" );
    pars.xdot = rescale_parameter( data, "xdot" );
    pars.edot = rescale_parameter( data, "edot" );
   
    pars.s = rescale_parameter( data, "s" );
    pars.dr = rescale_parameter( data, "dr" );
    pars.dth = rescale_parameter( data, "dth" );
    pars.a0 = rescale_parameter( data, "a0" );
    pars.b0 = rescale_parameter( data, "b0" ); 

    pars.M = rescale_parameter( data, "M" );
    pars.m2 = rescale_parameter( data, "m2" );
  }

  /* now get pulsar model */
  pulsar_model( pars, data );
    
}


/** \brief Rescale parameter back to its true value
 * 
 * This function will rescale a parameter to its true value using the scale
 * factor and minimum scale value.
 * 
 * \param data [in] data structure containing parameter information
 * \param parname [in] name of the parameter requiring rescaling
 * 
 * \return Rescaled parameter value
 */
REAL8 rescale_parameter( LALInferenceIFOData *data, const CHAR *parname ){
  REAL8 par = 0., scale = 0., offset = 0.;
  CHAR scaleName[VARNAME_MAX] = "";
  CHAR offsetName[VARNAME_MAX] = "";
  
  sprintf(scaleName, "%s_scale", parname);
  sprintf(offsetName, "%s_scale_min", parname);
  
  scale = *(REAL8*)LALInferenceGetVariable( data->dataParams, scaleName );
  offset = *(REAL8*)LALInferenceGetVariable( data->dataParams, offsetName );
  
  par = *(REAL8*)LALInferenceGetVariable( data->modelParams, parname );
  
  par = par*scale + offset;
  
  return par;
}


/** \brief Generate the model of the neutron star signal
 *
 * The function requires that the pulsar model is set using the \c model-type
 * command line argument (this is set in \c main, and if not specified defaults
 * to a \c triaxial model). Currently the model can be \c triaxial for
 * quadrupole emission from a triaxial star at twice the rotation freqeuncy, or
 * \c pinsf for a two component emission model with emission at the rotation
 * frequency <i>and</i> twice the rotation frequency. Depending on the specified
 * model the function calls the appropriate model function.
 * 
 * Firstly the time varying amplitude of the signal will be calculated based on 
 * the antenna pattern and amplitude parameters. Then, if searching over phase 
 * parameters, the phase evolution of the signal will be calculated. The
 * difference between the new phase model, \f$\phi(t)_n\f$, and that used to
 * heterodyne the data, \f$\phi(t)_h\f$, (stored in \c data->timeData->data)
 * will be calculated and the complex signal model, \f$M\f$, modified
 * accordingly:
 * \f[
 * M'(t) = M(t)\exp{i(-(\phi(t)_n - \phi(t)_h))}. 
 * \f]
 * 
 * \param params [in] A \c BinaryPulsarParams structure containing the model
 * parameters
 * \param data [in] The data structure containing the detector data and
 * additional info
 * 
 * \sa get_triaxial_amplitude_model
 * \sa get_pinsf_amplitude_model
 * \sa get_phase_model
 */
void pulsar_model( BinaryPulsarParams params, 
                   LALInferenceIFOData *data ){
  INT4 i = 0, length = 0;
  UINT4 j = 0;
  CHAR *modeltype = NULL;
  
  /* check model type to get amplitude model */
  modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );
  
  if ( !strcmp( modeltype, "triaxial" ) ){
    get_triaxial_amplitude_model( params, data );
  }
  else if ( !strcmp( modeltype, "pinsf" ) ){
    get_pinsf_amplitude_model( params, data );
  }
  /* ADD NEW MODELS HERE */
  else{
    fprintf(stderr, "Error... model '%s' is not defined!\n", modeltype);
    exit(0);
  }
   
  /* get difference in phase for f component and perform extra heterodyne */
  REAL8Vector *freqFactors = NULL;
  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams,
                                                          "freqfactors" );
  
  for( j = 0; j < freqFactors->length; j++ ){
    REAL8Vector *dphi = NULL;
    
    /* move data pointer along one as one iteration of the model is held over
    j data structures, moved this from bottom of loop so actioned for 2nd run through the loop.*/
    if ( j > 0 ){ 
      data = data->next;
    }
    
    length = data->compModelData->data->length;
    /* the timeData vector within the LALIFOData structure contains the
     phase calculated using the initial (heterodyne) values of the phase
     parameters */
    if ( varyphase ){
      if ( (dphi = get_phase_model( params, data, 
        freqFactors->data[j] )) != NULL ){
        for( i=0; i<length; i++ ){
          COMPLEX16 M;
          REAL8 dphit;
          REAL4 sp, cp;
    
          dphit = -fmod(dphi->data[i] - data->timeData->data->data[i], 1.);
    
          sin_cos_2PI_LUT( &sp, &cp, dphit );
    
          M.re = data->compModelData->data->data[i].re;
          M.im = data->compModelData->data->data[i].im;
    
          /* heterodyne */
          data->compModelData->data->data[i].re = M.re*cp - M.im*sp;
          data->compModelData->data->data[i].im = M.im*cp + M.re*sp;
        }
      }
    }
    XLALDestroyREAL8Vector( dphi );
  }
}


/** \brief The phase evolution of a source
 *
 * This function will calculate the phase evolution of a source at a particular
 * sky location as observed at Earth. The phase evolution is described by a 
 * Taylor expansion:
 * \f[
 * \phi(T) = \sum_{k=1}^n \frac{f^{(k-1)}{k!} T^k,
 * \f]
 * where \f$f^x\f$ is the xth time derivative of the gravitational wave
 * frequency, and \f$T\f$ is the pulsar proper time. Frequency time derivatives
 * are currently allowed up to the fifth derivative. The pulsar proper time is 
 * calculated by correcting the time of arrival at Earth, \f$t\f$ to the solar
 * system barycentre and if necessary the binary system barycenter, so \f$T =
 * t + \delta{}t_{\rm SSB} + \delta{}t_{\rm BSB}\f$.
 * 
 * In this function the time delay caused needed to correct to the solar system
 * barycenter is only calculated if required i.e. if it's not been previously
 * calculated and an update is required due to a change in the sky position. The
 * same is true for the binary system time delay, which is only calculated if
 * it has not previously been obtained or needs updating due to a change in the
 * binary system parameters.
 * 
 * The solar system barycentre delay does not have to be explicitly computed
 * for every time stamp passed to it, but instead will just use linear
 * interpolation within a time range set by \c interptime.
 * 
 * \param params [in] A set of pulsar parameters
 * \param data [in] The data structure containing the detector data and
 * additional info
 * \param freqFactor [in] the multiplicative factor on the pulsar frequency for
 * a particular model
 * 
 * \return A vector of rotational phase values
 * 
 * \sa get_ssb_delay
 * \sa get_bsb_delay
 */
REAL8Vector *get_phase_model( BinaryPulsarParams params, 
                              LALInferenceIFOData *data,
                              REAL8 freqFactor ){
  INT4 i = 0, length = 0;

  REAL8 T0 = 0., DT = 0., deltat = 0., deltat2 = 0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */
  
  REAL8Vector *phis = NULL, *dts = NULL, *bdts = NULL;
  
  

  /* if edat is NULL then return a NULL poniter */
  if( data->ephem == NULL )
    return NULL;

  length = data->dataTimes->length;
  
  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );

  /* get time delays */ 
  if( (dts = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams,
      "ssb_delays" )) == NULL || varyskypos == 1 ){
    /* get time delays with an interpolation of interptime (30 mins) */
    dts = get_ssb_delay( params, data->dataTimes, data->ephem, data->detector,
                         interptime );
  }
  
  if( (bdts = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams,
      "bsb_delays" )) == NULL || varybinary == 1 ){
    /* get binary system time delays */
    bdts = get_bsb_delay( params, data->dataTimes, dts );
  }
  
  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8( &data->dataTimes->data[i] );
    
    T0 = params.pepoch;

    DT = realT - T0;

    if ( params.model != NULL )
      deltat = DT + dts->data[i] + bdts->data[i];
    else
      deltat = DT + dts->data[i];
    
    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = freqFactor*deltat*(params.f0 + 
      inv_fact[2]*params.f1*deltat +
      inv_fact[3]*params.f2*deltat2 +
      inv_fact[4]*params.f3*deltat*deltat2 +
      inv_fact[5]*params.f4*deltat2*deltat2 +
      inv_fact[6]*params.f5*deltat2*deltat2*deltat);
  }
  return phis;
}


/** \brief Computes the delay between a GPS time at Earth and the solar system 
 * barycentre
 *
 * This function calculate the time delay between a GPS time at a specific 
 * location (e.g. a gravitational wave detector) on Earth and the solar system
 * barycentre. The delay consists of three components: the geometric time delay
 * (Roemer delay) \f$t_R = \mathbf{r}(t)\hat{n}/c\f$ (where \f$\mathbf{r}(t)\f$
 * is the detector's position vector at time \f$t\f$), the special relativistic
 * Einstein delay \f$t_E\f$, and the general relativistic Shapiro delay
 * \f$t_S\f$.
 * 
 * Rather than computing the time delay at every time stamp passed to the
 * function it is instead (if requested) able to perform linear interpolation
 * to a point within a range given by \c interptime. 
 *  
 * \param pars [in] A set of pulsar parameters
 * \param datatimes [in] A vector of GPS times at Earth
 * \param ephem [in] Information on the solar system ephemeris
 * \param detector [in] Information on the detector position on the Earth
 * \param interptime [in] The time (in seconds) between explicit recalculations
 * of the time delay
 * 
 * \return A vector of time delays in seconds
 *
 * \sa LALBarycenter
 * \sa LALBarycenterEarth
 */
REAL8Vector *get_ssb_delay( BinaryPulsarParams pars, 
                            LIGOTimeGPSVector *datatimes,
                            EphemerisData *ephem,
                            LALDetector *detector,
                            REAL8 interptime ){
  static LALStatus status;

  INT4 i = 0, length = 0;

  REAL8 T0 = 0., DT = 0., DTplus = 0.;

  EarthState earth, earth2;
  EmissionTime emit, emit2;

  BarycenterInput *bary = NULL;
  
  REAL8Vector *dts = NULL;

  /* if edat is NULL then return a NULL poniter */
  if( ephem == NULL )
    return NULL;
  
  /* copy barycenter and ephemeris data */
  bary = (BarycenterInput*)XLALCalloc( 1, sizeof(BarycenterInput) );
  memcpy( &bary->site, detector, sizeof(LALDetector) );
  
  bary->alpha = pars.ra;
  bary->delta = pars.dec;
  
   /* set the position and frequency epochs if not already set */
  if( pars.pepoch == 0. && pars.posepoch != 0.)
    pars.pepoch = pars.posepoch;
  else if( pars.posepoch == 0. && pars.pepoch != 0. )
    pars.posepoch = pars.pepoch;

  length = datatimes->length;
  
  /* allocate memory for times delays */
  dts = XLALCreateREAL8Vector( length );
 
  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( pars.px != 0. )
    bary->dInv = pars.px*1e-3*LAL_C_SI/LAL_PC_SI;
  else if( pars.dist != 0. )
    bary->dInv = LAL_C_SI/(pars.dist*1e3*LAL_PC_SI);
  else
    bary->dInv = 0.;
  
  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8( &datatimes->data[i] );
    
    T0 = pars.pepoch;

    DT = realT - T0;

    /* only do call to the barycentring routines once every interptime (unless
       interptime == 0), otherwise just linearly interpolate between them */
    if( i == 0 || DT > DTplus || interptime == 0 ){
      bary->tgps = datatimes->data[i];

      bary->delta = pars.dec + (realT-pars.posepoch) * pars.pmdec;
      bary->alpha = pars.ra + (realT-pars.posepoch) *
         pars.pmra/cos(bary->delta);
     
      /* call barycentring routines */
      LAL_CALL( LALBarycenterEarth( &status, &earth, &bary->tgps, ephem ),
                &status );
      
      LAL_CALL( LALBarycenter( &status, &emit, bary, &earth ), &status );

      /* add interptime to the time */
      if ( interptime > 0 ){
        DTplus = DT + interptime;
        XLALGPSAdd( &bary->tgps, interptime );

        /* No point in updating the positions as difference will be tiny */
        LAL_CALL( LALBarycenterEarth( &status, &earth2, &bary->tgps, ephem ),
                  &status );
        LAL_CALL( LALBarycenter( &status, &emit2, bary, &earth2), &status );
      }
    }

    /* linearly interpolate to get emitdt */
    if( interptime > 0. ){
      dts->data[i] = emit.deltaT + (DT - (DTplus - interptime)) *
        (emit2.deltaT - emit.deltaT)/interptime;
    }
    else
      dts->data[i] = emit.deltaT;
  }
  
  XLALFree( bary );
  
  return dts;
}


/** \brief Computes the delay between a pulsar in a binary system and the
 * barycentre of the system
 *
 * This function uses \c XLALBinaryPulsarDeltaT to calculate the time delay
 * between for a pulsar in a binary system between the time at the pulsar and
 * the time at the barycentre of the system. This includes Roemer delays and
 * relativistic delays. The orbit may be described by different models and can
 * be purely Keplarian or include various relativistic corrections.
 *
 * \param pars [in] A set of pulsar parameters
 * \param datatimes [in] A vector of GPS times
 * \param dts [in] A vector of solar system barycentre time delays
 * 
 * \return A vector of time delays in seconds
 * 
 * \sa XLALBinaryPulsarDeltaT
 */
REAL8Vector *get_bsb_delay( BinaryPulsarParams pars,
                            LIGOTimeGPSVector *datatimes,
                            REAL8Vector *dts ){
  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;
  REAL8Vector *bdts = NULL;
  
  INT4 i = 0, length = datatimes->length;
  
  bdts = XLALCreateREAL8Vector( length );
  
  for ( i = 0; i < length; i++ ){
    binput.tb = XLALGPSGetREAL8( &datatimes->data[i] ) + dts->data[i];
  
    XLALBinaryPulsarDeltaT( &boutput, &binput, &pars );
    
    bdts->data[i] = boutput.deltaT;
  }
  
  return bdts;
}


/** \brief The amplitude model of a complex heterodyned traxial neutron star
 * 
 * This function calculates the complex heterodyned time series model for a 
 * triaxial neutron star (see [\ref DupuisWoan2005]). It is defined as:
 * \f{eqnarray*}{
 * y(t) & = & \frac{h_0}{2} \left( \frac{1}{2}F_+(t,\psi)
 * (1+\cos^2\iota)\cos{\phi_0} + F_{\times}(t,\psi)\cos{\iota}\sin{\phi_0}
 * \right) + \\
 *  & & i\frac{h_0}{2}\left( \frac{1}{2}F_+(t,\psi)
 * (1+\cos^2\iota)\sin{\phi_0} - F_{\times}(t,\psi)\cos{\iota}\cos{\phi_0}
 * \right),
 * \f}
 * where \f$F_+\f$ and \f$F_{\times}\f$ are the antenna response functions for
 * the plus and cross polarisations.
 * 
 * The antenna pattern functions are contained in a 2D lookup table, so within
 * this function the correct value for the given time and \f$\psi\f$ are
 * interpolated from this lookup table using bilinear interpolation (e.g.):
 * \f{eqnarray*}{
 * F_+(\psi, t) = F_+(\psi_i, t_j)(1-\psi)(1-t) + F_+(\psi_{i+1}, t_j)\psi(1-t)
 * + F_+(\psi_i, t_{j+1})(1-\psi)t + F_+(\psi_{i+1}, t_{j+1})\psi{}t,
 * \f}
 * where \f$\psi\f$ and \f$t\f$ have been scaled to be within a unit square,
 * and \f$\psi_i\f$ and \f$t_j\f$ are the closest points within the lookup
 * table to the required values.
 * 
 * \param pars [in] A set of pulsar parameters
 * \param data [in] The data parameters giving information on the data and
 * detector
 * 
 */
void get_triaxial_amplitude_model( BinaryPulsarParams pars, 
                                   LALInferenceIFOData *data ){
  INT4 i = 0, length;
  
  REAL8 psteps, tsteps, psv, tsv;
  INT4 psibinMin, psibinMax, timebinMin, timebinMax;
  REAL8 tstart;
  REAL8 plus, cross;
  REAL8 plus00, plus01, plus10, plus11, cross00, cross01, cross10, cross11;
  REAL8 psiScaled, timeScaled;
  REAL8 psiMin, psiMax, timeMin, timeMax;
  REAL8 T;
  REAL8 Xplus, Xcross;
  REAL8 Xpcosphi, Xccosphi, Xpsinphi, Xcsinphi;
  REAL4 sinphi, cosphi;
  
  gsl_matrix *LU_Fplus, *LU_Fcross;
  
  length = data->dataTimes->length;
  
  /* set lookup table parameters */
  psteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "timeSteps" );
  
  LU_Fplus = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, 
                                                     "LU_Fplus");
  LU_Fcross = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, 
                                                      "LU_Fcross");
  
  sin_cos_LUT( &sinphi, &cosphi, pars.phi0 );
  
  /************************* CREATE MODEL *************************************/
  /* This model is a complex heterodyned time series for a triaxial neutron
     star emitting at twice its rotation frequency (as defined in Dupuis and
     Woan, PRD, 2005):
       real = (h0/2) * ((1/2)*F+*(1+cos(iota)^2)*cos(phi0) 
         + Fx*cos(iota)*sin(phi0))
       imag = (h0/2) * ((1/2)*F+*(1+cos(iota)^2)*sin(phi0)
         - Fx*cos(iota)*cos(phi0))
   ****************************************************************************/
  
  Xplus = 0.25*(1.+pars.cosiota*pars.cosiota)*pars.h0;
  Xcross = 0.5*pars.cosiota*pars.h0;
  Xpsinphi = Xplus*sinphi;
  Xcsinphi = Xcross*sinphi;
  Xpcosphi = Xplus*cosphi;
  Xccosphi = Xcross*cosphi;
 
  tstart = XLALGPSGetREAL8( &data->dataTimes->data[0] ); /*time of first B_k*/
  
  /* set the psi bin for the lookup table */
  psv = LAL_PI_2 / ( psteps - 1. );
  psibinMin = (INT4)floor( ( pars.psi + LAL_PI/4. )/psv );
  psiMin = -(LAL_PI/4.) + psibinMin*psv;
  psibinMax = psibinMin + 1;
  psiMax = psiMin + psv;
  
  /* rescale psi for bilinear interpolation on a unit square */
  psiScaled = (pars.psi - psiMin)/(psiMax - psiMin);
  
  tsv = LAL_DAYSID_SI / tsteps;
  
  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = fmod( XLALGPSGetREAL8(&data->dataTimes->data[i]) - tstart,
              LAL_DAYSID_SI );
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;
    
    /* get values of matrix for bilinear interpolation */
    plus00 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMin );
    plus01 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMax );
    plus10 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMin );
    plus11 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMax );
    
    cross00 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMin );
    cross01 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMax );
    cross10 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMin );
    cross11 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMax );
    
    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);
    
    plus = plus00*(1. - psiScaled)*(1. - timeScaled) + 
      plus10*psiScaled*(1. - timeScaled) + plus01*(1. - psiScaled)*timeScaled +
      plus11*psiScaled*timeScaled;
    cross = cross00*(1. - psiScaled)*(1. - timeScaled) + 
      cross10*psiScaled*(1. - timeScaled) + cross01*(1. - psiScaled)*timeScaled
      + cross11*psiScaled*timeScaled;
    
    /* create the complex signal amplitude model */
    data->compModelData->data->data[i].re = plus*Xpcosphi + cross*Xcsinphi;
    data->compModelData->data->data[i].im = plus*Xpsinphi - cross*Xccosphi;
  }
}


void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data ){
  INT4 i = 0, length;
  
  REAL8 psteps, tsteps, psv, tsv;
  INT4 psibinMin, psibinMax, timebinMin, timebinMax;
  REAL8 plus00, plus01, plus10, plus11, cross00, cross01, cross10, cross11;
  REAL8 psiScaled, timeScaled;
  REAL8 psiMin, psiMax, timeMin, timeMax;
  REAL8 tstart;
  REAL8 plus, cross;
  REAL8 T;
  REAL8 Xplusf, Xcrossf, Xplus2f, Xcross2f;
  REAL8 A1, A2, B1, B2;
  REAL4 sinphi, cosphi, sin2phi, cos2phi;
  
  gsl_matrix *LU_Fplus, *LU_Fcross;
  
  
  /* set lookup table parameters */
  psteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "timeSteps" );
  
  LU_Fplus = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fplus");
  LU_Fcross = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fcross");
  
  sin_cos_LUT( &sinphi, &cosphi, 0.5*pars.phi0 );
  sin_cos_LUT( &sin2phi, &cos2phi, pars.phi0 );
  
  /************************* CREATE MODEL *************************************/
  /* This model is a complex heterodyned time series for a pinned superfluid neutron
     star emitting at its roation frequency and twice its rotation frequency 
     (as defined in Jones 2009):

   ****************************************************************************/
  
  Xplusf = 0.125*sin(acos(pars.cosiota))*pars.cosiota*pars.h0;
  Xcrossf = 0.125*sin(acos(pars.cosiota))*pars.h0;
  Xplus2f = 0.25*(1.+pars.cosiota*pars.cosiota)*pars.h0;
  Xcross2f = 0.5*pars.cosiota*pars.h0;
  A1=( (cos(pars.lambda)*cos(pars.lambda)) - pars.h1 ) * (sin( (2*pars.theta) ));
  A2=sin(2*pars.lambda)*sin(pars.theta);
  B1=( (cos(pars.lambda)*cos(pars.lambda))*(cos(pars.theta)*cos(pars.theta)) ) - (sin(pars.lambda)*sin(pars.lambda)) 
    + ( pars.h1*(sin(pars.theta)*sin(pars.theta)) );
  B2=sin(2*pars.lambda)*cos(pars.theta);

  tstart = XLALGPSGetREAL8( &data->dataTimes->data[0] );
  
  /* set the psi bin for the lookup table */
  psv = LAL_PI_2 / ( psteps - 1. );
  psibinMin = (INT4)floor( ( pars.psi + LAL_PI/4. )/psv );
  psiMin = -(LAL_PI/4.) + psibinMin*psv;
  psibinMax = psibinMin + 1;
  psiMax = psiMin + psv;
  
  /* rescale psi for bilinear interpolation on a unit square */
  psiScaled = (pars.psi - psiMin)/(psiMax - psiMin);
  
  tsv = LAL_DAYSID_SI / tsteps;
  
  /* set model for 1f component */
  length = data->dataTimes->length;
  tstart = XLALGPSGetREAL8( &data->dataTimes->data[0] ); /*time of first B_k*/
  
  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/    
    T = fmod( XLALGPSGetREAL8(&data->dataTimes->data[i]) - tstart,
              LAL_DAYSID_SI );
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;
    
    /* get values of matrix for bilinear interpolation */
    plus00 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMin );
    plus01 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMax );
    plus10 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMin );
    plus11 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMax );
    
    cross00 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMin );
    cross01 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMax );
    cross10 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMin );
    cross11 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMax );
    
    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);
    
    plus = plus00*(1. - psiScaled)*(1. - timeScaled) + 
      plus10*psiScaled*(1. - timeScaled) + plus01*(1. - psiScaled)*timeScaled +
      plus11*psiScaled*timeScaled;
    cross = cross00*(1. - psiScaled)*(1. - timeScaled) + 
      cross10*psiScaled*(1. - timeScaled) + cross01*(1. - psiScaled)*timeScaled
      + cross11*psiScaled*timeScaled;
    
    /* create the complex signal amplitude model */
    /*at f*/
    data->compModelData->data->data[i].re = plus*Xplusf*((A1*cosphi)-(A2*sinphi)) + 
    ( cross*Xcrossf*((A2*cosphi)-(A1*sinphi)) );
    
    data->compModelData->data->data[i].im = plus*Xplusf*((A2*cosphi)+(A1*sinphi)) + 
    ( cross*Xcrossf*((A2*sinphi)-(A1*cosphi)) );

  }
  
  /* set model for 2f component */
  length = data->next->dataTimes->length;
  tstart = XLALGPSGetREAL8( &data->next->dataTimes->data[0] );
  
  for( i=0; i<length; i++ ){
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/    
    T = fmod( XLALGPSGetREAL8(&data->next->dataTimes->data[i]) - tstart,
              LAL_DAYSID_SI );
    timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
    timeMin = timebinMin*tsv;
    timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
    timeMax = timeMin + tsv;
    
    /* get values of matrix for bilinear interpolation */
    plus00 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMin );
    plus01 = gsl_matrix_get( LU_Fplus, psibinMin, timebinMax );
    plus10 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMin );
    plus11 = gsl_matrix_get( LU_Fplus, psibinMax, timebinMax );
    
    cross00 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMin );
    cross01 = gsl_matrix_get( LU_Fcross, psibinMin, timebinMax );
    cross10 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMin );
    cross11 = gsl_matrix_get( LU_Fcross, psibinMax, timebinMax );
    
    /* rescale time for bilinear interpolation on a unit square */
    timeScaled = (T - timeMin)/(timeMax - timeMin);
    
    plus = plus00*(1. - psiScaled)*(1. - timeScaled) + 
      plus10*psiScaled*(1. - timeScaled) + plus01*(1. - psiScaled)*timeScaled +
      plus11*psiScaled*timeScaled;
    cross = cross00*(1. - psiScaled)*(1. - timeScaled) + 
      cross10*psiScaled*(1. - timeScaled) + cross01*(1. - psiScaled)*timeScaled
      + cross11*psiScaled*timeScaled;
    
    /* create the complex signal amplitude model at 2f*/
    data->next->compModelData->data->data[i].re =
      plus*Xplus2f*((B1*cos2phi)-(B2*sin2phi)) +
      cross*Xcross2f*((B2*cos2phi)+(B1*sin2phi));
    
    data->next->compModelData->data->data[i].im =
      plus*Xplus2f*((B2*cos2phi)+(B1*sin2phi)) -
      cross*Xcross2f*((B1*cos2phi)-(B2*sin2phi));
  }
  
}


/** \brief Calculate the natural logarithm of the evidence that the data
 * consists of only Gaussian noise
 * 
 * The function will calculate the natural logarithm of the evidence that the
 * data (from one or more detectors) consists of stationary segments/chunks 
 * describe by a Gaussian with zero mean and unknown variance.
 * 
 * The evidence is obtained from the joint likelihood given in \c
 * pulsar_log_likelihood with the model term \f$y\f$ set to zero.
 * 
 * \param data [in] Structure containing detector data
 * 
 * \return The natural logarithm of the noise only evidence
 */
REAL8 noise_only_model( LALInferenceIFOData *data ){
  LALInferenceIFOData *datatemp = data;
  
  REAL8 logL = 0.0;
  UINT4 i = 0;
  
  while ( datatemp ){
    UINT4Vector *chunkLengths = NULL;
    REAL8Vector *sumDat = NULL;
  
    REAL8 chunkLength = 0.;
  
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams, 
                                                             "chunkLength" );
    sumDat = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams,
                                                        "sumData" );
  
    for (i=0; i<chunkLengths->length; i++){
      chunkLength = (REAL8)chunkLengths->data[i];
   
      logL -= chunkLength * log(sumDat->data[i]);
    }
  
    datatemp = datatemp->next;
  }
  
  return logL;
}

/*------------------------ END OF MODEL FUNCTIONS ----------------------------*/


/******************************************************************************/
/*                       SOFTWARE INJECTION FUNCTIONS                         */
/******************************************************************************/

/** \brief Inject a simulated signal into the data 
 *
 * This function will create an simulated signal (of the required model) to
 * inject into the data from multiple detectors. The parameters of the signal
 * to be injected must be specified in a TEMPO-stype .par file given with the
 * \c inject-file command line argument. The parameters do not have to be the 
 * same as those in the .par file controlling the analysis (although should
 * ideally contain a signal within the bandwidth of the data).
 * 
 * If a signal of a specific signal-to-noise ratio is required then the \c
 * scale-snr command line argument can be used to give the multi-detector SNR to
 * which the signal needs to be scaled.
 * 
 * The injected signal can be output if \c inject-output is set. Two
 * files will be output: one containing the signal only, and one containing the
 * signal plus noise. These will both be in the format of a standard data input 
 * file. The files will have names given by the \c inject-output value, with a
 * prefix of the detector name, and a suffix of of \c _signal_only,
 * respectively. 
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
  INT4 ndets = 0, j = 1, k = 0, numSNRs = 1;
  
  REAL8Vector *freqFactors = NULL;
  REAL8 *snrmulti = NULL, *snrscale = NULL;
  
  CHAR *modeltype = NULL;/*need to check model type in this function*/
  
  modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );
  
  ppt = LALInferenceGetProcParamVal( commandLine, "--inject-file" );
  if( ppt ){

    injectfile = XLALStringDuplicate( ppt->value );
    
    /* check that the file exists */
    if ( access(injectfile, F_OK) != 0 ){
      fprintf(stderr, "Error... Injection specified, but the injection \
parameter file %s is wrong.\n", injectfile);
      exit(3);
    }
    
    /* read in injection parameter file */
    XLALReadTEMPOParFile( &injpars, injectfile );
  }
  else{
    return;
  }
 
  freqFactors = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, 
                                                          "freqfactors" );
 
  snrscale = XLALCalloc(sizeof(REAL8), numSNRs);
  
  ppt = LALInferenceGetProcParamVal( commandLine, "--scale-snr" );
  if ( ppt ){
    /* if there are more than one data streams the SNRs for the individual
       (multi-detector) streams can be set, rather than having a combined SNR.
       The SNR vales are set as comma separated values to --snr-scale. If only
       one value is given, but the data has multiple streams then the combined
       multi-stream SNR will still be used. */
    CHAR *snrscales = NULL, *tmpsnrs = NULL, *tmpsnr = NULL, snrval[256];
    
    snrscales = XLALStringDuplicate( ppt->value );
      
    tmpsnrs = XLALStringDuplicate( snrscales );
      
    /* count the number of SNRs (comma seperated values) */
    numSNRs = count_csv( snrscales );
    
    if( numSNRs != 1 || numSNRs != (INT4)freqFactors->length ){
      fprintf(stderr, "Error... number of SNR values must either be 1, or equal\
 to the number of data streams required for your model!\n");
      exit(0);
    }
    
    snrscale = XLALRealloc(snrscale, sizeof(REAL8)*numSNRs);
    
    for( k = 0; k < numSNRs; k++ ){
      tmpsnr = strsep( &tmpsnrs, "," );
      XLALStringCopy( snrval, tmpsnr, strlen(tmpsnr)+1 );
      snrscale[k] = atof(snrval);
    }
  }
 
  snrmulti = XLALCalloc(sizeof(REAL8), numSNRs);
 
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
        
  /* create signal to inject */
  while( data ){
    /* for injection always attempt to include the signal phase model even if
       the search is not going to be over phase */
    UINT4 varyphasetmp = varyphase;
    varyphase = 1;
    
    pulsar_model( injpars, data );
    
    /* reset varyphase to its original value */
    varyphase = varyphasetmp;
    
    data = data->next;
    
    /* If modeltype uses more than one data stream need to advance data on to
       next, so this loop only runs once if there is only 1 det*/
    for ( k = 1; k < (INT4)freqFactors->length; k++ ) data = data->next;
  }
  
  /* reset data to head */
  data = runState->data;
  
  /* calculate SNRs */
  while ( data ){
    for ( k = 0; k < numSNRs; k++ ){
      REAL8 snrval = calculate_time_domain_snr( data );
   
      snrmulti[k] += SQUARE(snrval);
    
      if ( snrscale[k] == 0 ) fprintf(fpsnr, "%le\t", snrval);
                             
      data = data->next;
    }
    
    ndets++;
  }
  
  /* get overall multi-detector SNR */
  for ( k = 0; k < numSNRs; k++ ) snrmulti[k] = sqrt( snrmulti[k] );
  
  /* only need to print out multi-detector snr if the were multiple detectors */
  if( numSNRs == 1 && snrscale[0] == 0 ){
    if ( ndets > 1 ) fprintf(fpsnr, "%le\n", snrmulti[0]);
    else fprintf(fpsnr, "\n");
  }
  else{
    /* rescale the signal and calculate the SNRs */
    data = runState->data;
   
    for ( k = 0; k < numSNRs; k++ ) snrscale[k] /= snrmulti[k];
    
    /* rescale the h0 (and other amplitude factors for other models) */
    if ( !strcmp(modeltype, "triaxial") || !strcmp(modeltype, "pinsf") ){
      for ( k = 0; k < numSNRs; k++ ){
        if ( freqFactors->data[k] == 2. ) injpars.h0 *= snrscale[k];
        if ( freqFactors->data[k] == 1. ) injpars.h1 *= snrscale[k];
      }
    }
    
    /* recreate the signal with scale amplitude */
    while( data ){
      UINT4 varyphasetmp = varyphase;
      varyphase = 1;
      
      pulsar_model( injpars, data );
      
      /* reset varyphase to its original value */
      varyphase = varyphasetmp;
      
      for ( k = 1; k < (INT4)freqFactors->length; k++ ) data = data->next;
    }
    
    data = runState->data;
    
    /* get new snrs */
    for ( k = 0; k < numSNRs; k++ ) snrmulti[k] = 0;
    
    while( data ){
      for ( k = 0; k < numSNRs; k++ ){
        REAL8 snrval = 0.;

        /* recalculate the SNR */
        snrval = calculate_time_domain_snr( data );
      
        snrmulti[k] += SQUARE(snrval);
      
        fprintf(fpsnr, "%le\t", snrval);
      
        data = data->next;
      }
    }
    
    for ( k = 0; k < numSNRs; k++ ) snrmulti[k] = sqrt( snrmulti[k] );
    
    if( ndets > 1 ){ 
      for ( k = 0; k < numSNRs; k++ ) fprintf(fpsnr, "%le\t", snrmulti[k]);
      fprintf(fpsnr, "\n");
      /*fprintf(stderr,"snr multi: %f\n",snrmulti); */
    }
    else fprintf(fpsnr, "\n");

  }
  
  fclose( fpsnr );
   
  /* reset data to head */
  data = runState->data;
  
  /* add signal to data */
  while( data ){
    FILE *fp = NULL, *fpso = NULL;
    ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal( commandLine,
                                                            "--inject-output" );
    INT4 i = 0, length = data->dataTimes->length;
                        
    /* check whether to output the data */
    if ( ppt2 ){
      /* add the site prefix to the start of the output name */
      CHAR *outfile = NULL;
      CHAR *signalonly = NULL; /* file containing only signal and no noise */
      
      outfile = XLALStringDuplicate( data->detector->frDetector.prefix );
      
      outfile = XLALStringAppend( outfile, ppt2->value );
      
      if ( ( !strcmp( modeltype, "pinsf" ) ) && (fmod(j,2)==0) ){
        outfile = XLALStringAppend( outfile, "_2f");
      }
      else if ( !strcmp( modeltype, "pinsf" )){
        outfile = XLALStringAppend( outfile, "_1f");
      }

      if ( (fp = fopen(outfile, "w")) == NULL ){
        fprintf(stderr, "Non-fatal error... unable to open file %s to output \
injection\n", outfile);
      }
      
      signalonly = XLALStringDuplicate( outfile );
      signalonly = XLALStringAppend( signalonly, "_signal_only" );
      
      if ( (fpso = fopen(signalonly, "w")) == NULL ){
        fprintf(stderr, "Non-fatal error... unable to open file %s to output \
injection\n", signalonly);
      }
    }
    
    /* add the signal to the data */
    for ( i = 0; i < length; i++ ){
      data->compTimeData->data->data[i].re +=
        data->compModelData->data->data[i].re;
      data->compTimeData->data->data[i].im +=
        data->compModelData->data->data[i].im;
        
      /* write out injection to file */
      if( fp != NULL && fpso != NULL ){
        /* print out data - time stamp, real and imaginary parts of data
           (injected signal + noise) */
        fprintf(fp, "%.5lf\t%le\t%le\n", 
                XLALGPSGetREAL8( &data->dataTimes->data[i] ),
                data->compTimeData->data->data[i].re, 
                data->compTimeData->data->data[i].im );
        
        /* print signal only data - time stamp, real and imaginary parts of
           signal */
        fprintf(fpso, "%.5lf\t%le\t%le\n", 
                XLALGPSGetREAL8( &data->dataTimes->data[i] ),
                data->compModelData->data->data[i].re, 
                data->compModelData->data->data[i].im );
      }
    }

    if ( fp != NULL ) fclose( fp );
    if ( fpso != NULL ) fclose( fpso );
    
    data = data->next;
    j++;
  }
  
  XLALFree(snrmulti);
  XLALFree(snrscale);
}

/*-------------------- END OF SOFTWARE INJECTION FUNCTIONS -------------------*/


/******************************************************************************/
/*                            HELPER FUNCTIONS                                */
/******************************************************************************/

/** \brief Split the data into segments
 * 
 * This function is deprecated to \c chop_n_merge, but gives the functionality 
 * of the old code.
 * 
 * It cuts the data into as many contiguous segments of data as possible of 
 * length \c chunkMax. Where contiguous is defined as containing consecutive
 * point within 180 seconds of each other. The length of segments that do not
 * fit into a \c chunkMax length are also included. 
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
    
    /* if consecutive points are within 180 seconds of each other count as in
       the same chunk */
    if( t2 - t1 > 180. || count == chunkMax ){
      chunkLengths->data[j] = count;
      count = 0; /* reset counter */

      j++;
    }
  }

  chunkLengths = XLALResizeUINT4Vector( chunkLengths, j );
  
  return chunkLengths;
}


/* function to use change point analysis to chop up and remerge the data to
   find stationary chunks (i.e. lengths of data which look like they have the
   same statistics e.g. the same standard deviation) */
/** \brief Chops and remerges data into stationary segments
 * 
 * This function finds segments of data that appear to be stationary (have the
 * same standard deviation).
 * 
 * The function first attempts to chop up the data into as many stationary 
 * segments as possible. The splitting may not be optimal, so it then tries 
 * remerging consecutive segments to see if the merged segments show more
 * evidence of stationarity. <b>[NOTE: Remerging is currently turned off and
 * will make very little difference to the algorithm]</b>. It then, if
 * necessary, chops the segments again to make sure there are none greater
 * than the required \c chunkMax. The default \c chunkMax is 0, so this
 * rechopping will not normally happen.
 * 
 * This is all performed on data that has had a running median subtracted, to 
 * try and removed any underlying trends in the data (e.g. those caused by a 
 * strong signal), which might affect the calculations (which assume the data is
 * Gaussian with zero mean).
 * 
 * if the \c verbose flag is set then a list of the segments will be output to
 * a file called \c data_segment_list.txt, with a prefix of the detector name.
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
UINT4Vector *chop_n_merge( LALInferenceIFOData *data, INT4 chunkMin, 
                           INT4 chunkMax ){
  UINT4 j = 0;
  
  UINT4Vector *chunkLengths = NULL;
  UINT4Vector *chunkIndex = NULL;
  
  COMPLEX16Vector *meddata = NULL;
  
  /* subtract a running median value from the data to remove any underlying
     trends (e.g. caused by a string signal) that might affect the chunk
     calculations (which can assume the data is Gaussian with zero mean). */
  meddata = subtract_running_median( data->compTimeData->data );
  
  chunkIndex = chop_data( meddata, chunkMin );
  
  /* DON'T BOTHER WITH THE MERGING AS IT WILL MAKE VERY LITTLE DIFFERENCE */
  /* merge_data( meddata, chunkIndex ); */
  
  /* if a maximum chunk length is defined then rechop up the data, to segment
     any chunks longer than this value */
  if ( chunkMax > chunkMin )
    rechop_data( chunkIndex, chunkMax, chunkMin );
  
  chunkLengths = XLALCreateUINT4Vector( chunkIndex->length );
  
  /* go through segments and turn into vector of chunk lengths */
  for ( j = 0; j < chunkIndex->length; j++ ){
    if ( j == 0 )
      chunkLengths->data[j] = chunkIndex->data[j];
    else
      chunkLengths->data[j] = chunkIndex->data[j] - chunkIndex->data[j-1];
  }
  
  /* if verbose print out the segment end indices to a file */
  if ( verbose ){
    FILE *fpsegs = NULL;
    
    CHAR *outfile = NULL;
    
    /* set detector name as prefix */
    outfile = XLALStringDuplicate( data->detector->frDetector.prefix );
      
    outfile = XLALStringAppend( outfile, "data_segment_list.txt" );
    
    if ( (fpsegs = fopen(outfile, "w")) == NULL ){
      fprintf(stderr, "Non-fatal error open file to output segment list.\n");
      return chunkLengths;
    }
    
    for ( j = 0; j < chunkIndex->length; j++ )
      fprintf(fpsegs, "%u\n", chunkIndex->data[j]);
    
    /* add space at the end so that you can seperate lists from different
       detector data streams */
    fprintf(fpsegs, "\n");
    
    fclose( fpsegs );
  }
  
  return chunkLengths;
}


/** \brief Subtract the running median from complex data
 * 
 * This function uses \c gsl_stats_median_from_sorted_data to subtract a running
 * median, calculated from the 30 consecutive point around a set point, from the
 * data. At the start of the data running median is calculated from 30-15+(i-1)
 * points, and at the end it is calculated from 15+(N-i) points, where i is the
 * point index and N is the total number of data points.
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
    double *dre = NULL;
    double *dim = NULL;
    
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
    
    dre = XLALCalloc( n, sizeof(double) );
    dim = XLALCalloc( n, sizeof(double) );
    
    for ( j = 0; j < n; j++ ){
      dre[j] = data->data[j+sidx].re;
      dim[j] = data->data[j+sidx].im;
    }
    
    /* sort data */
    gsl_sort( dre, 1, n );
    gsl_sort( dim, 1, n );
    
    /* get median and subtract from data*/
    submed->data[i-1].re = data->data[i-1].re
      - gsl_stats_median_from_sorted_data( dre, 1, n );
    submed->data[i-1].im = data->data[i-1].im
      - gsl_stats_median_from_sorted_data( dim, 1, n );
      
    XLALFree( dre );
    XLALFree( dim );
  }
  
  return submed;
}


/** \brief Chops the data into stationary segments based on Bayesian change
 * point analysis
 * 
 * This function splits data into two (and recursively runs on those two
 * segments) if it is found that the odds ratio for them being from two 
 * independent Gaussian distributions is greater than a certain threshold.
 * 
 * The threshold is for the natural logarithm of the odds ratio is empirically
 * set to be:
 * \f[
 * T = 0.57\ln{N} + 2.71,
 * \f]
 * where \f$N\f$ is the length of the data set. This comes from a fit to the 
 * threshold value required to give a 1% chance of splitting actual Gaussian
 * data (drawn from one distribution) for data of various lengths. The first
 * two terms come from a fit to odds ratios for a Monte Carlo
 * of Gaussian noise (with real and imaginary components) of various lengths,
 * and the final term comes from an offset to give the 1% false alarm rate.
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
  
  /* set threshold based on empirical tests that only give a 1% chance of
     splitting Gaussian data of various lengths. The relation is approximately:
     T = 0.57*ln(length) + 2.71
     where this comes from a fit to odds ratios for a Monte Carlo
     of Gaussian noise (with real and imaginary components) of various lengths,
     with an offset to give the 1% false alarm rate. */
  threshold = 0.57*log(length) + 2.71;
  
  if ( logodds > threshold ){
    UINT4Vector *cp1 = NULL;
    UINT4Vector *cp2 = NULL;
    
    COMPLEX16Vector *data1 = XLALCreateCOMPLEX16Vector( changepoint );
    COMPLEX16Vector *data2 = XLALCreateCOMPLEX16Vector( length - changepoint );
  
    UINT4 i = 0, l = 0;
    
    /* fill in data */
    for (i = 0; i < changepoint; i++)
      data1->data[i] = data->data[i];

    for (i = 0; i < length - changepoint; i++)
      data2->data[i] = data->data[i+changepoint];
    
    cp1 = chop_data( data1, chunkMin );
    cp2 = chop_data( data2, chunkMin );
    
    l = cp1->length + cp2->length;
    
    chunkIndex = XLALResizeUINT4Vector( chunkIndex, l );
    
    /* combine new chunks */
    for (i = 0; i < cp1->length; i++)
      chunkIndex->data[i] = cp1->data[i];
    
    for (i = 0; i < cp2->length; i++)
      chunkIndex->data[i+cp1->length] = cp2->data[i] + changepoint;
  
    /* free memory */
    XLALDestroyCOMPLEX16Vector( data1 );
    XLALDestroyCOMPLEX16Vector( data2 );
  }
  else{
    chunkIndex->data[0] = length;
  }
  
  return chunkIndex;
}


/** \brief Find a change point in complex data
 * 
 * This function is based in the Bayesian Blocks algorithm of [\ref Scargle1998]
 * that finds "change points" in data - points at which the statistics of the 
 * data change. It is based on calculating evidence, or odds, ratios. The
 * function first computes the marginal likelihood (or evidence) that the whole
 * of the data is described by a single Gaussian (with mean of zero). This comes
 * from taking a Gaussian likelihood function and analytically marginalising
 * over the standard deviation (using a prior on the standard deviation of
 * \f$1/\sigma\f$), giving (see [\ref DupuisWoan2005]) a Students-t
 * distribution (see <a
href="https://wiki.ligo.org/foswiki/pub/CW/
PulsarParameterEstimationNestedSampling/studentst.pdf">here</a>).
 * Following this the data is split into two segments (with lengths greater
 * than, or equal to the minimum chunk length) for all possible combinations,
 * and the joint evidence for each of the two segments consisting of independent
 * Gaussian (basically multiplying the above equation calculated for each
 * segment separately) is calculated and the split point recorded. However, the
 * value required for comparing to that for the whole data set, to give the odds
 * ratio, is the evidence that having any split is better than having no split,
 * so the individual split data evidences need to be added incoherently to give
 * the total evidence for a split. The index at which the evidence for a single
 * split is maximum (i.e. the most favoured split point) that which is returned.
 * 
 * \param data [in] a complex data vector
 * \param logodds [in] a pointer to return the natural logarithm of the odds
 * ratio/Bayes factor
 * \param minlength [in] the minimum chunk length
 * 
 * \return The position of the change point
 */
UINT4 find_change_point( COMPLEX16Vector *data, REAL8 *logodds, 
                         INT4 minlength ){
  UINT4 changepoint = 0, i = 0;
  UINT4 length = data->length, lsum = 0;
  
  REAL8 datasum = 0.;
  
  REAL8 logsingle = 0., logtot = -LAL_REAL8_MAX;
  REAL8 logdouble = 0., logdouble_min = -LAL_REAL8_MAX;
  REAL8 logratio = 0.;
 
  REAL8Vector *sumforward = NULL, *sumback = NULL;
  
  /* check that data is at least twice the minimum length, if not return an
     odds ratio of zero (log odds = -inf [or close to that!]) */
  if ( length < (UINT4)(2*minlength) ){
    logratio = -LAL_REAL8_MAX;
    memcpy(logodds, &logratio, sizeof(REAL8));
    return 0;
  }
  
  /* calculate the sum of the data squared */
  for (i = 0; i < length; i++){
    datasum += SQUARE( data->data[i].re );
    datasum += SQUARE( data->data[i].im );
  }
  
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
      sumforward->data[0] += SQUARE( data->data[i].re );
      sumforward->data[0] += SQUARE( data->data[i].im );
      
      sumback->data[0] += SQUARE( data->data[length-(i+1)].re );
      sumback->data[0] += SQUARE( data->data[length-(i+1)].im );
    }
    else{
      sumforward->data[i+1-minlength] = sumforward->data[i-minlength] + 
        SQUARE( data->data[i].re );
      sumforward->data[i+1-minlength] += SQUARE( data->data[i].im );
      
      sumback->data[i+1-minlength] = sumback->data[i-minlength] +
        SQUARE( data->data[length-(i+1)].re );
      sumback->data[i+1-minlength] += SQUARE( data->data[length-(i+1)].im );
    }
  }
  
  /* go through each possible change point and calculate the evidence for the
     data consisting of two independent Gaussian's either side of the change
     point. Also calculate the total evidence for any change point */
  /* Don't allow single points, so start at the second data point. */
  for (i = 0; i < lsum; i++){ 
    UINT4 ln1 = i+minlength, ln2 = (length-i-minlength);
   
    REAL8 log_1 = 0., log_2 = 0.;
    
    /* get log evidences for the individual segments */
    log_1 = -2 + logfactorial[ln1-1] -
      (REAL8)ln1 * log( sumforward->data[i] );
    log_2 = -2 + logfactorial[ln2-1] -
      (REAL8)ln2 * log( sumback->data[lsum-i-1] );
      
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
 * This function chops any chunks that are greater than \c chunkMax into chunks
 * smaller than, or equal to \c chunkMax, and greater than \c chunkMin. On some
 * occasions this might result in a segment smaller than \c chunkMin, but these
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
  
  /* chop any chunks that are greater than chunkMax into chunks smaller than,
     or equal to chunkMax, and greater than chunkMin */
  for ( i = 0; i < length; i++ ){
    if ( i == 0 ) startindex = 0;
    else startindex = chunkIndex->data[i-1]+1;
    
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
          /* split the last two cells into one that is chunkMin long and one
             that is (chunkMax+remainder)-chunkMin long - this may leave a cell
             shorter than chunkMin, but we'll have to live with that! */
          INT4 n1 = (chunkMax + remain) - chunkMin;
       
          /* reset second to last value two values */
          newindex->data[count-1] = newindex->data[count] - chunkMin;
        
          if ( n1 < chunkMin && verbose ){
            fprintf(stderr, "Non-fatal error... segment no. %d is %d long, \
which is less than chunkMin = %d.\n", count, n1, chunkMin);
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

  for ( i = 0; i < count; i++ ) chunkIndex->data[i] = newindex->data[i];
  
  XLALDestroyUINT4Vector( newindex );  
}


/** \brief Merge adjacent segments
 * 
 * This function will attempt to remerge adjacent segments if statistically 
 * favourable (as calculated by the odds ratio). For each pair of adjacent
 * segments the joint likelihood of them being from two independent
 * distributions is compared to the likelihood that combined they are from one
 * distribution. If the likelihood is highest for the combined segments they are
 * merged.
 * 
 * \param data [in] A complex data vector
 * \param segs [in] A vector of split segment indexes
 */
void merge_data( COMPLEX16Vector *data, UINT4Vector *segs ){
  UINT4 j = 0;
  REAL8 threshold = 0.; /* may need to be passed to function in the future, or
                           defined globally */
                           
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
      if( j == 1 ) cellstarts1 = 0;
      else cellstarts1 = segs->data[j-2];
        
      cellends1 = segs->data[j-1];
      
      cellstarts2 = segs->data[j-1];
      cellends2 = segs->data[j];
      
      n1 = cellends1 - cellstarts1;
      n2 = cellends2 - cellstarts2;
      nm = cellends2 - cellstarts1;
      
      for( i = cellstarts1; i < cellends1; i++ ){
        sum1 += SQUARE( data->data[i].re );
        sum1 += SQUARE( data->data[i].im );
      }
      
      for( i = cellstarts2; i < cellends2; i++ ){
        sum2 += SQUARE( data->data[i].re );
        sum2 += SQUARE( data->data[i].im );
      }
      
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
    if ( minl < threshold ) break;
    else{ /* merge cells */
      /* remove the cell end value between the two being merged and shift */   
      for( UINT4 i=0; i < ncells-(mergepoint+1); i++ )
        segs->data[mergepoint+i] = segs->data[mergepoint+i+1];
        
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
 * for each stationary segment given in the \c chunkLength vector. These value
 * are used in the likelihood calculation in \c pulsar_log_likelihood and are 
 * precomputed here to speed that calculation up. The vector of value is output
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
  
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams, 
                                                             "chunkLength" );
    
    length = data->dataTimes->length + 1 -
      chunkLengths->data[chunkLengths->length - 1];
    
    sumdat = XLALCreateREAL8Vector( chunkLengths->length );
    
    for( i = 0 ; i < length ; i+= chunkLength ){
      chunkLength = chunkLengths->data[count];
    
      sumdat->data[count] = 0.;
    
      for( j = i ; j < i + chunkLength ; j++){
        B.re = data->compTimeData->data->data[j].re;
        B.im = data->compTimeData->data->data[j].im;

        /* sum up the data */
        sumdat->data[count] += (B.re*B.re + B.im*B.im);
      }

      count++;
    }
    
    LALInferenceAddVariable( data->dataParams, "sumData", &sumdat, 
                             LALINFERENCE_REAL8Vector_t,
                             LALINFERENCE_PARAM_FIXED);
  
    data = data->next;
  }
  
  return;
}


/** \brief Creates a lookup table of the detector antenna pattern
 *
 * This function creates a 2D lookup table of the 'plus' and 'cross' antenna 
 * patterns for a given detector orientation and source sky position. The 
 * lookup table spans one sidereal day in time (this being the period over which
 * the antenna pattern changes) and goes between \f$\pm\pi/4\f$ radians in
 * \f$\psi\f$.
 * 
 * Note: the may want to be converted into an XLAL function and moved into LAL 
 * at some point.
 * 
 * \param t0 [in] initial GPS time of the data
 * \param detNSource [in] structure containing the detector and source
 * orientations and locations
 * \param timeSteps [in] the number of grid bins to use in time
 * \param psiSteps [in] the number of grid bins to use in polarisation angle
 * \f$\psi\f$
 * \param LUfplus [in] a matrix into which the 'plus' antenna pattern lookup
 * table will be output
 * \param LUfcross [in] a matrix into which the 'cross' antenna pattern lookup
 * table will be output 
 */
void response_lookup_table( REAL8 t0, LALDetAndSource detNSource, 
                            INT4 timeSteps, INT4 psiSteps, gsl_matrix *LUfplus,
                            gsl_matrix *LUfcross ){ 
  LIGOTimeGPS gps;
  REAL8 T = 0;

  REAL8 fplus = 0., fcross = 0.;
  REAL8 psteps = (REAL8)psiSteps;
  REAL8 tsteps = (REAL8)timeSteps;

  INT4 i = 0, j = 0;
  
  for( i = 0 ; i < psiSteps ; i++ ){
    detNSource.pSource->orientation = -(LAL_PI/4.) +
        (REAL8)i*(LAL_PI/2.) / ( psteps - 1. );

    for( j = 0 ; j < timeSteps ; j++ ){
      T = t0 + (REAL8)j*LAL_DAYSID_SI / tsteps;

      XLALGPSSetREAL8(&gps, T);

      XLALComputeDetAMResponse( &fplus, &fcross,
                                detNSource.pDetector->response,
                                detNSource.pSource->equatorialCoords.longitude,
                                detNSource.pSource->equatorialCoords.latitude,
                                detNSource.pSource->orientation,
                                XLALGreenwichMeanSiderealTime( &gps ) );
                                
      gsl_matrix_set( LUfplus, i, j, fplus );
      gsl_matrix_set( LUfcross, i, j, fcross );
    }
  }
}


/** \brief Rescale the value output by the Nested Sampling algorithm
 * 
 * This function reads in the file of samples output from the Nested Sampling
 * algorithm (in the file specified by \c outfile) and scales them back to 
 * their true values. It removes the string variable "model" from the output
 * and shifts the logPrior and logLikelihood values to the end of the parameter
 * list.
 * 
 * Note: The output may soon be in an XML format, so this function will need to
 * be amended. 
 * 
 * \param runState [in] The analysis information structure
 */
void rescaleOutput( LALInferenceRunState *runState ){
  /* Open original output output file */
  CHAR *outfile, outfiletmp[256] = "";
  CHAR outfilepars[256] = "", outfileparstmp[256] = "";
  FILE *fp = NULL, *fptemp = NULL, *fppars = NULL, *fpparstmp = NULL;
  
  LALStringVector *paramsStr = NULL;
  
  ProcessParamsTable *ppt = LALInferenceGetProcParamVal( runState->commandLine,
                                                         "--outfile" );
  if( !ppt ){
    fprintf(stderr,"Must specify --outfile <filename.dat>\n");
    exit(1);
  }
  outfile = ppt->value;
  
  /* set temporary file for re-writing out samples */
  sprintf(outfiletmp, "%s_tmp", outfile);
  
  /* open output file */
  if( (fp = fopen(outfile, "r")) == NULL ){
    fprintf(stderr, "Error... cannot open output file %s.\n", outfile);
    exit(3);
  }
  
  /* open temporary output file for reading */
  if( (fptemp = fopen(outfiletmp, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open temporary output file %s.\n",
            outfile);
    exit(3);
  }
 
  /* open file for printing out list of parameter names - this should already 
     exist */
  sprintf(outfilepars, "%s_params.txt", outfile);
  if( (fppars = fopen(outfilepars, "r")) == NULL ){
    fprintf(stderr, "Error... cannot open parameter name output file %s.\n",
            outfilepars);
    exit(3);
  }
  /* read in the parameter names and remove the "model" value */
  sprintf(outfileparstmp, "%s_params.txt_tmp", outfile);
  if( (fpparstmp = fopen(outfileparstmp, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open parameter name output file %s.\n",
            outfileparstmp);
    exit(3);
  }
  
  CHAR v[128] = "";
  while( fscanf(fppars, "%s", v) != EOF ){
      paramsStr = XLALAppendString2Vector( paramsStr, v );
    
    /* re-output everything but the "model" value to a temporary file */
    if( strcmp(v, "model") != 0 || strcmp(v, "logL")!=0 
      || strcmp(v, "logPrior") )
      fprintf(fpparstmp, "%s\t", v);
  }
  
  /* we will put the logPrior and logLikelihood at the end of the lines */
  fprintf(fpparstmp, "logPrior\tlogL\n");
  
  fclose(fppars);
  fclose(fpparstmp);
  
  /* move the temporary file name to the standard outfile_param name */
  rename( outfileparstmp, outfilepars );
  
  while ( 1 ){
    UINT4 i = 0;
    
    REAL8 logPrior = 0., logL = 0.;
    
    /* scan through line, get value and reprint out scaled value to temporary
       file */
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
      
        scalefac = 
          *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, 
                                             scalename );
        scalemin = 
          *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, 
                                              scaleminname );
      
        fprintf(fptemp, "%.12le", atof(value)*scalefac + scalemin);
      }
      else if( !strcmp(paramsStr->data[i], "logL") )
        logL = atof(value);
      else if( !strcmp(paramsStr->data[i], "logPrior") )
        logPrior = atof(value);
        
      fprintf(fptemp, "\t");
    }
    
    if( feof(fp) ) break;
    
    /* print out the last two items to be the logPrior and logLikelihood */
    fprintf(fptemp, "%lf\t%lf\n", logPrior, logL);
  }
  
  fclose(fp);
  fclose(fptemp);
  
  /* move the temporary file name to the standard outfile name */
  rename( outfiletmp, outfile );
  
  return;
}


/** \brief Counts the number of comma separated values in a string
 *
 * This function counts the number of comma separated values in a given input
 * string.
 * 
 * \param csvline [in] Any string
 * 
 * \return The number of comma separated value in the input string
 */
INT4 count_csv( CHAR *csvline ){
  CHAR *inputstr = NULL;
  CHAR *tempstr = NULL;
  INT4 count = 0;
    
  inputstr = XLALStringDuplicate( csvline );

  /* count number of commas */
  while(1){
    tempstr = strsep(&inputstr, ",");
      
    if ( inputstr == NULL ) break;
      
    count++;
  }
    
  return count+1;
}


/** \brief Checks if a given parameter is recognised
 * 
 * This function checks whether a given parameter is one of the defined 
 * amplitude (\c amppars), frequency (\c freqpars), sky location (\c skypars)
 * or binary system (\c binpars) parameters given in the header file.
 * 
 * \param parname [in] The name of a parameter
 * 
 * \return true (1) if the parameter is recognised and false (0) if not
 */
INT4 recognised_parameter( CHAR *parname ){
  INT4 i = 0;
  
  for( i = 0; i < NUMAMPPARS; i++ )
    if (!strcmp(parname, amppars[i])) return 1;
    
  for( i = 0; i < NUMFREQPARS; i++ )
    if (!strcmp(parname, freqpars[i])) return 1;
    
  for( i = 0; i < NUMSKYPARS; i++ )
    if (!strcmp(parname, skypars[i])) return 1;
    
  for( i = 0; i < NUMBINPARS; i++ )
    if (!strcmp(parname, binpars[i])) return 1;
    
  return 0;
}


/** \brief Calculates the optimal matched filter signal-to-noise ratio for a
 * given signal 
 * 
 * This function calculates the optimal matched filter signal-to-noise ratio
 * (SNR) of a given signal model for a set of detector data via:
 * \f[
 * \rho = \sqrt{\sum_{i=1}^N \frac{d_i^2}{\sigma^2}},
 * \f]
 * where \f$\{d\}\f$ is a time series of data, and \f$\sigma^2\f$ is its
 * variance. As the data and model used here are complex the real and imaginary
 * SNRs are added in quadrature to give the total SNR.
 * 
 * The data variance \f$\sigma^2\f$ is calculated on data that has had the
 * running median subtracted in order to remove any underlying trends (e.g.
 * caused by a string signal). The variance is assumed constant over segments
 * given in the \c chunkLength vector and the SNR from each segment is added in
 * quadrature.
 * 
 * \param data [in] A data pointer containing detector data and the signal model
 * 
 * \return The optimal matched filter signal-to-noise ratio
 */
REAL8 calculate_time_domain_snr( LALInferenceIFOData *data ){
  REAL8 snrval = 0., chunkLength;
  COMPLEX16 snrc = {0., 0.}, vari = {0., 0.};
                                                   
  INT4 i = 0, j = 0, length = 0, cl = 0;
  
  COMPLEX16Vector *meddata = NULL;
  INT4 chunkMin = 0, count = 0;
  
  /* subtract a running median value from the data to remove any underlying
     trends (e.g. caused by a string signal) */
  meddata = subtract_running_median( data->compTimeData->data );
  
  UINT4Vector *chunkLengths = NULL;
  chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams,
                                                             "chunkLength" );
  chunkMin = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMin" );
  
  length = data->compTimeData->data->length;
  
  /* add the signal to the data */
  for ( i = 0; i < length; i+=chunkLength ){
    chunkLength = (REAL8)chunkLengths->data[count];
    
    /* skip section of data if its length is less than the minimum allowed
      chunk length */
    if( chunkLength < chunkMin ){
      count++;
      continue;
    }
    
    cl = i + (INT4)chunkLength;
    
    for( j = i ; j < cl ; j++ ){
      vari.re += SQUARE(meddata->data[j].re);
      vari.im += SQUARE(meddata->data[j].im);
      
      /* calculate optimal signal power */
      snrc.re += SQUARE(data->compModelData->data->data[j].re);
      snrc.im += SQUARE(data->compModelData->data->data[j].im);
    }
    
    vari.re /= (chunkLength - 1.);
    vari.im /= (chunkLength - 1.);
    
    /* add SNRs for each chunk in quadrature */
    snrval += (snrc.re/vari.re) + (snrc.im/vari.im);
    
    count++;
  }
  
  snrval = sqrt( snrval );
   
  return snrval;
}


/** \brief Get the signal-to-noise ratio of the maximum likelihood signal
 * 
 * The function uses the signal with the highest likelihood (which will be the 
 * final point in the live points array) and calculates the optimal 
 * signal-to-noise ratio (SNR) for it. This is output to a file based on the 
 * \c outfile value, but with \c _SNR appended to it. For multiple detector, 
 * and/or models with multiple data sets, the individual detector/data set SNR
 * values will be output, with the final value being the multi-detector SNR. If
 * a fake signal has been injected into the data this file will already
 * contain the optimal SNR of the true signal. 
 * 
 * \param runState [in] The analysis information structure
 * 
 * \sa calculate_time_domain_snr
 */
void get_loudest_snr( LALInferenceRunState *runState ){
  REAL8 lMax = 0.;
  INT4 ndets = 0;
  INT4 Nlive = *(INT4 *)LALInferenceGetVariable( runState->algorithmParams,
                                                 "Nlive" );
  REAL8 snrmulti = 0.;
  
  CHAR *snrfile = NULL;
  FILE *fpsnr = NULL;
  
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  
  LALInferenceVariables *loudestParams = NULL;
  LALInferenceIFOData *data;
  
  loudestParams = calloc( 1, sizeof(LALInferenceVariables) );
 
  /* max likelihood point should have been sorted to be the final value */
  LALInferenceCopyVariables( runState->livePoints[Nlive-1], loudestParams );
  
  lMax = runState->likelihood( loudestParams, runState->data,
                               runState->template );
 
  LALInferenceDestroyVariables( loudestParams );
  
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
  
  while( data ){
    REAL8 snrval = calculate_time_domain_snr( data );
    
    ndets++;
    
    snrmulti += SQUARE( snrval );
    
    /* print out SNR value */
    fprintf(fpsnr, "%le\t", snrval);
    
    data = data->next;
  }
 
  /* print out multi-detector SNR value */
  if ( ndets > 1 ) fprintf(fpsnr, "%le\n", sqrt(snrmulti));
  else fprintf(fpsnr, "\n");
  
  fclose( fpsnr );
}


/** \brief Automatically set the solar system ephemeris file based on
 * environment variables and data time span
 * 
 * This function will attempt to find Earth and Sun ephemeris files based on
 * LAL environment variables (as set up by <code> lalpulsar-user-env.(c)sh
 * </code>) and a given start and end GPS time (presumably taken from the data
 * that is to be analysed). It requires \c LALPULSAR is installed and the \c
 * LALPULSAR_PREFIX variable is set, which should mean that ephemeris files are
 * installed in the directory \c ${LALPULSAR_PREFIX}/share/lalpulsar/.
 *
 * Within this directory the should be ephemeris files of the format 
 * earthXX-YY.dat, or earthXX.dat, where XX are two digit years e.g.
 * earth11.dat and represent the year or years for which the ephemeris is
 * valid. The function will find the appropriate ephemeris files that cover the 
 * time range (in GPS seconds) required. If no files exist errors will be
 * returned.
 * 
 * The function uses requires <code>dirent,h</code> and <code>time.h</code>. 
 * 
 * NOTE: This may want to be moved into LAL at some point.
 * 
 * \param efile [in] a string that will return the Earth ephemeris file
 * \param sfile [in] a string that will return the Sun ephemeris file 
 * \param gpsstart [in] the GPS time of the start of the data
 * \param gpsend [in] the GPS time of the end of the data
 * 
 * \return Zero will be return on successful completion
 */
INT4 XLALAutoSetEphemerisFiles( CHAR *efile, CHAR *sfile, 
                                INT4 gpsstart, INT4 gpsend ){
  struct tm utcstart,  utcend;
  CHAR *eftmp = NULL, *sftmp = NULL;
  INT4 buf = 3;
  CHAR yearstart[buf], yearend[buf]; /* the two digit year i.e. 11 for 2011*/
  INT4 ys = 0, ye = 0;
  CHAR *lalpath = NULL, *lalpulsarpath = NULL;
  
  INT4 yr1 = 0, yr2 = 0, i = 0;
    
  CHAR tmpyr[6];
    
  struct dirent *entry;
  DIR *dp;
  
  /* first check that the path to the Ephemeris files is available in the
     environment variables */
  if((lalpath = getenv("LALPULSAR_PREFIX")) == NULL){
    XLALPrintError("LALPULSAR_PREFIX environment variable not set. Cannot \
automatically generate ephemeris files!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }
  
  /* get the utc times of the start and end GPS times */
  if( !XLALGPSToUTC( &utcstart, gpsstart ) )
    XLAL_ERROR(XLAL_EFUNC);
  if( !XLALGPSToUTC( &utcend, gpsend ) )
    XLAL_ERROR(XLAL_EFUNC);
  
  /* get years */
  if( strftime(yearstart, buf, "%y", &utcstart) != 2 )
    XLAL_ERROR(XLAL_EFUNC);
  
  if( strftime(yearend, buf, "%y", &utcend) != 2 )
    XLAL_ERROR(XLAL_EFUNC);
  
  ys = atoi(yearstart);
  ye = atoi(yearend);
  
  lalpulsarpath = XLALStringDuplicate( lalpath );
  
  if ( (lalpulsarpath = XLALStringAppend(lalpulsarpath, "/share/lalpulsar/")) ==
NULL )
    XLAL_ERROR(XLAL_EFUNC);
  
  eftmp = XLALStringDuplicate(lalpulsarpath);
  sftmp = XLALStringDuplicate(lalpulsarpath);
  
  eftmp = XLALStringAppend(eftmp, "earth");
  sftmp = XLALStringAppend(sftmp, "sun");
 
  /* find the ephemeris file that bounds the required range */
  if ( ( dp = opendir(lalpulsarpath) ) == NULL ){
    XLALPrintError("Error... cannot open directory path %s!\n", lalpulsarpath);
    XLAL_ERROR(XLAL_EFUNC);
  }
    
  while( (entry = readdir(dp) ) ){
    /* just use "earth" files rather than doubling up */
    if ( strstr(entry->d_name, "earth") != NULL ){
      /* get current set of ranges - filenames are of the format
         earthXX-YY.dat, so find the '-' and extract two characters either
         side */
      if ( strchr(entry->d_name, '-') ){
        sscanf(entry->d_name, "earth%d-%d.dat", &yr1, &yr2);
        
        if ( ys >= yr1 && ye <= yr2 ) {
          snprintf(tmpyr, 6, "%02d-%02d", yr1, yr2);
          
          eftmp = XLALStringAppend(eftmp, tmpyr);
          sftmp = XLALStringAppend(sftmp, tmpyr);
          i++;
          break;
        }
      }
      /* get the single year value ephemeris file e.g. earthXX.dat */
      else{
        sscanf(entry->d_name, "earth%d.dat", &yr1);
        
        if ( ys == yr1 && ye == yr1 ){
          snprintf(tmpyr, 3, "%02d", yr1);
          eftmp = XLALStringAppend(eftmp, tmpyr);
          sftmp = XLALStringAppend(sftmp, tmpyr);
          i++;
          break;
        }
      }
    }
  }
    
  closedir(dp);
    
  if( i == 0 ){
    XLALPrintError("No ephemeris files in the time range %02d-%02d found!\n",
                   ys, ye);
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* add .dat extension */
  eftmp = XLALStringAppend(eftmp, ".dat");
  sftmp = XLALStringAppend(sftmp, ".dat");
  
  if ( eftmp == NULL || sftmp == NULL )
    XLAL_ERROR(XLAL_EFUNC);

  XLALStringCopy( efile, eftmp, 1024 );
  XLALStringCopy( sfile, sftmp, 1024 );
  
  /* double check that the files exist */
  if( access(sfile, F_OK) != 0 || access(efile, F_OK) != 0 ){
    XLALPrintError("Error... ephemeris files not, or incorrectly, defined!\n");
    XLAL_ERROR(XLAL_EFUNC);
  }
  
  return 0;
}


/*----------------------- END OF HELPER FUNCTIONS ----------------------------*/

/******************************************************************************/
/*                          TESTING FUNCTIONS                                 */
/******************************************************************************/

/** \brief A test function to calculate a 1D posterior on a grid
 * 
 * This function is only to be used as a check/test of the code and will be run
 * if the \c grid command line argument is present. It will calculate the
 * posterior for one parameter (given by \c gridpar), between the ranges given 
 * by \c gridmin and \c gridmax (which default to 0 and 1) at a number of points
 * given by \c gridsteps (which default to 100).
 * 
 * \param runState [in] The analysis information structure
 */
void gridOutput( LALInferenceRunState *runState ){
  REAL8 h0min = 0.;
  REAL8 h0max = 0.;
  REAL8 h0range = 0, h0step = 0;
  INT4 h0steps = 0, i = 0;
 
  ProcessParamsTable *ppt;
  REAL8 scaleval = 1., minval = 0., tmpscale = 0., tmpmin = 0., tmpgridval = 0.;
  
  ProcessParamsTable *commandLine = runState->commandLine;
  
  FILE *fp = NULL;
  REAL8 minL = LAL_REAL8_MAX;
  REAL8 sumPost = 0.;
  
  REAL8Vector *logL = NULL;
  
  CHAR *parname = NULL, parscale[256], parmin[256], outputgrid[256];
  
  /*------------------------------------------------------------*/
  /* test output on a h0 grid */
  ppt = LALInferenceGetProcParamVal( commandLine, "--grid" );
  if ( ppt ){
    ProcessParamsTable *ppt2;
    
    /* parameters over which to perform the grid search */
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridpar" );
    
    if( ppt2 ){
      parname = XLALStringDuplicate( ppt2->value );
        
      if( !recognised_parameter( parname ) ){
        fprintf(stderr, "Error... parameter %s not recognised\n", parname );
        exit(0);
      }
        
      sprintf(parscale, "%s_scale", parname);
      sprintf(parmin, "%s_scale_min", parname);
    }
    else{
      fprintf(stderr, USAGEGRID, commandLine->program);
      exit(0);
    }
    
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridmin" );
    
    if( ppt2 ) h0min = atof( ppt2->value );
    else h0min = 0.; /* default to zero */
    
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridmax" );
    
    if( ppt2 ) h0max = atof( ppt2->value );  
    else h0max = 1.; /* default to 1 */
    
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridsteps" );
    
    if( ppt2 ) h0steps = atoi( ppt2->value );
    else h0steps = 100; /* default to 100 steps */
  }
  else{
    return;
  }
  
  if ( verbose ){
    fprintf(stderr, "Calculating posterior on %s over a grid from:\n", parname);
    fprintf(stderr, "\t%le --> %le in %d steps.\n", h0min, h0max, h0steps);
  }
  
  h0range = h0max - h0min;
  h0step = h0range / (REAL8)(h0steps-1.);
  
  logL = XLALCreateREAL8Vector( h0steps );
  
  /* reset rescale value for h0 */
  tmpscale = *(REAL8*)LALInferenceGetVariable( runState->data->dataParams,
                                               parscale );
  tmpmin = *(REAL8*)LALInferenceGetVariable( runState->data->dataParams,
                                             parmin );
  LALInferenceRemoveVariable( runState->data->dataParams, parscale );
  LALInferenceAddVariable( runState->data->dataParams, parscale, &scaleval,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceRemoveVariable( runState->data->dataParams, parmin );
  LALInferenceAddVariable( runState->data->dataParams, parmin, &minval,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  
  tmpgridval = *(REAL8*)LALInferenceGetVariable( runState->currentParams,
                                                 parname );
  
  sprintf(outputgrid, "%s_grid_posterior.txt", parname);
  
  if ( (fp = fopen(outputgrid, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open grid posterior file %s.\n",
            outputgrid);
    exit(0);
  }
  
  for( i = 0; i < h0steps; i++ ){
    REAL8 h0val = h0min + i*h0step;
    
    LALInferenceSetVariable( runState->currentParams, parname, &h0val );
    
    logL->data[i] = runState->likelihood( runState->currentParams,
                                          runState->data, runState->template );
    
    if ( logL->data[i] < minL ) minL = logL->data[i];
  }
  
  /* integrate area under posterior - trapezium rule */
  for( i = 0; i < h0steps-1; i++ ){
    sumPost += ( exp(logL->data[i] - minL) + exp(logL->data[i+1] - minL) ) *
      h0step / 2.;
  }
  
  /* output posterior */
  for( i = 0; i < h0steps; i++ ){
    REAL8 h0val = h0min + i*h0step;
    fprintf(fp, "%le\t%le\n", h0val, exp( logL->data[i] - minL ) / sumPost);
  }
  
  fclose(fp);
  
  XLALDestroyREAL8Vector( logL );
  
  /* reset scale value and parameter value in currentParams */
  LALInferenceRemoveVariable( runState->data->dataParams, parscale );
  LALInferenceAddVariable( runState->data->dataParams, parscale, &tmpscale,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  
  LALInferenceRemoveVariable( runState->data->dataParams, parmin );
  LALInferenceAddVariable( runState->data->dataParams, parmin, &tmpmin,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
                           
  LALInferenceSetVariable( runState->currentParams, parname, &tmpgridval );
}


/** \brief Test the sampler using a Gaussian likelihood
 * 
 * This is a testing function that can be substituted for the standard 
 * likelihood function. It calculates only the \c h0 parameter posterior based
 * on a Gaussian likelihood with mean of 0.5 and standard deviation of 0.025 - 
 * these values can be changed if required. It is just to be used to test the
 * sampling routine (e.g. Nested Sampling) with a well defined likelihood
 * function.
 * 
 * \param vars [in] A set of pulsar parameters
 * \param data [in] A data structure
 * 
 * \return Natural logarithm of the likelihood
 */
REAL8 test_gaussian_log_likelihood( LALInferenceVariables *vars,
                                    LALInferenceIFOData *data,
                                    LALInferenceTemplateFunction UNUSED
                                      *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  
  REAL8 like_mean = 0.5;
  REAL8 like_sigma = 0.025;
  REAL8 h0 = *(REAL8 *)LALInferenceGetVariable( vars, "h0" );
  REAL8 h0scale = *(REAL8 *)LALInferenceGetVariable( data->dataParams,
                                                     "h0_scale" );
  REAL8 h0min = *(REAL8 *)LALInferenceGetVariable( data->dataParams,
                                                   "h0_scale_min" );
  
  get_model = NULL;
                                                   
  h0 = h0*h0scale + h0min;

  /* search over a simple 1D Gaussian with x defined by the h0 variable */
  loglike = -log(sqrt(2.*LAL_PI)*like_sigma);
  loglike -= (h0-like_mean)*(h0-like_mean) / (2.*like_sigma*like_sigma);
  
  return loglike;
}


/*----------------------- END OF TESTING FUNCTIONS ---------------------------*/
