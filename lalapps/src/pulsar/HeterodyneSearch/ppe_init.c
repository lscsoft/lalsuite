/*
*  Copyright (C) 2014 Matthew Pitkin
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

/******************************************************************************/
/*                      INITIALISATION FUNCTIONS                              */
/******************************************************************************/

#include "ppe_init.h"

/** Array for conversion from lowercase to uppercase */
static const CHAR a2A[256] = {
  ['a'] = 'A', ['b'] = 'B', ['c'] = 'C', ['d'] = 'D', ['e'] = 'E', ['f'] = 'F', ['g'] = 'G', ['h'] = 'H',
  ['i'] = 'I', ['j'] = 'J', ['k'] = 'K', ['l'] = 'L', ['m'] = 'M', ['n'] = 'N', ['o'] = 'O', ['p'] = 'P',
  ['q'] = 'Q', ['r'] = 'R', ['s'] = 'S', ['t'] = 'T', ['u'] = 'U', ['v'] = 'V', ['w'] = 'W', ['x'] = 'X',
  ['y'] = 'Y', ['z'] = 'Z' };


/** \brief Convert string to uppercase */
static void strtoupper(CHAR *s) {
  /* convert everything to uppercase */
  CHAR c;
  for ( ; *s; ++s ) {
    if ( (c = a2A[(int)(*s)]) ) { *s = c; }
  }
}


/**
 * \brief A wrapper around \c LALInferenceNestedSamplingAlgorithm
 *
 * This function just calls \c LALInferenceNestedSamplingAlgorithm, but will time the algorithm
 * if required.
 *
 * \param runState [] A pointer to the \c LALInferenceRunState
 */
void nested_sampling_algorithm_wrapper( LALInferenceRunState *runState ){
  /* timing values */
  struct timeval time1, time2;
  REAL8 tottime;

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){ gettimeofday(&time1, NULL); }

  LALInferenceNestedSamplingAlgorithm(runState);

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){
    gettimeofday(&time2, NULL);

    FILE *timefile = *(FILE **)LALInferenceGetVariable( runState->algorithmParams, "timefile" );
    UINT4 timenum = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams, "timenum" );
    tottime = (REAL8)((time2.tv_sec + time2.tv_usec*1.e-6) - (time1.tv_sec + time1.tv_usec*1.e-6));
    fprintf(timefile, "[%d] %s: %.9le secs\n", timenum, __func__, tottime);
    timenum++;
    check_and_add_fixed_variable( runState->algorithmParams, "timenum", &timenum, LALINFERENCE_UINT4_t );
  }
}


/**
 * \brief A wrapper around \c LALInferenceSetupLivePointsArray
 *
 * This function just calls \c LALInferenceSetupLivePointsArray, but will time the algorithm
 * if required.
 *
 * \param runState [] A pointer to the \c LALInferenceRunState
 */
void setup_live_points_array_wrapper( LALInferenceRunState *runState ){
  /* timing values */
  struct timeval time1, time2;
  REAL8 tottime;

  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){ gettimeofday(&time1, NULL); }

  LALInferenceSetupLivePointsArray( runState );

  /* note that this time divided by the number of live points will give the approximate time per likelihood evaluation */
  if ( LALInferenceCheckVariable( runState->algorithmParams, "timefile" ) ){
    gettimeofday(&time2, NULL);

    FILE *timefile = *(FILE **)LALInferenceGetVariable( runState->algorithmParams, "timefile" );
    UINT4 timenum = *(UINT4 *)LALInferenceGetVariable( runState->algorithmParams, "timenum" );
    tottime = (REAL8)((time2.tv_sec + time2.tv_usec*1.e-6) - (time1.tv_sec + time1.tv_usec*1.e-6));
    fprintf(timefile, "[%d] %s: %.9le secs\n", timenum, __func__, tottime);
    timenum++;
    check_and_add_fixed_variable( runState->algorithmParams, "timenum", &timenum, LALINFERENCE_UINT4_t );
  }
}


/**
 * \brief Initialises the nested sampling algorithm control
 *
 * Memory is allocated for the parameters, priors and proposals. The nested sampling control parameters are set: the
 * number of live points \c Nlive, the number of points for each MCMC \c Nmcmc, the number of independent runs within
 * the algorithm \c Nruns, and the stopping criterion \c tolerance.
 *
 * The random number generator is initialise (the GSL Mersenne Twister algorithm \c gsl_rng_mt19937) using either a user
 * defined seed \c randomseed, the system defined \c /dev/random file, or the system clock time.
 *
 * \param runState [in] A pointer to the \c LALInferenceRunState
 */
void initialise_algorithm( LALInferenceRunState *runState )
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
    LALInferenceAddVariable( runState->algorithmParams, "tolerance", &tmp, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
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

  /* check if we want to time the program */
  ppt = LALInferenceGetProcParamVal( commandLine, "--time-it" );
  if ( ppt != NULL ){
    FILE *timefile = NULL;
    UINT4 timenum = 1;
    ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );

    if ( !ppt ){ XLAL_ERROR_VOID(XLAL_EFUNC, "Error... no output file is specified!"); }

    CHAR outtimefile[256] = "";
    sprintf(outtimefile, "%s_timings", ppt->value);

    if ( ( timefile = fopen(outtimefile, "w") ) == NULL ){
      fprintf(stderr, "Warning... cannot create a timing file, so proceeding without timings\n");
    }
    else{
      LALInferenceAddVariable( runState->algorithmParams, "timefile", &timefile, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( runState->algorithmParams, "timenum", &timenum, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
    }
  }

  gsl_rng_set( runState->GSLrandom, randomseed );

  /* log samples */
#ifdef HAVE_LIBLALXML
  runState->logsample = LogSampleToArray;
#else
  runState->logsample = LogSampleToFile;
#endif

  return;
}


/**
 * \brief Sets the time angle antenna response lookup table
 *
 * This function sets up an antenna repsonse lookup table in time for each detector from which data
 * exists (either real or fake). The time ranges over one sidereal day. The number of bins for the grid
 * over time can be specified on the command line via \c time-bins, but if this are not given then default
 * values are used. The data times as a fraction of a sidereal day from the start time will also be
 * calculated.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 * \param source [in] A pointer to a LALSource variable containing the source location
 */
void setup_lookup_tables( LALInferenceRunState *runState, LALSource *source ){
  /* Set up lookup tables */
  /* Using psi bins, time bins */
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;

  LALInferenceIFOData *data = runState->data;
  LALInferenceIFOModel *ifo_model = runState->model->ifo;

  REAL8 t0;
  LALDetAndSource detAndSource;

  ppt = LALInferenceGetProcParamVal( commandLine, "--time-bins" );
  INT4 timeBins;
  if( ppt ) { timeBins = atoi( ppt->value ); }
  else { timeBins = TIMEBINS; } /* default time bins */

  while( ifo_model ){
    REAL8Vector *arespT = NULL, *brespT = NULL;
    REAL8Vector *arespV = NULL, *brespV = NULL;
    REAL8Vector *arespS = NULL, *brespS = NULL;

    REAL8Vector *sidDayFrac = NULL;
    REAL8 dt = 0;
    UINT4 i = 0;

    LALInferenceAddVariable( ifo_model->params, "timeSteps", &timeBins, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED );

    t0 = XLALGPSGetREAL8( &data->compTimeData->epoch );

    sidDayFrac = XLALCreateREAL8Vector( ifo_model->times->length );

    /* set the time in sidereal days since the first data point (mod 1 sidereal day) */
    for( i = 0; i < ifo_model->times->length; i++ ){
      sidDayFrac->data[i] = fmod( XLALGPSGetREAL8( &ifo_model->times->data[i] ) - t0, LAL_DAYSID_SI );
    }

    LALInferenceAddVariable( ifo_model->params, "siderealDay", &sidDayFrac, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

    detAndSource.pDetector = data->detector;
    detAndSource.pSource = source;

    arespT = XLALCreateREAL8Vector( timeBins );
    brespT = XLALCreateREAL8Vector( timeBins );
    arespV = XLALCreateREAL8Vector( timeBins );
    brespV = XLALCreateREAL8Vector( timeBins );
    arespS = XLALCreateREAL8Vector( timeBins );
    brespS = XLALCreateREAL8Vector( timeBins );

    dt = LALInferenceGetREAL8Variable( ifo_model->params, "dt" );

    response_lookup_table( t0, detAndSource, timeBins, dt, arespT, brespT, arespV, brespV, arespS, brespS );

    LALInferenceAddVariable( ifo_model->params, "a_response_tensor", &arespT, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( ifo_model->params, "b_response_tensor", &brespT, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

    if ( LALInferenceCheckVariable( ifo_model->params, "nonGR" ) ){
      LALInferenceAddVariable( ifo_model->params, "a_response_vector", &arespV, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( ifo_model->params, "b_response_vector", &brespV, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( ifo_model->params, "a_response_scalar", &arespS, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
      LALInferenceAddVariable( ifo_model->params, "b_response_scalar", &brespS, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    }
    else{
      XLALDestroyREAL8Vector( arespV );
      XLALDestroyREAL8Vector( brespV );
      XLALDestroyREAL8Vector( arespS );
      XLALDestroyREAL8Vector( brespS );
    }

    data = data->next;
    ifo_model = ifo_model->next;
  }

  return;
}


/**
 * \brief Set up all the allowed variables for a known pulsar search
 * This functions sets up all possible variables that are possible in a known pulsar search. Parameter values read in
 * from a .par file and passed in via the \c pars variable will be set.
 *
 * \param ini [in] A pointer to a \c LALInferenceVariables type that will be filled in with pulsar parameters
 * \param pars [in] A \c BinaryPulsarParams type containing pulsar parameters read in from a TEMPO-style .par file
 */
void add_initial_variables( LALInferenceVariables *ini, PulsarParameters *pars ){
  /* amplitude model parameters for l=m=2 harmonic emission */
  add_variable_parameter( pars, ini, "H0", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PHI0", LALINFERENCE_PARAM_FIXED ); /* note that this is rotational phase */
  add_variable_parameter( pars, ini, "COSIOTA", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "IOTA", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PSI", LALINFERENCE_PARAM_FIXED );

  /* amplitude model parameters for l=2, m=1 and 2 harmonic emission from Jones (2010) */
  add_variable_parameter( pars, ini, "I21", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "I31", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "LAMBDA", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "COSTHETA", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "THETA", LALINFERENCE_PARAM_FIXED );

  /* amplitude model parameters in phase and amplitude form */
  add_variable_parameter( pars, ini, "C22", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "C21", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PHI22", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PHI21", LALINFERENCE_PARAM_FIXED );

  /***** phase model parameters ******/
  /* frequency */
  add_variable_parameter( pars, ini, "F0", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "F1", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "F2", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "F3", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "F4", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "F5", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PEPOCH", LALINFERENCE_PARAM_FIXED );

  /* add non-GR parameters */
  add_variable_parameter( pars, ini, "CGW", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "HPLUS", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "HCROSS", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PSITENSOR", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PHI0TENSOR", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "HSCALARB", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "HSCALARL", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PSISCALAR", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PHI0SCALAR", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "HVECTORX", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "HVECTORY", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PSIVECTOR", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PHI0VECTOR", LALINFERENCE_PARAM_FIXED );

  /* sky position */
  REAL8 ra = 0.;
  if ( PulsarCheckParam( pars, "RA" ) ) { ra = PulsarGetREAL8Param( pars, "RA" ); }
  else if ( PulsarCheckParam( pars, "RAJ" ) ) { ra = PulsarGetREAL8Param( pars, "RAJ" ); }
  else {
    XLALPrintError ("%s: No source right ascension specified!", __func__ );
    XLAL_ERROR_VOID( XLAL_EINVAL );
  }
  REAL8 dec = 0.;
  if ( PulsarCheckParam( pars, "DEC" ) ) { dec = PulsarGetREAL8Param( pars, "DEC" ); }
  else if ( PulsarCheckParam( pars, "DECJ" ) ) { dec = PulsarGetREAL8Param( pars, "DECJ" ); }
  else {
    XLALPrintError ("%s: No source declination specified!", __func__ );
    XLAL_ERROR_VOID( XLAL_EINVAL );
  }
  LALInferenceAddVariable( ini, "RA", &ra, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( ini, "DEC", &dec, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  add_variable_parameter( pars, ini, "PMRA", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "PMDEC", LALINFERENCE_PARAM_FIXED );
  add_variable_parameter( pars, ini, "POSEPOCH", LALINFERENCE_PARAM_FIXED );

  /* only add binary system parameters if required */
  if ( PulsarCheckParam( pars, "BINARY" ) ){
    CHAR *binary = XLALStringDuplicate(PulsarGetStringParam(pars, "BINARY"));
    LALInferenceAddVariable( ini, "BINARY", &binary, LALINFERENCE_string_t, LALINFERENCE_PARAM_FIXED );

    add_variable_parameter( pars, ini, "PB", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "ECC", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "EPS1", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "EPS2", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "T0", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "TASC", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "A1", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "OM", LALINFERENCE_PARAM_FIXED );

    add_variable_parameter( pars, ini, "PB_2", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "ECC_2", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "T0_2", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "A1_2", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "OM_2", LALINFERENCE_PARAM_FIXED );

    add_variable_parameter( pars, ini, "PB_3", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "ECC_3", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "T0_3", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "A1_3", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "OM_3", LALINFERENCE_PARAM_FIXED );

    add_variable_parameter( pars, ini, "XPBDOT", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "EPS1DOT", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "EPS2DOT", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "OMDOT", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "GAMMA", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "PBDOT", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "XDOT", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "EDOT", LALINFERENCE_PARAM_FIXED );

    add_variable_parameter( pars, ini, "SINI", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "DR", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "DTHETA", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "A0", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "B0", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "MTOT", LALINFERENCE_PARAM_FIXED );
    add_variable_parameter( pars, ini, "M2", LALINFERENCE_PARAM_FIXED );

    if ( PulsarCheckParam(pars, "FB") ){
      REAL8Vector *fb = NULL;
      fb = PulsarGetREAL8VectorParam( pars, "FB" );
      LALInferenceAddVariable( ini, "FB", &fb, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    }
  }
}


/**
 * \brief Sets up the parameters to be searched over and their prior ranges
 *
 * This function sets up any parameters that you require the code to search over and specifies the prior range and type
 * for each. This information is contained in a prior file specified by the command line argument \c prior-file. This
 * file should contain four columns: the first has the name of a parameter to be searched over; the second has the prior
 * type (e.g. "uniform" for a prior that is flat over the given range, or "gaussian" with a certain mean and standard
 * deviation, or "predefined", which means that the prior for that variable is already hardcoded into the prior function);
 * the third has the lower limit, or mean, of the prior, for "uniform"/"predefined" and "gaussian" priors respectively;
 * and the fourth has the upper limit, or standard deviation, for "uniform"/"predefined" and "gaussian" priors respectively.
 * E.g.
 * \code
 * H0 uniform 0 1e-21
 * PHI0 uniform 0 3.14159265359
 * COSIOTA uniform -1 1
 * PSI uniform -0.785398163397448 0.785398163397448
 * \endcode
 * or
 * \code
 * H0 uniform 0 1e-21
 * PHI0 uniform 0 3.14159265359
 * IOTA predefined 0 3.14159265359
 * PSI uniform -0.785398163397448 0.785398163397448
 * \endcode
 *
 * Any parameter specified in the file will have its vary type set to \c LALINFERENCE_PARAM_LINEAR.
 *
 * If a parameter correlation matrix is given by the \c cor-file command then this is used to construct a multi-variate
 * Gaussian prior for the given parameters (it is assumed that this file is created using TEMPO and the parameters it
 * contains are the same as those for which a standard deviation is defined in the par file). This overrules the
 * Gaussian priors that will have been set for these parameters.
 *
 * \param runState [in] A pointer to the LALInferenceRunState
 */
void initialise_prior( LALInferenceRunState *runState )
{
  CHAR *propfile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  FILE *fp=NULL;

  CHAR tempPar[VARNAME_MAX] = "", tempPrior[VARNAME_MAX] = "";
  REAL8 low, high;

  LALInferenceIFOModel *ifo = runState->model->ifo;

  /* parameters in correlation matrix */
  LALStringVector *corParams = NULL;
  REAL8Array *corMat = NULL;

  INT4 varyphase = 0, varyskypos = 0, varybinary = 0;

  /* check if non-GR parameters are going to be used */
  UINT4 nonGR = 0;
  if ( LALInferenceGetProcParamVal( commandLine, "--nonGR" ) ){ nonGR = 1; }

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
    INT4 isthere = 0, i = 0;

    LALInferenceParamVaryType varyType;

    /* convert tempPar to all uppercase letters */
    strtoupper( tempPar );

    if ( !strcmp(tempPrior, "uniform") || !strcmp(tempPrior, "predefined") ){
      if( high < low ){
        XLALPrintError("Error... In %s the %s parameters ranges are wrongly set.\n", propfile, tempPar);
        XLAL_ERROR_VOID( XLAL_EINVAL );
      }
    }

    if ( strcmp(tempPrior, "uniform") && strcmp(tempPrior, "predefined") && strcmp(tempPrior, "gaussian") ){
      XLALPrintError("Error... prior type '%s' not recognised\n", tempPrior);
      XLAL_ERROR_VOID( XLAL_EINVAL );
    }

    /* set variable type to LINEAR (as they are initialised as FIXED) */
    varyType = LALINFERENCE_PARAM_LINEAR;
    LALInferenceSetParamVaryType( runState->currentParams, tempPar, varyType );

    /* Add the prior variables */
    if ( !strcmp(tempPrior, "uniform") || !strcmp(tempPrior, "predefined") ){
      LALInferenceAddMinMaxPrior( runState->priorArgs, tempPar, &low, &high, LALINFERENCE_REAL8_t );
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

  LALInferenceIFOModel *ifotemp = ifo;
  while( ifotemp ){
    /* add in variables to say whether phase, sky position and binary parameter are varying */
    if( varyphase ) { LALInferenceAddVariable( ifotemp->params, "varyphase", &varyphase, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED ); }
    if( varyskypos ) { LALInferenceAddVariable( ifotemp->params, "varyskypos", &varyskypos, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED ); }
    if( varybinary ) { LALInferenceAddVariable( ifotemp->params, "varybinary", &varybinary, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED ); }

    ifotemp = ifotemp->next;
  }

  REAL8Vector *freqFactors = *(REAL8Vector **)LALInferenceGetVariable( ifo->params, "freqfactors" );

  if ( nonGR && freqFactors->length != 1 && freqFactors->data[0] != 2. ){
    fprintf(stderr, "Error... currently can only run with non-GR parameters for l=m=2 harmonic!\n");
    exit(3);
  }

  /* now check for a parameter correlation coefficient matrix file */
  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--cor-file" );
  if( ppt ){
    CHAR *corFile = XLALStringDuplicate( ppt->value );
    UINT4Vector *dims = XLALCreateUINT4Vector( 2 );
    dims->data[0] = 1;
    dims->data[1] = 1;

    corMat = XLALCreateREAL8Array( dims );

    corParams = XLALReadTEMPOCorFile( corMat, corFile );

    /* if the correlation matrix is given then add it as the prior for values with Gaussian errors specified in the par file */
    add_correlation_matrix( runState->currentParams, runState->priorArgs, corMat, corParams );

    XLALDestroyUINT4Vector( dims );
  }

  /* check if using a previous nested sampling file as a prior */
  samples_prior( runState );

  return;
}


/**
 * \brief Initialise the MCMC proposal distribution for sampling new points
 *
 * There are various proposal distributions that can be used to sample new live points via an MCMC. A combination of
 * different ones can be used to help efficiency for awkward posterior distributions. Here the proposals that can be
 * used are:
 * \c diffev Drawing a new point by differential evolution of two randomly chosen live points. All parameters are
 * evolved during a single draw.
 * \c freqBinJump Jumps that are the size of the Fourier frequency bins (can be used if searching over frequency).
 * \c ensembleStretch Ensemble stretch moves (WARNING: These can lead to long autocorrelation lengths).
 * \c ensembleWalk Ensemble walk moves. These are used as the default proposal.
 *
 * This function sets up the relative weights with which each of above distributions is used.
 *
 * \param runState [in] A pointer to the run state
 */
void initialise_proposal( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL;
  UINT4 defrac = 0, freqfrac = 0, esfrac = 0, ewfrac = 0;
  REAL8 temperature = 0.;
  const CHAR *defaultPropName = NULL;
  defaultPropName = XLALStringDuplicate( "none" );

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--diffev" );
  if( ppt ) { defrac = atoi( ppt->value ); }
  else { defrac = 0; } /* default value */

  ppt = LALInferenceGetProcParamVal( runState->commandLine, "--freqBinJump" );
  if( ppt ) { freqfrac = atoi( ppt->value ); }

  ppt = LALInferenceGetProcParamVal(runState->commandLine, "--ensembleStretch" );
  if ( ppt ) { esfrac = atoi( ppt->value ); }
  else { esfrac = 0; }

  ppt = LALInferenceGetProcParamVal(runState->commandLine, "--ensembleWalk" );
  if ( ppt ) { ewfrac = atoi( ppt->value ); }
  else { ewfrac = 1; }

  if( !defrac && !freqfrac && !ewfrac && !esfrac ){
    XLALPrintError("All proposal weights are zero!\n");
    XLAL_ERROR_VOID(XLAL_EFAILED);
  }

  runState->proposalStats = NULL;
  if(!runState->proposalStats) runState->proposalStats = XLALCalloc(1,sizeof(LALInferenceVariables));

  /* add proposals */
  if( defrac ){
    LALInferenceAddProposalToCycle( runState, differentialEvolutionFullName, &LALInferenceDifferentialEvolutionFull,
                                    defrac );
  }

  if ( freqfrac ){
    LALInferenceAddProposalToCycle( runState, frequencyBinJumpName, &LALInferenceFrequencyBinJump, freqfrac );
  }

  /* Use ensemble moves */
  if ( esfrac ){
    LALInferenceAddProposalToCycle( runState, ensembleStretchFullName, &LALInferenceEnsembleStretchFull, esfrac );
  }

  if ( ewfrac ){
    LALInferenceAddProposalToCycle( runState, ensembleWalkFullName, &LALInferenceEnsembleWalkFull, ewfrac );
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


/**
 * \brief Adds a correlation matrix for a multi-variate Gaussian prior
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
        if( !XLALStringCaseCompare(parMat->data[i], checkPrior->name) ){
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
        if( !XLALStringCaseCompare(parMat->data[i], checkPrior->name) ){
          /* replace it with the correlation matrix as a gsl_matrix */
          LALInferenceAddCorrelatedPrior( priors, checkPrior->name, &corMatg, &i );

          /* NOTE: the Gaussian prior will not be removed as the mean and standard deviation values are still
           * required when calculating the prior (see ppe_likelihood.c) */

          break;
        }
      }
    }
  }
}


/**
 * \brief Calculates the sum of the square of the data and model terms
 *
 * This function calculates the sum of the square of the data and model terms:
 * \f[
 * \sum_i^N \Re{d_i}^2 + \Im{d_i}^2, \sum_i^N \Re{h_i}^2, \sum_i^N \Im{h_i}^2, \sum_i^N \Re{d_i}\Re{h_i}, \sum_i^N Im{d_i}\Im{h_i}
 * \f]
 * for each stationary segment given in the \c chunkLength vector. These value are used in the likelihood calculation in
 * \c pulsar_log_likelihood and are precomputed here to speed that calculation up.
 *
 * \param runState [in] The analysis information structure
 */
void sum_data( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;
  LALInferenceIFOModel *ifomodel = runState->model->ifo;

  UINT4 gaussianLike = 0, roq = 0, nonGR = 0;

  if ( LALInferenceGetProcParamVal( runState->commandLine, "--gaussian-like" ) ){ gaussianLike = 1; }
  if ( LALInferenceGetProcParamVal( runState->commandLine, "--roq" ) ){ roq = 1; }
  if ( LALInferenceGetProcParamVal( runState->commandLine, "--nonGR" ) ){ nonGR = 1; }

  while( data ){
    REAL8Vector *sumdat = NULL;

    /* sums of the antenna pattern functions with themeselves and the data.
     * These won't be needed if searching over phase parameters, but there's
     * no harm in computing them anyway. */

    /* also have "explicitly" whitened (i.e. divided by variance) versions of these for use
     * by the function to calculate the signal-to-noise ratio (if using the Gaussian
     * likelihood the standard versions of these vectors will also be whitened too) */
    REAL8Vector *sumP = NULL; /* sum of tensor antenna pattern function a(t)^2 */
    REAL8Vector *sumC = NULL; /* sum of tensor antenna pattern function b(t)^2 */
    REAL8Vector *sumPWhite = NULL; /* sum of antenna pattern function a(t)^2 */
    REAL8Vector *sumCWhite = NULL; /* sum of antenna pattern function b(t)^2 */

    /* non-GR values */
    REAL8Vector *sumX = NULL; /* sum of vector antenna pattern function a(t)^2 */
    REAL8Vector *sumY = NULL; /* sum of vector antenna pattern function b(t)^2 */
    REAL8Vector *sumXWhite = NULL; /* whitened version */
    REAL8Vector *sumYWhite = NULL; /* whitened version */

    REAL8Vector *sumB = NULL; /* sum of scalar antenna pattern function a(t)^2 */
    REAL8Vector *sumL = NULL; /* sum of scalar antenna pattern function b(t)^2 */
    REAL8Vector *sumBWhite = NULL; /* whitened version */
    REAL8Vector *sumLWhite = NULL; /* whitened version */

    COMPLEX16Vector *sumDataP = NULL; /* sum of the data * a(t) */
    COMPLEX16Vector *sumDataC = NULL; /* sum of the data * b(t) */
    COMPLEX16Vector *sumDataX = NULL; /* sum of the data * a(t) */
    COMPLEX16Vector *sumDataY = NULL; /* sum of the data * b(t) */
    COMPLEX16Vector *sumDataB = NULL; /* sum of the data * a(t) */
    COMPLEX16Vector *sumDataL = NULL; /* sum of the data * b(t) */

    /* cross terms */
    REAL8Vector *sumPC = NULL, *sumPX = NULL, *sumPY = NULL, *sumPB = NULL, *sumPL = NULL;
    REAL8Vector *sumCX = NULL, *sumCY = NULL, *sumCB = NULL, *sumCL = NULL;
    REAL8Vector *sumXY = NULL, *sumXB = NULL, *sumXL = NULL;
    REAL8Vector *sumYB = NULL, *sumYL = NULL;
    REAL8Vector *sumBL = NULL;

    /* whitened versions cross terms */
    REAL8Vector *sumPCWhite = NULL, *sumPXWhite = NULL, *sumPYWhite = NULL, *sumPBWhite = NULL, *sumPLWhite = NULL;
    REAL8Vector *sumCXWhite = NULL, *sumCYWhite = NULL, *sumCBWhite = NULL, *sumCLWhite = NULL;
    REAL8Vector *sumXYWhite = NULL, *sumXBWhite = NULL, *sumXLWhite = NULL;
    REAL8Vector *sumYBWhite = NULL, *sumYLWhite = NULL;
    REAL8Vector *sumBLWhite = NULL;

    /* get antenna patterns */
    REAL8Vector *arespT = *(REAL8Vector **)LALInferenceGetVariable( ifomodel->params, "a_response_tensor" );
    REAL8Vector *brespT = *(REAL8Vector **)LALInferenceGetVariable( ifomodel->params, "b_response_tensor" );
    REAL8Vector *arespV = NULL, *brespV = NULL, *arespS = NULL, *brespS = NULL;

    if ( nonGR ){
      arespV = *(REAL8Vector **)LALInferenceGetVariable( ifomodel->params, "a_response_vector" );
      brespV = *(REAL8Vector **)LALInferenceGetVariable( ifomodel->params, "b_response_vector" );
      arespS = *(REAL8Vector **)LALInferenceGetVariable( ifomodel->params, "a_response_scalar" );
      brespS = *(REAL8Vector **)LALInferenceGetVariable( ifomodel->params, "b_response_scalar" );
    }

    INT4 tsteps = *(INT4 *)LALInferenceGetVariable( ifomodel->params, "timeSteps" );

    INT4 chunkLength = 0, length = 0, i = 0, j = 0, count = 0;
    COMPLEX16 B;
    REAL8 aT = 0., bT = 0., aV = 0., bV = 0., aS = 0., bS = 0.;

    UINT4Vector *chunkLengths;

    REAL8Vector *sidDayFrac = *(REAL8Vector**)LALInferenceGetVariable( ifomodel->params, "siderealDay" );

    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( ifomodel->params, "chunkLength" );

    length = ifomodel->times->length + 1 - chunkLengths->data[chunkLengths->length - 1];

    sumdat = XLALCreateREAL8Vector( chunkLengths->length );

    if ( !roq ){
      /* allocate memory */
      sumP = XLALCreateREAL8Vector( chunkLengths->length );
      sumC = XLALCreateREAL8Vector( chunkLengths->length );
      sumPC = XLALCreateREAL8Vector( chunkLengths->length );
      sumDataP = XLALCreateCOMPLEX16Vector( chunkLengths->length );
      sumDataC = XLALCreateCOMPLEX16Vector( chunkLengths->length );

      sumPWhite = XLALCreateREAL8Vector( chunkLengths->length );
      sumCWhite = XLALCreateREAL8Vector( chunkLengths->length );
      sumPCWhite = XLALCreateREAL8Vector( chunkLengths->length );

      if ( nonGR ){
        sumX = XLALCreateREAL8Vector( chunkLengths->length );
        sumY = XLALCreateREAL8Vector( chunkLengths->length );
        sumB = XLALCreateREAL8Vector( chunkLengths->length );
        sumL = XLALCreateREAL8Vector( chunkLengths->length );

        sumPX = XLALCreateREAL8Vector( chunkLengths->length );
        sumPY = XLALCreateREAL8Vector( chunkLengths->length );
        sumPB = XLALCreateREAL8Vector( chunkLengths->length );
        sumPL = XLALCreateREAL8Vector( chunkLengths->length );
        sumCX = XLALCreateREAL8Vector( chunkLengths->length );
        sumCY = XLALCreateREAL8Vector( chunkLengths->length );
        sumCB = XLALCreateREAL8Vector( chunkLengths->length );
        sumCL = XLALCreateREAL8Vector( chunkLengths->length );
        sumXY = XLALCreateREAL8Vector( chunkLengths->length );
        sumXB = XLALCreateREAL8Vector( chunkLengths->length );
        sumXL = XLALCreateREAL8Vector( chunkLengths->length );
        sumYB = XLALCreateREAL8Vector( chunkLengths->length );
        sumYL = XLALCreateREAL8Vector( chunkLengths->length );
        sumBL = XLALCreateREAL8Vector( chunkLengths->length );

        sumXWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumYWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumBWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumLWhite = XLALCreateREAL8Vector( chunkLengths->length );

        sumPXWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumPYWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumPBWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumPLWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumCXWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumCYWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumCBWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumCLWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumXYWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumXBWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumXLWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumYBWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumYLWhite = XLALCreateREAL8Vector( chunkLengths->length );
        sumBLWhite = XLALCreateREAL8Vector( chunkLengths->length );

        sumDataX = XLALCreateCOMPLEX16Vector( chunkLengths->length );
        sumDataY = XLALCreateCOMPLEX16Vector( chunkLengths->length );
        sumDataB = XLALCreateCOMPLEX16Vector( chunkLengths->length );
        sumDataL = XLALCreateCOMPLEX16Vector( chunkLengths->length );
      }
    }

    REAL8 tsv = LAL_DAYSID_SI / tsteps, T = 0., timeMin = 0., timeMax = 0.;
    REAL8 logGaussianNorm = 0.; /* normalisation constant for Gaussian distribution */

    for( i = 0, count = 0 ; i < length ; i+= chunkLength, count++ ){
      chunkLength = chunkLengths->data[count];

      sumdat->data[count] = 0.;

      if ( !roq ){
        sumP->data[count] = 0.;
        sumC->data[count] = 0.;
        sumPC->data[count] = 0.;
        sumDataP->data[count] = 0.;
        sumDataC->data[count] = 0.;

        sumPWhite->data[count] = 0.;
        sumCWhite->data[count] = 0.;
        sumPCWhite->data[count] = 0.;

        if ( nonGR ){
          sumX->data[count] = 0.;
          sumY->data[count] = 0.;
          sumB->data[count] = 0.;
          sumL->data[count] = 0.;

          sumPX->data[count] = 0.;
          sumPY->data[count] = 0.;
          sumPB->data[count] = 0.;
          sumPL->data[count] = 0.;
          sumCX->data[count] = 0.;
          sumCY->data[count] = 0.;
          sumCB->data[count] = 0.;
          sumCL->data[count] = 0.;
          sumXY->data[count] = 0.;
          sumXB->data[count] = 0.;
          sumXL->data[count] = 0.;
          sumYB->data[count] = 0.;
          sumYL->data[count] = 0.;
          sumBL->data[count] = 0.;

          sumXWhite->data[count] = 0.;
          sumYWhite->data[count] = 0.;
          sumBWhite->data[count] = 0.;
          sumLWhite->data[count] = 0.;

          sumPXWhite->data[count] = 0.;
          sumPYWhite->data[count] = 0.;
          sumPBWhite->data[count] = 0.;
          sumPLWhite->data[count] = 0.;
          sumCXWhite->data[count] = 0.;
          sumCYWhite->data[count] = 0.;
          sumCBWhite->data[count] = 0.;
          sumCLWhite->data[count] = 0.;
          sumXYWhite->data[count] = 0.;
          sumXBWhite->data[count] = 0.;
          sumXLWhite->data[count] = 0.;
          sumYBWhite->data[count] = 0.;
          sumYLWhite->data[count] = 0.;
          sumBLWhite->data[count] = 0.;

          sumDataX->data[count] = 0.;
          sumDataY->data[count] = 0.;
          sumDataB->data[count] = 0.;
          sumDataL->data[count] = 0.;
        }
      }

      for( j = i ; j < i + chunkLength ; j++){
        REAL8 vari = 1., a0 = 0., a1 = 0., b0 = 0., b1 = 0., timeScaled = 0.;
        INT4 timebinMin = 0, timebinMax = 0;

        B = data->compTimeData->data->data[j];

        /* if using a Gaussian likelihood divide all these values by the variance */
        if ( gaussianLike ) {
          vari = data->varTimeData->data->data[j];
          logGaussianNorm -= log(LAL_TWOPI*vari);
        }

        /* sum up the data */
        sumdat->data[count] += (creal(B)*creal(B) + cimag(B)*cimag(B))/vari;

        if ( !roq ){
          /* set the time bin for the lookup table and interpolate between bins */
          T = sidDayFrac->data[j];
          timebinMin = (INT4)fmod( floor(T / tsv), tsteps );
          timeMin = timebinMin*tsv;
          timebinMax = (INT4)fmod( timebinMin + 1, tsteps );
          timeMax = timeMin + tsv;

          /* get values of vector for linear interpolation */
          a0 = arespT->data[timebinMin];
          a1 = arespT->data[timebinMax];
          b0 = brespT->data[timebinMin];
          b1 = brespT->data[timebinMax];

          /* rescale time for linear interpolation on a unit square */
          timeScaled = (T - timeMin)/(timeMax - timeMin);

          aT = a0 + (a1-a0)*timeScaled;
          bT = b0 + (b1-b0)*timeScaled;

          /* sum up the other terms */
          sumP->data[count] += aT*aT/vari;
          sumC->data[count] += bT*bT/vari;
          sumPC->data[count] += aT*bT/vari;
          sumDataP->data[count] += B*aT/vari;
          sumDataC->data[count] += B*bT/vari;

          /* non-GR values */
          if ( nonGR ){
            a0 = arespV->data[timebinMin];
            a1 = arespV->data[timebinMax];
            b0 = brespV->data[timebinMin];
            b1 = brespV->data[timebinMax];

            aV = a0 + (a1-a0)*timeScaled;
            bV = b0 + (b1-b0)*timeScaled;

            a0 = arespS->data[timebinMin];
            a1 = arespS->data[timebinMax];
            b0 = brespS->data[timebinMin];
            b1 = brespS->data[timebinMax];

            aS = a0 + (a1-a0)*timeScaled;
            bS = b0 + (b1-b0)*timeScaled;

            sumX->data[count] += aV*aV/vari;
            sumY->data[count] += bV*bV/vari;
            sumB->data[count] += aS*aS/vari;
            sumL->data[count] += bS*bS/vari;

            sumPX->data[count] += aT*aV/vari;
            sumPY->data[count] += aT*bV/vari;
            sumPB->data[count] += aT*aS/vari;
            sumPL->data[count] += aT*bS/vari;
            sumCX->data[count] += bT*aV/vari;
            sumCY->data[count] += bT*bV/vari;
            sumCB->data[count] += bT*aS/vari;
            sumCL->data[count] += bT*bS/vari;
            sumXY->data[count] += aV*bV/vari;
            sumXB->data[count] += aV*aS/vari;
            sumXL->data[count] += aV*bS/vari;
            sumYB->data[count] += bV*aS/vari;
            sumYL->data[count] += bV*bS/vari;
            sumBL->data[count] += aS*bS/vari;

            sumDataX->data[count] += B*aV/vari;
            sumDataY->data[count] += B*bV/vari;
            sumDataB->data[count] += B*aS/vari;
            sumDataL->data[count] += B*bS/vari;
          }

          /* get "explicitly whitened" versions, i.e. for use in signal-to-noise ratio
           * calculations even when not using a Gaussian likelihood */
          vari = data->varTimeData->data->data[j];
          sumPWhite->data[count] += aT*aT/vari;
          sumCWhite->data[count] += bT*bT/vari;
          sumPCWhite->data[count] += aT*bT/vari;

          /* non-GR values */
          if ( nonGR ){
            sumXWhite->data[count] += aV*aV/vari;
            sumYWhite->data[count] += bV*bV/vari;
            sumBWhite->data[count] += aS*aS/vari;
            sumLWhite->data[count] += bS*bS/vari;

            sumPXWhite->data[count] += aT*aV/vari;
            sumPYWhite->data[count] += aT*bV/vari;
            sumPBWhite->data[count] += aT*aS/vari;
            sumPLWhite->data[count] += aT*bS/vari;
            sumCXWhite->data[count] += bT*aV/vari;
            sumCYWhite->data[count] += bT*bV/vari;
            sumCBWhite->data[count] += bT*aS/vari;
            sumCLWhite->data[count] += bT*bS/vari;
            sumXYWhite->data[count] += aV*bV/vari;
            sumXBWhite->data[count] += aV*aS/vari;
            sumXLWhite->data[count] += aV*bS/vari;
            sumYBWhite->data[count] += bV*aS/vari;
            sumYLWhite->data[count] += bV*bS/vari;
            sumBLWhite->data[count] += aS*bS/vari;
          }
        }
      }
    }

    /* add all the summed data values - remove if already there, so that sum_data can be called more
     * than once if required e.g. if needed in the injection functions */

    check_and_add_fixed_variable( ifomodel->params, "sumData", &sumdat, LALINFERENCE_REAL8Vector_t );

    if ( !roq ){
      check_and_add_fixed_variable( ifomodel->params, "sumP", &sumP, LALINFERENCE_REAL8Vector_t );
      check_and_add_fixed_variable( ifomodel->params, "sumC", &sumC, LALINFERENCE_REAL8Vector_t );
      check_and_add_fixed_variable( ifomodel->params, "sumPC", &sumPC, LALINFERENCE_REAL8Vector_t );
      check_and_add_fixed_variable( ifomodel->params, "sumDataP", &sumDataP, LALINFERENCE_COMPLEX16Vector_t );
      check_and_add_fixed_variable( ifomodel->params, "sumDataC", &sumDataC, LALINFERENCE_COMPLEX16Vector_t );
      check_and_add_fixed_variable( ifomodel->params, "sumPWhite", &sumPWhite, LALINFERENCE_REAL8Vector_t );
      check_and_add_fixed_variable( ifomodel->params, "sumCWhite", &sumCWhite, LALINFERENCE_REAL8Vector_t );
      check_and_add_fixed_variable( ifomodel->params, "sumPCWhite", &sumPCWhite, LALINFERENCE_REAL8Vector_t );

      if ( nonGR ){
        check_and_add_fixed_variable( ifomodel->params, "sumX", &sumX, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumY", &sumY, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumB", &sumB, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumL", &sumL, LALINFERENCE_REAL8Vector_t );

        check_and_add_fixed_variable( ifomodel->params, "sumPX", &sumPX, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumPY", &sumPY, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumPB", &sumPB, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumPL", &sumPL, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCX", &sumCX, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCY", &sumCY, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCB", &sumCB, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCL", &sumCL, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumXY", &sumXY, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumXB", &sumXB, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumXL", &sumXL, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumYB", &sumYB, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumYL", &sumYL, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumBL", &sumBL, LALINFERENCE_REAL8Vector_t );

        check_and_add_fixed_variable( ifomodel->params, "sumDataX", &sumDataX, LALINFERENCE_COMPLEX16Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumDataY", &sumDataY, LALINFERENCE_COMPLEX16Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumDataB", &sumDataB, LALINFERENCE_COMPLEX16Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumDataL", &sumDataL, LALINFERENCE_COMPLEX16Vector_t );

        check_and_add_fixed_variable( ifomodel->params, "sumXWhite", &sumXWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumYWhite", &sumYWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumBWhite", &sumBWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumLWhite", &sumLWhite, LALINFERENCE_REAL8Vector_t );

        check_and_add_fixed_variable( ifomodel->params, "sumPXWhite", &sumPXWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumPYWhite", &sumPYWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumPBWhite", &sumPBWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumPLWhite", &sumPLWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCXWhite", &sumCXWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCYWhite", &sumCYWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCBWhite", &sumCBWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumCLWhite", &sumCLWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumXYWhite", &sumXYWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumXBWhite", &sumXBWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumXLWhite", &sumXLWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumYBWhite", &sumYBWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumYLWhite", &sumYLWhite, LALINFERENCE_REAL8Vector_t );
        check_and_add_fixed_variable( ifomodel->params, "sumBLWhite", &sumBLWhite, LALINFERENCE_REAL8Vector_t );
      }
    }
    else{ /* add parameter defining the usage of RQO here (as this is after any injection generation, which
           * would fail if this was set */
      LALInferenceAddVariable( ifomodel->params, "roq", &roq, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
    }

    LALInferenceAddVariable( ifomodel->params, "logGaussianNorm", &logGaussianNorm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

    data = data->next;
    ifomodel = ifomodel->next;
  }

  return;
}


/**
 * \brief Print any non-fixed sample to file (based on \c LALInferencePrintSampleNonFixed)
 */
static void PrintNonFixedSample(FILE *fp, LALInferenceVariables *sample){
  UINT4 i;
  UINT4Vector *v=NULL;

  if(sample==NULL) { return; }

  LALInferenceVariableItem *ptr=sample->head;
  if(fp==NULL) { return; }

  while(ptr!=NULL) {
    if (LALInferenceGetVariableVaryType(sample, ptr->name) != LALINFERENCE_PARAM_FIXED && ptr->type != LALINFERENCE_gslMatrix_t ) {
      switch (ptr->type) {
        case LALINFERENCE_INT4_t:
          fprintf(fp, "%"LAL_INT4_FORMAT, *(INT4 *) ptr->value);
          break;
        case LALINFERENCE_INT8_t:
          fprintf(fp, "%"LAL_INT8_FORMAT, *(INT8 *) ptr->value);
          break;
        case LALINFERENCE_UINT4_t:
          fprintf(fp, "%"LAL_UINT4_FORMAT, *(UINT4 *) ptr->value);
          break;
        case LALINFERENCE_REAL4_t:
          fprintf(fp, "%9.20e", *(REAL4 *) ptr->value);
          break;
        case LALINFERENCE_REAL8_t:
          fprintf(fp, "%9.20le", *(REAL8 *) ptr->value);
          break;
        case LALINFERENCE_COMPLEX8_t:
          fprintf(fp, "%e + i*%e", (REAL4) crealf(*(COMPLEX8 *) ptr->value), (REAL4) cimagf(*(COMPLEX8 *) ptr->value));
          break;
        case LALINFERENCE_COMPLEX16_t:
          fprintf(fp, "%e + i*%e", (REAL8) creal(*(COMPLEX16 *) ptr->value), (REAL8) cimag(*(COMPLEX16 *) ptr->value));
          break;
        case LALINFERENCE_UINT4Vector_t:
          v = *((UINT4Vector **)ptr->value);
          for(i=0;i<v->length;i++){
            fprintf(fp,"%11.7f",(REAL8)v->data[i]);
            if( i!=(UINT4)(v->length-1) ) { fprintf(fp,"\t"); }
          }
          break;
        default:
          fprintf(stdout, "<can't print>");
      }
      fprintf(fp,"\t");
    }
    ptr=ptr->next;
  }
  return;
}


/**
 * \brief Print out only the variable (i.e. non-fixed) parameters to the file
 *
 * If the command line argument \c --output-all-params is given then this function
 * will output all parameters that have been stored, but by default it will only
 * output variable parameters to the nested samples to a file.
 */
void LogSampleToFile(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  FILE *outfile = NULL;
  if( LALInferenceCheckVariable(state->algorithmParams,"outfile") ){
    outfile = *(FILE **)LALInferenceGetVariable( state->algorithmParams, "outfile" );
  }
  /* Write out old sample */
  if( outfile == NULL ) { return; }
  LALInferenceSortVariablesByName(vars);

  /* only write out non-fixed samples if required */
  if ( LALInferenceGetProcParamVal( state->commandLine, "--output-all-params" ) ){
    LALInferencePrintSample(outfile, vars);
  }
  else{
    /* default to writing out only the non-fixed (i.e. variable) parameters */
    PrintNonFixedSample(outfile, vars);
  }
  fprintf(outfile,"\n");

  return;
}


/**
 * \brief Print out only the variable (i.e. non-fixed) parameters to the file whilst also creating
 * an array for all parameters
 *
 * This function (which is used if the XML library is present) will only output variable
 * parameters.
 */
void LogSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  LALInferenceVariables **output_array = NULL;
  UINT4 N_output_array = 0;
  LALInferenceSortVariablesByName(vars);

  LogSampleToFile(state, vars);

  /* Set up the array if it is not already allocated */
  if(LALInferenceCheckVariable(state->algorithmParams,"outputarray")){
    output_array = *(LALInferenceVariables ***)LALInferenceGetVariable(state->algorithmParams,"outputarray");
  }
  else{
    LALInferenceAddVariable(state->algorithmParams,"outputarray",&output_array,LALINFERENCE_void_ptr_t,LALINFERENCE_PARAM_OUTPUT);
  }

  if(LALInferenceCheckVariable(state->algorithmParams,"N_outputarray")){
    N_output_array = *(INT4 *)LALInferenceGetVariable(state->algorithmParams,"N_outputarray");
  }
  else{
    LALInferenceAddVariable(state->algorithmParams,"N_outputarray",&N_output_array,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
  }

  /* Expand the array for new sample */
  output_array = XLALRealloc(output_array, (N_output_array+1) *sizeof(LALInferenceVariables *));
  if( !output_array ){
    XLAL_ERROR_VOID(XLAL_EFAULT, "Unable to allocate array for samples.");
  }
  else{
    /* Save sample and update */
    output_array[N_output_array] = XLALCalloc(1,sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(vars, output_array[N_output_array]);
    N_output_array++;

    LALInferenceSetVariable(state->algorithmParams,"outputarray",&output_array);
    LALInferenceSetVariable(state->algorithmParams,"N_outputarray",&N_output_array);
  }

  return;
}

