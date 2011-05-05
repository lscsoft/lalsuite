/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

#include "pulsar_parameter_estimation_nested.h"


#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0),
(1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* maximum number of different detectors */
#define MAXDETS 6

#include <lal/Units.h>
#include <sys/time.h>
#include <lal/XLALError.h>

RCSID("$Id$");

/* global variable */
INT4 verbose=0;


/* Usage format string */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --detectors         all IFOs with data to be analysed e.g. H1,H2\n\
                     (delimited by commas)\n"\
" --pulsar            name of pulsar e.g. J0534+2200\n"\
" --par-file          pulsar parameter (.par) file (full path) \n"\
" --input-files       full paths and file names for the data for each\n\
                     detector in the list (must be in the same order)\n\
                     delimited by commas\n"\
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
" Nested sampling parameters:\n"\
" --Nlive             (INT4) no. of live points for nested sampling\n"\
" --Nmcmc             (INT4) no. of for MCMC used to find new live points\n"\
" --Nruns             (INT4) no. of parallel runs\n"\
" --tolerance         (REAL8) tolerance of nested sampling integrator\n"\
" --randomseed        seed for random number generator\n"\
"\n"


INT4 main( INT4 argc, CHAR *argv[] ){
  ProcessParamsTable *param_table;
  LALInferenceRunState runState;
  REAL8 logZnoise = 0.;
  
  REAL8Vector *logLikelihoods = NULL;
  
  /* Get ProcParamsTable from input arguments */
  param_table = parseCommandLine( argc, argv );
  runState.commandLine = param_table;
  
  /* Initialise data structures from the command line arguments */
  /* runState.data = readPulsarData(argc,argv); */
  
  /* Initialise the algorithm structures from the command line arguments */
  /* Include setting up random number generator etc */
  initialiseAlgorithm( &runState );
  printf("Done initialiseAlgorithm\n");
  
  /* read in data using command line arguments */
  readPulsarData( &runState );
  printf("Done readPulsarData\n");
  
  /* set algorithm to use Nested Sampling */
  runState.algorithm = &NestedSamplingAlgorithm;
  runState.evolve = &NestedSamplingOneStep;
  
  /* set likelihood function */
  runState.likelihood = &pulsar_log_likelihood;
  
  /* set prior function */
  runState.prior = &priorFunction;
  
  /* set signal model/template */
  runState.template = get_pulsar_model;
  
  /* Generate the lookup tables and read parameters from par file */
  setupFromParFile( &runState );
  printf("Done setupFromParFile\n");
  
  /* Initialise the prior distribution given the command line arguments */
  /* Initialise the proposal distribution given the command line arguments */
  initialiseProposal( &runState );
  printf("Done initialiseProposal\n");
 
  runState.proposal = LALInferenceProposalPulsarNS;
  
  /* get noise likelihood and add as variable to runState */
  logZnoise = noise_only_model( runState.data );
  addVariable( runState.algorithmParams, "logZnoise", &logZnoise, REAL8_t, 
               PARAM_FIXED );
  printf("Done setting noise likelihood\n");
  
  /* Create live points array and fill initial parameters */
  setupLivePointsArray( &runState );
  printf("Done setupLivePointsArray\n");
  
  logLikelihoods = *(REAL8Vector **)getVariable( runState.algorithmParams,
                                                 "logLikelihoods" );

  /* for(i=0; i< logLikelihoods->length; i++)
    fprintf(stderr, "logL = %le\n", logLikelihoods->data[i]); */
  
  /* Call the nested sampling algorithm */
  runState.algorithm( &runState );
  printf("Done nested sampling algorithm\n");

  /* re-read in output samples and rescale appropriately */
  rescaleOutput( &runState );
  
  return 0;
}

/******************************************************************************/
/*                      INITIALISATION FUNCTIONS                              */
/******************************************************************************/

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
  ppt = getProcParamVal( commandLine, "--help" );
  if(ppt){
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* Initialise parameters structure */
  runState->algorithmParams = XLALCalloc( 1, sizeof(LALVariables) );
  runState->priorArgs = XLALCalloc( 1, sizeof(LALVariables) );
  runState->proposalArgs = XLALCalloc( 1, sizeof(LALVariables) );
 
  ppt = getProcParamVal( commandLine, "--verbose" );
  if( ppt ) {
    verbose = 1;
    addVariable( runState->algorithmParams, "verbose", &verbose , INT4_t,
                 PARAM_FIXED);
  }

  /* Number of live points */
  tmpi = atoi( getProcParamVal(commandLine, "--Nlive")->value );
  addVariable( runState->algorithmParams,"Nlive", &tmpi, INT4_t, PARAM_FIXED );
        
  /* Number of points in MCMC chain */
  tmpi = atoi( getProcParamVal(commandLine, "--Nmcmc")->value );
  addVariable( runState->algorithmParams, "Nmcmc", &tmpi, INT4_t, PARAM_FIXED );

  /* Optionally specify number of parallel runs */
  ppt = getProcParamVal( commandLine, "--Nruns" );
  if(ppt) {
    tmpi = atoi( ppt->value );
    addVariable( runState->algorithmParams, "Nruns", &tmpi, INT4_t,
                 PARAM_FIXED );
  }
        
  /* Tolerance of the Nested sampling integrator */
  ppt = getProcParamVal( commandLine, "--tolerance" );
  if( ppt ){
    tmp = strtod( ppt->value, (char **)NULL );
    addVariable( runState->algorithmParams, "tolerance", &tmp, REAL8_t,
                 PARAM_FIXED );
  }
        
  /* Set up the random number generator */
  gsl_rng_env_setup();
  runState->GSLrandom = gsl_rng_alloc( gsl_rng_mt19937 );
        
  /* (try to) get random seed from command line: */
  ppt = getProcParamVal( commandLine, "--randomseed" );
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


void readPulsarData( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;
  
  CHAR *detectors = NULL;
  CHAR *outfile = NULL;
  CHAR *inputfile = NULL;
  
  CHAR *filestr = NULL;
  
  CHAR *efile = NULL;
  CHAR *sfile = NULL;
  
  CHAR *tempdets=NULL;
  CHAR *tempdet=NULL;
  
  CHAR dets[MAXDETS][256];
  INT4 numDets = 0, i = 0;
 
  LALIFOData *ifodata = NULL;
  LALIFOData *prev = NULL;
  
  runState->data = NULL;
  
  /* get the detectors */
  ppt = getProcParamVal( commandLine, "--detectors" );
  if( ppt ){
    detectors = 
      XLALStringDuplicate( getProcParamVal(commandLine,"--detectors")->value );
  }
  else{
    fprintf(stderr, "Error... --detectors needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  ppt = getProcParamVal( commandLine,"--input-files" );
  if( ppt ){
    inputfile = 
      XLALStringDuplicate(getProcParamVal(commandLine,"--input-files")->value);
  }
  
  if ( inputfile == NULL ){
    fprintf(stderr, "Error... either an input file needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* get the output directory */
  ppt = getProcParamVal( commandLine, "--outfile" );
  if( ppt ){
    outfile = 
      XLALStringDuplicate( getProcParamVal( commandLine, "--outfile" )->value );
  }
  else{
    fprintf(stderr, "Error... --outfile needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* get ephemeris files */
  ppt = getProcParamVal( commandLine, "--ephem-earth" );
  if( ppt ){
  efile = 
    XLALStringDuplicate( getProcParamVal(commandLine,"--ephem-earth")->value );
  }
  ppt = getProcParamVal( commandLine, "--ephem-sun" );
  if( ppt ){
    sfile = 
      XLALStringDuplicate( getProcParamVal(commandLine,"--ephem-sun")->value );
  }

  /* count the number of detectors from command line argument of comma separated
     vales and set their names */
  tempdets = XLALStringDuplicate( detectors );
    
  while(1){
    tempdet = strsep( &tempdets, "," );
    XLALStringCopy( dets[numDets], tempdet, strlen(tempdet)+1 );
    
    numDets++;
    
    if( tempdets == NULL ) break;
  }

  /* check ephemeris files exist and if not output an error message */
  if( access(sfile, F_OK) != 0 || access(efile, F_OK) != 0 ){
    fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
    exit(3);
  }
  
  /* count the number of input files (by counting commas) and check it's equal
     to the number of detectors */
  {
    CHAR *inputstr = NULL;
    CHAR *tempstr = NULL;
    INT4 count = 0;
    
    inputstr = XLALStringDuplicate( inputfile );

    /* count number of commas */
    while(1){
      tempstr = strsep(&inputstr, ",");
      
      if ( inputstr == NULL ) break;
      
      count++;
    }
    
    if ( count+1 != numDets ){
      fprintf(stderr, "Error... Number of input files given is not equal to the\
 number of detectors given!\n");
      exit(3);
    }
  }
  
  /* reset filestr */
  filestr = XLALStringDuplicate( inputfile );
  
  /* read in data */
  for( i = 0, prev=NULL ; i < numDets ; i++, prev=ifodata ){
    CHAR *datafile = NULL;
    REAL8 times = 0;
    LIGOTimeGPS gpstime;
    COMPLEX16 dataVals;
    REAL8Vector *temptimes = NULL;
    INT4 j = 0, k = 0, count = 0;
    
    FILE *fp = NULL;
    
    ifodata = XLALCalloc( 1, sizeof(LALIFOData) );
    ifodata->modelParams = XLALCalloc( 1, sizeof(LALVariables) );
    ifodata->modelDomain = timeDomain;
    ifodata->next = NULL;
    ifodata->dataParams = XLALCalloc( 1, sizeof(LALVariables) );
    
    if( i == 0 ) runState->data = ifodata;
    if( i > 0 ) prev->next = ifodata;

    /* set detector */
    ifodata->detector = XLALGetSiteInfo( dets[i] );

    /*============================ GET DATA ==================================*/
    /* get i'th filename from the comma separated list */
    while ( count <= i ){
      datafile = strsep(&filestr, ",");
      
      count++;
    }
   
    /* open data file */
    if( (fp = fopen(datafile, "r")) == NULL ){
      fprintf(stderr, "Error... can't open data file %s!\n", datafile);
      exit(0);
    }

    j=0;

    /* read in data */
    temptimes = XLALCreateREAL8Vector( MAXLENGTH );

    /* read in data */
    while(fscanf(fp, "%lf%lf%lf", &times, &dataVals.re, &dataVals.im) != EOF){
      /* check that size of data file is not to large */
      if ( j == 0 ){
        XLALGPSSetREAL8( &gpstime, times );

        ifodata->compTimeData = NULL;
        ifodata->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                               1.,
                                                               &lalSecondUnit,
                                                               MAXLENGTH );
      }
      
      if( j == MAXLENGTH ){
        fprintf(stderr, "Error... size of MAXLENGTH not large enough.\n");
        exit(3);
      }

      temptimes->data[j] = times;
      ifodata->compTimeData->data->data[j].re = dataVals.re;
      ifodata->compTimeData->data->data[j].im = dataVals.im;
      
      j++;
    }
    
    fclose(fp);
    
    /* resize the data */
    ifodata->compTimeData =
      XLALResizeCOMPLEX16TimeSeries( ifodata->compTimeData, 0, j );

    ifodata->compModelData = NULL;
    ifodata->compModelData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                            1., &lalSecondUnit,
                                                            j );
    
    /* fill in time stamps as LIGO Time GPS Vector */
    ifodata->dataTimes = NULL;
    ifodata->dataTimes = XLALCreateTimestampVector( j );
          
    for ( k = 0; k<j; k++ )
      XLALGPSSetREAL8( &ifodata->dataTimes->data[k], temptimes->data[k] );
      
    XLALDestroyREAL8Vector( temptimes );

    /* set ephemeris data */
    ifodata->ephem = XLALMalloc( sizeof(EphemerisData) );
    
    /* set up ephemeris information */
    ifodata->ephem = XLALInitBarycenter( efile, sfile );
  }
}


void setupFromParFile( LALInferenceRunState *runState )
/* Read the PAR file of pulsar parameters and setup the code using them */
/* Generates lookup tables also */
{
  LALSource psr;
  BinaryPulsarParams pulsar;
  REAL8Vector *phase_vector;
  LALIFOData *data = runState->data;
  LALVariables *scaletemp;
  ProcessParamsTable *ppt = NULL;
  UINT4 withphase=0;
  
  ppt = getProcParamVal( runState->commandLine, "--par-file" );
  if( ppt == NULL ) { fprintf(stderr,"Must specify --par-file!\n"); exit(1); }
  char *parFile = ppt->value;
        
  /* get the pulsar parameters */
  XLALReadTEMPOParFile( &pulsar, parFile );
  psr.equatorialCoords.longitude = pulsar.ra;
  psr.equatorialCoords.latitude = pulsar.dec;
  psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* Setup lookup tables for amplitudes */
  setupLookupTables( runState, &psr );
  printf("Setup lookup tables.\n");
  
  runState->currentParams = XLALCalloc( 1, sizeof(LALVariables) );
  
  scaletemp = XLALCalloc( 1, sizeof(LALVariables) );
  
  /* if no binary model set the value to "None" */
  if ( pulsar.model == NULL )
    pulsar.model = XLALStringDuplicate("None");  
  
  /* Add initial (unchanging) variables for the model, initial (unity) scale
     factors, and any Gaussian priors defined from the par file */
  withphase = add_initial_variables( runState->currentParams, scaletemp,
                                     runState->priorArgs, pulsar );
  
  /* Setup initial phase */
  while( data ){
    UINT4 i = 0;
    LALVariableItem *scaleitem = scaletemp->head;
    
    phase_vector = get_phase_model( pulsar, data );
    
    data->timeData = NULL;
    data->timeData = XLALCreateREAL8TimeSeries( "", &data->dataTimes->data[0], 
                                                0., 1., &lalSecondUnit,
                                                phase_vector->length );
    for ( i=0; i<phase_vector->length; i++ )
      data->timeData->data->data[i] = phase_vector->data[i];
  
    /* add the withphase variable to the data->dataParams variable, which
       defines whether the phase model needs to be calculated in the
       likelihood, or not (i.e. are any phase parameters required to be
       searched over) */
    addVariable( data->dataParams, "withphase", &withphase, UINT4_t,
                 PARAM_FIXED);
      
    /* add the scale factors from scaletemp into the data->dataParams
       structure */
    for( ; scaleitem; scaleitem = scaleitem->next ){
      addVariable( data->dataParams, scaleitem->name, scaleitem->value,
                   scaleitem->type, scaleitem->vary );
    }
    
    data = data->next;
  }
                                     
  return;
}


void setupLookupTables( LALInferenceRunState *runState, LALSource *source ){    
  /* Set up lookup tables */
  /* Using psi bins, time bins */
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  REAL8Vector *exclamation = NULL;
  
  INT4 chunkMin, chunkMax, count = 0;

  /* Get chunk min and chunk max */
  ppt = getProcParamVal( commandLine, "--chunk-min" );
  if( ppt ) chunkMin = atoi( ppt->value );
  else chunkMin = CHUNKMIN; /* default minimum chunk length */
    
  ppt = getProcParamVal( commandLine, "--chunk-max" );
  if( ppt ) chunkMax = atoi( ppt->value );
  else chunkMax = CHUNKMAX; /* default maximum chunk length */
   
  LALIFOData *data = runState->data;

  gsl_matrix *LUfplus = NULL;
  gsl_matrix *LUfcross = NULL;

  REAL8 t0;
  LALDetAndSource detAndSource;
        
  ppt = getProcParamVal( commandLine, "--psi-bins" );
  INT4 psiBins;
  if( ppt ) psiBins = atoi( ppt->value );
  else psiBins = PSIBINS; /* default psi bins */
                
  ppt = getProcParamVal( commandLine, "--time-bins" );
  INT4 timeBins;
  if( ppt ) timeBins = atoi( ppt->value );
  else timeBins = TIMEBINS; /* default time bins */  
        
  if(verbose) fprintf(stdout,"psi-bins = %i, time-bins =%i\n",psiBins,timeBins);
        
  while(data){
    REAL8Vector *sumData = NULL;
    UINT4Vector *chunkLength = NULL;
    INT4 i=0;
    
    addVariable( data->dataParams, "psiSteps", &psiBins, INT4_t, PARAM_FIXED ); 
    addVariable( data->dataParams, "timeSteps", &timeBins, INT4_t, 
                 PARAM_FIXED );
    
    t0 = XLALGPSGetREAL8( &data->dataTimes->data[0] );
    detAndSource.pDetector = data->detector;
    detAndSource.pSource = source;
                
    LUfplus = gsl_matrix_alloc( psiBins, timeBins );

    LUfcross = gsl_matrix_alloc( psiBins, timeBins );
    response_lookup_table( t0, detAndSource, timeBins, psiBins, LUfplus, 
                           LUfcross );
    addVariable( data->dataParams, "LU_Fplus", &LUfplus, gslMatrix_t,
                 PARAM_FIXED );
    addVariable( data->dataParams, "LU_Fcross", &LUfcross, gslMatrix_t,
                 PARAM_FIXED );

    addVariable( data->dataParams, "chunkMin", &chunkMin, INT4_t, PARAM_FIXED );
    addVariable( data->dataParams, "chunkMax", &chunkMax, INT4_t, PARAM_FIXED );
    
    /* get chunk lengths of data */
    
    /* chunkLength = get_chunk_lengths( data, chunkMax ); */
    chunkLength = chop_n_merge( data, chunkMin );
    
    addVariable( data->dataParams, "chunkLength", &chunkLength, UINT4Vector_t,
                 PARAM_FIXED );

    if ( count == 0 ){
      UINT4 k = 0;
      UINT4 maxcl = 0;
      
      /* get max chunklength */
      for ( k=0; k < chunkLength->length; k++ ){
        if ( chunkLength->data[k] > maxcl ) 
          maxcl = chunkLength->data[k] > maxcl;
      }
      
      exclamation = XLALCreateREAL8Vector( maxcl + 1 );

      /* to save time get all log factorials up to chunkMax */
      for( i = 0 ; i < (INT4)maxcl+1 ; i++ )
        exclamation->data[i] = log_factorial(i);
    
      count++;
    }
    
    addVariable( data->dataParams, "logFactorial", &exclamation, REAL8Vector_t,
                 PARAM_FIXED );
                 
    /* get sum of data for each chunk */
    sumData = sum_data( data );
    addVariable( data->dataParams, "sumData", &sumData, REAL8Vector_t,
                 PARAM_FIXED);

    data = data->next;
  }
        
  return;
}


UINT4 add_initial_variables( LALVariables *ini, LALVariables *scaleFac,
                             LALVariables *priorArgs, 
                             BinaryPulsarParams pars ){ 
  /* include a scale factor of 1 scaling values and if the parameter file
     contains an uncertainty then set the prior to be Gaussian with the
     uncertainty as the standard deviation */
  
  UINT4 withphase=0, dummy;
  
  /* amplitude model parameters */
  dummy = add_variable_scale_prior( ini, scaleFac, priorArgs, "h0", pars.h0,
                                    pars.h0Err );
  dummy = add_variable_scale_prior( ini, scaleFac, priorArgs, "phi0", pars.phi0,
                                    pars.phi0Err );
  dummy = add_variable_scale_prior( ini, scaleFac, priorArgs, "cosiota",
                                    pars.cosiota, pars.cosiotaErr );
  dummy = add_variable_scale_prior( ini, scaleFac, priorArgs, "psi", pars.psi,
                                    pars.psiErr );
    
  /* phase model parameters */
  
  /* frequency */
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "f0",
                                        pars.f0, pars.f0Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "f1",
                                        pars.f1, pars.f1Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "f2",
                                        pars.f2, pars.f2Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "f3",
                                        pars.f3, pars.f3Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "f4",
                                        pars.f4, pars.f4Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "f5",
                                        pars.f5, pars.f5Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "pepoch",
                                        pars.pepoch, pars.pepochErr );
  
  /* sky position */
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "ra",
                                        pars.ra, pars.raErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "pmra",
                                        pars.pmra, pars.pmraErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "dec",
                                        pars.dec, pars.decErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "pmdec",
                                        pars.pmdec, pars.pmdecErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "posepoch",
                                        pars.posepoch, pars.posepochErr );
  
  /* binary system parameters */
  addVariable( ini, "model", &pars.model, string_t, PARAM_FIXED );
  
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "Pb",
                                        pars.Pb, pars.PbErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "e",
                                        pars.e, pars.eErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "eps1",
                                        pars.eps1, pars.eps1Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "eps2",
                                        pars.eps2, pars.eps2Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "T0",
                                        pars.T0, pars.T0Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "Tasc",
                                        pars.Tasc, pars.TascErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "x",
                                        pars.x, pars.xErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "w0",
                                        pars.w0, pars.w0Err );

  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "Pb2",
                                        pars.Pb2, pars.Pb2Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "e2",
                                        pars.e2, pars.e2Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "T02",
                                        pars.T02, pars.T02Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "x2",
                                        pars.x2, pars.x2Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "w02",
                                        pars.w02, pars.w02Err );

  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "Pb3",
                                        pars.Pb3, pars.Pb3Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "e3",
                                        pars.e3, pars.e3Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "T03",
                                        pars.T03, pars.T03Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "x3",
                                        pars.x3, pars.x3Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "w03",
                                        pars.w03, pars.w03Err );

  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "xpbdot",
                                        pars.xpbdot, pars.xpbdotErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "eps1dot",
                                        pars.eps1dot, pars.eps1dotErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "eps2dot",
                                        pars.eps2dot, pars.eps2dotErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "wdot",
                                        pars.wdot, pars.wdotErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "gamma",
                                        pars.gamma, pars.gammaErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "Pbdot",
                                        pars.Pbdot, pars.PbdotErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "xdot",
                                        pars.xdot, pars.xdotErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "edot",
                                        pars.edot, pars.edotErr );
 
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "s",
                                        pars.s, pars.sErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "dr",
                                        pars.dr, pars.drErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "dth",
                                        pars.dth, pars.dthErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "a0",
                                        pars.a0, pars.a0Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "b0",
                                        pars.b0, pars.b0Err );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "M",
                                        pars.M, pars.MErr );
  withphase = add_variable_scale_prior( ini, scaleFac, priorArgs, "m2",
                                        pars.m2, pars.m2Err );

  return withphase;
}


UINT4 add_variable_scale_prior( LALVariables *var, LALVariables *scale, 
                                LALVariables *prior, const char *name, 
                                REAL8 value, REAL8 sigma ){
  REAL8 scaleVal = 1.;
  ParamVaryType vary;
  CHAR scaleName[VARNAME_MAX] = "";
 
  UINT4 nonzero = 0;
  
  /* if the sigma is non-zero then set a Gaussian prior */
  if ( sigma != 0. ){
    vary = PARAM_LINEAR;
    
    /* set the prior to a Gaussian prior with mean value and sigma */
    addGaussianPrior( prior, name, (void *)&value, (void *)&sigma, REAL8_t );
    
    nonzero = 1;
  }
  else vary = PARAM_FIXED;
  
  /* add the variable */
  addVariable( var, name, &value, REAL8_t, vary );
  
  /* add the initial scale factor of 1 */
  sprintf( scaleName, "%s_scale", name );
  addVariable( scale, scaleName, &scaleVal, REAL8_t, PARAM_FIXED );
              
  return nonzero;
}


void initialiseProposal( LALInferenceRunState *runState )
/* sets up the parameters to be varied, and the extent of the proposal
   distributions from a proposal distribution file */
{
  CHAR *propfile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  FILE *fp=NULL;
  
  CHAR tempPar[VARNAME_MAX] = "";
  REAL8 low, high;
  
  LALIFOData *data = runState->data;
  
  ppt = getProcParamVal( commandLine, "--prop-file" );
  if( ppt ){
  propfile =
    XLALStringDuplicate( getProcParamVal(commandLine,"--prop-file")->value );
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
    VariableType type;
    
    REAL8 scale = 0.;
    VariableType scaleType;
    CHAR tempParScale[VARNAME_MAX] = "";
    CHAR tempParPrior[VARNAME_MAX] = "";
 
    LALIFOData *datatemp = data;
    
    ParamVaryType varyType;
    
    if( high < low ){
      fprintf(stderr, "Error... In %s the %s parameters ranges are wrongly \
set.\n", propfile, tempPar);
      exit(3);
    }
    
    sprintf(tempParScale, "%s_scale", tempPar);
    sprintf(tempParPrior, "%s_gaussian_mean", tempPar);

    tempVar = *(REAL8*)getVariable( runState->currentParams, tempPar );
    type = getVariableType( runState->currentParams, tempPar );
    varyType = getVariableVaryType( runState->currentParams, tempPar );
    
    /* remove variable value */
    removeVariable( runState->currentParams, tempPar );
    
    /* if a Gaussian prior has already been defined (i.e. from the par file)
       remove this and overwrite with values from the propfile */
    if ( checkVariable(runState->priorArgs, tempParPrior) )
      removeGaussianPrior( runState->priorArgs, tempPar );
    
    scale = high;
    
    /* set the scale factor to be the high value of the prior */
    while( datatemp ){
      scaleType = getVariableType( datatemp->dataParams, tempParScale );
      removeVariable( datatemp->dataParams, tempParScale );
    
      addVariable( datatemp->dataParams, tempParScale, &scale, scaleType,
        PARAM_FIXED );
    
      datatemp = datatemp->next;
    }
      
    /* scale variable and priors */
    tempVar /= scale;
    low /= scale;
    high /= scale;
    
    /* re-add variable */
    if( !strcmp(tempPar, "phi0") ){
      addVariable( runState->currentParams, tempPar, &tempVar, type,
                   PARAM_CIRCULAR );
    }
    else{
      addVariable( runState->currentParams, tempPar, &tempVar, type,
                   PARAM_LINEAR );
    }
    /* Add the prior variables */
    addMinMaxPrior( runState->priorArgs, tempPar, (void *)&low, (void *)&high,
                    type );
    
    /* if there is a phase parameter defined in the proposal then set withphase
       to 1 */
    if ( !strcmp(tempPar, "phi0") || !strcmp(tempPar, "h0") || 
         !strcmp(tempPar, "cosiota") || !strcmp(tempPar, "psi") ){
      UINT4 withphase = 1;
      setVariable( data->dataParams, "withphase", &withphase );
    }
  }
  
  /* check for any parameters with Gaussian priors and rescale to mean value */
  LALVariableItem *checkPrior = runState->currentParams->head;
  for( ; checkPrior ; checkPrior = checkPrior->next ){
    REAL8 scale = 0.;
    REAL8 tempVar;
    
    REAL8 mu, sigma;
    
    LALIFOData *datatemp = data;
    
    VariableType scaleType;
    CHAR tempParScale[VARNAME_MAX] = "";
    CHAR tempParPrior[VARNAME_MAX] = "";

    sprintf(tempParPrior,"%s_gaussian_mean",checkPrior->name);

    if ( checkVariable(runState->priorArgs, tempParPrior) ){
      tempVar = *(REAL8 *)checkPrior->value;
      
      /* get the mean and standard deviation of the Gaussian prior */
      getGaussianPrior( runState->priorArgs, checkPrior->name, (void *)&mu,
                        (void *)&sigma );
      
      /* set the scale factor to be the mean value */
      scale = mu;
      tempVar /= scale;
      
      /* scale the parameter value and reset it */
      memcpy( checkPrior->value, &tempVar, typeSize[checkPrior->type] );
      
      mu /= scale;
      sigma /= scale;
      
      /* remove the Gaussian prior values and reset as scaled values */
      removeGaussianPrior( runState->priorArgs, checkPrior->name );
      addGaussianPrior( runState->priorArgs, checkPrior->name, 
                        (void *)&mu, (void *)&sigma, checkPrior->type );
      
      sprintf(tempParScale, "%s_scale", checkPrior->name);
        
      /* set scale factor in data structure */
      while( datatemp ){
        scaleType = getVariableType( datatemp->dataParams, tempParScale );
        removeVariable( datatemp->dataParams, tempParScale );
    
        addVariable( datatemp->dataParams, tempParScale, &scale, scaleType,
                     PARAM_FIXED );
      
        datatemp = datatemp->next;
      }
    }
  }
    
  return;
}


void setupLivePointsArray( LALInferenceRunState *runState ){
/* Set up initial basket of live points, drawn from prior,
   by copying runState->currentParams to all entries in the array*/
  UINT4 Nlive = (UINT4)*(INT4 *)getVariable(runState->algorithmParams,"Nlive");
  UINT4 i;
  REAL8Vector *logLs;
        
  LALVariableItem *current;

  /* Allocate the array */
  runState->livePoints = XLALCalloc( Nlive, sizeof(LALVariables *) );
  
  if( runState->livePoints == NULL ){
    fprintf(stderr,"Unable to allocate memory for %i live points\n",Nlive);
    exit(1);
  }

  logLs = XLALCreateREAL8Vector( Nlive );
     
  addVariable( runState->algorithmParams, "logLikelihoods",
               &logLs, REAL8Vector_t, PARAM_FIXED);
               
  fprintf(stdout, "Sprinkling %i live points, may take some time\n", Nlive);
  
  for( i=0; i<Nlive; i++){
    runState->livePoints[i] = XLALCalloc( 1, sizeof(LALVariables) );
                
    /* Copy the param structure */
    copyVariables( runState->currentParams, runState->livePoints[i] );
    
    /* Sprinkle the varying points among prior */
    do{
      for( current=runState->livePoints[i]->head; current!=NULL;
           current=current->next){
        CHAR tempParPrior[VARNAME_MAX] = "";
        UINT4 gp = 0;
      
        sprintf(tempParPrior,"%s_gaussian_mean",current->name);
      
        if( checkVariable( runState->priorArgs, tempParPrior ) ) gp = 1;
        
        if( current->vary==PARAM_CIRCULAR || current->vary==PARAM_LINEAR )
        {
          switch (current->type){
            case REAL4_t:
            {
              REAL4 tmp;
              REAL4 min, max, mu, sigma;
                                                       
              if( gp ){
                getGaussianPrior( runState->priorArgs, current->name, 
                                  (void *)&mu, (void *)&sigma );
                tmp = mu + gsl_ran_gaussian(runState->GSLrandom, (double)sigma);
              }
              else{
                getMinMaxPrior( runState->priorArgs, current->name, 
                                (void *)&min, (void *)&max );
                tmp = min + (max-min)*gsl_rng_uniform( runState->GSLrandom );
              }
                                                       
              setVariable( runState->livePoints[i], current->name, &tmp );
              break;
            }
            case REAL8_t:
            {
              REAL8 tmp;
              REAL8 min, max, mu, sigma;
                                                       
              if( gp ){
                getGaussianPrior( runState->priorArgs, current->name, 
                                  (void *)&mu, (void *)&sigma );
                tmp = mu + gsl_ran_gaussian(runState->GSLrandom, (double)sigma);
              }
              else{
                getMinMaxPrior( runState->priorArgs, current->name, 
                                (void *)&min, (void *)&max );
                tmp = min + (max-min)*gsl_rng_uniform( runState->GSLrandom );
              }
                                                       
              setVariable( runState->livePoints[i], current->name, &tmp );
              break;
            }
            case INT4_t:
            {
              INT4 tmp;
              INT4 min,max;
                                                       
              getMinMaxPrior( runState->priorArgs, current->name, (void *)&min,
                              (void *)&max );
                                                       
              tmp = min + (max-min)*gsl_rng_uniform(runState->GSLrandom);
                                                       
              setVariable( runState->livePoints[i], current->name, &tmp );
              break;
            }
            case INT8_t:
            {
              INT8 tmp;
              INT8 min, max;
                                                       
              getMinMaxPrior( runState->priorArgs, current->name, (void *)&min,
                              (void *)&max );
                                                       
              tmp = min + (max-min)*gsl_rng_uniform(runState->GSLrandom);
                                                       
              setVariable( runState->livePoints[i], current->name, &tmp);
              break;
            }
            default:
              fprintf(stderr,"Trying to randomise a non-numeric parameter!");
          }
        }
      }
               
    }while( runState->prior( runState,runState->livePoints[i] ) == -DBL_MAX );
    
    /* Populate log likelihood */           
    logLs->data[i] = runState->likelihood( runState->livePoints[i],
                                           runState->data, runState->template);
  }
        
}

/*------------------- END INITIALISATION FUNCTIONS ---------------------------*/


/******************************************************************************/
/*                     LIKELIHOOD AND PRIOR FUNCTIONS                         */
/******************************************************************************/

REAL8 pulsar_log_likelihood( LALVariables *vars, LALIFOData *data,
                             LALTemplateFunction *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  UINT4 i = 0;
  
  while ( data ){
    UINT4 j = 0, count = 0, cl = 0;
    UINT4 length = 0, chunkMin, chunkMax;
    REAL8 chunkLength = 0.;
    REAL8 logliketmp = 0.;

    REAL8 sumModel = 0., sumDataModel = 0.;
    REAL8 chiSquare = 0.;
    COMPLEX16 B, M;

    REAL8Vector *exclamation = NULL; /* all factorials up to chunkMax */

    INT4 first = 0, through = 0;
  
    REAL8Vector *sumData = NULL;
    UINT4Vector *chunkLengths = NULL;

    sumData = *(REAL8Vector **)getVariable( data->dataParams, "sumData" );
    chunkLengths = *(UINT4Vector **)getVariable( data->dataParams,
                                                 "chunkLength" );
    chunkMin = *(INT4*)getVariable( data->dataParams, "chunkMin" );
    chunkMax = *(INT4*)getVariable( data->dataParams, "chunkMax" );
  
    exclamation = *(REAL8Vector **)getVariable( data->dataParams,
                                                "logFactorial" );
  
    /* copy model parameters to data parameters */
    copyVariables( vars, data->modelParams );
  
    /* get pulsar model */
    get_model( data );
  
    length = data->compTimeData->data->length;
  
    for( i = 0 ; i < length ; i += chunkLength ){
      chunkLength = (REAL8)chunkLengths->data[count];
    
      /* skip section of data if its length is less than the minimum allowed
        chunk length */
      if( chunkLength < chunkMin ){
        count++;

        if( through == 0 ) first = 0;

        continue;
      }

      through = 1;

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
      }
 
      chiSquare = sumData->data[count];
      chiSquare -= 2.*sumDataModel;
      chiSquare += sumModel;
      
      if( first == 0 ){
        logliketmp = 2. * (chunkLength - 1.)*LAL_LN2;
        first++;
      }
      else logliketmp += 2. * (chunkLength - 1.)*LAL_LN2;

      logliketmp += 2. * exclamation->data[(INT4)chunkLength];
      logliketmp -= chunkLength*log(chiSquare);
    
      count++;
    }
  
    loglike += logliketmp;
  
    data = data->next;
  }
  
  return loglike;
}


REAL8 priorFunction( LALInferenceRunState *runState, LALVariables *params ){
  LALIFOData *data = runState->data;
  (void)runState;
  LALVariableItem *item = params->head;
  REAL8 min, max, mu, sigma, prior = 0, value = 0.;
  
  for(; item; item = item->next ){
    /* get scale factor */
    CHAR scalePar[VARNAME_MAX] = "";
    CHAR priorPar[VARNAME_MAX] = "";
    REAL8 scale;
    
    if( item->vary == PARAM_FIXED || item->vary == PARAM_OUTPUT ){ continue; }
    
    sprintf(scalePar, "%s_scale", item->name);
    scale = *(REAL8 *)getVariable( data->dataParams, scalePar );
    
    if( item->vary == PARAM_LINEAR || item->vary == PARAM_CIRCULAR ){
      sprintf(priorPar, "%s_gaussian_mean", item->name);
      /* Check for a gaussian */
      if ( checkVariable(runState->priorArgs, priorPar) ){
        getGaussianPrior( runState->priorArgs, item->name, (void *)&mu, 
                          (void *)&sigma );
      
       value = (*(REAL8 *)item->value) * scale;
       mu *= scale;
       sigma *= scale;
       prior -= log(sqrt(2.*LAL_PI)*sigma);
       prior -= (value - mu)*(value - mu) / (2.*sigma*sigma);
      }
      /* Otherwise use a flat prior */
      else{
	getMinMaxPrior( runState->priorArgs, item->name, (void *)&min, 
                      (void *)&max );

        if( (*(REAL8 *) item->value)*scale < min*scale || 
          (*(REAL8 *)item->value)*scale > max*scale ){
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
void get_pulsar_model( LALIFOData *data ){
  INT4 i = 0, length = 0;
  
  BinaryPulsarParams pars;
  
  REAL8Vector *dphi = NULL;
  
  REAL8 rescale = 1.;
  
  UINT4 withphase = 0; /* set this to 1 if the phase model needs calculating */
  
  /* set model parameters (including rescaling) */
  rescale = *(REAL8*)getVariable( data->dataParams, "h0_scale" );
  pars.h0 = *(REAL8*)getVariable( data->modelParams, "h0" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "cosiota_scale" );
  pars.cosiota = *(REAL8*)getVariable( data->modelParams, "cosiota" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "psi_scale" );
  pars.psi = *(REAL8*)getVariable( data->modelParams, "psi" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "phi0_scale" );
  pars.phi0 = *(REAL8*)getVariable( data->modelParams, "phi0" ) * rescale;
  
  /* set the potentially variable parameters */
  rescale = *(REAL8*)getVariable( data->dataParams, "pepoch_scale" );
  pars.pepoch = *(REAL8*)getVariable( data->modelParams, "pepoch" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "posepoch_scale" );
  pars.posepoch = *(REAL8*)getVariable( data->modelParams, "posepoch" ) *
    rescale;
  
  rescale = *(REAL8*)getVariable( data->dataParams, "ra_scale" );
  pars.ra = *(REAL8*)getVariable( data->modelParams, "ra" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "pmra_scale" );
  pars.pmra = *(REAL8*)getVariable( data->modelParams, "pmra" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "dec_scale" );
  pars.dec = *(REAL8*)getVariable( data->modelParams, "dec" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "pmdec_scale" );
  pars.pmdec = *(REAL8*)getVariable( data->modelParams, "pmdec" ) * rescale;
  
  rescale = *(REAL8*)getVariable( data->dataParams, "f0_scale" );
  pars.f0 = *(REAL8*)getVariable( data->modelParams, "f0" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "f1_scale" );
  pars.f1 = *(REAL8*)getVariable( data->modelParams, "f1" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "f2_scale" );
  pars.f2 = *(REAL8*)getVariable( data->modelParams, "f2" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "f3_scale" );
  pars.f3 = *(REAL8*)getVariable( data->modelParams, "f3" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "f4_scale" );
  pars.f4 = *(REAL8*)getVariable( data->modelParams, "f4" ) * rescale;
  rescale = *(REAL8*)getVariable( data->dataParams, "f5_scale" );
  pars.f5 = *(REAL8*)getVariable( data->modelParams, "f5" ) * rescale;
  
  pars.model = *(CHAR**)getVariable( data->modelParams, "model" );

  /* binary parameters */
  if( pars.model != NULL ){
    rescale = *(REAL8*)getVariable( data->dataParams, "e_scale" );
    pars.e = *(REAL8*)getVariable( data->modelParams, "e" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "w0_scale" );
    pars.w0 = *(REAL8*)getVariable( data->modelParams, "w0" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "Pb_scale" );
    pars.Pb = *(REAL8*)getVariable( data->modelParams, "Pb" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "x_scale" );
    pars.x = *(REAL8*)getVariable( data->modelParams, "x" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "T0_scale" );
    pars.T0 = *(REAL8*)getVariable( data->modelParams, "T0" ) * rescale;
    
    rescale = *(REAL8*)getVariable( data->dataParams, "e2_scale" );
    pars.e2 = *(REAL8*)getVariable( data->modelParams, "e2" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "w02_scale" );
    pars.w02 = *(REAL8*)getVariable( data->modelParams, "w02" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "Pb2_scale" );
    pars.Pb2 = *(REAL8*)getVariable( data->modelParams, "Pb2" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "x2_scale" );
    pars.x2 = *(REAL8*)getVariable( data->modelParams, "x2" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "T02_scale" );
    pars.T02 = *(REAL8*)getVariable( data->modelParams, "T02" ) * rescale;
    
    rescale = *(REAL8*)getVariable( data->dataParams, "e3_scale" );
    pars.e3 = *(REAL8*)getVariable( data->modelParams, "e3" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "w03_scale" );
    pars.w03 = *(REAL8*)getVariable( data->modelParams, "w03" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "Pb3_scale" );
    pars.Pb3 = *(REAL8*)getVariable( data->modelParams, "Pb3" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "x3_scale" );
    pars.x3 = *(REAL8*)getVariable( data->modelParams, "x3" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "T03_scale" );
    pars.T03 = *(REAL8*)getVariable( data->modelParams, "T03" ) * rescale;
    
    rescale = *(REAL8*)getVariable( data->dataParams, "xpbdot_scale" );
    pars.xpbdot = *(REAL8*)getVariable( data->modelParams, "xpbdot" ) * rescale;
    
    rescale = *(REAL8*)getVariable( data->dataParams, "eps1_scale" );
    pars.eps1 = *(REAL8*)getVariable( data->modelParams, "eps1" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "eps2_scale" );
    pars.eps2 = *(REAL8*)getVariable( data->modelParams, "eps2" ) * rescale;   
     rescale = *(REAL8*)getVariable( data->dataParams, "eps1dot_scale" );
    pars.eps1dot = *(REAL8*)getVariable( data->modelParams, "eps1dot" ) *
      rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "eps2dot_scale" );
    pars.eps2dot = *(REAL8*)getVariable( data->modelParams, "eps2dot" ) *
      rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "Tasc_scale" );
    pars.Tasc = *(REAL8*)getVariable( data->modelParams, "Tasc" ) * rescale;
    
    rescale = *(REAL8*)getVariable( data->dataParams, "wdot_scale" );
    pars.wdot = *(REAL8*)getVariable( data->modelParams, "wdot" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "gamma_scale" );
    pars.gamma = *(REAL8*)getVariable( data->modelParams, "gamma" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "Pbdot_scale" );
    pars.Pbdot = *(REAL8*)getVariable( data->modelParams, "Pbdot" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "xdot_scale" );
    pars.xdot = *(REAL8*)getVariable( data->modelParams, "xdot" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "edot_scale" );
    pars.edot = *(REAL8*)getVariable( data->modelParams, "edot" ) * rescale;
    
    rescale = *(REAL8*)getVariable( data->dataParams, "s_scale" );
    pars.s = *(REAL8*)getVariable( data->modelParams, "s" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "dr_scale" );
    pars.dr = *(REAL8*)getVariable( data->modelParams, "dr" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "dth_scale" );
    pars.dth = *(REAL8*)getVariable( data->modelParams, "dth" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "a0_scale" );
    pars.a0 = *(REAL8*)getVariable( data->modelParams, "a0" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "b0_scale" );
    pars.b0 = *(REAL8*)getVariable( data->modelParams, "b0" ) * rescale; 

    rescale = *(REAL8*)getVariable( data->dataParams, "M_scale" );
    pars.M = *(REAL8*)getVariable( data->modelParams, "M" ) * rescale;
    rescale = *(REAL8*)getVariable( data->dataParams, "m2_scale" );
    pars.m2 = *(REAL8*)getVariable( data->modelParams, "m2" ) * rescale;
  }

  get_amplitude_model( pars, data );
  length = data->compModelData->data->length;
  
  /* assume that timeData vector within the LALIFOData structure contains the
     phase calculated using the initial (heterodyne) values of the phase
     parameters */
  
  withphase = *(UINT4 *)getVariable( data->dataParams, "withphase" );
  
  /* get difference in phase and perform extra heterodyne with it */ 
  if ( withphase && (dphi = get_phase_model( pars, data )) != NULL ){
    
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
  
  XLALDestroyREAL8Vector( dphi );
}


REAL8Vector *get_phase_model( BinaryPulsarParams params, LALIFOData *data ){
  static LALStatus status;

  INT4 i = 0, length = 0;

  REAL8 T0 = 0., DT = 0., DTplus = 0., deltat = 0., deltat2 = 0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */

  EarthState earth, earth2;
  EmissionTime emit, emit2;
  REAL8 emitdt = 0.;

  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;

  BarycenterInput *bary = NULL;
  
  REAL8Vector *phis = NULL;

  /* if edat is NULL then return a NULL poniter */
  if( data->ephem == NULL )
    return NULL;
  
  /* copy barycenter and ephemeris data */
  bary = (BarycenterInput*)XLALCalloc( 1, sizeof(BarycenterInput) );
  memcpy( &bary->site, data->detector, sizeof(LALDetector) );
  
  bary->alpha = params.ra;
  bary->delta = params.dec;
  
   /* set the position and frequency epochs if not already set */
  if( params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if( params.posepoch == 0. && params.pepoch != 0. )
    params.posepoch = params.pepoch;

  length = data->dataTimes->length;
  
  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );
 
  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( params.px != 0. )
    bary->dInv = params.px*1e-3*LAL_C_SI/LAL_PC_SI;
  else if( params.dist != 0. )
    bary->dInv = LAL_C_SI/(params.dist*1e3*LAL_PC_SI);
  else
    bary->dInv = 0.;

  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8( &data->dataTimes->data[i] );
    
    T0 = params.pepoch;

    DT = realT - T0;

    /* only do call the barycentring routines every 30 minutes, otherwise just
       linearly interpolate between them */
    if( i==0 || DT > DTplus ){
      bary->tgps = data->dataTimes->data[i];

      bary->delta = params.dec + (realT-params.posepoch) * params.pmdec;
      bary->alpha = params.ra + (realT-params.posepoch) *
         params.pmra/cos(bary->delta);
     
      /* call barycentring routines */
      LAL_CALL( LALBarycenterEarth( &status, &earth, &bary->tgps, data->ephem ),
                &status );
      
      LAL_CALL( LALBarycenter( &status, &emit, bary, &earth ), &status );

      /* add interptime to the time */
      DTplus = DT + interptime;
      XLALGPSAdd( &bary->tgps, interptime );

      /* No point in updating the positions as difference will be tiny */
      LAL_CALL( LALBarycenterEarth( &status, &earth2, &bary->tgps, 
                                    data->ephem ), &status );
      LAL_CALL( LALBarycenter( &status, &emit2, bary, &earth2), &status );
    }

    /* linearly interpolate to get emitdt */
    emitdt = emit.deltaT + (DT - (DTplus - interptime)) *
      (emit2.deltaT - emit.deltaT)/interptime;

    /* check if need to perform binary barycentring */
    if( params.model != NULL ){
      binput.tb = realT + emitdt;

      XLALBinaryPulsarDeltaT( &boutput, &binput, &params );

      deltat = DT + emitdt + boutput.deltaT;
    }
    else
      deltat = DT + emitdt;

    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = 2.*deltat*(params.f0 + 
      inv_fact[2]*params.f1*deltat +
      inv_fact[3]*params.f2*deltat2 +
      inv_fact[4]*params.f3*deltat*deltat2 +
      inv_fact[5]*params.f4*deltat2*deltat2 +
      inv_fact[6]*params.f5*deltat2*deltat2*deltat);
  }

  XLALFree( bary );

  return phis;
}


void get_amplitude_model( BinaryPulsarParams pars, LALIFOData *data ){
  INT4 i = 0, length;
  
  REAL8 psteps, tsteps;
  INT4 psibin, timebin;
  REAL8 tstart;
  REAL8 plus, cross;
  REAL8 T;
  REAL8 Xplus, Xcross;
  REAL8 Xpcosphi, Xccosphi, Xpsinphi, Xcsinphi;
  REAL4 sinphi, cosphi;
  
  gsl_matrix *LU_Fplus, *LU_Fcross;
  
  length = data->dataTimes->length;
  
  /* set lookup table parameters */
  psteps = *(INT4*)getVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)getVariable( data->dataParams, "timeSteps" );
  
  LU_Fplus = *(gsl_matrix**)getVariable( data->dataParams, "LU_Fplus");
  LU_Fcross = *(gsl_matrix**)getVariable( data->dataParams, "LU_Fcross");
  
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
  
  for( i=0; i<length; i++ ){
    /* set the psi bin for the lookup table */
    psibin = (INT4)ROUND( ( pars.psi + LAL_PI/4. ) * ( psteps-1. )/LAL_PI_2 );

    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = fmod( XLALGPSGetREAL8(&data->dataTimes->data[i]) - tstart,
              LAL_DAYSID_SI );
    timebin = (INT4)fmod( ROUND(T*tsteps/LAL_DAYSID_SI), tsteps );

    plus = gsl_matrix_get( LU_Fplus, psibin, timebin );
    cross = gsl_matrix_get( LU_Fcross, psibin, timebin );
    
    /* create the complex signal amplitude model */
    data->compModelData->data->data[i].re = plus*Xpcosphi + cross*Xcsinphi;
    data->compModelData->data->data[i].im = plus*Xpsinphi - cross*Xccosphi;
  }
  
}


/* calculate the likelihood for the data just being noise - we still have a
students-t likelihood, so this basically involves taking that likelihood and
setting the signal to zero */
REAL8 noise_only_model( LALIFOData *data ){
  REAL8 logL = 0.;
  UINT4 i = 0;
  
  UINT4Vector *chunkLengths = NULL;
  REAL8Vector *sumData = NULL;
  
  INT4 chunkMin = 0, chunkMax = 0;
  REAL8 chunkLength = 0.;
  
  chunkMin = *(INT4*)getVariable( data->dataParams, "chunkMin" );
  chunkMax = *(INT4*)getVariable( data->dataParams, "chunkMax" );
  
  chunkLengths = *(UINT4Vector **)getVariable( data->dataParams, 
                                               "chunkLength" );
  sumData = *(REAL8Vector **)getVariable( data->dataParams, "sumData" );
  
  for (i=0; i<chunkLengths->length; i++){
    chunkLength = (REAL8)chunkLengths->data[i];
    
    logL += 2. * (chunkLength - 1.) * LAL_LN2;
    logL += 2. * log_factorial((INT4)chunkLength);
    logL -= chunkLength * log(sumData->data[i]);
  }
  
  return logL;
}

/*------------------------ END OF MODEL FUNCTIONS ----------------------------*/


/******************************************************************************/
/*                            HELPER FUNCTIONS                                */
/******************************************************************************/

/* function to get the lengths of consecutive chunks of data */
UINT4Vector *get_chunk_lengths( LALIFOData *data, INT4 chunkMax ){
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
UINT4Vector *chop_n_merge( LALIFOData *data, INT4 chunkMin ){
  UINT4 j = 0;
  
  UINT4Vector *chunkLengths = NULL;
  UINT4Vector *chunkIndex = NULL;
  
  chunkIndex = chop_data( data->compTimeData->data, chunkMin );
  
  merge_data( data->compTimeData->data, chunkIndex );
  
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
    
    if ( (fpsegs = fopen("data_segment_list.txt", "a")) == NULL ){
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


/* function to find change points and chop up the data */
UINT4Vector *chop_data( COMPLEX16Vector *data, INT4 chunkMin ){
  UINT4Vector *chunkIndex = NULL;
  
  UINT4 length = data->length;
  
  REAL8 logodds = 0;
  UINT4 changepoint = 0;
  
  REAL8 threshold = 0.; /* may need tuning or setting globally */
  
  chunkIndex = XLALCreateUINT4Vector( 1 );
  
  changepoint = find_change_point( data, &logodds );
  
  if ( logodds > threshold && changepoint > (UINT4)chunkMin 
    && length - changepoint > (UINT4)chunkMin ){
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


/* this function will go through the data and compare the "evidence" for it
   being Gaussian with a single standard deviation, or it being two Gaussian's
   with independent standard deviations. It will return the index of the index
   of the point at with the data statistics change, and the log odds ratio of
   the two hypotheses at that point. */ 
UINT4 find_change_point( COMPLEX16Vector *data, REAL8 *logodds ){
  UINT4 changepoint = 0, i = 0;
  UINT4 length = data->length;
  
  REAL8 datameanre = 0., datameanim = 0., datasum = 0.;
  
  REAL8 logsingle = 0., logtot = -LAL_REAL8_MAX;
  REAL8 logdouble = 0., logdouble_min = -LAL_REAL8_MAX;
  REAL8 logratio = 0.;
  
  /* calculate the mean of the data (put real and imaginary parts together) */
  for (i = 0; i < length; i++){
    datameanre += data->data[i].re;
    datameanim += data->data[i].im;
  }
    
  datameanre /= (REAL8)length;
  datameanim /= (REAL8)length;
  
  /* calculate the sum of the data minus the mean squared */
  for (i = 0; i < length; i++){
    datasum += (data->data[i].re - datameanre) * (data->data[i].re -
      datameanre);
    datasum += (data->data[i].im - datameanim) * (data->data[i].im -
      datameanim);
  }
  
  /* calculate the evidence that the data consists of a Gaussian data with a
     single standard deviation */
  logsingle = 2. * ( ((REAL8)length - 1.)*LAL_LN2 + log_factorial( length ) ) -
    (REAL8)length * log( datasum );
  
  /* go through each possible change point and calculate the evidence for the
     data consisting of two independent Gaussian's either side of the change
     point. Also calculate the total evidence for any change point */
  /* Don't allow single points, so start at the second data point. */
  for (i = 2; i < length-1; i++){
    UINT4 ln1 = i, ln2 = (length - i), j = 0;
    
    REAL8 mu1re = 0., mu1im = 0., mu2re = 0., mu2im = 0.;
    REAL8 sum1 = 0, sum2 = 0.;
    REAL8 log_1 = 0., log_2 = 0.;
    
    /* get means of the two segments */
    for (j = 0; j < ln1; j++){
      mu1re += data->data[j].re;
      mu1im += data->data[j].im;
    }
    
    mu1re /= (REAL8)ln1;
    mu1im /= (REAL8)ln1;
    
    for (j = ln1; j < length; j++){
      mu2im += data->data[j].re;
      mu2im += data->data[j].im;
    }
    
    mu2re /= (REAL8)ln2;
    mu2im /= (REAL8)ln2;
    
    /* get the sum data squared */
    for (j = 0; j < ln1; j++){
      sum1 += (data->data[j].re - mu1re) * (data->data[j].re - mu1re);
      sum1 += (data->data[j].im - mu1im) * (data->data[j].im - mu1im);
    }
    
    for (j = ln1; j < length; j++){
      sum2 += (data->data[j].re - mu2re) * (data->data[j].re - mu2re);
      sum2 += (data->data[j].im - mu2im) * (data->data[j].im - mu2im);
    }
    
    /* get log evidences for the individual segments */
    log_1 = 2. * ( ((REAL8)ln1 - 1.)*LAL_LN2 + log_factorial( ln1 ) ) -
      (REAL8)ln1 * log( sum1 );
    log_2 = 2. * ( ((REAL8)ln2 - 1.)*LAL_LN2 + log_factorial( ln2 ) ) -
      (REAL8)ln2 * log( sum2 );
      
    /* get evidence for the two segments */
    logdouble = log_1 + log_2;
    
    /* add to total evidence for a change point */
    logtot += LOGPLUS(logtot, logdouble);
    
    /* find maximum value of logdouble and record that as the change point */
    if ( logdouble > logdouble_min ){
      changepoint = i;
      logdouble_min = logdouble;
    }
  }
  
  /* get the log odds ratio of segmented versus non-segmented model */
  logratio = logtot - logsingle;
  logodds = &logratio;
  
  return changepoint;
}


/* function to merge chunks */
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
      REAL8 meanmergedre = 0., meanmergedim = 0.; 
      REAL8 mu1re = 0., mu1im = 0., mu2re = 0., mu2im = 0.;
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
      
      n1 = cellends1 - cellstarts1 - 1;
      n2 = cellends2 - cellstarts2 - 1;
      nm = cellends2 - cellstarts1 - 1;
      
      for( i = cellstarts1; i < cellends1; i++ ){
        mu1re += data->data[i].re;
        mu1im += data->data[i].im;
      }
        
      for( i = cellstarts2; i < cellends2; i++ ){
        mu2re += data->data[i].re;
        mu2im += data->data[i].im;
      }
      
      meanmergedre = (mu1re + mu2re) / (REAL8)nm;
      meanmergedim = (mu1im + mu2im) / (REAL8)nm;
      mu1re /= (REAL8)n1;
      mu1im /= (REAL8)n1;
      mu2re /= (REAL8)n2;
      mu2im /= (REAL8)n2;
      
      for( i = cellstarts1; i < cellends1; i++ ){
        sum1 += (data->data[i].re - mu1re) * (data->data[i].re - mu1re);
        sum1 += (data->data[i].im - mu1im) * (data->data[i].im - mu1im);
      }
      
      for( i = cellstarts2; i < cellends2; i++ ){
        sum2 += (data->data[i].re - mu2re) * (data->data[i].re - mu2re);
        sum2 += (data->data[i].im - mu2im) * (data->data[i].im - mu2im);
      }
      
      for( i = cellstarts1; i < cellends2; i++ ){
        summerged += (data->data[i].re - meanmergedre) * (data->data[i].re -
          meanmergedre);
        summerged += (data->data[i].im - meanmergedim) * (data->data[i].im -
          meanmergedim);
      }
      
      /* calculated evidences */
      log_merged = 2. * ( ((REAL8)nm - 1.)*LAL_LN2 + log_factorial( nm ) ) -
        (REAL8)nm * summerged;
        
      log_individual = 2. * ( ((REAL8)n1 - 1)*LAL_LN2 + log_factorial( n1 ) ) -
        (REAL8)n1 * sum1;
      log_individual += 2. * ( ((REAL8)n2 - 1)*LAL_LN2 + log_factorial( n2 ) ) -
        (REAL8)n2 * sum2;
      
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
      memcpy( &segs->data[mergepoint], &segs->data[mergepoint+1], (ncells -
              mergepoint - 1) * sizeof(UINT4) );
      
      segs = XLALResizeUINT4Vector( segs, ncells - 1 );
    }
  }
}


/* a function to sum over the data */
REAL8Vector *sum_data( LALIFOData *data ){
  INT4 chunkLength = 0, length = 0, i = 0, j = 0, count = 0;
  COMPLEX16 B;
  REAL8Vector *sumData = NULL; 

  UINT4Vector *chunkLengths;
  
  chunkLengths = *(UINT4Vector **)getVariable( data->dataParams, 
                                               "chunkLength" );
  
  length = data->dataTimes->length + 1 -
    chunkLengths->data[chunkLengths->length - 1];

  sumData = XLALCreateREAL8Vector( chunkLengths->length );
  
  for( i = 0 ; i < length ; i+= chunkLength ){
    chunkLength = chunkLengths->data[count];
    sumData->data[count] = 0.;

    for( j = i ; j < i + chunkLength ; j++){
      B.re = data->compTimeData->data->data[j].re;
      B.im = data->compTimeData->data->data[j].im;

      /* sum up the data */
      sumData->data[count] += (B.re*B.re + B.im*B.im);
    }

    count++;
  }
  
  return sumData;
}


/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a sidereal day from the start time
(t0) and from -pi/4 to pi/4 in psi */
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


REAL8 log_factorial( UINT4 num ){
  UINT4 i = 0;
  REAL8 logFac = 0.;
  
  for( i=2 ; i <= num ; i++ ) logFac += log((REAL8)i);
  
  return logFac;
}


void rescaleOutput( LALInferenceRunState *runState ){
  /* Open orginal output output file */
  CHAR *outfile, outfiletmp[256] = "";
  FILE *fp = NULL, *fptemp = NULL;
  
  LALVariables *current=XLALCalloc(1,sizeof(LALVariables));
  
  ProcessParamsTable *ppt = getProcParamVal(runState->commandLine,"--outfile");
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
  
  /* open tempoary output file for reading */
  if( (fptemp = fopen(outfiletmp, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open temporary output file %s.\n",
            outfile);
    exit(3);
  }
 
  /* copy variables from runState to current (this seems to switch the order,
     which is required! */
  copyVariables(runState->currentParams, current);
  
  do{
    LALVariableItem *item = current->head;
    
    CHAR line[1000];
    CHAR value[128] = "";
    UINT4 i=0;
    
    /* read in one line of the file */
    if( fgets(line, 1000*sizeof(CHAR), fp) == NULL && !feof(fp) ){
      fprintf(stderr, "Error... cannot read line from file %s.\n", outfile);
      exit(3);
    }
    
    if( feof(fp) ) break;
    
    /* scan through line, get value and reprint out scaled value to temporary
       file */
    while (item != NULL){
      CHAR scalename[VARNAME_MAX] = "";
      REAL8 scalefac = 1.;
      
      if( i == 0 ){
        sprintf(value, "%s", strtok(line, "'\t'"));
        i++;
      }
      else
        sprintf(value, "%s", strtok(NULL, "'\t'"));
      
      sprintf(scalename, "%s_scale", item->name);
      
      if( strcmp(item->name, "model") ){
        scalefac = *(REAL8 *)getVariable( runState->data->dataParams, 
                                          scalename );
      }
      
      switch (item->type) {
        case INT4_t:
          fprintf(fptemp, "%d", (INT4)(atoi(value)*scalefac));
          break;
        case REAL4_t:
          fprintf(fptemp, "%e", atof(value)*scalefac);
          break;
        case REAL8_t:
          fprintf(fptemp, "%le", atof(value)*scalefac);
          break;
        case string_t:
          /* don't reprint out any string values */
          break;
        default:
          fprintf(stderr, "No type specified for %s.\n", item->name);
      }
      
      fprintf(fptemp, "\t");
      
      item = item->next;
    }
    
    /* last item in the line should be the logLikelihood (which is not in the
       currentParams structure) */
    sprintf(value, "%s", strtok(NULL, "'\t'"));
    fprintf(fptemp, "%lf\n", atof(value));
    
  }while( !feof(fp) );
  
  fclose(fp);
  fclose(fptemp);
  
  /* move the temporary file name to the standard outfile name */
  rename( outfiletmp, outfile );
  
  return;
}

/*----------------------- END OF HELPER FUNCTIONS ----------------------------*/
