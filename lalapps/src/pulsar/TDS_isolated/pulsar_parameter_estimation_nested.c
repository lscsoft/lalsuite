/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

#include "pulsar_parameter_estimation_nested.h"


#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0),
(1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* maximum number of different detectors */
#define MAXDETS 6

RCSID("$Id$");

/* global variables */
INT4 verbose = 0;

/* array to contain the log of factorials up to a certain number */
REAL8 *logfactorial = NULL;

/* set if phase parameters are being searched over and therefore the pulsar
   model requires phase evolution to be re-calculated */ 
UINT4 varyphase = 0; 

/* set if the sky position will be searched over, and therefore whether the
   solar system barycentring needs recalculating */
UINT4 varyskypos = 0; 

/* set if the binary system parameters will be searched over, and therefore
   whether the binary system barycentring needs recalculating */
UINT4 varybinary = 0; 


/* Usage format string */
#define USAGE \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --detectors         all IFOs with data to be analysed e.g. H1,H2\n\
                     (delimited by commas) (if generating fake data these\n\
                     should not be set)\n"\
" --pulsar            name of pulsar e.g. J0534+2200\n"\
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
"\n"


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
  ProcessParamsTable *param_table, *ppt;
  LALInferenceRunState runState;
  REAL8 logZnoise = 0.;
  CHAR *modeltype = NULL;
  
  REAL8Vector *logLikelihoods = NULL;
  
  /* Get ProcParamsTable from input arguments */
  param_table = LALInferenceParseCommandLine( argc, argv );
  runState.commandLine = param_table;
  
  /* Initialise the algorithm structures from the command line arguments */
  /* Include setting up random number generator etc */
  initialiseAlgorithm( &runState );
  
  /*read in data, 2 datastreams if pinsf, one if triaxial*/
  ppt = LALInferenceGetProcParamVal( param_table, "--model-type" );
  if( ppt ){
    modeltype  = XLALStringDuplicate( 
    LALInferenceGetProcParamVal( param_table, "--model-type" )->value );
    if ( !strcmp( modeltype, "triaxial" ) ){
      fprintf(stderr,"Model type set as triaxial\n");
      readPulsarData( &runState );
    }
    else if ( !strcmp( modeltype, "pinsf" ) ){
      fprintf(stderr,"Model type set as pinsf\n");
      readDoublePulsarData( &runState );
    }
    else{
      fprintf(stderr, "Error... model '%s' is not defined!\n", modeltype);
      exit(0);
    }
  }
  else readPulsarData( &runState );/*defaults to read in one stream of data per ifo*/
  
  /* read in data using command line arguments */
  /*readPulsarData( &runState );*/
  
  /* set the pulsar model type */
  setSignalModelType( &runState );
  
  /* set algorithm to use Nested Sampling */
  runState.algorithm = &LALInferenceNestedSamplingAlgorithm;
  runState.evolve = &LALInferenceNestedSamplingOneStep;
  
  if ( !strcmp( modeltype, "triaxial" ) ){
      runState.likelihood = &pulsar_log_likelihood;
    }
    else if ( !strcmp( modeltype, "pinsf" ) ){
      runState.likelihood = &pulsar_double_log_likelihood;
    }
    else{
      fprintf(stderr, "Error... model '%s' is not defined!\n", modeltype);
      exit(0);
    }
  /* set likelihood function */
  /*runState.likelihood = &pulsar_log_likelihood;*/
  
  /* set prior function */
  runState.prior = &priorFunction;
  
  /* set signal model/template */
  runState.template = get_pulsar_model;
  
  /* Generate the lookup tables and read parameters from par file */
  setupFromParFile( &runState );
  
  /* add injections if requested */
  injectSignal( &runState );
  
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
  setupLivePointsArray( &runState );
  
  logLikelihoods = *(REAL8Vector **)
    LALInferenceGetVariable( runState.algorithmParams, "logLikelihoods" );
  
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


void readPulsarData( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL, *ppt2 = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;
  
  CHAR *detectors = NULL;
  CHAR *outfile = NULL;
  CHAR *inputfile = NULL;
  
  CHAR *filestr = NULL;
  
  CHAR *efile = NULL;
  CHAR *sfile = NULL;
  
  CHAR *tempdets=NULL;
  CHAR *tempdet=NULL;
 
  REAL8 fpsds[MAXDETS];
  CHAR *fakestarts = NULL, *fakelengths = NULL, *fakedt = NULL;
  REAL8 fstarts[MAXDETS], flengths[MAXDETS], fdt[MAXDETS];
  
  CHAR dets[MAXDETS][256];
  INT4 numDets = 0, i = 0, numPsds = 0, numLengths = 0, numStarts = 0;
  INT4 numDt = 0, count = 0;
  UINT4 maxlen = 0;
  
  LALInferenceIFOData *ifodata = NULL;
  LALInferenceIFOData *prev = NULL;
  
  UINT4 seed = 0; /* seed for data generation */
  RandomParams *randomParams = NULL;
  
  runState->data = NULL;
  
  /* get the detectors - must */
  ppt = LALInferenceGetProcParamVal( commandLine, "--detectors" );
  ppt2 = LALInferenceGetProcParamVal( commandLine, "--fake-data" );
  if( ppt && !ppt2 ){
    detectors = XLALStringDuplicate( 
      LALInferenceGetProcParamVal(commandLine,"--detectors")->value );
      
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
  else if( ppt2 && !ppt ){
    detectors = XLALStringDuplicate( 
      LALInferenceGetProcParamVal(commandLine,"--fake-data")->value );
   
    if ( (numDets = count_csv( detectors )) > MAXDETS ){
      fprintf(stderr, "Error... too many detectors specified. Increase MAXDETS\
 to be greater than %d if necessary.\n", MAXDETS);
      exit(0);
    }
      
    /* get noise psds if specified */
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-psd");
    if( ppt ){
      CHAR *psds = NULL, *tmppsds = NULL, *tmppsd = NULL, psdval[256];
      
      psds = XLALStringDuplicate(
        LALInferenceGetProcParamVal(commandLine,"--fake-psd")->value );
      
      tmppsds = XLALStringDuplicate( psds );
      tempdets = XLALStringDuplicate( detectors );
      
      /* count the number of PSDs (comma seperated values) to compare to number
         of detectors */
      if( (numPsds = count_csv( psds )) != numDets ){
        fprintf(stderr, "Error... number of PSDs for fake data must be equal to\
 the number of detectors specified = %d!\n", numDets);
        exit(0);
      }
      
      for( i = 0; i < numDets; i++ ){
        CHAR *tmpstr = NULL;
        
        tmppsd = strsep( &tmppsds, "," );
        XLALStringCopy( psdval, tmppsd, strlen(tmppsd)+1 );
        fpsds[i] = atof(psdval);
        
        /* set detector */
        tempdet = strsep( &tempdets, "," );
        
        if( (tmpstr = strstr(tempdet, "A")) != NULL ) /* have advanced */
          XLALStringCopy( dets[i], tmpstr+1, strlen(tmpstr)+1 );
        else
          XLALStringCopy( dets[i], tempdet, strlen(tempdet)+1 );
      }
    }
    else{ /* get PSDs from model functions and set detectors */
      CHAR *parFile = NULL;
      BinaryPulsarParams pulsar;
      REAL8 pfreq = 0.;
      
      ppt = LALInferenceGetProcParamVal( runState->commandLine, "--par-file" );
      if( ppt == NULL ) { 
        fprintf(stderr,"Must specify --par-file!\n"); exit(1);
      }
      parFile = ppt->value;
        
      /* get the pulsar parameters to give a value of f */
      XLALReadTEMPOParFile( &pulsar, parFile );
      
      /* putting in pulsar frequency at 2*f here - this will need to be
         check when making the code more flexible for new models */
      pfreq = 2.0 * pulsar.f0;
      
      tempdets = XLALStringDuplicate( detectors );
      
      for( i = 0; i < numDets; i++ ){
        LALStatus status;
        
        CHAR *tmpstr = NULL;
        REAL8 psdvalf = 0.;
        
        numPsds++;
       
        tempdet = strsep( &tempdets, "," );
        
        if( (tmpstr = strstr(tempdet, "A")) != NULL ){ /* have Advanced */
          XLALStringCopy( dets[i], tmpstr+1, strlen(tmpstr)+1 );
          
          if( !strcmp(dets[i], "H1") || !strcmp(dets[i], "L1") ||  
              !strcmp(dets[i], "H2") ){ /* ALIGO */
            LALAdvLIGOPsd( &status, &psdvalf, pfreq );
            psdvalf *= 1.e-49; /* scale factor in LALAdvLIGOPsd.c */
          }
          else if( !strcmp(dets[i], "V1") ){ /* AVirgo */
            LALEGOPsd( &status, &psdvalf, pfreq ); 
          }
          else{
            fprintf(stderr, "Error... trying to use Advanced detector that is\
 not available!\n");
            exit(0);
          }
        }
        else{ /* initial detector */
          XLALStringCopy( dets[i], tempdet, strlen(tempdet)+1 );
          
          if( !strcmp(dets[i], "H1") || !strcmp(dets[i], "L1") ||  
              !strcmp(dets[i], "H2") ){ /* Initial LIGO */
            psdvalf = XLALLIGOIPsd( pfreq );
          
            /* divide H2 psd by 2 */
            if( !strcmp(dets[i], "H2") ) psdvalf /= 2.; 
          }
          else if( !strcmp(dets[i], "V1") ){ /* Initial Virgo */
            LALVIRGOPsd( &status, &psdvalf, pfreq ); 
          }
          else if( !strcmp(dets[i], "G1") ){ /* GEO 600 */
            LALGEOPsd( &status, &psdvalf, pfreq );
            psdvalf *= 1.e-46; /* scale factor in LALGEOPsd.c */
          }
          else if( !strcmp(dets[i], "T1") ){ /* TAMA300 */
            LALTAMAPsd( &status, &psdvalf, pfreq );
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
    
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-starts");
    if( ppt ){
      CHAR *tmpstarts = NULL, *tmpstart = NULL, startval[256];
      
      fakestarts = XLALStringDuplicate( 
        LALInferenceGetProcParamVal(commandLine,"--fake-starts")->value );
      
      if( (numStarts = count_csv( fakestarts )) != numDets ){
        fprintf(stderr, "Error... number of start times for fake data must be\
 equal to the number of detectors specified = %d!\n", numDets);
        exit(0);
      }
      
      tmpstarts = XLALStringDuplicate( fakestarts );
      
      for( i = 0; i < numDets; i++ ){
        tmpstart = strsep( &tmpstarts, "," );
        XLALStringCopy( startval, tmpstart, strlen(tmpstart)+1 );
        fstarts[i] = atof(startval);
      }
    }
    else /* set default GPS 900000000 - 13th July 2008 at 15:59:46 */
      for(i = 0; i < numDets; i++) fstarts[i] = 900000000.;
      
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-lengths");
    if( ppt ){
      CHAR *tmplengths = NULL, *tmplength = NULL, lengthval[256];
      
      fakelengths = XLALStringDuplicate( 
        LALInferenceGetProcParamVal(commandLine,"--fake-lengths")->value );
      
      if( (numLengths = count_csv( fakelengths )) != numDets ){
        fprintf(stderr, "Error... number of lengths for fake data must be\
 equal to the number of detectors specified = %d!\n", numDets);
        exit(0);
      }
      
      tmplengths = XLALStringDuplicate( fakelengths );
      
      for( i = 0; i < numDets; i++ ){
        tmplength = strsep( &tmplengths, "," );
        XLALStringCopy( lengthval, tmplength, strlen(tmplength)+1 );
        flengths[i] = atof(lengthval);
      }
    }
    else /* set default (86400 seconds or 1 day) */
      for(i = 0; i < numDets; i++) flengths[i] = 86400.;
      
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-dt");
    if( ppt ){
      CHAR *tmpdts = NULL, *tmpdt = NULL, dtval[256];
      
      fakedt = XLALStringDuplicate( 
        LALInferenceGetProcParamVal(commandLine,"--fake-dt")->value );
      
      if( (numDt = count_csv( fakedt )) != numDets ){
        fprintf(stderr, "Error... number of sample rates for fake data must be\
 equal to the number of detectors specified = %d!\n", numDets);
        exit(0);
      }
      
      tmpdts = XLALStringDuplicate( fakedt );
      
      for( i = 0; i < numDets; i++ ){
        tmpdt = strsep( &tmpdts, "," );
        XLALStringCopy( dtval, tmpdt, strlen(tmpdt)+1 );
        fdt[i] = atof(dtval);
      }
    }
    else /* set default (60 sesonds) */
      for(i = 0; i < numDets; i++) fdt[i] = 60.;
      
  }
  else{
    fprintf(stderr, "Error... --detectors OR --fake-data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  ppt = LALInferenceGetProcParamVal( commandLine,"--input-files" );
  if( ppt ){
    inputfile = 
      XLALStringDuplicate(
        LALInferenceGetProcParamVal(commandLine,"--input-files")->value );
  }
  
  if ( inputfile == NULL && !ppt2 ){
    fprintf(stderr, "Error... an input file or fake data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* get the output directory */
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( ppt ){
    outfile = 
      XLALStringDuplicate( LALInferenceGetProcParamVal( commandLine, "--outfile" )->value );
  }
  else{
    fprintf(stderr, "Error... --outfile needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* get ephemeris files */
  ppt = LALInferenceGetProcParamVal( commandLine, "--ephem-earth" );
  if( ppt ){
  efile = 
    XLALStringDuplicate( LALInferenceGetProcParamVal(commandLine,"--ephem-earth")->value );
  }
  ppt = LALInferenceGetProcParamVal( commandLine, "--ephem-sun" );
  if( ppt ){
    sfile = 
      XLALStringDuplicate( LALInferenceGetProcParamVal(commandLine,"--ephem-sun")->value );
  }

  /* check ephemeris files exist and if not output an error message */
  if( access(sfile, F_OK) != 0 || access(efile, F_OK) != 0 ){
    fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
    exit(3);
  }
  
  /* count the number of input files (by counting commas) and check it's equal
     to the number of detectors */
  if ( !ppt2 ){ /* if using real data */
    count = count_csv( inputfile );
    
    if ( count != numDets ){
      fprintf(stderr, "Error... Number of input files given is not equal to the\
 number of detectors given!\n");
      exit(0);
    }
  }
  
  /* set random number generator in case when that fake data is used */
  ppt = LALInferenceGetProcParamVal( commandLine, "--randomseed" );
  if ( ppt != NULL ) seed = atoi( ppt->value );
  else seed = 0; /* will be set from system clock */
      
  /* reset filestr if using real data (i.e. not fake) */
  if ( !ppt2 )
    filestr = XLALStringDuplicate( inputfile );
  
  /* read in data */
  for( i = 0, prev=NULL ; i < numDets ; i++, prev=ifodata ){
    CHAR *datafile = NULL;
    REAL8 times = 0;
    LIGOTimeGPS gpstime;
    COMPLEX16 dataVals;
    REAL8Vector *temptimes = NULL;
    INT4 j = 0, k = 0;
    
    FILE *fp = NULL;
    
    count = 0;
    
    /* initialise random number generator */
    /* Moved into det loop so same random seed can be used with */
    /* Different detector combos and still get same noise realisation */
    randomParams = XLALCreateRandomParams( seed +i );
  
    ifodata = XLALCalloc( 1, sizeof(LALInferenceIFOData) );
    ifodata->modelParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
    ifodata->modelDomain = LALINFERENCE_DOMAIN_TIME;
    ifodata->next = NULL;
    ifodata->dataParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
    
    if( i == 0 ) runState->data = ifodata;
    if( i > 0 ) prev->next = ifodata;

    /* set detector */
    ifodata->detector = XLALGetSiteInfo( dets[i] );
    
    /* set dummy initial time */
    gpstime.gpsSeconds = 0;
    gpstime.gpsNanoSeconds = 0;
    
    /* allocate time domain data */
    ifodata->compTimeData = NULL;
    ifodata->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                            1., &lalSecondUnit,
                                                            MAXLENGTH );
    
    /* allocate data time stamps */
    ifodata->dataTimes = NULL;
    ifodata->dataTimes = XLALCreateTimestampVector( MAXLENGTH );
    
    /* allocate time domain model */
    ifodata->compModelData = NULL;
    ifodata->compModelData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                            1., &lalSecondUnit,
                                                            MAXLENGTH );
    
    /*============================ GET DATA ==================================*/
    /* get i'th filename from the comma separated list */
    if ( !ppt2 ){ /* if using real data read in from the file */
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

      ifodata->compModelData = 
        XLALResizeCOMPLEX16TimeSeries( ifodata->compModelData, 0, j );
    
      /* fill in time stamps as LIGO Time GPS Vector */
      for ( k = 0; k<j; k++ )
        XLALGPSSetREAL8( &ifodata->dataTimes->data[k], temptimes->data[k] );
      
      ifodata->dataTimes->data = XLALRealloc( ifodata->dataTimes->data, 
                                              j * sizeof(LIGOTimeGPS) );
      ifodata->dataTimes->length = (UINT4)j;
      
      ifodata->compTimeData->epoch = ifodata->dataTimes->data[0];
      ifodata->compModelData->epoch = ifodata->dataTimes->data[0];
      
      XLALDestroyREAL8Vector( temptimes );
    }
    else{ /* set up fake data */
      INT4 datalength = flengths[i] / fdt[i];
      
      /* temporary real and imaginary data vectors */
      REAL4Vector *realdata = NULL;
      REAL4Vector *imagdata = NULL;
      
      REAL8 psdscale = 0.;
      
      /* resize data time stamps */
      ifodata->dataTimes->data = XLALRealloc( ifodata->dataTimes->data, 
                                              datalength*sizeof(LIGOTimeGPS) );
      ifodata->dataTimes->length = (UINT4)datalength;
      
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
}

void readDoublePulsarData( LALInferenceRunState *runState ){
  ProcessParamsTable *ppt = NULL, *ppt2 = NULL;
  ProcessParamsTable *commandLine = runState->commandLine;
  
  CHAR *detectors = NULL;
  CHAR *outfile = NULL;
  CHAR *inputfile = NULL;
  
  CHAR *filestr = NULL;
  
  CHAR *efile = NULL;
  CHAR *sfile = NULL;
  
  CHAR *tempdets=NULL;
  CHAR *tempdet=NULL;
 
  REAL8 fpsds[MAXDETS];
  CHAR *fakestarts = NULL, *fakelengths = NULL, *fakedt = NULL;
  REAL8 fstarts[MAXDETS], flengths[MAXDETS], fdt[MAXDETS];
  
  CHAR dets[MAXDETS][256];
  INT4 numDets = 0, i = 0, numPsds = 0, numLengths = 0, numStarts = 0, i2=0;
  INT4 numDt = 0, count = 0;
  UINT4 maxlen = 0;
  
  LALInferenceIFOData *ifodata = NULL;
  LALInferenceIFOData *prev = NULL;
  
  UINT4 seed = 0; /* seed for data generation */
  RandomParams *randomParams = NULL;
  
  runState->data = NULL;
  
  /* get the detectors - must */
  ppt = LALInferenceGetProcParamVal( commandLine, "--detectors" );
  ppt2 = LALInferenceGetProcParamVal( commandLine, "--fake-data" );
  if( ppt && !ppt2 ){
    detectors = XLALStringDuplicate( 
      LALInferenceGetProcParamVal(commandLine,"--detectors")->value );
      
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
    detectors = XLALStringDuplicate( 
      LALInferenceGetProcParamVal(commandLine,"--fake-data")->value );
   
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
      
      psds = XLALStringDuplicate(
        LALInferenceGetProcParamVal(commandLine,"--fake-psd")->value );
      
      tmppsds = XLALStringDuplicate( psds );
      tempdets = XLALStringDuplicate( detectors );
      
      /* count the number of PSDs (comma seperated values) to compare to number
         of detectors */
      if( (numPsds = count_csv( psds )) != 2*numDets ){
        fprintf(stderr, "Error... number of PSDs for fake data must be equal to\
 twice the number of detectors specified for model type pinsf = %d!\n", numDets);
        exit(0);
      }
      
      for( i = 0; i < 2*numDets; i+=2 ){
        CHAR *tmpstr = NULL;
        
        tmppsd = strsep( &tmppsds, "," );
        XLALStringCopy( psdval, tmppsd, strlen(tmppsd)+1 );
        fpsds[i] = atof(psdval);
	fprintf(stderr,"psdval1: %e,\t",fpsds[i]);
	tmppsd = strsep( &tmppsds, "," );
	XLALStringCopy( psdval, tmppsd, strlen(tmppsd)+1 );
	fpsds[i+1] = atof(psdval);
	fprintf(stderr,"psdval2: %e,\n",fpsds[i+1]);
	
        
        /* set detector */
        tempdet = strsep( &tempdets, "," );
        
        if( (tmpstr = strstr(tempdet, "A")) != NULL ) /* have advanced */
          XLALStringCopy( dets[i], tmpstr+1, strlen(tmpstr)+1 );
        else
          XLALStringCopy( dets[i], tempdet, strlen(tempdet)+1 );
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
      parFile = ppt->value;
        
      /* get the pulsar parameters to give a value of f */
      XLALReadTEMPOParFile( &pulsar, parFile );
      
      /* putting in pulsar frequency at f here - need to multiply by two for 2f datatstream */
      pfreq = pulsar.f0;
      
      tempdets = XLALStringDuplicate( detectors );
      
      for( i = 0; i < 2*numDets; i+=2 ){
        LALStatus status;
	
        CHAR *tmpstr = NULL;
        REAL8 psdvalf = 0., psdval2f =0.;
        
        numPsds++;
       
        tempdet = strsep( &tempdets, "," );
        
        if( (tmpstr = strstr(tempdet, "A")) != NULL ){ /* have Advanced */
          XLALStringCopy( dets[i2], tmpstr+1, strlen(tmpstr)+1 );
          
          if( !strcmp(dets[i2], "H1") || !strcmp(dets[i2], "L1") ||  
              !strcmp(dets[i2], "H2") ){ /* ALIGO */
            LALAdvLIGOPsd( &status, &psdvalf, pfreq );
	    LALAdvLIGOPsd( &status, &psdval2f, 2*pfreq );
            psdvalf *= 1.e-49; /* scale factor in LALAdvLIGOPsd.c */
	    psdval2f *= 1.e-49; /* scale factor in LALAdvLIGOPsd.c */
          }
          else if( !strcmp(dets[i2], "V1") ){ /* AVirgo */
            LALEGOPsd( &status, &psdvalf, pfreq ); 
	    LALEGOPsd( &status, &psdval2f, 2*pfreq ); 
          }
          else{
            fprintf(stderr, "Error... trying to use Advanced detector that is\
 not available!\n");
            exit(0);
          }
        }
        else{ /* initial detector */
          XLALStringCopy( dets[i2], tempdet, strlen(tempdet)+1 );
          
          if( !strcmp(dets[i2], "H1") || !strcmp(dets[i2], "L1") ||  
              !strcmp(dets[i2], "H2") ){ /* Initial LIGO */
            psdvalf = XLALLIGOIPsd( pfreq );
	    psdval2f = XLALLIGOIPsd( 2*pfreq );
          
            /* divide H2 psds by 2 */
            if( !strcmp(dets[i2], "H2") ){
		psdvalf /= 2.; 
		psdval2f /= 2.; 
	    }
          }
          else if( !strcmp(dets[i2], "V1") ){ /* Initial Virgo */
            LALVIRGOPsd( &status, &psdvalf, pfreq ); 
	    LALVIRGOPsd( &status, &psdval2f, 2*pfreq );
          }
          else if( !strcmp(dets[i2], "G1") ){ /* GEO 600 */
            LALGEOPsd( &status, &psdvalf, pfreq );
	    LALGEOPsd( &status, &psdval2f, 2*pfreq );
            psdvalf *= 1.e-46; /* scale factor in LALGEOPsd.c */
	    psdval2f *= 1.e-46; /* scale factor in LALGEOPsd.c */
          }
          else if( !strcmp(dets[i2], "T1") ){ /* TAMA300 */
            LALTAMAPsd( &status, &psdvalf, pfreq );
	    LALTAMAPsd( &status, &psdval2f, 2*pfreq );
            psdvalf *= 75. * 1.e-46; /* scale factor in LALTAMAPsd.c */
	    psdval2f *= 75. * 1.e-46; /* scale factor in LALTAMAPsd.c */
          }
          else{
            fprintf(stderr, "Error... trying to use detector that is\
 not available!\n");
            exit(0);
          }
        }
        
        fpsds[i] = psdvalf;
	fpsds[i+1] = psdval2f;
	fprintf(stderr,"psds1: %e, psd2: %e\n",psdvalf, psdval2f);
	i2++;
      }
    }
    /*generate the fake data timestamps.*/
    /*====================================================================*/
    
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-starts");
    if( ppt ){
      CHAR *tmpstarts = NULL, *tmpstart = NULL, startval[256];
      
      fakestarts = XLALStringDuplicate( 
        LALInferenceGetProcParamVal(commandLine,"--fake-starts")->value );
      
      if( (numStarts = count_csv( fakestarts )) != numDets ){
        fprintf(stderr, "Error... number of start times for fake data must be\
 equal to the number of detectors specified = %d!\n", numDets);
        exit(0);
      }
      
      tmpstarts = XLALStringDuplicate( fakestarts );
      for( i = 0; i < 2*numDets; i+=2 ){
	tmpstart = strsep( &tmpstarts, "," );
	XLALStringCopy( startval, tmpstart, strlen(tmpstart)+1 );
	
        fstarts[i] = atof(startval);
	fstarts[i+1] = atof(startval);
      }
    }
    else{ /* set default GPS 900000000 - 13th July 2008 at 15:59:46 */
      for(i = 0; i < 2*numDets; i+=2){
	  fstarts[i] = 900000000.;
	  fstarts[i+1] = 900000000.;
      }
    }
      
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-lengths");
    if( ppt ){
      CHAR *tmplengths = NULL, *tmplength = NULL, lengthval[256];
      
      fakelengths = XLALStringDuplicate( 
        LALInferenceGetProcParamVal(commandLine,"--fake-lengths")->value );
      
      if( (numLengths = count_csv( fakelengths )) != numDets ){
        fprintf(stderr, "Error... number of lengths for fake data must be\
 equal to the number of detectors specified = %d!\n", numDets);
        exit(0);
      }
      
      tmplengths = XLALStringDuplicate( fakelengths );
      for( i = 0; i < 2*numDets; i+=2 ){
        tmplength = strsep( &tmplengths, "," );
        XLALStringCopy( lengthval, tmplength, strlen(tmplength)+1 );
        flengths[i] = atof(lengthval);
	flengths[i+1] = atof(lengthval);
      }
    }
    else{ /* set default (86400 seconds or 1 day) */
      for(i = 0; i < 2*numDets; i+=2){
	  flengths[i] = 86400.;
	  flengths[i+1] = 86400.;
      }
    }
      
    ppt = LALInferenceGetProcParamVal(commandLine,"--fake-dt");
    if( ppt ){
      CHAR *tmpdts = NULL, *tmpdt = NULL, dtval[256];
      
      fakedt = XLALStringDuplicate( 
        LALInferenceGetProcParamVal(commandLine,"--fake-dt")->value );
      
      if( (numDt = count_csv( fakedt )) != numDets ){
        fprintf(stderr, "Error... number of sample rates for fake data must be\
 equal to the number of detectors specified = %d!\n", numDets);
        exit(0);
      }
      
      tmpdts = XLALStringDuplicate( fakedt );
      
      for( i = 0; i < 2*numDets; i+=2 ){
        tmpdt = strsep( &tmpdts, "," );
        XLALStringCopy( dtval, tmpdt, strlen(tmpdt)+1 );
        fdt[i] = atof(dtval);
	fdt[i+1] = atof(dtval);
      }
    }
    else{ /* set default (60 sesonds) */
      for(i = 0; i < 2*numDets; i+=2){
	fdt[i] = 60.;
	fdt[i+1] = 60.;
      }
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
  if( ppt ){
    inputfile = 
      XLALStringDuplicate(
        LALInferenceGetProcParamVal(commandLine,"--input-files")->value );
  }
  
  if ( inputfile == NULL && !ppt2 ){
    fprintf(stderr, "Error... an input file or fake data needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* get the output directory */
  ppt = LALInferenceGetProcParamVal( commandLine, "--outfile" );
  if( ppt ){
    outfile = 
      XLALStringDuplicate( LALInferenceGetProcParamVal( commandLine, "--outfile" )->value );
  }
  else{
    fprintf(stderr, "Error... --outfile needs to be set.\n");
    fprintf(stderr, USAGE, commandLine->program);
    exit(0);
  }
  
  /* get ephemeris files */
  ppt = LALInferenceGetProcParamVal( commandLine, "--ephem-earth" );
  if( ppt ){
  efile = 
    XLALStringDuplicate( LALInferenceGetProcParamVal(commandLine,"--ephem-earth")->value );
  }
  ppt = LALInferenceGetProcParamVal( commandLine, "--ephem-sun" );
  if( ppt ){
    sfile = 
      XLALStringDuplicate( LALInferenceGetProcParamVal(commandLine,"--ephem-sun")->value );
  }

  /* check ephemeris files exist and if not output an error message */
  if( access(sfile, F_OK) != 0 || access(efile, F_OK) != 0 ){
    fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
    exit(3);
  }
  
  /* count the number of input files (by counting commas) and check it's equal
     to twice the number of detectors */
  if ( !ppt2 ){ /* if using real data */
    count = count_csv( inputfile );
    
    if ( count != 2*numDets ){
      fprintf(stderr, "Error... Number of input files given is not equal to twice the\
 number of detectors given!\n");
      exit(0);
    }
  }
  
  /* set random number generator in case when that fake data is used */
  ppt = LALInferenceGetProcParamVal( commandLine, "--randomseed" );
  if ( ppt != NULL ) seed = atoi( ppt->value );
  else seed = 0; /* will be set from system clock */
      
  /* reset filestr if using real data (i.e. not fake) */
  if ( !ppt2 )
    filestr = XLALStringDuplicate( inputfile );
  
  i2=0;
  /* read in data, needs to read in two sets of data for each ifo! */
  for( i = 0, prev=NULL ; i < 2*numDets ; i++, prev=ifodata ){
    CHAR *datafile = NULL;
    REAL8 times = 0;
    LIGOTimeGPS gpstime;
    COMPLEX16 dataVals;
    REAL8Vector *temptimes = NULL;
    INT4 j = 0, k = 0;
    
    FILE *fp = NULL;
    
    count = 0;
    fprintf(stderr,"reading in data file\n");
    
    /* initialise random number generator */
    /* Moved into det loop so same random seed can be used with */
    /* Different detector combos and still get same noise realisation */
    randomParams = XLALCreateRandomParams( seed +i );
  
    ifodata = XLALCalloc( 1, sizeof(LALInferenceIFOData) );
    ifodata->modelParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
    ifodata->modelDomain = LALINFERENCE_DOMAIN_TIME;
    ifodata->next = NULL;
    ifodata->dataParams = XLALCalloc( 1, sizeof(LALInferenceVariables) );
    
    if( i == 0 ) runState->data = ifodata;
    if( i > 0 ) prev->next = ifodata;

    /* set detector */
    ifodata->detector = XLALGetSiteInfo( dets[i2] );
    
    /* set dummy initial time */
    gpstime.gpsSeconds = 0;
    gpstime.gpsNanoSeconds = 0;
    
    /* allocate time domain data */
    ifodata->compTimeData = NULL;
    ifodata->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                            1., &lalSecondUnit,
                                                            MAXLENGTH );
    
    /* allocate data time stamps */
    ifodata->dataTimes = NULL;
    ifodata->dataTimes = XLALCreateTimestampVector( MAXLENGTH );
    
    /* allocate time domain model */
    ifodata->compModelData = NULL;
    ifodata->compModelData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
                                                            1., &lalSecondUnit,
                                                            MAXLENGTH );
    
    /*============================ GET DATA ==================================*/
    /* get i'th filename from the comma separated list */
    if ( !ppt2 ){ /* if using real data read in from the file */
      /*while ( count <= i ){
        datafile = strsep(&filestr, ",");
        count++;
	fprintf(stderr,"data file count: %d\n",count);
	fprintf(stderr,"data file: %s\n",datafile);
      }*/
      datafile = strsep(&filestr, ",");
      fprintf(stderr,"data file: %s\n",datafile);
   
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

      ifodata->compModelData = 
        XLALResizeCOMPLEX16TimeSeries( ifodata->compModelData, 0, j );
    
      /* fill in time stamps as LIGO Time GPS Vector */
      for ( k = 0; k<j; k++ )
        XLALGPSSetREAL8( &ifodata->dataTimes->data[k], temptimes->data[k] );
      
      ifodata->dataTimes->data = XLALRealloc( ifodata->dataTimes->data, 
                                              j * sizeof(LIGOTimeGPS) );
      ifodata->dataTimes->length = (UINT4)j;
      
      ifodata->compTimeData->epoch = ifodata->dataTimes->data[0];
      ifodata->compModelData->epoch = ifodata->dataTimes->data[0];
      
      XLALDestroyREAL8Vector( temptimes );
    }
    else{ /* set up fake data */
      INT4 datalength = flengths[i] / fdt[i];
      
      /* temporary real and imaginary data vectors */
      REAL4Vector *realdata = NULL;
      REAL4Vector *imagdata = NULL;
      
      REAL8 psdscale = 0.;
      
      /* resize data time stamps */
      ifodata->dataTimes->data = XLALRealloc( ifodata->dataTimes->data, 
                                              datalength*sizeof(LIGOTimeGPS) );
      ifodata->dataTimes->length = (UINT4)datalength;
      
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
    
    /* set up ephemeris information */
    ifodata->ephem = XLALInitBarycenter( efile, sfile );
    XLALDestroyRandomParams( randomParams );
    
    /* get maximum data length */
    if ( ifodata->compTimeData->data->length > maxlen )
      maxlen = ifodata->compTimeData->data->length;
    
    i2+=(i%2);/*adds only onto i2 when i is odd, so it does not increment when on same det but 2nd file*/
  }
  
  /* set global variable logfactorial */
  logfactorial = XLALCalloc( maxlen+1, sizeof(REAL8) );
  for ( i = 2; i < (INT4)(maxlen+1); i++ )
    logfactorial[i] = logfactorial[i-1] + log((REAL8)i);
}


/* function to set the model of the pulsar that is being used */
void setSignalModelType( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;
  
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  
  CHAR *modeltype = NULL;
  
  ppt = LALInferenceGetProcParamVal( commandLine, "--model-type" );
  if( ppt ){
    modeltype  = XLALStringDuplicate( 
      LALInferenceGetProcParamVal( commandLine, "--model-type" )->value );
  }
  else{ /* set default model to triaxial */
    modeltype = XLALStringDuplicate( "triaxial" );
    fprintf(stderr,"Signal model set to triaxial as default\n");
  }
  
  while( data ){
    LALInferenceAddVariable( data->dataParams, "modeltype", &modeltype,
                             LALINFERENCE_string_t, LALINFERENCE_PARAM_FIXED );
    data = data->next;
  }
}


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
    UINT4 i = 0;
    LALInferenceVariableItem *scaleitem = scaletemp->head;
    REAL8Vector *dts = NULL, *bdts = NULL;
    
    if ( !varyskypos ){
      dts = get_ssb_delay( pulsar, data->dataTimes, data->ephem, data->detector,
                           0. );
    
      if ( !varybinary )
        bdts = get_bsb_delay( pulsar, data->dataTimes, dts );
    }
    
    LALInferenceAddVariable( data->dataParams, "ssb_delays", &dts,
                             LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    
    LALInferenceAddVariable( data->dataParams, "bsb_delays", &bdts,
                             LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
    
    phase_vector = get_phase_model( pulsar, data );
    
    data->timeData = NULL;
    data->timeData = XLALCreateREAL8TimeSeries( "", &data->dataTimes->data[0], 
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
                                     
  return;
}


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
    REAL8Vector *sumData = NULL;
    UINT4Vector *chunkLength = NULL;
    
    LALInferenceAddVariable( data->dataParams, "psiSteps", &psiBins, LALINFERENCE_INT4_t, 
                             LALINFERENCE_PARAM_FIXED ); 
    LALInferenceAddVariable( data->dataParams, "timeSteps", &timeBins, LALINFERENCE_INT4_t, 
                             LALINFERENCE_PARAM_FIXED );
    
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

    LALInferenceAddVariable( data->dataParams, "chunkMin", &chunkMin, LALINFERENCE_INT4_t, 
                             LALINFERENCE_PARAM_FIXED );
    LALInferenceAddVariable( data->dataParams, "chunkMax", &chunkMax, LALINFERENCE_INT4_t, 
                             LALINFERENCE_PARAM_FIXED );
    
    /* get chunk lengths of data */
    
    chunkLength = get_chunk_lengths( data, chunkMax );
    
    /* chunkLength = chop_n_merge( data, chunkMin, chunkMax ); */
    
    LALInferenceAddVariable( data->dataParams, "chunkLength", &chunkLength, 
                             LALINFERENCE_UINT4Vector_t, LALINFERENCE_PARAM_FIXED );

    /* get sum of data for each chunk */
    sumData = sum_data( data );
    
    LALInferenceAddVariable( data->dataParams, "sumData", &sumData, 
                             LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);

    data = data->next;
  }
  fprintf(stderr,"-----------------------------\n");
        
  return;
}


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
  LALInferenceAddVariable( ini, "model", &pars.model, LALINFERENCE_string_t, LALINFERENCE_PARAM_FIXED );
  
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


void add_variable_scale_prior( LALInferenceVariables *var, 
                               LALInferenceVariables *scale, 
                               LALInferenceVariables *prior, const CHAR *name, 
                               REAL8 value, REAL8 sigma ){
  REAL8 scaleVal = 1.;
  LALInferenceParamVaryType vary;
  CHAR scaleName[VARNAME_MAX] = "";
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
  LALInferenceAddVariable( scale, scaleName, &scaleVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );
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
    
    REAL8 scale = 0.;
    LALInferenceVariableType scaleType;
    CHAR tempParScale[VARNAME_MAX] = "";
    CHAR tempParPrior[VARNAME_MAX] = "";
 
    LALInferenceIFOData *datatemp = data;
    
    LALInferenceParamVaryType varyType;
    
    if( high < low ){
      fprintf(stderr, "Error... In %s the %s parameters ranges are wrongly \
set.\n", propfile, tempPar);
      exit(3);
    }
    
    sprintf(tempParScale, "%s_scale", tempPar);
    sprintf(tempParPrior, "%s_gaussian_mean", tempPar);

    tempVar = *(REAL8*)LALInferenceGetVariable( runState->currentParams, tempPar );
    type = LALInferenceGetVariableType( runState->currentParams, tempPar );
    varyType = LALInferenceGetVariableVaryType( runState->currentParams, tempPar );
    
    /* remove variable value */
    LALInferenceRemoveVariable( runState->currentParams, tempPar );
    
    /* if a Gaussian prior has already been defined (i.e. from the par file)
       remove this and overwrite with values from the propfile */
    if ( LALInferenceCheckVariable(runState->priorArgs, tempParPrior) )
      LALInferenceRemoveGaussianPrior( runState->priorArgs, tempPar );
    
    scale = high;
    
    /* set the scale factor to be the high value of the prior */
    while( datatemp ){
      scaleType = LALInferenceGetVariableType( datatemp->dataParams, tempParScale );
      LALInferenceRemoveVariable( datatemp->dataParams, tempParScale );
    
      LALInferenceAddVariable( datatemp->dataParams, tempParScale, &scale, scaleType,
        LALINFERENCE_PARAM_FIXED );
    
      datatemp = datatemp->next;
    }
      
    /* scale variable and priors */
    tempVar /= scale;
    low /= scale;
    high /= scale;
    
    /* re-add variable */
    if( !strcmp(tempPar, "phi0") ){
      LALInferenceAddVariable( runState->currentParams, tempPar, &tempVar, type,
                   LALINFERENCE_PARAM_CIRCULAR );
    }
    else{
      LALInferenceAddVariable( runState->currentParams, tempPar, &tempVar, type,
                   LALINFERENCE_PARAM_LINEAR );
    }
    /* Add the prior variables */
    LALInferenceAddMinMaxPrior( runState->priorArgs, tempPar, (void *)&low, 
                                (void *)&high, type );
    
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
    REAL8 scale = 0.;
    REAL8 tempVar;
    
    REAL8 mu, sigma;
    
    LALInferenceIFOData *datatemp = data;
    
    LALInferenceVariableType scaleType;
    CHAR tempParScale[VARNAME_MAX] = "";
    CHAR tempParPrior[VARNAME_MAX] = "";

    sprintf(tempParPrior,"%s_gaussian_mean",checkPrior->name);

    if ( LALInferenceCheckVariable(runState->priorArgs, tempParPrior) ){
      tempVar = *(REAL8 *)checkPrior->value;
      
      /* get the mean and standard deviation of the Gaussian prior */
      LALInferenceGetGaussianPrior( runState->priorArgs, checkPrior->name, (void *)&mu,
                        (void *)&sigma );
      
      /* set the scale factor to be the mean value */
      scale = mu;
      tempVar /= scale;
      
      /* scale the parameter value and reset it */
      memcpy( checkPrior->value, &tempVar, LALInferenceTypeSize[checkPrior->type] );
      
      mu /= scale;
      sigma /= scale;
      
      /* remove the Gaussian prior values and reset as scaled values */
      LALInferenceRemoveGaussianPrior( runState->priorArgs, checkPrior->name );
      LALInferenceAddGaussianPrior( runState->priorArgs, checkPrior->name, 
                        (void *)&mu, (void *)&sigma, checkPrior->type );
      
      sprintf(tempParScale, "%s_scale", checkPrior->name);
        
      /* set scale factor in data structure */
      while( datatemp ){
        scaleType = LALInferenceGetVariableType( datatemp->dataParams, tempParScale );
        LALInferenceRemoveVariable( datatemp->dataParams, tempParScale );
    
        LALInferenceAddVariable( datatemp->dataParams, tempParScale, &scale, scaleType,
                     LALINFERENCE_PARAM_FIXED );
      
        datatemp = datatemp->next;
      }
    }
  }
    
  return;
}


void setupLivePointsArray( LALInferenceRunState *runState ){
/* Set up initial basket of live points, drawn from prior,
   by copying runState->currentParams to all entries in the array*/
  UINT4 Nlive = (UINT4)*(INT4 *)LALInferenceGetVariable(runState->algorithmParams,"Nlive");
  UINT4 i;
  REAL8Vector *logLs;
        
  LALInferenceVariableItem *current;
  
  /* Allocate the array */
  runState->livePoints = XLALCalloc( Nlive, sizeof(LALInferenceVariables *) );
  
  if( runState->livePoints == NULL ){
    fprintf(stderr,"Unable to allocate memory for %i live points\n",Nlive);
    exit(1);
  }

  logLs = XLALCreateREAL8Vector( Nlive );
     
  LALInferenceAddVariable( runState->algorithmParams, "logLikelihoods",
               &logLs, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
               
  fprintf(stdout, "Sprinkling %i live points, may take some time\n", Nlive);
  
  for( i=0; i<Nlive; i++){
    runState->livePoints[i] = XLALCalloc( 1, sizeof(LALInferenceVariables) );
                
    /* Copy the param structure */
    LALInferenceCopyVariables( runState->currentParams, runState->livePoints[i] );
    
    /* Sprinkle the varying points among prior */
    do{
      for( current=runState->livePoints[i]->head; current!=NULL;
           current=current->next){
        CHAR tempParPrior[VARNAME_MAX] = "";
        UINT4 gp = 0;
      
        sprintf(tempParPrior,"%s_gaussian_mean",current->name);
      
        if( LALInferenceCheckVariable( runState->priorArgs, tempParPrior ) ) gp = 1;
        
        if( current->vary==LALINFERENCE_PARAM_CIRCULAR || current->vary==LALINFERENCE_PARAM_LINEAR )
        {
          switch (current->type){
            case LALINFERENCE_REAL4_t:
            {
              REAL4 tmp;
              REAL4 min, max, mu, sigma;
                                                       
              if( gp ){
                LALInferenceGetGaussianPrior( runState->priorArgs, current->name, 
                                  (void *)&mu, (void *)&sigma );
                tmp = mu + gsl_ran_gaussian(runState->GSLrandom, (double)sigma);
              }
              else{
                LALInferenceGetMinMaxPrior( runState->priorArgs, current->name, 
                                (void *)&min, (void *)&max );
                tmp = min + (max-min)*gsl_rng_uniform( runState->GSLrandom );
              }
                                                       
              LALInferenceSetVariable( runState->livePoints[i], current->name, &tmp );
              break;
            }
            case LALINFERENCE_REAL8_t:
            {
              REAL8 tmp;
              REAL8 min, max, mu, sigma;
                                                       
              if( gp ){
                LALInferenceGetGaussianPrior( runState->priorArgs, current->name, 
                                  (void *)&mu, (void *)&sigma );
                tmp = mu + gsl_ran_gaussian(runState->GSLrandom, (double)sigma);
              }
              else{
                LALInferenceGetMinMaxPrior( runState->priorArgs, current->name, 
                                (void *)&min, (void *)&max );
                tmp = min + (max-min)*gsl_rng_uniform( runState->GSLrandom );
              }
                                                       
              LALInferenceSetVariable( runState->livePoints[i], current->name, &tmp );
              break;
            }
            case LALINFERENCE_INT4_t:
            {
              INT4 tmp;
              INT4 min,max;
                                                       
              LALInferenceGetMinMaxPrior( runState->priorArgs, current->name, (void *)&min,
                              (void *)&max );
                                                       
              tmp = min + (max-min)*gsl_rng_uniform(runState->GSLrandom);
                                                       
              LALInferenceSetVariable( runState->livePoints[i], current->name, &tmp );
              break;
            }
            case LALINFERENCE_INT8_t:
            {
              INT8 tmp;
              INT8 min, max;
                                                       
              LALInferenceGetMinMaxPrior( runState->priorArgs, current->name, (void *)&min,
                              (void *)&max );
                                                       
              tmp = min + (max-min)*gsl_rng_uniform(runState->GSLrandom);
                                                       
              LALInferenceSetVariable( runState->livePoints[i], current->name, &tmp);
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
                                           runState->data, runState->template );
  }
        
}

/*------------------- END INITIALISATION FUNCTIONS ---------------------------*/


/******************************************************************************/
/*                     LIKELIHOOD AND PRIOR FUNCTIONS                         */
/******************************************************************************/

REAL8 pulsar_log_likelihood( LALInferenceVariables *vars, LALInferenceIFOData *data,
                             LALInferenceTemplateFunction *get_model ){
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
  
    REAL8Vector *sumData = NULL;
    UINT4Vector *chunkLengths = NULL;

    sumData = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "sumData" );
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams,
                                                 "chunkLength" );
    chunkMin = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMin" );
    chunkMax = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMax" );
  
    /* copy model parameters to data parameters */
    LALInferenceCopyVariables( vars, data->modelParams );
      
    /* get pulsar model */
    get_model( data );
  
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
      }
 
      chiSquare = sumData->data[count];
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

REAL8 pulsar_double_log_likelihood( LALInferenceVariables *vars, LALInferenceIFOData *data,
                             LALInferenceTemplateFunction *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  UINT4 i = 0, counter=0;
  
  while ( data ){
    counter++;
    /*fprintf(stderr,"Log like calc counter: %d\n",counter);*/
    UINT4 j = 0, count = 0, cl = 0;
    UINT4 length = 0, chunkMin, chunkMax;
    REAL8 chunkLength = 0.;
    REAL8 logliketmp = 0.;

    REAL8 sumModel = 0., sumDataModel = 0.;
    REAL8 sumModel2 = 0., sumDataModel2 = 0.;
    REAL8 sumdata=0., sumdata2=0.;
    REAL8 chiSquare = 0., chiSquare2 = 0.;
    REAL8 test;
    COMPLEX16 B, M;

    INT4 first = 0, through = 0;
  
    REAL8Vector *sumData = NULL, *sumData2 = NULL;
    UINT4Vector *chunkLengths = NULL;

    sumData = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams, "sumData" );
    sumData2 = *(REAL8Vector **)LALInferenceGetVariable( data->next->dataParams, "sumData" );
    
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams,
                                                 "chunkLength" );
    chunkMin = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMin" );
    chunkMax = *(INT4*)LALInferenceGetVariable( data->dataParams, "chunkMax" );
  
    /* copy model parameters to data parameters */
    LALInferenceCopyVariables( vars, data->modelParams );
  
    /* get pulsar model */
    get_model( data );
  
    length = data->compTimeData->data->length;
  
    for( i = 0 ; i < length ; i += chunkLength ){
      chunkLength = (REAL8)chunkLengths->data[count];
      /*fprintf(stderr,"chunkLength: %e\n", chunkLength);*/
    
      /* skip section of data if its length is less than the minimum allowed
        chunk length */
      if( chunkLength < chunkMin ){
        count++;
        if( through == 0 ) first = 0;
	fprintf(stderr,"Chunk length less than allowed\n");
        continue;
      }
      through = 1;

      sumdata=0.;
      sumdata=0.;
      sumModel = 0.;
      sumDataModel = 0.;
      sumModel2 = 0.;
      sumDataModel2 = 0.;

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
	sumdata+= (B.re-M.re)*(B.re-M.re) + (B.im-M.im)*(B.im-M.im);
	
	B.re = data->next->compTimeData->data->data[j].re;
        B.im = data->next->compTimeData->data->data[j].im;

        M.re = data->next->compModelData->data->data[j].re;
        M.im = data->next->compModelData->data->data[j].im;
	
	/* sum over the model */
        sumModel2 += M.re*M.re + M.im*M.im;
        /* sum over that data and model */
        sumDataModel2 += B.re*M.re + B.im*M.im;
	sumdata2+= (B.re-M.re)*(B.re-M.re) + (B.im-M.im)*(B.im-M.im);
      }
 
      chiSquare = sumdata;
      chiSquare -= 2.*sumDataModel;
      chiSquare += sumModel;
      chiSquare2 = sumdata2;
      chiSquare2 -= 2.*sumDataModel2;
      chiSquare2 += sumModel2;
      
      if( first == 0 ){
        logliketmp = 0.;
        first++;
      }
      /*fprintf(stderr,"logliketmp: %e, chunkLength: %e, chisquare: %e\n", logliketmp, chunkLength, chiSquare);*/
      test=(chiSquare+chiSquare2);	
      if(test<=0){
	  /*fprintf(stderr,"NAN warning, sumData2: %e, sumData %e, sumModel2: %e, sumDataModel2: %e\n", sumData2->data[count],sumData->data[count], sumModel2, 2*sumDataModel2);*/
	  /*fprintf(stderr,"B.re: %e, B.im: %e, M.re: %e, M.im: %e\n", B.re, B.im, M.re, M.im);*/
	  /*fprintf(stderr,"----------------------------------------------------------------\n");*/
      }
      logliketmp -= chunkLength*log(chiSquare);
      logliketmp -= chunkLength*log(chiSquare2);
      if(isnan(logliketmp)){
          /*fprintf(stderr,"nan....logliketmp: %e, chunkLength: %e, chisquare: %e, chisquare2: %e\n", logliketmp, chunkLength, chiSquare, chiSquare2);
	  fprintf(stderr,"model: %e, model2: %e, data: %e, data2: %e\n", sumDataModel,
	  sumDataModel2, sumModel, sumModel2);*/
	  /*fprintf(stderr,"sumData: %e, check: %e\n",sumData->data[count], sumdata);
	  fprintf(stderr,"----------------------------------\n");*/
      }
    
      count++;
    }
  
    loglike += logliketmp;
  
    data = data->next;
    data = data->next;
  }
  return loglike;
}

REAL8 priorFunction( LALInferenceRunState *runState, LALInferenceVariables *params ){
  LALInferenceIFOData *data = runState->data;
  (void)runState;
  LALInferenceVariableItem *item = params->head;
  REAL8 min, max, mu, sigma, prior = 0, value = 0.;
  
  for(; item; item = item->next ){
    /* get scale factor */
    CHAR scalePar[VARNAME_MAX] = "";
    CHAR priorPar[VARNAME_MAX] = "";
    REAL8 scale;
    
    if( item->vary == LALINFERENCE_PARAM_FIXED || item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }
    
    sprintf(scalePar, "%s_scale", item->name);
    scale = *(REAL8 *)LALInferenceGetVariable( data->dataParams, scalePar );
    
    if( item->vary == LALINFERENCE_PARAM_LINEAR || item->vary == LALINFERENCE_PARAM_CIRCULAR ){
      sprintf(priorPar, "%s_gaussian_mean", item->name);
      /* Check for a gaussian */
      if ( LALInferenceCheckVariable(runState->priorArgs, priorPar) ){
        LALInferenceGetGaussianPrior( runState->priorArgs, item->name, (void *)&mu, 
                          (void *)&sigma );
      
       value = (*(REAL8 *)item->value) * scale;
       mu *= scale;
       sigma *= scale;
       prior -= log(sqrt(2.*LAL_PI)*sigma);
       prior -= (value - mu)*(value - mu) / (2.*sigma*sigma);
      }
      /* Otherwise use a flat prior */
      else{
	LALInferenceGetMinMaxPrior( runState->priorArgs, item->name, (void *)&min, 
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

void get_pulsar_model( LALInferenceIFOData *data ){
  BinaryPulsarParams pars;
  
  REAL8 rescale = 1.;
  
  CHAR *modeltype = NULL;
  
  /* set model parameters (including rescaling) */
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "h0_scale" );
  pars.h0 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "h0" ) *
    rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, 
                                              "cosiota_scale" );
  pars.cosiota = *(REAL8*)LALInferenceGetVariable( data->modelParams, 
                                                   "cosiota" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "psi_scale" );
  pars.psi = *(REAL8*)LALInferenceGetVariable( data->modelParams, "psi" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "phi0_scale" );
  pars.phi0 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "phi0" ) * rescale;
  /*pinned superf;uid parameters*/
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "h1_scale" );
  pars.h1 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "h1" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "lambda_scale" );
  pars.lambda = *(REAL8*)LALInferenceGetVariable( data->modelParams, "lambda" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "theta_scale" );
  pars.theta = *(REAL8*)LALInferenceGetVariable( data->modelParams, "theta" ) * rescale;
  
  /* set the potentially variable parameters */
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "pepoch_scale" );
  pars.pepoch = *(REAL8*)LALInferenceGetVariable( data->modelParams, "pepoch" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "posepoch_scale" );
  pars.posepoch = *(REAL8*)LALInferenceGetVariable( data->modelParams, "posepoch" ) *
    rescale;
  
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "ra_scale" );
  pars.ra = *(REAL8*)LALInferenceGetVariable( data->modelParams, "ra" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "pmra_scale" );
  pars.pmra = *(REAL8*)LALInferenceGetVariable( data->modelParams, "pmra" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "dec_scale" );
  pars.dec = *(REAL8*)LALInferenceGetVariable( data->modelParams, "dec" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "pmdec_scale" );
  pars.pmdec = *(REAL8*)LALInferenceGetVariable( data->modelParams, "pmdec" ) * rescale;
  
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "f0_scale" );
  pars.f0 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "f0" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "f1_scale" );
  pars.f1 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "f1" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "f2_scale" );
  pars.f2 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "f2" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "f3_scale" );
  pars.f3 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "f3" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "f4_scale" );
  pars.f4 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "f4" ) * rescale;
  rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "f5_scale" );
  pars.f5 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "f5" ) * rescale;
  
  pars.model = *(CHAR**)LALInferenceGetVariable( data->modelParams, "model" );

  /* binary parameters */
  if( pars.model != NULL ){
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "e_scale" );
    pars.e = *(REAL8*)LALInferenceGetVariable( data->modelParams, "e" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "w0_scale" );
    pars.w0 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "w0" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "Pb_scale" );
    pars.Pb = *(REAL8*)LALInferenceGetVariable( data->modelParams, "Pb" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "x_scale" );
    pars.x = *(REAL8*)LALInferenceGetVariable( data->modelParams, "x" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "T0_scale" );
    pars.T0 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "T0" ) * rescale;
    
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "e2_scale" );
    pars.e2 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "e2" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "w02_scale" );
    pars.w02 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "w02" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "Pb2_scale" );
    pars.Pb2 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "Pb2" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "x2_scale" );
    pars.x2 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "x2" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "T02_scale" );
    pars.T02 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "T02" ) * rescale;
    
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "e3_scale" );
    pars.e3 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "e3" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "w03_scale" );
    pars.w03 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "w03" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "Pb3_scale" );
    pars.Pb3 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "Pb3" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "x3_scale" );
    pars.x3 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "x3" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "T03_scale" );
    pars.T03 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "T03" ) * rescale;
    
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "xpbdot_scale" );
    pars.xpbdot = *(REAL8*)LALInferenceGetVariable( data->modelParams, "xpbdot" ) * rescale;
    
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "eps1_scale" );
    pars.eps1 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "eps1" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "eps2_scale" );
    pars.eps2 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "eps2" ) * rescale;   
     rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "eps1dot_scale" );
    pars.eps1dot = *(REAL8*)LALInferenceGetVariable( data->modelParams, "eps1dot" ) *
      rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "eps2dot_scale" );
    pars.eps2dot = *(REAL8*)LALInferenceGetVariable( data->modelParams, "eps2dot" ) *
      rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "Tasc_scale" );
    pars.Tasc = *(REAL8*)LALInferenceGetVariable( data->modelParams, "Tasc" ) * rescale;
    
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "wdot_scale" );
    pars.wdot = *(REAL8*)LALInferenceGetVariable( data->modelParams, "wdot" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "gamma_scale" );
    pars.gamma = *(REAL8*)LALInferenceGetVariable( data->modelParams, "gamma" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "Pbdot_scale" );
    pars.Pbdot = *(REAL8*)LALInferenceGetVariable( data->modelParams, "Pbdot" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "xdot_scale" );
    pars.xdot = *(REAL8*)LALInferenceGetVariable( data->modelParams, "xdot" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "edot_scale" );
    pars.edot = *(REAL8*)LALInferenceGetVariable( data->modelParams, "edot" ) * rescale;
    
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "s_scale" );
    pars.s = *(REAL8*)LALInferenceGetVariable( data->modelParams, "s" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "dr_scale" );
    pars.dr = *(REAL8*)LALInferenceGetVariable( data->modelParams, "dr" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "dth_scale" );
    pars.dth = *(REAL8*)LALInferenceGetVariable( data->modelParams, "dth" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "a0_scale" );
    pars.a0 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "a0" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "b0_scale" );
    pars.b0 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "b0" ) * rescale; 

    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "M_scale" );
    pars.M = *(REAL8*)LALInferenceGetVariable( data->modelParams, "M" ) * rescale;
    rescale = *(REAL8*)LALInferenceGetVariable( data->dataParams, "m2_scale" );
    pars.m2 = *(REAL8*)LALInferenceGetVariable( data->modelParams, "m2" ) * rescale;
  }

  modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );

  /* model specific for a triaxial pulsar emitting at twice the rotation
     frequency - other models can be added later */
  if ( !strcmp( modeltype, "triaxial" ) ){
    get_triaxial_pulsar_model( pars, data );
  }
  else if ( !strcmp( modeltype, "pinsf" ) ){
    get_pinsf_pulsar_model( pars, data );
  }
  else{
    fprintf(stderr, "Error... model '%s' is not defined!\n", modeltype);
    exit(0);
  }
    
}


void get_triaxial_pulsar_model( BinaryPulsarParams params, 
                                LALInferenceIFOData *data ){
  REAL8Vector *dphi = NULL;
  INT4 i = 0, length = 0;
  
  get_amplitude_model( params, data );
  length = data->compModelData->data->length;
  
  /* the timeData vector within the LALIFOData structure contains the
     phase calculated using the initial (heterodyne) values of the phase
     parameters */
  
  /* get difference in phase and perform extra heterodyne with it */
  if ( varyphase ){
    if ( (dphi = get_phase_model( params, data )) != NULL ){
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

void get_pinsf_pulsar_model( BinaryPulsarParams params, 
                                LALInferenceIFOData *data ){
  REAL8Vector *dphi = NULL;
  INT4 i = 0, length = 0;
  
  get_pinsf_amplitude_model( params, data );
  length = data->compModelData->data->length;
  
  /* the timeData vector within the LALIFOData structure contains the
     phase calculated using the initial (heterodyne) values of the phase
     parameters */
  
  /* get difference in phase and perform extra heterodyne with it */ 
  if ( varyphase ){
    if ( (dphi = get_phase_model( params, data )) != NULL ){
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


REAL8Vector *get_phase_model( BinaryPulsarParams params, 
                              LALInferenceIFOData *data ){
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
      "ssb_delays" )) == NULL ){
    /* get time delays with an interpolation of interptime (30 mins) */
    dts = get_ssb_delay( params, data->dataTimes, data->ephem, data->detector,
                         interptime );
  }
  
  if( (bdts = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams,
      "bsb_delays" )) == NULL ){
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
    phis->data[i] = 2.*deltat*(params.f0 + 
      inv_fact[2]*params.f1*deltat +
      inv_fact[3]*params.f2*deltat2 +
      inv_fact[4]*params.f3*deltat*deltat2 +
      inv_fact[5]*params.f4*deltat2*deltat2 +
      inv_fact[6]*params.f5*deltat2*deltat2*deltat);
  }

  return phis;
}


/* function to get the solar system barycentring time delay */
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


/* function to get the binary system barycentring time delay */
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


void get_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data ){
    
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
  psteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "timeSteps" );
  
  LU_Fplus = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fplus");
  LU_Fcross = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fcross");
  
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

void get_pinsf_amplitude_model( BinaryPulsarParams pars, LALInferenceIFOData *data ){
  INT4 i = 0, length;
  
  REAL8 psteps, tsteps;
  INT4 psibin, timebin;
  REAL8 tstart;
  REAL8 plus, cross;
  REAL8 T;
  REAL8 Xplusf, Xcrossf, Xplus2f, Xcross2f;
  REAL8 A1, A2, B1, B2;
  REAL4 sinphi, cosphi, sin2phi, cos2phi;
  
  gsl_matrix *LU_Fplus, *LU_Fcross;
  
  length = data->dataTimes->length;
  
  /* set lookup table parameters */
  psteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "psiSteps" );
  tsteps = *(INT4*)LALInferenceGetVariable( data->dataParams, "timeSteps" );
  
  LU_Fplus = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fplus");
  LU_Fcross = *(gsl_matrix**)LALInferenceGetVariable( data->dataParams, "LU_Fcross");
  
  sin_cos_LUT( &sinphi, &cosphi, pars.phi0 );
  sin_cos_LUT( &sin2phi, &cos2phi, 2*pars.phi0 );
  
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
 
  tstart = XLALGPSGetREAL8( &data->dataTimes->data[0] ); /*time of first B_k*/
  
  /*fprintf(stderr,"lambda: %f, theta: %f, h1: %f\n",pars.lambda, pars.theta, pars.h1);*/
  
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
    /*at f*/
    data->compModelData->data->data[i].re = plus*Xplusf*((A1*cosphi)-(A2*sinphi)) + ( cross*Xcrossf*((A2*cosphi)-(A1*sinphi)) );
    
    data->compModelData->data->data[i].im = plus*Xplusf*((A2*cosphi)+(A1*sinphi)) + ( cross*Xcrossf*((A2*sinphi)-(A1*cosphi)) );
    
    /*at 2f*/
    data->next->compModelData->data->data[i].re = plus*Xplus2f*((B1*cos2phi)-(B2*sin2phi)) + cross*Xcross2f*((B2*cos2phi)+(B1*sin2phi));
    
    data->next->compModelData->data->data[i].im = plus*Xplus2f*((B2*cos2phi)+(B1*sin2phi)) - cross*Xcross2f*((B1*cos2phi)-(B2*sin2phi));
    
  }
  
}


/* calculate the likelihood for the data just being noise - we still have a
students-t likelihood, so this basically involves taking that likelihood and
setting the signal to zero */
REAL8 noise_only_model( LALInferenceIFOData *data ){
  LALInferenceIFOData *datatemp = data;
  
  REAL8 logL = 0.0;
  UINT4 i = 0;
  
  while ( datatemp ){
    UINT4Vector *chunkLengths = NULL;
    REAL8Vector *sumData = NULL;
  
    REAL8 chunkLength = 0.;
  
    chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams, 
                                                             "chunkLength" );
    sumData = *(REAL8Vector **)LALInferenceGetVariable( data->dataParams,
                                                        "sumData" );
  
    for (i=0; i<chunkLengths->length; i++){
      chunkLength = (REAL8)chunkLengths->data[i];
   
      logL -= chunkLength * log(sumData->data[i]);
    }
  
    datatemp = datatemp->next;
  }
  
  return logL;
}

/*------------------------ END OF MODEL FUNCTIONS ----------------------------*/


/******************************************************************************/
/*                       SOFTWARE INJECTION FUNCTIONS                         */
/******************************************************************************/


void injectSignal( LALInferenceRunState *runState ){
  LALInferenceIFOData *data = runState->data;
  
  CHAR *injectfile = NULL;
  ProcessParamsTable *ppt;
  ProcessParamsTable *commandLine = runState->commandLine;
  CHAR *modeltype = NULL;
  
  BinaryPulsarParams injpars;
  
  ppt = LALInferenceGetProcParamVal( commandLine, "--inject-file" );
  if( ppt ){

    injectfile = XLALStringDuplicate( 
      LALInferenceGetProcParamVal( commandLine, "--inject-file" )->value );
    
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
  
  /* create signal to inject and add to the data */
  while ( data ){
    FILE *fp = NULL, *fpso = NULL;
    ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal( commandLine,
                                                            "--inject-output" );
    
    /* check whether to output the data */
    if ( ppt2 ){
      /* add the site prefix to the start of the output name */
      CHAR *outfile = NULL;
      CHAR *signalonly = NULL; /* file containing only signal and no noise */
      
      outfile = XLALStringDuplicate( data->detector->frDetector.prefix );
      
      outfile = XLALStringAppend( outfile, LALInferenceGetProcParamVal(
        commandLine, "--inject-output" )->value );
      
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
                                                   
    INT4 i = 0, length = data->dataTimes->length;
   
    /* for injection always attempt to include the signal phase model even if
       the search is not going to be over phase */
    UINT4 varyphasetmp = varyphase;
    varyphase = 1;
                                                            
    /* create the signal */
    modeltype = *(CHAR**)LALInferenceGetVariable( data->dataParams, "modeltype" );
    if ( !strcmp( modeltype, "triaxial" ) ){
	fprintf(stderr,"injecting triaxial signal into data\n");
      get_triaxial_pulsar_model( injpars, data );
    }
    else if ( !strcmp( modeltype, "pinsf" ) ){
      fprintf(stderr,"injecting signal into data\n");
      get_pinsf_pulsar_model( injpars, data );

      for ( i = 0; i < length; i++ ){
	data->compTimeData->data->data[i].re +=
	    data->compModelData->data->data[i].re;
	data->compTimeData->data->data[i].im +=
	    data->compModelData->data->data[i].im;
      }
      data = data->next; 
    }
    else{
	exit(0);
	fprintf(stderr,"Model type incorrectly specified\n");
    }
    /*get_triaxial_pulsar_model( injpars, data );*/
    
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
    
    /* reset varyphase to its original value */
    varyphase = varyphasetmp;
    
    data = data->next; 
  }
}

/*-------------------- END OF SOFTWARE INJECTION FUNCTIONS -------------------*/


/******************************************************************************/
/*                            HELPER FUNCTIONS                                */
/******************************************************************************/

/* function to get the lengths of consecutive chunks of data */
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
  
  merge_data( meddata, chunkIndex );
  
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


/* function to find change points and chop up the data */
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
     T = 1.305*log10(length) + 0.254 + 2.531
     where the first two terms come from a fit to odds ratios for a Monte Carlo
     of Gaussian noise (with real and imaginary components) of various lengths,
     and the final term comes from an offset to give the 1% false alarm rate. */
  threshold = 1.305*log10(length) + 0.254 + 2.531;
  
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


/* this function will go through the data and compare the "evidence" for it
   being Gaussian with a single standard deviation, or it being two Gaussian's
   with independent standard deviations. It will return the index of the index
   of the point at with the data statistics change, and the log odds ratio of
   the two hypotheses at that point. */ 
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


/* a function to sum over the data */
REAL8Vector *sum_data( LALInferenceIFOData *data ){
  INT4 chunkLength = 0, length = 0, i = 0, j = 0, count = 0;
  COMPLEX16 B;
  REAL8Vector *sumData = NULL; 

  UINT4Vector *chunkLengths;
  
  chunkLengths = *(UINT4Vector **)LALInferenceGetVariable( data->dataParams, 
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


void rescaleOutput( LALInferenceRunState *runState ){
  /* Open original output output file */
  CHAR *outfile, outfiletmp[256] = "", outfilepars[256] = "";
  FILE *fp = NULL, *fptemp = NULL, *fppars = NULL;
  UINT4 j = 0;
  
  LALInferenceVariables *current=XLALCalloc(1,sizeof(LALInferenceVariables));
  
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
 
  /* open file for printing out list of parameter names */
  sprintf(outfilepars, "%s_parnames.txt", outfile);
  if( (fppars = fopen(outfilepars, "w")) == NULL ){
    fprintf(stderr, "Error... cannot open parameter name output file %s.\n",
            outfilepars);
    exit(3);
  }
  
  /* copy variables from runState to current (this seems to switch the order,
     which is required! */
  LALInferenceCopyVariables(runState->currentParams, current);
  
  do{
    LALInferenceVariableItem *item = current->head;
    
    CHAR line[1000];
    CHAR value[128] = "";
    UINT4 i = 0;
    
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
        scalefac = *(REAL8 *)LALInferenceGetVariable( runState->data->dataParams, 
                                          scalename );
      }
      
      switch (item->type) {
        case LALINFERENCE_INT4_t:
          fprintf(fptemp, "%d", (INT4)(atoi(value)*scalefac));
          if ( j == 0 ) fprintf(fppars, "%s\n", item->name);
          break;
        case LALINFERENCE_REAL4_t:
          fprintf(fptemp, "%e", atof(value)*scalefac);
          if ( j == 0 ) fprintf(fppars, "%s\n", item->name);
          break;
        case LALINFERENCE_REAL8_t:
          fprintf(fptemp, "%le", atof(value)*scalefac);
          if ( j == 0 ) fprintf(fppars, "%s\n", item->name);
          break;
        case LALINFERENCE_string_t:
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
    
    if ( j == 0 ) fprintf(fppars, "logLikelihood\n");
    
    j++;
  }while( !feof(fp) );
  
  fclose(fp);
  fclose(fptemp);
  fclose(fppars);
  
  /* move the temporary file name to the standard outfile name */
  rename( outfiletmp, outfile );
  
  return;
}


/* function to count number of comma separated values in a string */
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


/* function to check whether a given parameter is one of the recognised
   parameters listed in the header file */
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

/*----------------------- END OF HELPER FUNCTIONS ----------------------------*/

/******************************************************************************/
/*                          TESTING FUNCTIONS                                 */
/******************************************************************************/

/* this is a testing function to output the posterior of 1 parameter as
   calculated on a specified 1D grid, with all other values held fixed at there
   values from the input par file */
void gridOutput( LALInferenceRunState *runState ){
  REAL8 h0min = 0.;
  REAL8 h0max = 0.;
  REAL8 h0range = 0, h0step = 0;
  INT4 h0steps = 0, i = 0;
 
  ProcessParamsTable *ppt;
  REAL8 scaleval = 1., tmpscale = 0., tmpgridval = 0.;
  
  ProcessParamsTable *commandLine = runState->commandLine;
  
  FILE *fp = NULL;
  REAL8 minL = LAL_REAL8_MAX;
  REAL8 sumPost = 0.;
  
  REAL8Vector *logL = NULL;
  
  CHAR *parname = NULL, parscale[256], outputgrid[256];
  
  /*------------------------------------------------------------*/
  /* test output on a h0 grid */
  ppt = LALInferenceGetProcParamVal( commandLine, "--grid" );
  if ( ppt ){
    ProcessParamsTable *ppt2;
    
    /* parameters over which to perform the grid search */
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridpar" );
    
    if( ppt2 ){
      parname = XLALStringDuplicate( LALInferenceGetProcParamVal( commandLine,
        "--gridpar" )->value );
        
      if( !recognised_parameter( parname ) ){
        fprintf(stderr, "Error... parameter %s not recognised\n", parname );
        exit(0);
      }
        
      sprintf(parscale, "%s_scale", parname);
    }
    else{
      fprintf(stderr, USAGEGRID, commandLine->program);
      exit(0);
    }
    
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridmin" );
    
    if( ppt2 ){
      h0min = atof( LALInferenceGetProcParamVal(commandLine, 
                                                "--gridmin")->value );
    }
    else h0min = 0.; /* default to zero */
    
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridmax" );
    
    if( ppt2 ){
      h0max = atof( LALInferenceGetProcParamVal(commandLine, 
                                                "--gridmax")->value );
    }
    else h0max = 1.; /* default to 1 */
    
    ppt2 = LALInferenceGetProcParamVal( commandLine, "--gridsteps" );
    
    if( ppt2 ){
      h0steps = atoi( LALInferenceGetProcParamVal(commandLine, 
                                                "--gridsteps")->value );
    }
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
  LALInferenceRemoveVariable( runState->data->dataParams, parscale );
  LALInferenceAddVariable( runState->data->dataParams, parscale, &scaleval,
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
  LALInferenceSetVariable( runState->currentParams, parname, &tmpgridval );
}

REAL8 test_gaussian_log_likelihood( LALInferenceVariables *vars,
                                    LALInferenceIFOData *data,
                                    LALInferenceTemplateFunction *get_model ){
  REAL8 loglike = 0.; /* the log likelihood */
  
  REAL8 like_mean = 0.5;
  REAL8 like_sigma = 0.025;
  REAL8 h0 = *(REAL8 *)LALInferenceGetVariable( vars, "h0" );
  REAL8 h0scale = *(REAL8 *)LALInferenceGetVariable( data->dataParams,
                                                     "h0_scale" );
  
  h0 *= h0scale;
                                                     
  /* search over a simple 1D Gaussian with x defined by the h0 variable */
  loglike = -log(sqrt(2.*LAL_PI)*like_sigma);
  loglike -= (h0-like_mean)*(h0-like_mean) / (2.*like_sigma*like_sigma);
  
  return loglike;
}


/*----------------------- END OF TESTING FUNCTIONS ---------------------------*/
