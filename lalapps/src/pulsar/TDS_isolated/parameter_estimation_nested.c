/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

#include "pulsar_parameter_estimation.h"
#include "LALInference.h"


#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0),
(1.0/24.0), (1.0/120.0), (1.0/720.0) };


#include "pulsar_parameter_estimation.h"

RCSID("$Id$");

/* global variable */
INT4 verbose=0;


/* Usage format string */
#define USAGE1 \
"Usage: %s [options]\n\n"\
" --help              display this message\n"\
" --verbose           display all error messages\n"\
" --detectors         all IFOs with data to be analysed e.g. H1,H2\n\
                     (delimited by commas)\n"\
" --pulsar            name of pulsar e.g. J0534+2200\n"\
" --par-file          pulsar parameter (.par) file (full path) \n"\
" --input-dir         directory containing the input data files\n"\
" --output-dir        directory for output data files\n"\
" --chunk-min         (INT4) minimum stationary length of data to be used in\n\
                     the likelihood e.g. 5 mins\n"\
" --chunk-max         (INT4) maximum stationary length of data to be used in\n\
                     the likelihood e.g. 30 mins\n"\
" --psi-bins          (INT4) no. of psi bins in the time-psi lookup table\n"\
" --time-bins         (INT4) no. of time bins in the time-psi lookup table\n"\
" --prop-file         file containing the parameters to search over and their\n\
                     upper and lower ranges\n"\
" --nlive             (INT4) no. of live points for nested sampling\n"\
" --ephem-earth       Earth ephemeris file\n"\
" --ephem-sun         Sun ephemeris file\n"\
"\n"


LALIFOData *readPulsarData(int argc, char **argv)
{
  CHAR *detectors;
  CHAR *pulsar;
  CHAR *parfile;
  CHAR *inputdir;
  CHAR *outputdir;
  CHAR *propfile;
  
  CHAR *efile;
  CHAR *sfile;
  
  BinaryPulsarParams pars;
  
  struct option long_options[] =
  {
    { "help",           no_argument,       0, 'h' },
    { "verbose",        no_argument,    NULL, 'R' },
    { "detectors",      required_argument, 0, 'D' },
    { "pulsar",         required_argument, 0, 'p' },
    { "par-file",       required_argument, 0, 'P' },
    { "input-dir",      required_argument, 0, 'i' },
    { "output-dir",     required_argument, 0, 'o' },
    { "prop-file",      required_argument, 0, 'F' },
    { "ephem-earth",    required_argument, 0, 'J' },
    { "ephem-sun",      required_argument, 0, 'M' },
    { 0, 0, 0, 0 }
  };

  CHAR args[] = "hD:p:P:i:o:F:n:";
  CHAR *program = argv[0];

  /* parse input arguments */
  while( 1 ){
    INT4 option_index = 0;
    INT4 c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if( c == -1 ) /* end of options */
      break;

    switch( c ){
      case 0:
        if( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error passing option %s with argument %s\n",
            long_options[option_index].name, optarg);
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 'R': /* verbose */
        verbose = 1;
        break;
      case 'D': /* detectors */
        detectors = optarg;
        break;
      case 'p': /* pulsar name */
        pulsar = optarg;
        break;
      case 'P': /* pulsar parameter file */
        parfile = optarg;
        break;
      case 'i': /* input data file directory */
        inputdir = optarg;
        break;
      case 'o': /* output directory */
        outputdir = optarg;
        break;
      case 'J':
        efile = optarg;
        break;
      case 'M':
        sfile = optarg;
        break;
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n");
      default:
        break;
    }
  }

  /* get pulsar information from .par file */
	XLALReadTEMPOParFile(&pars, parfile);
  
  /* count the number of detectors from command line argument */
  
  /* read in data */
  for( i = 0 ; i < numDets ; i++ ){
    CHAR datafile[256];
    REAL8 times=0;
    REAL8Vector *temptimes=NULL;
    
    /*============================ GET DATA ==================================*/
    /* get detector B_ks data file in form finehet_JPSR_DET */
    sprintf(dataFile, "%s/data%s/finehet_%s_%s", inputs.inputDir, dets[i],
      inputs.pulsar, dets[i]);

    /* open data file */
    if((fp = fopen(dataFile, "r"))==NULL){
      fprintf(stderr, "Error... can't open data file %s!\n", dataFile);
      exit(0);
    }

    j=0;

    /* read in data */
    data->compTimeData = NULL;
    data->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &times, 0., 1.,
      &lalSecondUnit, MAXLENGTH );
    
    data->dataTimes = NULL;
    data->dataTimes = XLALCreateTimestampVector( MAXLENGTH );

    stdh0 = 0.;
    /* read in data */
    while(fscanf(fp, "%lf%lf%lf", &times, &dataVals.re, &dataVals.im) != EOF){
      /* check that size of data file is not to large */
      if( j == MAXLENGTH ){
        fprintf(stderr, "Error... size of MAXLENGTH not large enough.\n");
        exit(0);
      }

      /* exclude values smaller than 1e-28 as most are spurious points caused
         during a the heterodyne stage (e.g. when frame files were missing in
         the original S5 analysis) */
        temptimes->data[j] times);
        data->compTimeData->data->data[j] = dataVals;

        j++;
    }
    
    /* resize the data */
    data->compTimeData = XLALResizeCOMPLEX16TimeSeries(data->compTimeData,0,j);
    
    /* fill in time stamps as LIGO Time GPS Vector */
    data->dataTimes = NULL;
    data->dataTimes = XLALCreateTimestampVector( j );
    
    for ( k=0; k<j; k++ )
      XLALGPSSetREAL8(data->timeData->data[j], temptime);
    
    XLALDestroyREAL8Vector(temptimes);
  }
  
  /* set ephemeris data */
  data->edat = XLALMalloc(sizeof(*edat));

  data->(*edat).ephiles.earthEphemeris = efile;
  data->(*edat).ephiles.sunEphemeris = sfile;

  /* check files exist and if not output an error message */
  if( access(inputs.earthfile, F_OK) != 0 || 
      access(inputs.earthfile, F_OK) != 0 ){
    fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
    exit(3);
  }

  /* set up ephemeris information */
  LAL_CALL( LALInitBarycenter(&status, data->edat), &status );
}


void initialiseAlgorithm(ProcParamsTable *commandLine, LALInferenceRunState *runState)
/* Populates the structures for the algorithm control in runState, given the
 commandLine arguments. Includes setting up a random number generator.*/
{
	UINT4 verbose=0;
	ProcessParamsTable *ppt=NULL;
	
	ppt=getProcParamVal(commandLine,"--verbose");
	if(ppt) {
		verbose=1;
		addVariable(runState->algorithmParams,"verbose",1,INT4_t);
	}
		
	/* Number of live points */
	addVariable(runState->algorithmParams,"Nlive",atoi(getProcParamVal(commandLine,"--Nlive")->value),
				INT4_t);
	/* Number of points in MCMC chain */
	addVariable(runState->algorithmParams,"Nmcmc",atoi(getProcParamVal(commandLine,"--Nmcmc")->value),
				INT4_t);
	/* Optionally specify number of parallel runs */
	ppt=getProcParamVal(commandLine,"--Nruns");
	if(ppt) addVariable(runState->algorithmParams,"Nruns",atoi(ppt->value),INT4_t);

	/* Tolerance of the Nested sampling integrator */
	ppt=getProcParamVal(commandLine,"--tolerance");
	if(ppt) addVariable(runState->algorithmParams,"tolerance",atof(ppt->value),REAL4_t);
	
	/* Set up the random number generator */
	gsl_rng_env_setup();
	runState->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	/* (try to) get random seed from command line: */
	ppt = getProcParamVal(commandLine, "--randomseed");
	if (ppt != NULL)
		randomseed = atoi(ppt->value);
	else { /* otherwise generate "random" random seed: */
		if ((devrandom = fopen("/dev/random","r")) == NULL) {
			gettimeofday(&tv, 0);
			randomseed = tv.tv_sec + tv.tv_usec;
		} 
		else {
			fread(&randomseed, sizeof(randomseed), 1, devrandom);
			fclose(devrandom);
		}
	}
	fprintf(stdout, " initialize(): random seed: %lu\n", randomseed);
	gsl_rng_set(irs->GSLrandom, randomseed);
	
	/* Get chunk min and chunk max */
	ppt=getProcParamVal(commandLine,"--chunk-min");
	INT4 chunkMin;
	if(ppt) chunkMin=atoi(ppt-value)
		else chunkMin=5;
	addVariable(runState->algorithmParams,"chunk-min",chunkMin,INT4_t);
	
	ppt=getProcParamVal(commandLine,"--chunk-max");
	INT4 chunkMax;
	if(ppt) chunkMax=atoi(ppt-value)
		else chunkMax=30;
	addVariable(runState->algorithmParams,"chunk-max",chunkMax,INT4_t);
	
	if(verbose) fprintf(stdout,"Chunkmin = %i, chunkmax = %i\n",chunkMin,chunkMax);

	return;
}
	
void setupLookupTables(LALInferenceRunState runState, LALSource *source){	
	/* Set up lookup tables */
	/* Using psi bins, time bins */
	LALProcessParamsTable *ppt;
	LALProcessParamsTable *commandLine=runState->commandLine;
	
	ppt=getProcParamVal(commandLine,"--psi-bins");
	INT4 psiBins;
	if(ppt) psiBins=atoi(ppt-value)
		else psiBins=50;
	addVariable(runState->algorithmParams,"psi-bins",psiBins,INT4_t);
	
	ppt=getProcParamVal(commandLine,"--time-bins");
	INT4 timeBins;
	if(ppt) timeBins=atoi(ppt-value)
		else timeBins=1440;
	addVariable(runState->algorithmParams,"time-bins",timeBins,INT4_t);

	if(verbose) fprintf(stdout,"psi-bins = %i, time-bins = %i\n",psiBins,timeBins);

	LALIFOData *data=runState->data;

	gsl_matrix *LUfplus=NULL;
	gsl_matrix *LUfcross=NULL;
	
	REAL8 t0;
	LALDetAndSource detAndSource;
	
	while(data){
		t0=XLALGPSGetREAL8(data->dataTimes->data[0]);
		detAndSource.pDetector=data->detector;
		detAndSource.pSource=source;
		
		response_lookup_table(REAL8 t0, LALDetAndSource detAndSource,
						  timeBins, psiBins, LUfplus, LUfcross);
	
		addVariable(data->dataParams,"LU_Fplus",LUfplus,gslMatrix_t);
		addVariable(data->dataParams,"LU_Fcross",LUfcross,gslMatrix_t);
		data=data->next;
	}
	
	return;
}

void initialisePrior(ProcParamsTable *commandLine, LALInferenceRunState *runState)
/* Populates the priorArgs list in the runState using the command line arguments */
{
	
	return;
}


INT4 main(INT4 argc, CHAR *argv[]){
  static LALStatus status;

  REAL8 ****singleLike=NULL;
  REAL8 ****jointLike=NULL;

  InputParams inputs;
  BinaryPulsarParams pulsar;
  LALDetAndSource detAndSource;

  INT4 numDets=0; /* number of detectors */
  CHAR dets[5][3]; /* we'll have a max of five detectors */
  LALDetector detPos[5];

  INT4 i=0, j=0, k=0, n=0;

  DataStructure *data=NULL;
  REAL8 times=0.;
  COMPLEX16 dataVals;
  REAL8 stdh0=0.;       /* approximate h0 limit from data */

  FILE *fp=NULL;
  CHAR dataFile[256];
  CHAR outputFile[256];

  OutputParams output = empty_OutputParams;
  REAL8 maxPost=0.;
  REAL8 logNoiseEv[5]; /* log evidence for noise only (no signal) */
  Results results;
  REAL8 h0ul=0.;

  CHAR params[][10]={"h0", "phi", "psi", "ciota"};

  EphemerisData *edat=NULL;

	ProcParamsTable *param_table;
	LALInferenceRunState runState;
	
	
	/* Get ProcParamsTable from input arguments */
	param_table=parseCommandLine(argc,argv);
	runState.commandLine=param_table;
	
	/* Initialise data structures from the command line arguments */
	
	
	/* Initialise the algorithm structures from the command line arguments */
	/* Include setting up random number generator etc */
	/* Including generate lookup tables etc */
	
	
	/* Initialise the prior distribution given the command line arguments */
	
	
	
	/* Initialise the proposal distribution given the command line arguments */
	
	
	
	/* Call the nested sampling algorithm */
	runState.algorithm(&runState);
	
	/* Do any post-processing and output */
	
	
	/* End */
	
	
  /*===================== GET AND SET THE INPUT PARAMETER ====================*/
  get_input_args(&inputs, argc, argv);

  /* if we want to output in verbose mode set global variable */
  if(inputs.verbose) verbose = 1;

  /* get the pulsar parameters */
  XLALReadTEMPOParFile(&pulsar, inputs.parFile);
  inputs.psr.equatorialCoords.longitude = pulsar.ra;
  inputs.psr.equatorialCoords.latitude = pulsar.dec;
  inputs.psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* find the number of detectors being used */
  if( strstr(inputs.detectors, "H1") != NULL ){
    sprintf(dets[numDets], "H1");
    detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
    numDets++;
  }
  if( strstr(inputs.detectors, "H2") != NULL ){
     sprintf(dets[numDets], "H2");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }
  if( strstr(inputs.detectors, "L1") != NULL ){
     sprintf(dets[numDets], "L1");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }
  if( strstr(inputs.detectors, "G1") != NULL ){
     sprintf(dets[numDets], "G1");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }
  if( strstr(inputs.detectors, "V1") != NULL ){
     sprintf(dets[numDets], "V1");
     detPos[numDets] = *XLALGetSiteInfo( dets[numDets] );
     numDets++;
  }

  if( verbose ){
    fprintf(stderr, "Analysing data from %d detector(s):\n  ", numDets);
    for( i = 0 ; i < numDets ; i++ )
      fprintf(stderr, "%s ", dets[i]);
    fprintf(stderr, "\n");
  }
  /*==========================================================================*/

  /*====================== SET OUTPUT PARAMETERS =============================*/
  output.outputDir = inputs.outputDir;
  output.psr = inputs.pulsar;
  output.dob = inputs.dob; /* set degree of belief for UL */
  output.outPost = inputs.outputPost;
  sprintf(outputFile, "%s/evidence_%s", output.outputDir, output.psr);
  /*==========================================================================*/

  if( inputs.mcmc.doMCMC == 0 ){
    /* allocate likelihood memory */
    singleLike = allocate_likelihood_memory(inputs.mesh);

    /* if more than one detector create a joint likelihood */
    if( numDets > 1 )
      jointLike = allocate_likelihood_memory(inputs.mesh);

    if( verbose )
      fprintf(stderr, "I've allocated the memory for the likelihood.\n");

    /* if we're doing the grid search we only need to store one data set at a
       time */
    data = XLALMalloc(sizeof(DataStructure));
  }
  else{
    /* if we're doing the MCMC we need to store all the Bks */
    data = XLALMalloc(numDets*sizeof(DataStructure));

    /* if there's a covariance matrix file then set up the earth and sun
       ephemeris */
    if( inputs.matrixFile != NULL ){
      edat = XLALMalloc(sizeof(*edat));

      (*edat).ephiles.earthEphemeris = inputs.earthfile;
      (*edat).ephiles.sunEphemeris = inputs.sunfile;

      /* check files exist and if not output an error message */
      if( access(inputs.earthfile, F_OK) != 0 || 
          access(inputs.earthfile, F_OK) != 0 ){
        fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
        return 0;
      }

      LAL_CALL( LALInitBarycenter(&status, edat), &status );
    }
    else
      edat = NULL;
  }

  k = -1;

  /* read in data for each detector in turn an compute the likelihood */
  for( i = 0 ; i < numDets ; i++ ){
    /*============================ GET DATA ==================================*/
    /* get detector B_ks data file in form finehet_JPSR_DET */
    sprintf(dataFile, "%s/data%s/finehet_%s_%s", inputs.inputDir, dets[i],
      inputs.pulsar, dets[i]);

    /* open data file */
    if((fp = fopen(dataFile, "r"))==NULL){
      fprintf(stderr, "Error... can't open data file %s!\n", dataFile);
      exit(0);
    }

    j=0;

    /* only store one set of data for grid search (saves memory), but store
       them all for the MCMC */
    if( inputs.mcmc.doMCMC == 0 ) k = 0;
    else k++;

    /* read in data */
    data[k].data = NULL;
    data[k].data = XLALCreateCOMPLEX16Vector( MAXLENGTH );
    data[k].times = NULL;
    data[k].times = XLALCreateREAL8Vector( MAXLENGTH );

    /* set the minimum and maximum data segments lengths */
    data[k].chunkMin = inputs.chunkMin;
    data[k].chunkMax = inputs.chunkMax;

    stdh0 = 0.;
    /* read in data */
    while(fscanf(fp, "%lf%lf%lf", &times, &dataVals.re, &dataVals.im) != EOF){
      /* check that size of data file is not to large */
      if( j == MAXLENGTH ){
        fprintf(stderr, "Error... size of MAXLENGTH not large enough.\n");
        exit(0);
      }

      /* exclude values smaller than 1e-28 as most are spurious points caused
         during a the heterodyne stage (e.g. when frame files were missing in
         the original S5 analysis) */
      if( fabs(dataVals.re) > 1.e-28 && fabs(dataVals.im) > 1.e-28 ){
        data[k].times->data[j] = times;
        data[k].data->data[j] = dataVals;

        /* get the power from the time series */
        stdh0 += dataVals.re*dataVals.re + dataVals.im*dataVals.im;

        j++;
      }
    }

    fclose(fp);

    if( verbose )
      fprintf(stderr, "I've read in the data for %s.\n", dets[i]);

    data[k].data = XLALResizeCOMPLEX16Vector(data[k].data, j);
    data[k].times = XLALResizeREAL8Vector(data[k].times, j);

    /* if there is no input range for h0 then estimate it from the data */
    /* only do this once if performing grid search, but do for each seperate
       data set if doing MCMC */
    if( ( inputs.mesh.maxVals.h0 == 0 || inputs.mcmc.sigmas.h0 == 0 ) && 
        ( inputs.mcmc.doMCMC == 1 || i == 0 ) ){
      if( verbose ) fprintf(stderr, "Calculating h0 UL estimate: ");

      /* get the power spectral density power/bandwidth (1/60 Hz) */
      stdh0 = stdh0/((REAL8)j*(1./60.));

      /* upper limit estimate comes from ~ h0 = 10.8*sqrt(Sn/T) */
      stdh0 = 10.8*sqrt(stdh0/((REAL8)j*60.));

      /* set the MCMC h0 proposal step size at stdh0*scalefac */
      if( inputs.mcmc.doMCMC == 1 ){
        inputs.mcmc.sigmas.h0 = stdh0*inputs.mcmc.h0scale;
        
        if( inputs.mesh.maxVals.h0 == 0 )
          inputs.mesh.maxVals.h0 = stdh0;
      }

      /* set h0 max value for the grid at 5 times the expected ul */
      if( inputs.mesh.maxVals.h0 == 0 ){
        inputs.mesh.maxVals.h0 = 5.*stdh0;
        inputs.mesh.delta.h0 = (inputs.mesh.maxVals.h0 -
          inputs.mesh.minVals.h0)/(REAL8)(inputs.mesh.h0Steps - 1.);
      }

      if( verbose ) fprintf(stderr, "%le\n", stdh0);
    }

    /*========================================================================*/

    output.det = dets[i];

    /* create lookup table */
    data[k].lookupTable = NULL;
    data[k].lookupTable = XLALCalloc(1, sizeof(DetRespLookupTable));
    data[k].lookupTable->lookupTable=NULL;
    detAndSource.pSource = &inputs.psr;
    detAndSource.pDetector = &detPos[i];

    /* create memory for the lookup table */
    data[k].lookupTable->lookupTable = XLALCalloc(inputs.mesh.psiRangeSteps, 
      sizeof(LALDetAMResponse *));

    for( j = 0 ; j < inputs.mesh.psiRangeSteps ; j++ ){
      data[k].lookupTable->lookupTable[j] =
        XLALCalloc(inputs.mesh.timeRangeSteps, sizeof(LALDetAMResponse));
    }

    data[k].lookupTable->psiSteps = inputs.mesh.psiRangeSteps;
    data[k].lookupTable->timeSteps = inputs.mesh.timeRangeSteps;

    /* create lookup table */
    response_lookup_table(data[k].times->data[0], detAndSource,
      data[k].lookupTable);

    if( verbose ) fprintf(stderr, "Created look-up table.\n");

    /*======================== CALCULATE LIKELIHOOD ==========================*/
    if( inputs.mcmc.doMCMC == 0 ){
      logNoiseEv[i] = create_likelihood_grid(data[k], singleLike, inputs.mesh);

      if( verbose )
        fprintf(stderr, "I've calculated the likelihood for %s.\n", dets[i]);

      /* if there's more than one detector calculate the joint likelihood */
      if( numDets > 1 ){
        /* add the single detector log likelihood onto the joint likelihood */
        combine_likelihoods(singleLike, jointLike, inputs.mesh);
        output.outPost = 0; /* don't output individual detector posteriors */
      }
      /*======================================================================*/

      /*========== CREATE THE SINGLE DETECTOR POSTERIORS =====================*/
      maxPost = log_posterior(singleLike, inputs.priors, inputs.mesh, output);

      /* marginalise over each parameter and output the data */
      for( n = 0 ; n < 4 ; n++ ){
        output.margParam = params[n];
        if( verbose )
          fprintf(stderr, "Marginalising over %s.\n", output.margParam);

        results = marginalise_posterior(singleLike, inputs.mesh, output);

        if( output.dob != 0. && strcmp( output.margParam, "h0" ) == 0)
          h0ul = results.h0UpperLimit;
      }

      /* open file to output the evidence */
      if((fp = fopen(outputFile, "a"))==NULL){
        fprintf(stderr, "Error... can't open file %s!\n", outputFile);
        return 0;
      }

      /* output the log odds ratio and an UL if requested */
      fprintf(fp, "%s\t%le\n", output.det, results.evidence-logNoiseEv[i]);
      if( output.dob != 0. )
        fprintf(fp, "%.1lf%% h0 upper limit = %le\n", output.dob, h0ul);

      fclose(fp);

      XLALDestroyCOMPLEX16Vector(data[k].data);
      XLALDestroyREAL8Vector(data[k].times);
    }
    /*========================================================================*/

    /*================== PERFORM THE MCMC ====================================*/
    if( inputs.mcmc.doMCMC != 0 && inputs.onlyjoint == 0 )
      perform_mcmc(&data[k], inputs, 1, dets[i], &detPos[i], edat);

    /*========================================================================*/
  }

  /*=================== CREATE THE JOINT POSTERIOR IF REQUIRED ===============*/
  if( numDets > 1 ){
    output.det = joint_string;

    if( inputs.mcmc.doMCMC == 0 ){
      REAL8 totLogNoiseEv=0.;

      for( n = 0 ; n < numDets ; n++ ) totLogNoiseEv += logNoiseEv[n];

      output.outPost = inputs.outputPost; /* set for whether we want to output
                                            the full posterior */

      maxPost = log_posterior(jointLike, inputs.priors, inputs.mesh, output);
      if( verbose )
        fprintf(stderr, "I've calculated the joint posterior.\n");

      /* marginalise over each parameter and output the data */
      for( n = 0 ; n < 4 ; n++ ){
        output.margParam = params[n];
        if( verbose )
          fprintf(stderr, "Marginalising over %s.\n", output.margParam);

        results = marginalise_posterior(jointLike, inputs.mesh, output);
      }

      /* open file to output the evidence */
      if((fp = fopen(outputFile, "a"))==NULL){
        fprintf(stderr, "Error... can't open file %s!\n", outputFile);
        return 0;
      }

      /* output the evidence */
      fprintf(fp, "%s\t%le\n", output.det, results.evidence - totLogNoiseEv);
      if( output.dob != 0. ){
        fprintf(fp, "%.1lf%% h0 upper limit = %le\n", output.dob,
          results.h0UpperLimit);
      }
      fclose(fp);
    }
    /*========================================================================*/

    /*======================= PERFORM JOINT MCMC =============================*/
    if( inputs.mcmc.doMCMC == 1 ){
      perform_mcmc(data, inputs, numDets, output.det, detPos, edat);

      /* destroy data */
      for( i = 0 ; i < numDets ; i++ ){
        XLALDestroyCOMPLEX16Vector(data[i].data);
        XLALDestroyREAL8Vector(data[i].times);
      }
    }
  }

  /*====================== FREE THE LIKELIHOOD MEMORY ========================*/
  if( inputs.mcmc.doMCMC == 0 ){
    for( i = 0 ; i < inputs.mesh.phiSteps ; i++ ){
      for( j = 0 ; j < inputs.mesh.ciotaSteps ; j++ ){
        if( numDets > 1 )
          XLALFree(jointLike[i][j]);

        XLALFree(singleLike[i][j]);
      }
      if( numDets > 1 )
        XLALFree(jointLike[i]);

      XLALFree(singleLike[i]);
    }

    if( numDets > 1 )
      XLALFree(jointLike);

    XLALFree(singleLike);
  }
  /*=========================================================================*/ 

  return 0;
}


/* a function to sum over the data */
REAL8Vector * sum_data(COMPLEX16Vect *data, UINT4Vector *chunkLength){
  INT4 chunkLength=0, length=0, i=0, j=0, count=0;
  COMPLEX16 B;
  REAL8Vector *sumData=NULL; 
  
  length = data->length + 1 - chunkLengths->data[chunkLengths->length-1];

  sumData = XLALCreateVector( length );
  
  for( i = 0 ; i < length ; i+= chunkLength ){
    chunkLength = chunkLengths->data[count];
    sumData->data[count] = 0.;

    for( j = i ; j < i + chunkLength ; j++){
      B.re = data->data[j].re;
      B.im = data->data[j].im;

      /* sum up the data */
      sumData->data[count] += (B.re*B.re + B.im*B.im);
    }

    count++;
  }
  
  sumData = XLALResizeREAL8Vector( sumData, count );
  
  return sumData;
}

/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a sidereal day from the start time
(t0) and from -pi/4 to pi/4 in psi */
void response_lookup_table(REAL8 t0, LALDetAndSource detAndSource,
  INT4 timeSteps, INT4 psiSteps, gsl_matrix *LUfplus, gsl_matrix *LUfcross){ 
  LIGOTimeGPS gps;
  REAL8 T=0;

  REAL8 fplus=0., fcross=0.;
  REAL8 psteps = (REAL8)lookupTable->psiSteps;
  REAL8 tsteps = (REAL8)lookupTable->timeSteps;

  INT4 i=0, j=0;
  
  for( i = 0 ; i < psiSteps ; i++ ){
    detAndSource.pSource->orientation = -(LAL_PI/4.) +
        (REAL8)i*(LAL_PI/2.) / ( psteps - 1. );

    for( j = 0 ; j < timeSteps ; j++ ){
      T = t0 + (REAL8)j*LAL_DAYSID_SI / tsteps;

      gps = XLALGPSSetREAL8(&gps, T);

      XLALComputeDetAMResponse(&fplus, &fcross,
        detAndSource.pDetector->response,
        detAndSource.pSource->equatorialCoords.longitude,
        detAndSource.pSource->equatorialCoords.latitude,
        detAndSource.pSource->orientation,
        XLALGreenwichMeanSiderealTime(&gps));
        
      gsl_matrix_set(LUfplus, i, j, fplus);
      gsl_matrix_set(LUfcross, i, j, fcross);
    }
  }
}


/* function to calculate the log(likelihood) given some data and a set of
   particular pulsar parameters */
LALLikelihoodFunction pulsar_log_likelihood( LALVariables *vars, 
  LALIFOData *data, LALTemplateFunction *get_pulsar_model ){
  INT4 i=0, j=0, count=0, k=0, cl=0;
  INT4 length=0;
  REAL8 chunkLength=0.;

  REAL8 tstart=0., T=0.;

  COMPLEX16 model;
  INT4 psibin=0, timebin=0;

  REAL8 plus=0., cross=0.;
  REAL8 sumModel=0., sumDataModel=0.;
  REAL8 chiSquare=0.;
  COMPLEX16 B, M;

  REAL8 exclamation[data.chunkMax+1]; /* all factorials up to chunkMax */
  REAL8 logOf2=log(2.);

  REAL8 loglike=0.; /* the log likelihood */

  INT4 first=0, through=0;

  REAL8 phiIni=0., phi=0.;

  /* copy model parameters to data parameters */
  data.modelParams = vars;
  
  /* get pulsar model */
  get_pulsar_model( data );
  
  /* to save time get all log factorials up to chunkMax */
  for( i = 0 ; i < data.chunkMax+1 ; i++ )
    exclamation[i] = log_factorial(i);
  
  for( i = 0 ; i < length ; i += chunkLength ){
    chunkLength = (REAL8)data->chunkLengths->data[count];

    if( chunkLength < data.chunkMin ){
      count++;

      if( through == 0 ) first = 0;

      continue;
    }

    through = 1;

    sumModel = 0.;
    sumDataModel = 0.;

    cl = i + (INT4)chunkLength;
    
    for( j = i ; j < cl ; j++){
      B.re = data->compTimeData->data->data[j].re;
      B.im = data->compTimeData->data->data[j].im;

      M.re = data->compModelData->data->data[j].re;
      M.im = data->compModelData->data->data[j].im;
      
      /* sum over the model */
      sumModel += M.re*M.re + M.im*M.im;

      /* sum over that data and model */
      sumDataModel += B.re*M.re + B.im*M.im;
    }

    chiSquare = data->sumData->data[count];
    chiSquare -= 2.*sumDataModel;
    chiSquare += sumModel;

    if( first == 0 ){
      loglike = (chunkLength - 1.)*logOf2;
      first++;
    }
    else loglike += (chunkLength - 1.)*logOf2;

    loglike += exclamation[(INT4)chunkLength];
    loglike -= chunkLength*log(chiSquare);

    count++;
  }

  return loglike;
}

/* function to creates a vector of the phase differences between that used for
an initial set of phase parameters (i.e. those used to heterodyne the data),
and a new set of phase parameters */
void get_pulsar_model( LALIFOData *data ){
  INT4 i=0, length=0;
  
  BinaryPulsarParams pars;
  
  REAL8Vector *dphi=NULL;
  REAL8Vector *amp=NULL;
  
  /* set model parameters */
  pars.h0 = *(REAL8*)getVariable( data->modelParams, "h0");
  pars.cosiota = *(REAL8*)getVariable( data->modelParams, "cosiota");
  pars.psi = *(REAL8*)getVariable( data->modelParams, "psi");
  pars.phi0 = *(REAL8*)getVariable( data->modelParams, "phi0");
  
  /* set the potentially variable parameters */
  pars.pepoch = *(REAL8*)getVariable( data->modelParams, "pepoch");
  pars.posepoch = *(REAL8*)getVariable( data->modelParams, "posepoch");
  
  pars.ra = *(REAL8*)getVariable( data->modelParams, "ra");
  pars.pmra = *(REAL8*)getVariable( data->modelParams, "pmra");
  pars.dec = *(REAL8*)getVariable( data->modelParams, "dec");
  pars.pmdec = *(REAL8*)getVariable( data->modelParams, "pmdec");
  
  pars.f0 = *(REAL8*)getVariable( data->modelParams, "f0");
  pars.f1 = *(REAL8*)getVariable( data->modelParams, "f1");
  pars.f2 = *(REAL8*)getVariable( data->modelParams, "f2");
  pars.f3 = *(REAL8*)getVariable( data->modelParams, "f3");
  pars.f4 = *(REAL8*)getVariable( data->modelParams, "f4");
  pars.f5 = *(REAL8*)getVariable( data->modelParams, "f5");
  
  pars.model = *(CHAR*)getVariable( data->modelParams, "model");
  
  /* binary parameters */
  if( pars.model != NULL ){
    pars.e = *(REAL8*)getVariable( data->modelParams, "e");
    pars.w0 = *(REAL8*)getVariable( data->modelParams, "w0");
    pars.Pb = *(REAL8*)getVariable( data->modelParams, "Pb");
    pars.x = *(REAL8*)getVariable( data->modelParams, "x");
    pars.T0 = *(REAL8*)getVariable( data->modelParams, "T0");
    
    pars.e2 = *(REAL8*)getVariable( data->modelParams, "e2");
    pars.w02 = *(REAL8*)getVariable( data->modelParams, "w02");
    pars.Pb2 = *(REAL8*)getVariable( data->modelParams, "Pb2");
    pars.x2 = *(REAL8*)getVariable( data->modelParams, "x2");
    pars.T02 = *(REAL8*)getVariable( data->modelParams, "T02");
    
    pars.e3 = *(REAL8*)getVariable( data->modelParams, "e3");
    pars.w03 = *(REAL8*)getVariable( data->modelParams, "w03");
    pars.Pb3 = *(REAL8*)getVariable( data->modelParams, "Pb3");
    pars.x3 = *(REAL8*)getVariable( data->modelParams, "x3");
    pars.T03 = *(REAL8*)getVariable( data->modelParams, "T03");
    
    pars.xpbdot = *(REAL8*)getVariable( data->modelParams, "xpbdot");
    
    pars.eps1 = *(REAL8*)getVariable( data->modelParams, "eps1");    
    pars.eps2 = *(REAL8*)getVariable( data->modelParams, "eps2");       
    pars.eps1dot = *(REAL8*)getVariable( data->modelParams, "eps1dot");
    pars.eps2dot = *(REAL8*)getVariable( data->modelParams, "eps2dot");
    pars.Tasc = *(REAL8*)getVariable( data->modelParams, "Tasc");
    
    pars.wdot = *(REAL8*)getVariable( data->modelParams, "wdot"); 
    pars.gamma = *(REAL8*)getVariable( data->modelParams, "gamma");
    pars.Pbdot = *(REAL8*)getVariable( data->modelParams, "Pbdot");  
    pars.xdot = *(REAL8*)getVariable( data->modelParams, "xdot");   
    pars.edot = *(REAL8*)getVariable( data->modelParams, "edot");
    
    pars.s = *(REAL8*)getVariable( data->modelParams, "s"); 
    pars.dr = *(REAL8*)getVariable( data->modelParams, "dr");
    pars.dth = *(REAL8*)getVariable( data->modelParams, "dth");   
    pars.a0 = *(REAL8*)getVariable( data->modelParams, "a0");
    pars.b0 = *(REAL8*)getVariable( data->modelParams, "b0"); 

    pars.M = *(REAL8*)getVariable( data->modelParams, "M"); 
    pars.m2 = *(REAL8*)getVariable( data->modelParams, "m2");
  }

  amp = get_amplitude_model( pars, data );
  length = amp->length;
  
  /* assume that timeData vector within the LALIFOData structure contains the
     phase calculated using the initial (heterodyne) values of the phase
     parameters */
   
  /* get difference in phase and perform extra heterodyne with it */ 
  if ( (dphi = get_phase_model( pars, data )) != NULL){
    for( i=0; i<length; i++ ){
      COMPLEX16 M;
      REAL8 dphit;
      REAL4 sp, cp;
    
      dphit = -fmod(dphi->data[i] - data->timeData->data[i], 1.);
    
      sin_cos_2PI_LUT( &sp, &cp, dphit );
    
      M.re = data->compModelData->data->data[i].re;
      M.im = data->compModelData->data->data[i].im;
    
      /* heterodyne */
      data->compModelData->data->data[i].re = M.re*cp - M.im*sp;
      data->compModelData->data->data[i].im = M.im*cp + M.re*sp;
    }
  }
}

REAL8Vector *get_phase_model( BinaryPulsarParams params, LALIFOData *data ){
  static LALStatus status;

  INT4 i=0, length;

  REAL8 T0=0., DT=0., DTplus=0., deltat=0., deltat2=0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */

  EarthState earth, earth2;
  EmissionTime emit, emit2;
  REAL8 emitdt=0.;

  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;

  BarycenterInput bary;
  EphemerisData edat;
  
  REAL8Vector *phis=NULL;

  /* if edat is NULL then return a NULL poniter */
  if( data->ephem == NULL )
    return NULL;
  
  /* copy barycenter and ephemeris data */
  memcpy(bary, data->bary, sizeof(data->bary));
  memcpy(edat, data->ephem, sizeof(data->ephem));

   /* set the position and frequency epochs if not already set */
  if(params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if(params.posepoch == 0. && params.pepoch != 0.)
    params.posepoch = params.pepoch;

  length = data->dataTimes.length;
  
  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );

  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( params.px != 0. )
    bary.dInv = params.px*1e-3*LAL_C_SI/LAL_PC_SI;
  else if( params.dist != 0. )
    bary.dInv = LAL_C_SI/(params.dist*1e3*LAL_PC_SI);
  else
    bary.dInv = 0.;

  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8(data->dataTimes->data[i]);
    
    T0 = params.pepoch;

    DT = realT - T0;

    /* only do call the barycentring routines every 30 minutes, otherwise just
       linearly interpolate between them */
    if( i==0 || DT > DTplus ){
      bary.tgps.gpsSeconds = data->dataTimes->data[i].gpsSeconds;
      bary.tgps.gpsNanoSeconds = data->dataTimes->data[i].gpsNanoSeconds;

      bary.delta = params.dec + (realT-params.posepoch) * params.pmdec;
      bary.alpha = params.ra + (realT-params.posepoch) *
         params.pmra/cos(bary.delta);

      /* call barycentring routines */
      LAL_CALL( LALBarycenterEarth(&status, &earth, &bary.tgps, edat),
        &status );
      LAL_CALL( LALBarycenter(&status, &emit, &bary, &earth), &status );

      /* add interptime to the time */
      DTplus = DT + interptime;
      bary.tgps = XLALGPSAdd(bary.tgps, interptime);

      /* No point in updating the positions as difference will be tiny */
      LAL_CALL( LALBarycenterEarth(&status, &earth2, &bary.tgps, edat),
        &status );
      LAL_CALL( LALBarycenter(&status, &emit2, &bary, &earth2), &status );
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

  return phis;
}

REAL8Vector *get_amplitude_model( BinaryPulsarParams pars, LALIFOData *data ){
  INT4 i=0, length;
  
  REAL8Vector *amp=NULL;
  
  REAL8 psteps, tsteps;
  INT4 psibin, timebin;
  REAL8 tstart;
  REAL8 plus, cross;

  REAL8 Xplus, Xcross;
  REAL8 Xpcosphi_2, Xccosphi_2, Xpsinphi_2, Xcsinphi_2;
  REAL8 sinphi, cosphi;
  
  length = data->dataTimes.length;
  
  /* set lookup table parameters */
  psteps = *(INT4*)getVariable( data->modelParams, "psiSteps" );
  tsteps = *(INT4*)getVariable( data->modelParams, "timeSteps" );
  
  /* allocate memory for amplitudes */
  amp = XLALCreateREAL8Vector( length );
  
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
 
  for( i=0; i<length; i++ ){
    /* set the psi bin for the lookup table */
    psibin = (INT4)ROUND( ( pars.psi + LAL_PI/4. ) * ( psteps-1. )/LAL_PI_2 );

    tstart = XLALGPSGetREAL8( data->dataTimes->data[0] ); /*time of first B_k*/
  
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = fmod(XLALGPSGetREAL8(data->dataTimes->data[i]) - tstart,
      LAL_DAYSID_SI);
    timebin = (INT4)fmod( ROUND(T*tsteps/LAL_DAYSID_SI), tsteps );

    plus = gsl_matrix_get(data->lookupTablePlus, psibin, timebin);
    cross = gsl_matrix_get(data->lookupTableCross, psibin, timebin);
    
    /* create the complex signal amplitude model */
    data->compModelData->data->data[i].re = plus*Xpcosphi + cross*Xcsinphi;
    data->compModelData->data->data[i].im = plus*Xpsinphi - cross*Xccosphi;
  }
}

void add_initial_variables( LALVariables *ini, BinaryPulsarParams pars ){
  /* amplitude model parameters */
  addVariable(&ini, "h0", &pars.h0, REAL8_t);
  addVariable(&ini, "phi0", pars.phi0, REAL8_t);
  addVariable(&ini, "cosiota", pars.cosiota, REAL8_t);
  addVariable(&ini, "psi", pars.psi, REAL8_t);
  
  /* phase model parameters */
  
  /* frequency */
  addVariable(&ini, "f0", &pars.f0, REAL8_t);
  addVariable(&ini, "f1", &pars.f1, REAL8_t);
  addVariable(&ini, "f2", &pars.f2, REAL8_t);
  addVariable(&ini, "f3", &pars.f3, REAL8_t);
  addVariable(&ini, "f4", &pars.f4, REAL8_t);
  addVariable(&ini, "f5", &pars.f5, REAL8_t);
  addVariable(&ini, "pepoch", &pars.pepoch, REAL8_t);
  
  /* sky position */
  addVariable(&ini, "ra", &pars.ra, REAL8_t);
  addVariable(&ini, "pmra", &pars.pmra, REAL8_t);
  addVariable(&ini, "dec", &pars.dec, REAL8_t);
  addVariable(&ini, "pmdec", &pars.pmdec, REAL8_t);
  addVariable(&ini, "posepoch", &pars.posepoch, REAL8_t);
  
  /* binary system parameters */
  addVariable(&ini, "model", &pars.model, string);
  
  addVariable(&ini, "Pb", &pars.Pb, REAL8_t);
  addVariable(&ini, "e", &pars.e, REAL8_t);
  addVariable(&ini, "eps1", &pars.eps1, REAL8_t);
  addVariable(&ini, "eps2", &pars.eps2, REAL8_t);
  addVariable(&ini, "T0", &pars.T0, REAL8_t);
  addVariable(&ini, "Tasc", &pars.Tasc, REAL8_t);
  addVariable(&ini, "x", &pars.x, REAL8_t);
  addVariable(&ini, "w0", &pars.w0, REAL8_t);

  addVariable(&ini, "Pb2", &pars.Pb2, REAL8_t);
  addVariable(&ini, "e2", &pars.e2, REAL8_t);
  addVariable(&ini, "T02", &pars.T02, REAL8_t);
  addVariable(&ini, "x2", &pars.x2, REAL8_t);
  addVariable(&ini, "w02", &pars.w02, REAL8_t);
  
  addVariable(&ini, "Pb3", &pars.Pb3, REAL8_t);
  addVariable(&ini, "e3", &pars.e3, REAL8_t);
  addVariable(&ini, "T03", &pars.T03, REAL8_t);
  addVariable(&ini, "x3", &pars.x3, REAL8_t);
  addVariable(&ini, "w03", &pars.w03, REAL8_t);
  
  addVariable(&ini, "xpbdot", &pars.xpbdot, REAL8_t);
  addVariable(&ini, "eps1dot", &pars.eps1dot, REAL8_t);
  addVariable(&ini, "eps2dot", &pars.eps2dot, REAL8_t);
  addVariable(&ini, "wdot", &pars.wdot, REAL8_t);
  addVariable(&ini, "gamma", &pars.gamma, REAL8_t);
  addVariable(&ini, "Pbdot", &pars.Pbdot, REAL8_t);
  addVariable(&ini, "xdot", &pars.xdot, REAL8_t);
  addVariable(&ini, "edot", &pars.edot, REAL8_t);
  
  addVariable(&ini, "s", &pars.s, REAL8_t);
  addVariable(&ini, "dr", &pars.dr, REAL8_t);
  addVariable(&ini, "dth", &pars.dth, REAL8_t);
  addVariable(&ini, "a0", &pars.a0, REAL8_t);
  addVariable(&ini, "b0", &pars.b0, REAL8_t);
  addVariable(&ini, "M", &pars.M, REAL8_t);
  addVariable(&ini, "m2", &pars.m2, REAL8_t);
}

/* things I need to pass to the likelihood function via the LALInferenceRunState
  - 
  - a REAL8 value giving the initial time stamp for the above detector response 


/* FOR REFERENCE - using LIGOTimeGPSVector requires the
XLALCreateTimestampVector function in SFTutils */

/* functions to read in, or set, specific entries from an array */
/* REAL8 get_array_value( REAL8Array *array, UINT4Vector *entry ){
  UINT4 idx=0;
  UINT4 i=0; */
  
  /* check that the length of the vector containin the required array entry
and the array dimension length are the same */
  /* if ( array->dimLength->length != entry->length ){
    fprintf(stderr, "Error... entries not the same length!\n");
    exit(3);
  } */
  
  /* check that required entries are within the dimesion lengths */
  /* for( i=0; i<entry->length; i++ ){
    if (entry->data[i] > array->dimLength->data[i]){
      fprintf(stderr, "Error... entries greater than maximum value!\n");
      exit(4);
    }
  } */
  
  /* get the index of the entry in the array */
/*  for ( i=0; i<entry->length - 1; i++ )
    idx += entry->data[i] * array->dimLength->data[0];
  
  idx += entry->data[i];
  
  return array->data[idx];
} */

/* void set_array_value( REAL8Array *array, UINT4Vector *entry, REAL8 val ){
  UINT4 idx=0;
  UINT4 i=0; */
  
  /* check that the length of the vector containin the required array entry
and the array dimension length are the same */
/*  if ( array->dimLength->length != entry->length ){
    fprintf(stderr, "Error... entries not the same length!\n");
    exit(3);
  } */
  
  /* check that required entries are within the dimesion lengths */
  /*for( i=0; i<entry->length; i++ ){
    if (entry->data[i] > array->dimLength->data[i]){
      fprintf(stderr, "Error... entries greater than maximum value!\n");
      exit(4);
    }
  }*/
  
  /* get the index of the entry in the array */
  /*for ( i=0; i<entry->length - 1; i++ )
    idx += entry->data[i] * array->dimLength->data[0];
  
  idx += entry->data[i];
  
  array->data[idx] = val;
}*/