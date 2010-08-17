/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

#include "parameter_estimation_nested.h"


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
" --input-dir         directory containing the input data directories:\n\
                       within this the files should be in a directory\n\
                       strcuture with the format dataDETNAME, and have files\n\
                       with names finehet_PSRNAME_DETNAME. If not in this\n\
                       format use --force-file\n"\
" --force-file        full paths and file names for the data for each\n\
                     detector in the list (must be in the same order)\n\
                     delimited by commas\n"\
" --output-dir        directory for output data files\n"\
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


INT4 main(INT4 argc, CHAR *argv[]){
  ProcessParamsTable *param_table;
  LALInferenceRunState runState;
    
  /* Get ProcParamsTable from input arguments */
  param_table=parseCommandLine(argc,argv);
  runState.commandLine=param_table;
  
  /* Initialise data structures from the command line arguments */
  runState.data = readPulsarData(argc,argv);
  printf("Done readPulsarData\n");
  
  /* Initialise the algorithm structures from the command line arguments */
  /* Include setting up random number generator etc */
  initialiseAlgorithm(&runState);
  printf("Done initialiseAlgorithm\n");
  
  runState.algorithm=&NestedSamplingAlgorithm;
  runState.evolve=&NestedSamplingOneStep;
  
  runState.likelihood=&pulsar_log_likelihood;
  runState.prior=&priorFunction;
  runState.template=get_pulsar_model;
  
  /* Generate the lookup tables and read parameters from par file */
  setupFromParFile(&runState);
  printf("Done setupFromParFile\n");
  
  /* Initialise the prior distribution given the command line arguments */
  /* Initialise the proposal distribution given the command line arguments */
  initialiseProposal(&runState);
  printf("Done initialiseProposal\n");
 
  runState.proposal=LALInferenceProposalDifferentialEvolution;
  
  /* Create live points array and fill initial parameters */
  setupLivePointsArray(&runState);
  printf("Done setupLivePointsArray\n");
  
  /* Call the nested sampling algorithm */
  runState.algorithm(&runState);
  printf("Done nested sampling algorithm\n");
  
  /* Do any post-processing and output */
  
  
  /* End */
  

  /* if we want to output in verbose mode set global variable */


  return 0;
}

REAL8 priorFunction(LALInferenceRunState *runState, LALVariables *params)
{
	(void)runState;
	LALVariableItem *item=params->head;
	REAL8 min, max;
	for(;item;item=item->next)
	{
		if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR) continue;
		else
		{
			getMinMaxPrior(params, item->name, (void *)&min, (void *)&max);
			if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
		}
	}
	return (0);	
}


LALIFOData *readPulsarData(int argc, char *argv[])
{
  CHAR *detectors=NULL;
  CHAR *pulsar=NULL;
  CHAR *inputdir=NULL;
  CHAR *fname=NULL;
  CHAR *outputdir=NULL;
  CHAR *propfile=NULL;
  CHAR *forcefile=NULL;
  
  CHAR *filestr=NULL;
  
  CHAR *efile=NULL;
  CHAR *sfile=NULL;
  
  CHAR dets[MAXDETS][256];
  INT4 numDets=0, i=0;
 
  BinaryPulsarParams pars;
 
  LALIFOData *ifodata=NULL, *head=NULL;
  LALIFOData *prev=NULL;
  
  struct option long_options[] =
  {
    { "help",           no_argument,       0, 'h' },
    { "detectors",      required_argument, 0, 'D' },
    { "pulsar",         required_argument, 0, 'p' },
    { "input-dir",      required_argument, 0, 'i' },
    { "output-dir",     required_argument, 0, 'o' },
    { "force-file",     required_argument, 0, 'F' },
    { "ephem-earth",    required_argument, 0, 'J' },
    { "ephem-sun",      required_argument, 0, 'M' },
    { 0, 0, 0, 0 }
  };

  
  CHAR args[] = "hD:p:i:o:F:n:";
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
        detectors = XLALStringDuplicate(optarg);
        break;
      case 'p': /* pulsar name */
        pulsar = XLALStringDuplicate(optarg);
        break;
      case 'i': /* input data file directory */
        inputdir = XLALStringDuplicate(optarg);
        break;
      case 'o': /* output directory */
        outputdir = XLALStringDuplicate(optarg);
        break;
      case 'F': /* force file names to be specifed values */
        forcefile = XLALStringDuplicate(optarg);
        break;
      case 'J':
        efile = XLALStringDuplicate(optarg);
        break;
      case 'M':
        sfile = XLALStringDuplicate(optarg);
        break;
      //case '?':
        //fprintf(stderr, "Unknown error while parsing options\n");
      default:
        break;
    }
  }
  
  /* count the number of detectors from command line argument of comma separated
     vales and set their names */
  {
    CHAR *tempdets=NULL;
    CHAR *tempdet=NULL;
    
    tempdets = XLALStringDuplicate( detectors );
    
    while(1){
      tempdet = strsep(&tempdets, ",");
      XLALStringCopy(dets[numDets], tempdet, strlen(tempdet)+1);
      
      numDets++;
      
      if( tempdets == NULL ) break;
    }
  }
  
  /* check ephemeris files exist and if not output an error message */
  if( access(sfile, F_OK) != 0 || access(efile, F_OK) != 0 ){
    fprintf(stderr, "Error... ephemeris files not, or incorrectly, \
defined!\n");
    exit(3);
  }
  
  /* if file names are pre-defined (and seperated by comma, count the number of
     input files (by counting commas) and check it's equal to the number of
     detectors */
  if ( forcefile != NULL ){
    CHAR *inputstr=NULL;
    CHAR *tempstr=NULL;
    INT4 count=0;
    
    inputstr = XLALStringDuplicate( forcefile );

    /* count number of commas */
    while(1){
      tempstr = strsep(&inputstr, ",");
      
      if (inputstr == NULL) break;
      
      count++;
    }
    
    if ( count+1 != numDets ){
      fprintf(stderr, "Error... Number of input files given is not equal to the\
 number of detectors given!\n");
      exit(3);
    }
  }
  
  /* reset filestr */
  if ( forcefile != NULL ) filestr = XLALStringDuplicate(forcefile);

  /* read in data */
  for( i = 0,prev=NULL ; i < numDets ; i++,prev=ifodata ){
    CHAR *datafile=NULL;
    REAL8 times=0;
    LIGOTimeGPS gpstime;
    COMPLEX16 dataVals;
    REAL8Vector *temptimes=NULL;
    INT4 j=0, k=0;
    
    FILE *fp=NULL;
    
    ifodata=XLALCalloc(1,sizeof(LALIFOData));
    ifodata->next=NULL;
	  ifodata->dataParams=XLALCalloc(1,sizeof(LALVariables));
    if(i==0) head=ifodata;
    if(i>0) prev->next=ifodata;
	
    /* set detector */
    ifodata->detector = XLALGetSiteInfo( dets[i] );

    /*============================ GET DATA ==================================*/
    /* get detector B_ks data file in form finehet_JPSR_DET */
    if (forcefile == NULL){
      datafile = XLALMalloc( 256 ); /* allocate memory for file name */

      sprintf(datafile, "%s/data%s/finehet_%s_%s", inputdir, dets[i],
        pulsar, dets[i]);
    }
    else{ /* get i'th filename from the comma separated list */
      INT4 count=0;
      while (count <= i){
        datafile = strsep(&filestr, ",");
      
        count++;
      }
    }
	  printf("datafile name: %s\n",datafile);
    /* open data file */
    if((fp = fopen(datafile, "r"))==NULL){
      fprintf(stderr, "Error... can't open data file %s!\n", datafile);
      exit(0);
    }

    j=0;

    /* read in data */
	  temptimes = XLALCreateREAL8Vector( MAXLENGTH );

    /* read in data */
    while(fscanf(fp, "%lf%lf%lf", &times, &dataVals.re, &dataVals.im) != EOF){
      /* check that size of data file is not to large */
      if (j == 0){
        XLALGPSSetREAL8( &gpstime, times );

        ifodata->compTimeData = NULL;
        ifodata->compTimeData = XLALCreateCOMPLEX16TimeSeries( "", &gpstime, 0.,
1., &lalSecondUnit, MAXLENGTH );
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
      XLALResizeCOMPLEX16TimeSeries(ifodata->compTimeData,0,j);

    /* fill in time stamps as LIGO Time GPS Vector */
    ifodata->dataTimes = NULL;
    ifodata->dataTimes = XLALCreateTimestampVector( j );
	  
    for ( k=0; k<j; k++ )
      XLALGPSSetREAL8(&ifodata->dataTimes->data[k], temptimes->data[k]);

    XLALDestroyREAL8Vector(temptimes);

    /* set ephemeris data */
    ifodata->ephem = XLALMalloc(sizeof(EphemerisData));
    
    /* set up ephemeris information */
    ifodata->ephem = XLALInitBarycenter( efile, sfile );
  }

  return head;
}


void initialiseAlgorithm(LALInferenceRunState *runState)
/* Populates the structures for the algorithm control in runState, given the
 commandLine arguments. Includes setting up a random number generator.*/
{
	ProcessParamsTable *ppt=NULL;
	ProcessParamsTable *commandLine=runState->commandLine;
	REAL8 tmp;
	INT4 tmpi;
	INT4 randomseed;
	
  FILE *devrandom=NULL;
  struct timeval tv;
  
	ppt=getProcParamVal(commandLine,"--verbose");
	if(ppt) {
		verbose=1;
		addVariable(runState->algorithmParams,"verbose", &verbose , INT4_t,
PARAM_FIXED);
	}
	
	/* Initialise parameters structure */
	runState->algorithmParams=XLALCalloc(1,sizeof(LALVariables));
	runState->priorArgs=XLALCalloc(1,sizeof(LALVariables));
	runState->proposalArgs=XLALCalloc(1,sizeof(LALVariables));
	
	printf("set number of live points.\n");
	/* Number of live points */
	tmpi=atoi(getProcParamVal(commandLine,"--Nlive")->value);
	addVariable(runState->algorithmParams,"Nlive",&tmpi, INT4_t,PARAM_FIXED);
	
  printf("set number of MCMC points.\n");
  /* Number of points in MCMC chain */
	tmpi=atoi(getProcParamVal(commandLine,"--Nmcmc")->value);
	addVariable(runState->algorithmParams,"Nmcmc",&tmpi,
				INT4_t,PARAM_FIXED);
	
  printf("set number of parallel runs.\n");
  /* Optionally specify number of parallel runs */
	ppt=getProcParamVal(commandLine,"--Nruns");
	if(ppt) {
		tmpi=atoi(ppt->value);
		addVariable(runState->algorithmParams,"Nruns",&tmpi,INT4_t,PARAM_FIXED);
	}
	
  printf("set tolerance.\n");
  /* Tolerance of the Nested sampling integrator */
	ppt=getProcParamVal(commandLine,"--tolerance");
	if(ppt){
		tmp=strtod(ppt->value,(char **)NULL);
		addVariable(runState->algorithmParams,"tolerance",&tmp, REAL8_t,
PARAM_FIXED);
	}
	
  printf("set random seed.\n");
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
	fprintf(stdout, " initialize(): random seed: %u\n", randomseed);
	gsl_rng_set(runState->GSLrandom, randomseed);

	return;
}
	
void setupLookupTables(LALInferenceRunState *runState, LALSource *source){	
	/* Set up lookup tables */
	/* Using psi bins, time bins */
	ProcessParamsTable *ppt;
	ProcessParamsTable *commandLine=runState->commandLine;

  INT4 chunkMin, chunkMax;
  
	ppt=getProcParamVal(commandLine,"--psi-bins");
	INT4 psiBins;
	if(ppt) psiBins=atoi(ppt->value);
		else psiBins=50;
	addVariable(runState->algorithmParams,"psi-bins",&psiBins,INT4_t,PARAM_FIXED);
	
	ppt=getProcParamVal(commandLine,"--time-bins");
	INT4 timeBins;
	if(ppt) timeBins=atoi(ppt->value);
		else timeBins=1440;
	addVariable(runState->algorithmParams,"time-bins",&timeBins,INT4_t,
    PARAM_FIXED);

	if(verbose) fprintf(stdout,"psi-bins = %i, time-bins =%i\n",psiBins,timeBins);

  /* Get chunk min and chunk max */
  ppt=getProcParamVal(commandLine,"--chunk-min");
  if(ppt) chunkMin=atoi(ppt->value);
  else chunkMin=5;
    
  ppt=getProcParamVal(commandLine,"--chunk-max");
  if(ppt) chunkMax=atoi(ppt->value);
  else chunkMax=30; 
  
  if(verbose) fprintf(stdout,"Chunkmin = %i, chunkmax =%i\n", chunkMin,
    chunkMax);
  
	LALIFOData *data=runState->data;

	gsl_matrix *LUfplus=NULL;
	gsl_matrix *LUfcross=NULL;
	
	REAL8 t0;
	LALDetAndSource detAndSource;
	
	while(data){
		REAL8Vector *sumData=NULL;
    UINT4Vector *chunkLength=NULL;

    t0=XLALGPSGetREAL8(&data->dataTimes->data[0]);
		detAndSource.pDetector=data->detector;
		detAndSource.pSource=source;
		
		LUfplus = gsl_matrix_alloc(psiBins, timeBins);

		LUfcross = gsl_matrix_alloc(psiBins, timeBins);
	response_lookup_table(t0, detAndSource, timeBins, psiBins, LUfplus, LUfcross);
    addVariable(data->dataParams,"LU_Fplus",LUfplus,gslMatrix_t,PARAM_FIXED);
    addVariable(data->dataParams,"LU_Fcross",LUfcross,gslMatrix_t,PARAM_FIXED);

    addVariable(data->dataParams,"chunk-min",&chunkMin,INT4_t,PARAM_FIXED);
    addVariable(data->dataParams,"chunk-max",&chunkMax,INT4_t,PARAM_FIXED);

    /* get chunk lengths of data */
    chunkLength = get_chunk_lengths( data, chunkMax );
    addVariable(data->dataParams, "chunkLength", &chunkLength, UINT4Vector_t,
      PARAM_FIXED);

    /* get sum of data for each chunk */
    sumData = sum_data( data );
    addVariable(data->dataParams, "sumData", &sumData, REAL8Vector_t,
      PARAM_FIXED);

		data=data->next;
	}
	
	return;
}

void initialiseProposal(LALInferenceRunState *runState)
/* sets up the parameters to be varied, and the extent of the proposal
   distributions from a proposal distribution file */
{
  CHAR *propfile;
  ProcessParamsTable *commandLine=runState->commandLine;
  FILE *fp=NULL;
  
  CHAR tempPar[100]="";
  REAL8 low, high;
  
  BinaryPulsarParams pulsar;
  REAL8Vector *phase_vector;
  
  propfile = (CHAR *)getProcParamVal(commandLine,"--prop-file")->value;
  
  runState->priorArgs = XLALCalloc(1,sizeof(LALVariables));
  
  /* open file */
  fp = fopen(propfile, "r");
  
  while(fscanf(fp, "%s %lf %lf", tempPar, &low, &high) != EOF){
    REAL8 tempVar;
    VariableType type;
    REAL8 tempmin, tempmax;
    
    tempVar = *(REAL8*)getVariable(runState->currentParams, tempPar);
    type = getVariableType(runState->currentParams, tempPar);
    
    /* remove variable value */
    removeVariable(runState->currentParams, tempPar);
    
    /* re-add variable */
    if(strstr(tempPar, "ra") || strstr(tempPar, "phi0")){
      addVariable(runState->currentParams, tempPar, &tempVar, type,
        PARAM_CIRCULAR);
    }
    else{
      addVariable(runState->currentParams, tempPar, &tempVar, type,
        PARAM_LINEAR);
    }
	  
    /* Add the prior variables */
	  addMinMaxPrior(runState->priorArgs, tempPar, (void *)&low, (void *)&high,
type);
  }
  
	return;
}

void setupFromParFile(LALInferenceRunState *runState)
/* Read the PAR file of pulsar parameters and setup the code using them */
/* Generates lookup tables also */
{
	
	LALSource psr;
	BinaryPulsarParams pulsar;
	REAL8Vector *phase_vector;
	LALIFOData *data=runState->data;
	ProcessParamsTable *ppt=NULL;
	
	ppt=getProcParamVal(runState->commandLine,"--par-file");
	if(ppt==NULL) {fprintf(stderr,"Must specify --par-file!\n"); exit(1);}
	char *parFile=ppt->value;
	
	/* get the pulsar parameters */
	XLALReadTEMPOParFile(&pulsar, parFile);
	psr.equatorialCoords.longitude = pulsar.ra;
	psr.equatorialCoords.latitude = pulsar.dec;
	psr.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  printf("Read info from par file.\n");
 
	/* Setup lookup tables for amplitudes */
	setupLookupTables(runState, &psr);
	printf("Setup lookup tables.\n");
  
	/* Setup initial phase */
	while(data){
		phase_vector=get_phase_model(pulsar, data );
		data->timeData=XLALCalloc(1,sizeof(REAL8TimeSeries));
		data->timeData->data=phase_vector;
		data->timeData->epoch=data->dataTimes->data[0];
		data->timeData->sampleUnits=lalSecondUnit;
		data=data->next;
	}
  
  runState->currentParams = XLALCalloc(1,sizeof(LALVariables));
  
	/* Add initial (unchanging) variables for the model. */
	add_initial_variables( runState->currentParams, pulsar );

	return;
}


void setupLivePointsArray(LALInferenceRunState *runState){
/* Set up initial basket of live points, drawn from prior,
 by copying runState->currentParams to all entries in the array*/
	
	UINT4 Nlive=(UINT4)*(INT4 *)getVariable(runState->algorithmParams,"Nlive");
	UINT4 i;
	LALVariableItem *current;
	
  REAL8 temp, mini, maxi; 
  
	/* Allocate the array */
	/* runState->livePoints=XLALCalloc(Nlive,sizeof(LALVariables *)); */
  runState->livePoints=XLALCalloc(1,sizeof(LALVariables *));
	
  for(i=0;i<Nlive;i++)
	{
		runState->livePoints[i]=XLALCalloc(1,sizeof(LALVariables));
   
    /* Copy the param structure */
		copyVariables(runState->currentParams,runState->livePoints[i]);
		
    /* Sprinkle the varying points among prior */
    
    for(current=runState->livePoints[i]->head ;current!=NULL;
      current=current->next){
			if(current->vary==PARAM_CIRCULAR || current->vary==PARAM_LINEAR)
			{
				switch (current->type){
					case REAL4_t:
					case REAL8_t:
						{
							REAL8 tmp;
							REAL8 min,max;
							getMinMaxPrior(runState->priorArgs,current->name, 
                (void *)&min,(void *)&max);
							tmp=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
						}
						break;
					case INT4_t:
					case INT8_t:
						{
							INT4 tmp;
							INT4 min,max;
							getMinMaxPrior(runState->priorArgs,current->name,
                (void *)&min,(void *)&max);
							tmp=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
						}
						break;
					default:
						fprintf(stderr,"Trying to randomise a non-numeric parameter!");
				}
			}
		}
	}
}


/* function to get the lengths of consecutive chunks of data */
UINT4Vector *get_chunk_lengths( LALIFOData *data, INT4 chunkMax ){
  INT4 i=0, j=0, count=0;
  INT4 length;
  
  REAL8 t1, t2;
  
  UINT4Vector *chunkLengths=NULL;
  
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

    t1 = XLALGPSGetREAL8(&data->dataTimes->data[i-1]);
    t2 = XLALGPSGetREAL8(&data->dataTimes->data[i]);
    
    /* if consecutive points are within 180 seconds of each other count as in
       the same chunk */
    if( t2 - t1 > 180. || count == chunkMax ){
      chunkLengths->data[j] = count;
      count = 0; /* reset counter */

      j++;
    }
  }

  chunkLengths = XLALResizeUINT4Vector(chunkLengths, j);
  
  return chunkLengths;
}


/* a function to sum over the data */
REAL8Vector * sum_data( LALIFOData *data ){
  INT4 chunkLength=0, length=0, i=0, j=0, count=0;
  COMPLEX16 B;
  REAL8Vector *sumData=NULL; 

  UINT4Vector *chunkLengths;
  
  chunkLengths = *(UINT4Vector **)getVariable(data->dataParams, "chunkLength");
  
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
void response_lookup_table(REAL8 t0, LALDetAndSource detAndSource, INT4
timeSteps, INT4 psiSteps, gsl_matrix *LUfplus, gsl_matrix *LUfcross){ 
  LIGOTimeGPS gps;
  REAL8 T=0;

  REAL8 fplus=0., fcross=0.;
  REAL8 psteps = (REAL8)psiSteps;
  REAL8 tsteps = (REAL8)timeSteps;

  INT4 i=0, j=0;
  
  for( i = 0 ; i < psiSteps ; i++ ){
    detAndSource.pSource->orientation = -(LAL_PI/4.) +
        (REAL8)i*(LAL_PI/2.) / ( psteps - 1. );

    for( j = 0 ; j < timeSteps ; j++ ){
      T = t0 + (REAL8)j*LAL_DAYSID_SI / tsteps;

      XLALGPSSetREAL8(&gps, T);

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
REAL8 pulsar_log_likelihood( LALVariables *vars, LALIFOData *data,
LALTemplateFunction *get_model ){
  INT4 i=0, j=0, count=0, k=0, cl=0;
  INT4 length=0, chunkMin, chunkMax;
  REAL8 chunkLength=0.;

  REAL8 tstart=0., T=0.;

  COMPLEX16 model;
  INT4 psibin=0, timebin=0;

  REAL8 plus=0., cross=0.;
  REAL8 sumModel=0., sumDataModel=0.;
  REAL8 chiSquare=0.;
  COMPLEX16 B, M;

  REAL8 *exclamation=NULL; /* all factorials up to chunkMax */
  REAL8 logOf2=log(2.);

  REAL8 loglike=0.; /* the log likelihood */

  INT4 first=0, through=0;

  REAL8 phiIni=0., phi=0.;
  
  REAL8Vector *sumData=NULL;
  UINT4Vector *chunkLengths=NULL;

  sumData = (REAL8Vector*)getVariable(data->dataParams, "sumData");
  chunkLengths = (UINT4Vector*)getVariable(data->dataParams, "chunkLength");
  chunkMin = *(INT4*)getVariable(data->dataParams, "chunkMin");
  chunkMax = *(INT4*)getVariable(data->dataParams, "chunkMax");
  
  exclamation = XLALCalloc(chunkMax+1, sizeof(REAL8));
  
  /* copy model parameters to data parameters */
  copyVariables(data->modelParams, vars);
  
  /* get pulsar model */
  get_model( data );
  
  /* to save time get all log factorials up to chunkMax */
  for( i = 0 ; i < chunkMax+1 ; i++ )
    exclamation[i] = log_factorial(i);
  
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

    chiSquare = sumData->data[count];
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
  
  pars.model = (CHAR*)getVariable( data->modelParams, "model");
  
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

  get_amplitude_model( pars, data );
  length = data->compModelData->data->length;
  
  /* assume that timeData vector within the LALIFOData structure contains the
     phase calculated using the initial (heterodyne) values of the phase
     parameters */
   
  /* get difference in phase and perform extra heterodyne with it */ 
  if ( (dphi = get_phase_model( pars, data )) != NULL){
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

  BarycenterInput *bary=NULL;
  EphemerisData *edat=NULL;
  
  REAL8Vector *phis=NULL;

  /* if edat is NULL then return a NULL poniter */
  if( data->ephem == NULL )
    return NULL;
  
  /* copy barycenter and ephemeris data */
	bary = (BarycenterInput*)XLALCalloc(1,sizeof(data->bary));
  bary->site = *data->detector;
  bary->alpha = params.ra;
  bary->delta = params.dec;
  
  edat = (EphemerisData*)XLALCalloc(1,sizeof(data->ephem));  
  edat = data->ephem;
  
   /* set the position and frequency epochs if not already set */
  if(params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if(params.posepoch == 0. && params.pepoch != 0.)
    params.posepoch = params.pepoch;

  length = data->dataTimes->length;
  
  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );
  printf("allocated memory for phis\n"); 
  
  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( params.px != 0. )
    bary->dInv = params.px*1e-3*LAL_C_SI/LAL_PC_SI;
  else if( params.dist != 0. )
    bary->dInv = LAL_C_SI/(params.dist*1e3*LAL_PC_SI);
  else
    bary->dInv = 0.;

  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8(&data->dataTimes->data[i]);
    
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
      LAL_CALL( LALBarycenterEarth(&status, &earth, &bary->tgps, edat),
        &status );
      LAL_CALL( LALBarycenter(&status, &emit, bary, &earth), &status );

      /* add interptime to the time */
      DTplus = DT + interptime;
      XLALGPSAdd(&bary->tgps, interptime);

      /* No point in updating the positions as difference will be tiny */
      LAL_CALL( LALBarycenterEarth(&status, &earth2, &bary->tgps, edat),
        &status );
      LAL_CALL( LALBarycenter(&status, &emit2, bary, &earth2), &status );
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

void get_amplitude_model( BinaryPulsarParams pars, LALIFOData *data ){
  INT4 i=0, length;
  
  REAL8Vector *amp=NULL;
  
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
  
  LU_Fplus = (gsl_matrix*)getVariable( data->dataParams, "LU_Fplus");
  LU_Fcross = (gsl_matrix*)getVariable( data->dataParams, "LU_Fcross");
  
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

    tstart = XLALGPSGetREAL8( &data->dataTimes->data[0] ); /*time of first B_k*/
  
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = fmod(XLALGPSGetREAL8(&data->dataTimes->data[i]) - tstart,
      LAL_DAYSID_SI);
    timebin = (INT4)fmod( ROUND(T*tsteps/LAL_DAYSID_SI), tsteps );

    plus = gsl_matrix_get(LU_Fplus, psibin, timebin);
    cross = gsl_matrix_get(LU_Fcross, psibin, timebin);
    
    /* create the complex signal amplitude model */
    data->compModelData->data->data[i].re = plus*Xpcosphi + cross*Xcsinphi;
    data->compModelData->data->data[i].im = plus*Xpsinphi - cross*Xccosphi;
  }
}

void add_initial_variables( LALVariables *ini, BinaryPulsarParams pars ){
  /* amplitude model parameters */
  addVariable(ini, "h0", &pars.h0, REAL8_t,PARAM_FIXED);
  addVariable(ini, "phi0", &pars.phi0, REAL8_t,PARAM_FIXED);
  addVariable(ini, "cosiota", &pars.cosiota, REAL8_t,PARAM_FIXED);
  addVariable(ini, "psi", &pars.psi, REAL8_t,PARAM_FIXED);
  
  /* phase model parameters */
  
  /* frequency */
  addVariable(ini, "f0", &pars.f0, REAL8_t,PARAM_FIXED);
  addVariable(ini, "f1", &pars.f1, REAL8_t,PARAM_FIXED);
  addVariable(ini, "f2", &pars.f2, REAL8_t,PARAM_FIXED);
  addVariable(ini, "f3", &pars.f3, REAL8_t,PARAM_FIXED);
  addVariable(ini, "f4", &pars.f4, REAL8_t,PARAM_FIXED);
  addVariable(ini, "f5", &pars.f5, REAL8_t,PARAM_FIXED);
  addVariable(ini, "pepoch", &pars.pepoch, REAL8_t,PARAM_FIXED);
  
  /* sky position */
  addVariable(ini, "ra", &pars.ra, REAL8_t,PARAM_FIXED);
  addVariable(ini, "pmra", &pars.pmra, REAL8_t,PARAM_FIXED);
  addVariable(ini, "dec", &pars.dec, REAL8_t,PARAM_FIXED);
  addVariable(ini, "pmdec", &pars.pmdec, REAL8_t,PARAM_FIXED);
  addVariable(ini, "posepoch", &pars.posepoch, REAL8_t,PARAM_FIXED);
  
  /* binary system parameters */
  addVariable(ini, "model", &pars.model, string_t, PARAM_FIXED);
  
  addVariable(ini, "Pb", &pars.Pb, REAL8_t,PARAM_FIXED);
  addVariable(ini, "e", &pars.e, REAL8_t,PARAM_FIXED);
  addVariable(ini, "eps1", &pars.eps1, REAL8_t,PARAM_FIXED);
  addVariable(ini, "eps2", &pars.eps2, REAL8_t,PARAM_FIXED);
  addVariable(ini, "T0", &pars.T0, REAL8_t,PARAM_FIXED);
  addVariable(ini, "Tasc", &pars.Tasc, REAL8_t,PARAM_FIXED);
  addVariable(ini, "x", &pars.x, REAL8_t,PARAM_FIXED);
  addVariable(ini, "w0", &pars.w0, REAL8_t,PARAM_FIXED);

  addVariable(ini, "Pb2", &pars.Pb2, REAL8_t,PARAM_FIXED);
  addVariable(ini, "e2", &pars.e2, REAL8_t,PARAM_FIXED);
  addVariable(ini, "T02", &pars.T02, REAL8_t,PARAM_FIXED);
  addVariable(ini, "x2", &pars.x2, REAL8_t,PARAM_FIXED);
  addVariable(ini, "w02", &pars.w02, REAL8_t,PARAM_FIXED);
  
  addVariable(ini, "Pb3", &pars.Pb3, REAL8_t,PARAM_FIXED);
  addVariable(ini, "e3", &pars.e3, REAL8_t,PARAM_FIXED);
  addVariable(ini, "T03", &pars.T03, REAL8_t,PARAM_FIXED);
  addVariable(ini, "x3", &pars.x3, REAL8_t,PARAM_FIXED);
  addVariable(ini, "w03", &pars.w03, REAL8_t,PARAM_FIXED);
  
  addVariable(ini, "xpbdot", &pars.xpbdot, REAL8_t,PARAM_FIXED);
  addVariable(ini, "eps1dot", &pars.eps1dot, REAL8_t,PARAM_FIXED);
  addVariable(ini, "eps2dot", &pars.eps2dot, REAL8_t,PARAM_FIXED);
  addVariable(ini, "wdot", &pars.wdot, REAL8_t,PARAM_FIXED);
  addVariable(ini, "gamma", &pars.gamma, REAL8_t,PARAM_FIXED);
  addVariable(ini, "Pbdot", &pars.Pbdot, REAL8_t,PARAM_FIXED);
  addVariable(ini, "xdot", &pars.xdot, REAL8_t,PARAM_FIXED);
  addVariable(ini, "edot", &pars.edot, REAL8_t,PARAM_FIXED);
  
  addVariable(ini, "s", &pars.s, REAL8_t,PARAM_FIXED);
  addVariable(ini, "dr", &pars.dr, REAL8_t,PARAM_FIXED);
  addVariable(ini, "dth", &pars.dth, REAL8_t,PARAM_FIXED);
  addVariable(ini, "a0", &pars.a0, REAL8_t,PARAM_FIXED);
  addVariable(ini, "b0", &pars.b0, REAL8_t,PARAM_FIXED);
  addVariable(ini, "M", &pars.M, REAL8_t,PARAM_FIXED);
  addVariable(ini, "m2", &pars.m2, REAL8_t,PARAM_FIXED);
}

REAL8 log_factorial(INT4 num){
  INT4 i=0;
  REAL8 logFac=0.;
  
  for( i=2 ; i <= num ; i++ ) logFac += log((REAL8)i);
  
  return logFac;
}

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