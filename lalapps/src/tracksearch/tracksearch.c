/*
 * Author: Torres C. (Univ of TX at Brownsville)
 */

#include "tracksearch.h"
#include <lal/TSSearch.h>

#define PROGRAM_NAME "tracksearch"

typedef struct
{
  INT4    argc;
  CHAR**  argv;
}LALInitSearchParams;

/* 
 * Declare local subroutines to use 
 */

void initializeTSSearch(LALStatus*,int argc,char* argv[],TSSearchParams*,CHAR*,CHARVector*);
void Get_Data(LALStatus*,TSSearchParams*,REAL4TimeSeries*,CHARVector*,CHAR*);
void Get_Ascii_Data(LALStatus*,TSSearchParams*,REAL4TimeSeries*,CHARVector*);
void Do_TrackSearch(LALStatus*,REAL4Vector*,TSSearchParams,INT4);
void Dump_Search_Data(TSSearchParams,TrackSearchOut,CHAR*);
void QuickDump_Data(TSSearchParams,TrackSearchOut,CHAR*);
void fakeDataGeneration(LALStatus*,REAL4TimeSeries*,INT4,INT4);
void Create_Curve_DataSection(LALStatus*,Curve**);
void Destroy_Search_DataSection(LALStatus*,Curve**,INT4);


/* Code Identifying information */
NRCSID( TRACKSEARCHC, "tracksearch $Id$");
RCSID( "tracksearch $Id$");
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"





/* Non-Error Code Define Statements */
#define TRACKSEARCHC_NARGS 8
/* End N-E Defines */

/* Define Error Codes */
#define TRACKSEARCHC_ENORM              0
#define TRACKSEARCHC_ESUB               1
#define TRACKSEARCHC_EARGS              2
#define TRACKSEARCHC_EVAL               4
#define TRACKSEARCHC_EFILE              8
#define TRACKSEARCHC_EREAD              16
#define TRACKSEARCHC_EMEM               32
#define TRACKSEARCHC_EMISC              64

#define TRACKSEARCHC_MSGENORM          "Normal Exit"
#define TRACKSEARCHC_MSGESUB           "Subroutine Fail"
#define TRACKSEARCHC_MSGEARGS          "Arguement Parse Error"
#define TRACKSEARCHC_MSGEVAL           "Invalid Argument(s)"
#define TRACKSEARCHC_MSGEFILE          "Error Opening File"
#define TRACKSEARCHC_MSGEREAD          "Error Reading File"
#define TRACKSEARCHC_MSGEMEM           "Memory Error"
#define TRACKSEARCHC_MSGEMISC          "Unknown Error"

#define TRUE     1
#define FALSE    0

/* Usage format string. */
#define USAGE "Still in flux"

/* Code Body */
int main (int argc, char *argv[])
{
  /* Global Variable Declarations */
  /* Variables that will be in final produciton code */
  static LALStatus      status;/* Error containing structure */
  TSSearchParams       *params;/*Large struct for all params*/
  TSSegmentVector      *SegVec = NULL;/*Holds array of data sets*/
  TSCreateParams        SegVecParams;/*Hold info for memory allocation*/
  REAL4TimeSeries       dataset;/*Structure holding entire dataset to analyze*/
  REAL4Vector          *datavector = NULL;/*Vector holding copy of incoming data*/
  COMPLEX8Vector       *dummyresponse = NULL;/* Specifically not null */
  REAL8FrequencySeries *dummyspectra = NULL;/* Specifically not null */
  CHAR                 *cachefile = NULL;/* name of file with frame cache info */
  CHARVector           *dirpath=NULL;/*Void contain for path string */
  CHARVector           *dirname=NULL;/* name of directory with frames      */
  /* 
   * Variables that will be most likely removed after debugging 
   */
  UINT4             i;  /* Counter for data breaking */
  INT4              j;
  UINT4             productioncode = 1; /* Variable for ascii or frame */
                                        /* Set to zero to use ascii files */
  /*
   *Sleep for Attaching DDD 
   */
  unsigned int doze = 5;
  pid_t myPID;
  myPID = getpid( );
  fprintf( stdout, "pid %d sleeping for %d seconds\n", myPID, doze );
  fflush( stdout );
  sleep( doze );
  fprintf( stdout, "pid %d awake\n", myPID );
  fflush( stdout );

  /* End Global Variable Declarations */

  /* SET LAL DEBUG STUFF */
  set_debug_level("MEMDBG");

  /*Setup structure for path name temp*/
  LALCHARCreateVector(&status,&dirpath,1024);
  params=LALMalloc(sizeof(TSSearchParams));
  /* 
   * Parse incoming search arguments initialize variables 
   */
  initializeTSSearch(&status,argc,argv,params,cachefile,dirpath);
  /* 
   * Setup string for path to data frames
   */
  LALCHARCreateVector(&status,&dirname,strlen(dirpath->data)+1);
  strncpy(dirname->data,dirpath->data,dirname->length);
  LALCHARDestroyVector(&status,&dirpath);

  /* 
   * This is the data reading section of code
   * Sources: Frame file or Single Ascii file (1C)
   * All pathing is relative for this code
   * 0 - use frames
   * 1 - use Ascii file
   */
  if (productioncode == 1) /* Use production frame file inputs */
    {
      if (params->makenoise < 1 )
	{  
	  Get_Data(&status,params,&dataset,dirname,NULL);
	}
      else
	{
	  /* 
	   * Fake Data Generation routing generates UnitGaussian 
	   * random values 
	   *  the output array is desire length in a 
	   * time series structure 
	   */
	  fakeDataGeneration(&status,&dataset,
			     (params->NumSeg*params->SegLengthPoints),
			     params->makenoise);
	};
    }
  else
    {
      Get_Ascii_Data(&status,params,&dataset,dirname);
    }

  LALCHARDestroyVector(&status,&dirname);

  /* 
   * Prepare the data for the call to Tracksearch 
   */
  /* All functionality is still not implemented yet here */
  SegVecParams.dataSegmentPoints = dataset.data->length;
  SegVecParams.responseSegmentPoints = 0;
  SegVecParams.spectraSegmentPoints = 0;
  SegVecParams.numberDataSegments = params->NumSeg;
  LALCreateTSDataSegmentVector(&status,&SegVec,&SegVecParams);
  LALCreateVector(&status,&datavector,dataset.data->length);
  for(i = 0;i < dataset.data->length; i++)
    datavector->data[i] = dataset.data->data[i];
  /* 
   * TracksearchPrep creates a collection of overlapped segments
   * each of the segments will be used to make an individual TF map.
   * We also will apply whitening and condition the data via calls
   * inside this LAL routine.  The TSSearchParams struct has all
   * the information needed to accoplish this
   *
   * Function protype to change soon....
   */
  TrackSearchPrep(&status,&dataset,dummyresponse,dummyspectra,SegVec,*params);
  j=0;
  for(i = 0;i < params->NumSeg;i++)
    {
      printf(".");
      Do_TrackSearch(&status,SegVec->dataSeg[j].data,*(params),j);
      j++;
    };
  printf("\n");
  /* Free some of the memory used to do the analysis */
  LALDestroyVector(&status,&datavector);
  LALDestroyTSDataSegmentVector(&status,&SegVec);
  LALDestroyVector(&status,&(dataset.data));
  /* 
   * Free searchParams input arguement structure from InitializeRoutine
   * Free elements first 
   */
  if (params->channelName)
    LALFree(params->channelName);
  if (params->auxlabel)
    LALFree(params->auxlabel);
  if (params->dataSegVec)
    {
      LALDestroyTSDataSegmentVector(&status,&(params->dataSegVec));
    }
  if (params->numSlaves)
    LALFree(params->numSlaves);
  LALFree(params);   
  /* 
   * Done freeing for search params struct memory
   */
  LALCheckMemoryLeaks();
  return 0;
}
/* End Main Code Body */


/*
 * These are local scope subroutines mean for use 
 * by tracksearch only if it may be globally useful then we will place
 * it in the LAL libraries
 */

/*
 * Setup params structure by parsing the command line
 */
void initializeTSSearch(
			LALStatus     *status,
			int            argc,
			char          *argv[], 
			TSSearchParams  *params,
			CHAR*         cachefile,
			CHARVector   *dirname
			)
{

  /* Setup file option to read ascii types and not frames */
  /* getop arguments */
  struct option long_options[] =
    {
      {"gpsstart_seconds",    required_argument,  0,    'a'},
      {"total_time_points",   required_argument,  0,    'b'},
      {"map_type",            required_argument,  0,    'c'},
      {"line_width",          required_argument,  0,    'd'},
      {"start_threshold",     required_argument,  0,    'e'},
      {"member_threshold",    required_argument,  0,    'f'},
      {"length_threshold",    required_argument,  0,    'g'},
      {"power_threshold",     required_argument,  0,    'h'},
      {"channel_name",        required_argument,  0,    'i'},
      {"number_of_maps",      required_argument,  0,    'j'},
      {"number_of_time_bins", required_argument,  0,    'k'},
      {"overlap",             required_argument,  0,    'l'},
      {"window_type",         required_argument,  0,    'm'},
      {"number_of_freq_bins", required_argument,  0,    'n'},
      {"whiten_level",        required_argument,  0,    'o'},
      {"transform_type",      required_argument,  0,    'p'},
      {"segment_time_points", required_argument,  0,    'q'},
      {"cachefile",           required_argument,  0,    'v'},
      {"directory",           required_argument,  0,    's'},
      {"window_size",         required_argument,  0,    't'},
      {"fake",                no_argument,        0,    'u'},
      {"auxlabel",            required_argument,  0,    'r'},
      {"join_lines",          no_argument,        0,    'w'},
      {"PSD_framefile",       required_argument,  0,    'x'},
      {"spectrum_average_method", required_argument,0,  'y'},
      {"calibration_cache_file", required_argument,0,   'z'},
      {0,                     0,                  0,      0}
    };
  
  int              C;

  INITSTATUS (status, "getparams", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);

  /*
   * Set default values for optional arguments 
   * These options will be used if omitted from the command line 
   * Required arguements will be initialized as well to simply debugging
   */
  params->LineWidth = 0;
  params->TransformType = Spectrogram;
  params->NumSeg = 1;
  params->SegLengthPoints = 0; /* Set for debuging this with overlap # */
  params->overlapFlag = 0; /* Change this to an INT4 for overlapping */
  params->TimeLengthPoints = 0;
  params->discardTLP = 0;
  params->window = Hann; /*Welch*/
  params->whiten = 0;
  params->avgSpecMethod = -1; /*Will trigger error*/
  params->avgSpecWindow = Rectangular; /* Default no window */
  params->FreqBins = 0;
  params->TimeBins = 0;
  params->windowsize = 0;/*256*/ /*30*/
  params->StartThresh =-1; /* We choose invalid values to check incoming args*/
  params->LinePThresh =-1; /* We choose invalid values to check incoming args*/
  params->MinPower = 0;
  params->MinLength = 0;
  params->GPSstart.gpsSeconds = 0;
  params->GPSstart.gpsNanoSeconds = 0;
  /* 
   * In future...
   * Must include error checking for these values 
   */
  params->TimeLengthPoints = 0;
  params->SamplingRate = 1;
  params->makenoise = -1;
  params->calChannelType=SimDataChannel;/*See lsd*/
  params->channelName=NULL; /* Leave NULL to avoid frame read problems */
  params->dataDirPath=NULL;
  params->singleDataCache=NULL;
  params->detectorPSDCache=NULL;
  params->auxlabel=NULL;
  params->joinCurves=FALSE; /*default no curve joining */
  params->dataSegVec=NULL;
  params->numSlaves=0;
  params->calFrameCache=NULL;

  /* 
   * Parse the command line arguments 
   */
  if (argc < 8 ) /* Not enough arguments to run */
    {
      fprintf(stderr,TRACKSEARCHC_MSGEARGS);
      fprintf(stderr,"\n");
      exit(TRACKSEARCHC_EARGS);
    }
  
  while (TRUE)
    {
      int option_index=0;
      C = getopt_long_only(argc,
			   argv,
			   "a:b:c:d:e:f:h:i:j:k:l:m:o:p:q:r:s:t:u",
			   long_options, 
			   &option_index);
      /* The end of the arguments is C = -1 */
      if ( C == -1)
	{
	  break;
	}
      switch( C )
	{
	case 0:
	  /* if this option set a flag, do nothing else now */
	  if ( long_options[option_index].flag != 0 )
	    {
	      break;
	    }
	  else
	    {
	      fprintf( stderr, "error parsing option %s with argument %s\n",
		       long_options[option_index].name, optarg );
	      exit( 1 );
	    }
	  break;
	  
	case 'a':
	  /* Setting the GPS start time parameter */
	  {
	    UINT4 startTime = atoi(optarg);
	    params->GPSstart.gpsSeconds = startTime;
	    params->GPSstart.gpsNanoSeconds = 0;
	  }
	  break;
	  
	case 'b':
	  /* Setting the total length of data to analyze */
	  {
	    params->TimeLengthPoints = ((UINT4) atoi(optarg));
	  }
	  break;
	  
	case 'c':
	  /* Selecting Map Type to use */
	  {
	    /* Change the name field to be a integer instead of *char */
	    char *name = NULL;
	    name = (CHAR*) LALMalloc(strlen(optarg)+1);
	    if (!name)
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEMEM);
		exit(TRACKSEARCHC_EMEM);
	      }
	    /* We want to get the name as an enum value */
	    /*Add appropriate code HERE--Fix Later*/
	  }
	  break;
	  
	case 'd':
	  /* Setting Line Width parameter */
	  {
	    params->LineWidth = atoi(optarg);

	    /* We don't limit the max line width as an error the user
	     * is responsible to pick this with some care
	     */
	    if ((params->LineWidth) < 1)
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		exit(TRACKSEARCHC_EARGS);
	      }
	  }
	  break;
	  
	case 'e':
	  /* Setting Start Threshold parameter */
	  {
	    params->StartThresh = atof(optarg);
	    if (params->LinePThresh != -1)
	      {
		if (params->StartThresh < params->LinePThresh)
		  {
		    fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		    exit(TRACKSEARCHC_EARGS);
		  }
	      }
	  }
	  break;
	  
	case 'f':
	  /* Setting member point threshold paramter */
	  {
	    params->LinePThresh = atof(optarg);
	    if (params->StartThresh != -1)
	      {
		if (params->StartThresh < params->LinePThresh)
		  {
		    fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		    exit(TRACKSEARCHC_EARGS);
		  }
	      }
	  }
	  break;
	  
	case 'g':
	  /* Setting Thresholding on line length for search */
	  {
	    params->MinLength = atoi(optarg);
	    if (params->MinLength < 2)
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		exit(TRACKSEARCHC_EARGS);
	      }
	  }
	  break; 	  

	case 'h':
	  /* Setting Thresholding on integrated line power */
	  {
	    params->MinPower = atof(optarg);
	    if(params->MinPower < 0)
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		exit(TRACKSEARCHC_EARGS);
	      }
	  }
	  break;
	  
	case 'i':
	  /* Setting data source name to look for CHANNEL */
	  {
	    params->channelName = (CHAR*) LALMalloc(strlen(optarg)+1);
	    if (!(params->channelName))
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEMEM);
		exit(TRACKSEARCHC_EMEM);
	      };
	    strcpy(params->channelName,optarg);
	  }
	  break;
	  
	case 'j':
	  /* Setting number of maps to split requested data into */
	  {
	    fprintf(stderr,"Can not specify number of maps\n");
	    fprintf(stderr,TRACKSEARCHC_MSGEMEM);
	    exit(TRACKSEARCHC_EMEM);
	    /*  params->NumSeg = atoi(optarg);
		if (params->NumSeg < 1)
		{
		fprintf(stderr,TRACKSEARCHC_MSGEMEM);
		exit(TRACKSEARCHC_EMEM);
		} */
	  }
	  break;
	  
	case 'k':
	  /* Setting number of time bins for each TF map */
	  {
	    params->TimeBins = ((UINT4) atoi(optarg));
	  }
	  break;
	  
	case 'l':
	  {
	    params->overlapFlag = ((UINT4) atoi(optarg));
	  }
	  break;
	  
	case 'm':
	  /* Setup the designated windowing function for FFTing */
	  {
	    params->window = atoi(optarg);
	    /* Error check */
	    if ((params->window > NumberWindowTypes))
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		exit(TRACKSEARCHC_EARGS);
	      }
	  }
	  break;
	  
	case 'n':
	  /* Setup number of desire frequency bins to make map with */
	  {
	    params->FreqBins = atoi(optarg);
	  }
	  break;
	  
	case 'o':
	  /* Setup whiten level is specified */
	  {
	    params->whiten = atoi(optarg);
	    if (params->whiten > 2 )
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		exit(TRACKSEARCHC_EARGS);
	      }
	  }
	  break;

	case 'p':
	  /* Choose transform type */
	  {
	    /*Create string reader code*/
	    /* Implement Spectrogram type selection via einteger later */
	    /* For now we rely on numeric representation */
	    params->TransformType = atoi(optarg);
	    /*params->TransformType =  Spectrogram;*/
	  }
	  break;

	case 'q':
	  {/* Specify length of each time segment */
	    params->SegLengthPoints  = ((UINT4) atoi(optarg));
	    break;
	  }
	  break;

	case 'r':
	  {/* Auxlabel file name to read in */
	    params->auxlabel = (CHAR*) LALMalloc(strlen(optarg)+1);
	    if (!(params->auxlabel))
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEMEM);
		exit(TRACKSEARCHC_EMEM);
	      };
	    strncpy(params->auxlabel,optarg,(strlen(optarg)+1));
	  }
	  break;

	case 's':
	  {
	    INT4 len;
	    len= strlen(optarg) +1;
	    strncpy(dirname->data,optarg,len);
	  }
	  break;

	case 't':
	  { /* Looking further in code realize this isn't used !!!!*/
	    params->windowsize = atoi(optarg);
	  }
	  break;

	case 'u':
	  { /* Setting up the seed for random values */
	    params->makenoise = atoi(optarg);
	  }
	  break;

	case 'v':
	  { /* Setting up aux label if present */
	    /* Insert string copy code here */
	    INT4 len;
	    len = strlen(optarg) +1;
	    cachefile = (CHAR *) LALCalloc(len,sizeof(CHAR));
	    memcpy(cachefile,optarg,len);
	  }
	  break;

	case 'w':
	  {
	    params->joinCurves=TRUE;
	  }
	  break;

	case 'x':
	  {
	    INT4 len;
	    len=strlen(optarg)+1;
	    params->detectorPSDCache=(CHAR *)
	      LALCalloc(len,sizeof(CHAR));
	    memcpy(params->detectorPSDCache,optarg,len);
	  }

	case 'y':
	  {
	    if(!strcmp(optarg, "useMean"))
	      params->avgSpecMethod = useMean;
	    else if(!strcmp(optarg, "useMedian"))
	      params->avgSpecMethod = useMedian;
	    else if(!strcmp(optarg, "useUnity"))
	      params->avgSpecMethod = useUnity;
	    else 
	      {
		fprintf(stderr,TRACKSEARCHC_MSGEARGS);
		exit(TRACKSEARCHC_EARGS);
	      };
	    break;
	  }
	case 'z':
	  { /* Setting up aux label if present */
	    /* Insert string copy code here */
	    INT4 len;
	    len = strlen(optarg) +1;
	    params->calFrameCache = (CHAR *) LALCalloc(len,sizeof(CHAR));
	    memcpy(params->calFrameCache,optarg,len);
	  }
	  break;

	default :
	  {
	    fprintf(stderr,TRACKSEARCHC_MSGEMISC);
	    exit(TRACKSEARCHC_EMISC);
	  }
	  break;
	  
	};
    };
  /* 
   * Include here a simple section of code to santiy check all params 
   * structure arguements for reasonable values. This will involve 
   * moving error catching code out of the case statement 
   */

  /* 
   * Tabulate the number of segments to create for specified overlap
   * Check for valid overlap number which < segment length number 
   */
  if (params->overlapFlag >= params->SegLengthPoints)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEVAL);
      exit(TRACKSEARCHC_EVAL);
    };

  if (params->TimeLengthPoints < params->SegLengthPoints)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEVAL);
      exit(TRACKSEARCHC_EVAL);
    };

  if ( params->overlapFlag == 0)
    {
      params->NumSeg = floor(params->TimeLengthPoints/params->SegLengthPoints);
    }
  else
    {
      /* Determine the number of maps */
      params->NumSeg=floor((params->TimeLengthPoints-params->overlapFlag)
			   /(params->SegLengthPoints - params->overlapFlag));
      /* Determine number of points to throw away (not process) */
      params->discardTLP=params->TimeLengthPoints -
	(params->NumSeg*(params->SegLengthPoints -
			 params->overlapFlag)
	 /(params->TimeLengthPoints - params->SegLengthPoints));
      /* 
       * Reset points to process by N-discardTLP, this gives us
       * uniform segments to make maps with
       */
      params->TimeLengthPoints=params->TimeLengthPoints-params->discardTLP;
    };
  DETATCHSTATUSPTR(status);
  return;
}
/* End initialization subroutine for search */

/*
 * Local routine to actually carry out search
 */
void Do_TrackSearch(
		    LALStatus       *status,
		    REAL4Vector     *signal,
		    TSSearchParams   params,
		    INT4             CallNum
		    )
{
  CHAR                   strCallNum[16];
  CHAR                   tempchar[128];
  CreateTimeFreqIn       tfInputs;/*Input params for map making*/
  INT4                   j;
  LALWindowParams        windowParams;/*Needed to generate windowing funcs*/
  REAL4Window           *tempWindow = NULL;
  TimeFreqParam         *autoparams = NULL;/*SelfGenerated values for TFmap*/
  TimeFreqRep           *tfmap = NULL;/*TF map of dataset*/
  TrackSearchMapMarkingParams  mapMarkerParams;/*Struct of params for marking*/
  TrackSearchOut         outputCurves;/*Native LAL format for results*/
  TrackSearchOut         outputCurvesThreshold; /*Curve data thresholded */
  TrackSearchParams      inputs;/*Native LAL format for tracksearch module*/

  INITSTATUS (status, "Do_Tracksearch", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  tfInputs.type = params.TransformType;
  tfInputs.fRow = 2*params.FreqBins;
  tfInputs.tCol = params.TimeBins;
  /*  tfInputs.wlengthT = params.FreqBins+1;*/
  /*  tfInputs.wlengthF = params.FreqBins+1;*/
  tfInputs.wlengthF = params.windowsize;
  tfInputs.wlengthT = params.windowsize;

  LALCreateTimeFreqParam(status->statusPtr,&autoparams,&tfInputs);
  CHECKSTATUSPTR(status);

  LALCreateTimeFreqRep(status->statusPtr,&tfmap,&tfInputs);
  CHECKSTATUSPTR(status);

  for (j=0;j<tfmap->tCol;j++)
    {
      tfmap->timeInstant[j]=j;
    }
  windowParams.length = params.windowsize;
  windowParams.type = params.window;

  LALCreateREAL4Window(status->statusPtr,&(tempWindow),&(windowParams));
  CHECKSTATUSPTR(status);

  /* 
   * Case statment that look to do the various TF Reps 
   */
  switch ( tfInputs.type )
    {
    case Undefined:
      {
	fprintf(stderr,TRACKSEARCHC_MSGEVAL);
	exit(TRACKSEARCHC_EVAL);
      }
      break;
    case Spectrogram:
      {
	/* Required from deprication of LALWindow function */
	memcpy(autoparams->windowT->data,
	       tempWindow->data->data,
	       (windowParams.length * sizeof(REAL4)));
	LALTfrSp(status->statusPtr,signal,tfmap,autoparams);
      }
      break;
    case WignerVille:
      {

	LALTfrWv(status->statusPtr,signal,tfmap,autoparams);
      }
      break;
    case PSWignerVille:
      {
	/* Required from deprication of LALWindow function */
	memcpy(autoparams->windowT->data,
	       tempWindow->data->data,
	       (windowParams.length * sizeof(REAL4)));
	memcpy(autoparams->windowF->data,
	       tempWindow->data->data,
	       (windowParams.length * sizeof(REAL4)));
	LALTfrPswv(status->statusPtr,signal,tfmap,autoparams);
      }
      break;
    case RSpectrogram:
      {
	/* Required from deprication of LALWindow function */
	memcpy(autoparams->windowT->data,
	       tempWindow->data->data,
	       (windowParams.length * sizeof(REAL4)));
	LALTfrRsp(status->statusPtr,signal,tfmap,autoparams);
      }
      break;
    default:
      {
	fprintf(stderr,TRACKSEARCHC_MSGEMISC);
	exit(TRACKSEARCHC_EMISC);
      }
      break;
    };
  /*
   * Destroy window memory 
   * Destroy TF params also
   */
  LALDestroyREAL4Window(status->statusPtr,&tempWindow);
  CHECKSTATUSPTR(status);
  LALDestroyTimeFreqParam(status->statusPtr,&autoparams);
  CHECKSTATUSPTR(status);

  /*
   * Fill in LALSignalTrackSearch params structure via the use of
   * the TSSearch huge struct elements
   */
  inputs.sigma=params.LineWidth;
  inputs.high=params.StartThresh;
  inputs.low=params.LinePThresh;
  inputs.width=((tfmap->fRow/2)+1);
  inputs.height=tfmap->tCol;
  /* 
   * The LALTracksearch seems to map the 
   * width to tCol and the height to fRow
   */

  /* Flag to prepare structure inside LALTrackSearch */
  inputs.allocFlag = 1; 
  outputCurves.curves = NULL;
  /* Comment out WriteMap to reduce extra data file creation */
  WriteMap(*tfmap,*signal);

  /*
   * The power thresholding is not done in the LALroutine
   * This information is forwarded to this function for a post processing 
   * Should be default be givin minimum curve length parameter 
   * We want to apply thresholds in a seperate routine 
   */
  LALSignalTrackSearch(status->statusPtr,&outputCurves,tfmap,&inputs);
  CHECKSTATUSPTR(status);

  /* 
   * Call tracksearch again to free any temporary ram in 
   * variable outputCurves which is no longer required
   */
  inputs.allocFlag = 2;
  LALSignalTrackSearch(status->statusPtr,&outputCurves,tfmap,&inputs);
  CHECKSTATUSPTR(status);

  /*
   * Setup for call to function to do map marking
   * We mark maps is convert Tbin and Fbin to 
   * Hz and GPSseconds
   */
  mapMarkerParams.deltaT=1/(params.SamplingRate);
  mapMarkerParams.mapStartGPS.gpsSeconds=params.GPSstart.gpsSeconds;
  mapMarkerParams.mapStartGPS.gpsNanoSeconds=params.GPSstart.gpsNanoSeconds;
  mapMarkerParams.mapStopGPS.gpsSeconds=0;
  mapMarkerParams.mapStopGPS.gpsNanoSeconds=0;
  mapMarkerParams.mapTimeBins=inputs.height;/*Height -> time bins */
  mapMarkerParams.mapFreqBins=inputs.width; /*Width -> freq bins */
  LALTrackSearchInsertMarkers(status->statusPtr,
			      &outputCurves,
			      &mapMarkerParams); 
  CHECKSTATUSPTR(status);


  /* Call the connect curve routine if argument is specified */
  if (params.joinCurves)
    {
      LALTrackSearchConnectSigma(status->statusPtr,
				 &outputCurves,
				 *tfmap,
				 inputs);
      CHECKSTATUSPTR(status);
    };

  /* 
   * Diagnostic Code
   * Dump out this particular pgm map file for visual inspection 
   */
  sprintf(tempchar,"%s_%i",params.auxlabel,CallNum);
  fprintf("%s\n",tempchar);
  DumpTFImage(tfmap->map,tempchar,inputs.height,inputs.width,1);
  sprintf(strCallNum,"Pre_%i",CallNum);
  Dump_Search_Data(params,outputCurves,strCallNum);
  QuickDump_Data(params,outputCurves,strCallNum);

  /* 
   * Destroy the memory holding the actual TF map
   */
  LALDestroyTimeFreqRep(status->statusPtr,&tfmap);
  CHECKSTATUSPTR(status);
  
  /*
   * Apply the user requested thresholds 
   */
  outputCurvesThreshold.curves = NULL;
  outputCurvesThreshold.numberOfCurves = 0;
  LALTrackSearchApplyThreshold(status->statusPtr,
			       &outputCurves,
			       &outputCurvesThreshold,
			       params);
  CHECKSTATUSPTR(status);

  /* 
   * Dump out list of surving candidates
   * Long List 
   * Short List
   */
  sprintf(strCallNum,"%i",CallNum);
  Dump_Search_Data(params,outputCurvesThreshold,strCallNum);
  QuickDump_Data(params,outputCurvesThreshold,strCallNum);

  /*
   * General Memory Cleanup
   */
  Destroy_Search_DataSection(status->statusPtr,
			     &(outputCurves.curves),
			     outputCurves.numberOfCurves);
  CHECKSTATUSPTR(status);
  Destroy_Search_DataSection(status->statusPtr,
			     &(outputCurvesThreshold.curves),
			     outputCurvesThreshold.numberOfCurves);
  CHECKSTATUSPTR(status);
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
/* End Do_tracksearch function*/

/*
 * Diagnostic function to dump information to disk in Ascii form
 */
void Dump_Search_Data(
		      TSSearchParams     params,
		      TrackSearchOut     outCurve,
		      char*              strCallNum
      		      )
{
  INT4             i;
  INT4             j;
  FILE            *fp;
  FILE            *fp2;
  CHAR             appendname[128];
  CHAR             rawname[128];
  /* Intialize file IO naming scheme GPSstartTime-DateString.dat*/
  sprintf(appendname,"%s_%s.dat",params.auxlabel,strCallNum);
  fp = fopen(appendname,"w");
  sprintf(rawname,"%s_%s.raw",params.auxlabel,strCallNum);
  fp2= fopen(rawname,"w");
  /* Need more formal way to pack TrackSearchEvent Datatype for writing */
  /* to a sql of xml file for later post-processing */
  fprintf(fp,"Search Specific Information\n");
  fprintf(fp,"Contents of STRUCTURE: TSSearchParams\n\n");
  fprintf(fp,"GPS Start Time (seconds):      %d\n",params.GPSstart.gpsSeconds);
  fprintf(fp,"Total Data Length Points:      %d\n",params.TimeLengthPoints);
  fprintf(fp,"Segment lengths:               %d\n",params.SegLengthPoints);
  fprintf(fp,"Number of Segments Processed:  %s+1 of %d\n",strCallNum,params.NumSeg);
  fprintf(fp,"Data Sampling Rate:            %d\n",params.SamplingRate);
  fprintf(fp,"TF Transform Type used:        %i\n",params.TransformType);
  fprintf(fp,"Overlap Flag:                  %d\n",params.overlapFlag);
  /*  fprintf(fp,"Number of points in each FFT:  %d\n",params.fftPoints);*/
  fprintf(fp,"Line Width specified for map:  %d\n",params.LineWidth);
  fprintf(fp,"Ridge Start Threshold:         %e\n",params.StartThresh);
  fprintf(fp,"Ridge Member Threshold:        %e\n",params.LinePThresh);
  fprintf(fp,"Freq Bins for Maps:            %d\n",params.FreqBins);
  fprintf(fp,"Time Bins for Maps:            %d\n",params.TimeBins);
  fprintf(fp,"Window size for ffts:          %d\n",params.windowsize);
  fprintf(fp,"Window Type Used:              %i\n",params.window);
  fprintf(fp,"Channel:                       %s\n",params.channelName);
  fprintf(fp,"\n\n\n");

  fprintf(fp,"Total Number of Curves:%i\n",outCurve.numberOfCurves);
  fprintf(fp,"\n");
  for (i = 0;i < outCurve.numberOfCurves;i++)
    {
      fprintf(fp,"Curve      #:%i\n",i);
      fprintf(fp,"Points      :%i\n",outCurve.curves[i].n);
      fprintf(fp,"Junction    :%c\n",outCurve.curves[i].junction);
      fprintf(fp,"Total Power :%f\n",outCurve.curves[i].totalPower);
      fprintf(fp,"Freq      Time      Depth      Hz   gpsSec  gpsNanoSec \n");
      for (j = 0;j < outCurve.curves[i].n;j++)
	{
	  fprintf(fp,"%d    %d     %e     %f     %d   %d\n",
		  outCurve.curves[i].col[j],
		  outCurve.curves[i].row[j],
		  outCurve.curves[i].depth[j],
		  outCurve.curves[i].fBinHz[j],
		  outCurve.curves[i].gpsStamp[j].gpsSeconds,outCurve.curves[i].gpsStamp[j].gpsNanoSeconds);
	  fprintf(fp2,"%d    %d     %e\n",
		  outCurve.curves[i].col[j],
		  outCurve.curves[i].row[j],
		  outCurve.curves[i].depth[j]);
	}
      
      fprintf(fp,"\n\n");
    }
  fclose(fp);
  return;
}

/*
 * Dump candidates in short list Totals and Endpoints only
 */
void QuickDump_Data(
		    TSSearchParams     params,
		    TrackSearchOut     outCurve,
		    char*              strCallNum
		    )
{
  FILE            *fp;
  INT4             i;
  CHAR             appendname[128];

  sprintf(appendname,"Quick_%s_%s.dat",params.auxlabel,strCallNum);
  fp = fopen(appendname,"w");
  fprintf(fp,"%10s %10s %10s %10s %10s %10s %10s\n",
	  "CurveNum",
	  "Power",
	  "Length",
	  "StartF",
	  "FinishF",
	  "StartT",
	  "FinishT");
  for (i = 0;i < outCurve.numberOfCurves;i++)
    {
      fprintf(fp,"%10i %10.4f %10i %10i %10i %10i %10i\n",
	      i,
	      outCurve.curves[i].totalPower,
	      outCurve.curves[i].n,
	      outCurve.curves[i].col[0],
	      outCurve.curves[i].col[outCurve.curves[i].n - 1],
	      outCurve.curves[i].row[0],
	      outCurve.curves[i].row[outCurve.curves[i].n - 1]);
    }
  fclose(fp);
  return;
}

/*
 * This routine should be getting a seed to create the fake noise 
 * This is an inline routine not used anylonger
 */
void fakeDataGeneration(LALStatus              *status,
			REAL4TimeSeries        *fakeData,
			INT4                    pointnum,
			INT4                    seed
			)
{
  INT4               k;
  REAL4Vector       *fd = NULL;
  RandomParams      *RP = NULL;
  INITSTATUS (status, "makefakenoise", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  seed = 0;
  LALCreateVector(status->statusPtr,&fd,pointnum);
  LALCreateRandomParams(status->statusPtr,&RP,seed);
  LALNormalDeviates(status->statusPtr,fd,RP);
  LALDestroyRandomParams(status->statusPtr,&RP);
  fakeData->data = NULL;
  LALCreateVector(status->statusPtr,&(fakeData->data),pointnum);
  for (k = 0;k < ((INT4) fd->length);k++)
    {
      fakeData->data->data[k] = fd->data[k];
    }
  LALDestroyVector(status->statusPtr,&fd);
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/*
 * Local routine to allocate memory for the curve structure
 */
void Create_Curve_DataSection(LALStatus    *status,
			      Curve        **curvein)
{
  INT4    counter;
  Curve  *curve;
 
  INITSTATUS (status, "fixCurveStrucAllocation", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  *curvein = curve = LALMalloc(MAX_NUMBER_OF_CURVES * sizeof(Curve));
  for (counter = 0; counter < MAX_NUMBER_OF_CURVES;counter++)
    {
      curve[counter].n = 0;
      curve[counter].junction = 0;
      curve[counter].totalPower = 0;
      curve[counter].row = NULL;
      curve[counter].col = NULL;
      curve[counter].depth = NULL;
      curve[counter].row = LALMalloc(MAX_CURVE_LENGTH * sizeof(INT4));
      curve[counter].col = LALMalloc(MAX_CURVE_LENGTH * sizeof(INT4));
      curve[counter].depth = LALMalloc(MAX_CURVE_LENGTH * sizeof(REAL4));
      if (curve[counter].row == NULL)
	{
	  ABORT(status,TRACKSEARCHC_EMEM,TRACKSEARCHC_MSGEMEM);
	}
      if (curve[counter].col == NULL)
	{
	  ABORT(status,TRACKSEARCHC_EMEM,TRACKSEARCHC_MSGEMEM);
	}
      if (curve[counter].depth == NULL)
	{
	  ABORT(status,TRACKSEARCHC_EMEM,TRACKSEARCHC_MSGEMEM);
	}
      /* Need to gracdfully exit here by unallocating later */
    };
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*
 * Local routine to clean up memory once used to hold
 * collection of curve candidates
 */
void Destroy_Search_DataSection(LALStatus    *status,
				Curve        **curvein,
				INT4         numCurves)
{
  INT4    counter;
  Curve  *curve;

  INITSTATUS (status, "fixCurveStrucAllocation", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  curve = *curvein;

  for (counter = 0; counter < numCurves;counter++)
    {
      LALFree(curve[counter].fBinHz);
      LALFree(curve[counter].row);
      LALFree(curve[counter].col);
      LALFree(curve[counter].depth);
      LALFree(curve[counter].gpsStamp);
      /* Need to gracefully exit here by unallocating later */
      /* Just in case there is an error in deallocation */
    };
  if (numCurves > 0)
    {
      LALFree(curve);
    }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/* 
 * This is the frame reading routine taken from power.c
 * it has been modified very slightly
 */
void Get_Data(LALStatus*          status,
	      TSSearchParams*     params,
	      REAL4TimeSeries*    DataIn,
	      CHARVector*         dirname,
	      CHAR*               cachefile
	      )


{
  MetadataTable        procparams;
  MetadataTable procTable;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  ProcessParamsTable   *this_proc_param;
  MetadataTable         searchsumm;
  FrStream             *stream = NULL;
  FrCache              *frameCache = NULL;
  FrChanIn              channelIn;
  CHAR                 *comment=NULL;

  /* Set all variables from params structure here */
  channelIn.name = params->channelName;
  /* create the process and process params tables */
  procTable.processTable = (ProcessTable *) 
    LALCalloc( 1, sizeof(ProcessTable) );
  LAL_CALL( 
	   LALGPSTimeNow(status, 
			 &(procTable.processTable->start_time),
			 &accuracy ), 
	   status );
  LAL_CALL( 
	   populate_process_table( status, 
				   procTable.processTable, 
				   PROGRAM_NAME, 
				   CVS_REVISION, 
				   CVS_SOURCE,
				   CVS_DATE ),
	   status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    LALCalloc( 1, sizeof(ProcessParamsTable) );

  /* create the search summary table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    LALCalloc( 1, sizeof(SearchSummaryTable) );

  /* fill the comment, if a user has specified one, or leave it blank */
  if ( comment )
    {
      LALSnprintf( procTable.processTable->comment, 
		   LIGOMETA_COMMENT_MAX,
		   " " );
      LALSnprintf( searchsumm.searchSummaryTable->comment, 
		   LIGOMETA_COMMENT_MAX,
		   " " );    
    } 
  else 
    {
      LALSnprintf( procTable.processTable->comment, 
		   LIGOMETA_COMMENT_MAX,
		   "%s", 
		   comment );
      LALSnprintf( searchsumm.searchSummaryTable->comment,
		   LIGOMETA_COMMENT_MAX,
		   "%s",
		   comment );
    }

  /* the number of nodes for a standalone job is always 1 */
  searchsumm.searchSummaryTable->nnodes = 1;

  /* create and initialize the time series vector */
  DataIn->data = NULL;
  /*Bug HERE */
  LAL_CALL( 
	   LALCreateVector( status, 
			    &(DataIn->data), 
			    params->TimeLengthPoints), 
	   status);
  memset( DataIn->data->data, 0, DataIn->data->length*sizeof(REAL4) );
  /* Need to use NULL DataIN->name so we will skip this line of code
     Nullifying the DataIn->name pointer by hand hope it works */
  strcpy(DataIn->name, params->channelName);
  /*    DataIn->deltaT = 1.0/((REAL8) sampleRate); */
  DataIn->deltaT = 1.0;
  DataIn->f0 = 0.0;
  DataIn->sampleUnits = lalADCCountUnit;
  /* only try to load frame if name is specified */
  if (dirname->data || cachefile)
    {
      REAL8 tmpTime=0;
      if(dirname){
	/* Open frame stream */
	LAL_CALL( LALFrOpen( status, &stream, dirname->data, "*.gwf" ), status);
      }
      else if (cachefile){
	/* Open frame cache */
	LAL_CALL( LALFrCacheImport( status, &frameCache, cachefile ), status);
	LAL_CALL( LALFrCacheOpen( status, &stream, frameCache ), status);
	LAL_CALL( LALDestroyFrCache( status, &frameCache ), status );
      }
      /*
       * Determine information about the channel and seek to the
       * right place in the fram files 
       */
      DataIn->epoch.gpsSeconds     = params->GPSstart.gpsSeconds;
      DataIn->epoch.gpsNanoSeconds = params->GPSstart.gpsNanoSeconds;
      LAL_CALL( LALFrSeek(status, &(DataIn->epoch), stream), status);

      /* get the data */

      LAL_CALL( LALFrGetREAL4TimeSeries( status, 
					 DataIn, 
					 &channelIn, 
					 stream), 
		status);

      /* 
       * store the start and end time of the raw channel 
       * in the search summary 
       */
      searchsumm.searchSummaryTable->in_start_time = DataIn->epoch;
      LAL_CALL( LALGPStoFloat( status, &tmpTime, &(DataIn->epoch) ), 
		status );
      tmpTime += DataIn->deltaT * (REAL8) DataIn->data->length;
      LAL_CALL( LALFloatToGPS( status, 
			       &(searchsumm.searchSummaryTable->in_end_time), 
			       &tmpTime ), 
		status );

      /* close the frame stream */
      LAL_CALL( LALFrClose( status, &stream ), status);
    }
  /* Free memory allocation on procTable at start of routine */
  LALFree(procTable.processTable);
  LALFree(procparams.processParamsTable);
  LALFree(searchsumm.searchSummaryTable);
}

/*
 * Routine to allow use of acsii input rather than frames
 */
void Get_Ascii_Data(LALStatus*          status,
		    TSSearchParams*     params,
		    REAL4TimeSeries*    DataIn,
		    CHARVector*         dirname
		    )
{
  FILE                 *fp=NULL;
  INT4                  i;
  INITSTATUS (status, "readascii", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  /* File opening via an absolute path */
  fp = fopen(dirname->data,"r");
  if (!fp)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEREAD);
      exit(TRACKSEARCHC_EREAD);
    };
  DataIn->data=NULL;
  LALCreateVector(status->statusPtr,&(DataIn->data),params->TimeLengthPoints);
  for(i=0;i<(INT4)params->TimeLengthPoints;i++)
    {
      fscanf(fp,"%f\n",&(DataIn->data->data[i]));
    }
  fclose(fp);
  if (DataIn->data->length != params->TimeLengthPoints)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEREAD);
      exit(TRACKSEARCHC_EREAD);
      LALDestroyVector(status->statusPtr,&(DataIn->data));
    }
  DETATCHSTATUSPTR(status);
  return;
}



/*****************************************/
/* Scratch stuff */


/* WINDOW ENUM TYPES
   Rectangular,
   Hann,
   Welch,
   Bartlett,
   Parzen,
   Papoulis,
   Hamming,
   Kaiser,
   Creighton,
   The code will not know of any new window enum types!!!!!
*/
/* 
 * To Generate the fake incoming frame files under
 * ligotools dir is a utility called ascii2frame 
 */
