/*
 * Copyright (C) 2004, 2005 Cristina V. Torres
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
/*
 * Author: Torres Cristina 
 */

#include "tracksearch.h"
#include "lal/FrameStream.h"

#define PROGRAM_NAME "tracksearch"

typedef struct
{
  INT4    argc;
  CHAR**  argv;
}LALInitSearchParams;

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
#define TRACKSEARCHC_ENULL              128
#define TRACKSEARCHC_EDATA              256

#define TRACKSEARCHC_MSGENORM          "Normal Exit"
#define TRACKSEARCHC_MSGESUB           "Subroutine Fail"
#define TRACKSEARCHC_MSGEARGS          "Arguement Parse Error"
#define TRACKSEARCHC_MSGEVAL           "Invalid Argument(s)"
#define TRACKSEARCHC_MSGEFILE          "Error Opening File"
#define TRACKSEARCHC_MSGEREAD          "Error Reading File"
#define TRACKSEARCHC_MSGEMEM           "Memory Error"
#define TRACKSEARCHC_MSGEMISC          "Unknown Error"
#define TRACKSEARCHC_MSGENULL          "Unexpected NULL pointer"
#define TRACKSEARCHC_MSGEDATA          "Unknown data corruption error"

#define TRUE     1
#define FALSE    0

/* Usage format string. */
#define USAGE "Still in flux"

/* Code Body */
int main (int argc, char *argv[])
{
  /* Global Variable Declarations */
  /* Variables that will be in final produciton code */
  LALStatus           status;/* Error containing structure */
  TSSearchParams     *params;/*Large struct for all params*/
  TSappsInjectParams *injectParams;/*Struct for performing injects*/
  CHARVector         *cachefile=NULL;/* name of file with frame cache info */
  CHAR               *cachefileDataField=NULL;/*field of cachefile*/
  CHARVector         *dirpath=NULL;

  /*
   *Sleep for Attaching DDD 
   */
  unsigned int doze = 0;
  pid_t myPID;
  myPID = getpid( );
  fprintf( stdout, "pid %d sleeping for %d seconds\n", myPID, doze );
  fflush( stdout );
  sleep( doze );
  fprintf( stdout, "pid %d awake\n", myPID );
  fflush( stdout );

  /* End Global Variable Declarations */

  /* SET LAL DEBUG STUFF */
  /*set_debug_level("ERROR");*/
  /*set_debug_level("ERROR | WARNING | TRACE");*/
  /*set_debug_level("ERROR | WARNING | MEMDBG | TRACE");*/
  memset(&status, 0, sizeof(status));
  lal_errhandler = LAL_ERR_ABRT;
  lal_errhandler = LAL_ERR_DFLT;
  lal_errhandler = LAL_ERR_RTRN;

  /*
   * Initialize status structure 
   */
  
  params=XLALMalloc(sizeof(TSSearchParams));
  injectParams=XLALMalloc(sizeof(TSappsInjectParams));

  /* 
   * Parse incoming search arguments initialize variables 
   */
  LALappsTrackSearchInitialize(&status,
			       argc,
			       argv,
			       params,
			       injectParams,
			       &(cachefile),
			       &(dirpath));
  if (params->verbosity >= verbose)
    {
      fprintf(stdout,"Done initializing input structure.\n");
      fflush(stdout);
    }
  /*
   * Check params structure what analysis are we doing?
   * Do we look for tseries data files or map data files?
   */
  if (params->tSeriesAnalysis)
    {
      /*
       * We proceed as if this is a tSeries analysis
       */
      if (cachefile != NULL)
	cachefileDataField=cachefile->data;
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Doing a time series analysis.\n");
	  fflush(stdout);
	}
      LALappsDoTimeSeriesAnalysis(&status,
				  *params,
				  *injectParams,
				  cachefileDataField,
				  dirpath);
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Done performing a time series analysis.\n");
	  fflush(stdout);
	}
    }
  else
    {
      /*
       * Assume we have preformed maps to use and do a map
       * analysis instead
       */
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Performing a TF Map analysis.\n");
	  fflush(stdout);
	}
      LALappsDoTSAMapAnalysis(&status,*params);
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Done performing a TF Map analysis.\n");
	  fflush(stdout);
	}
    }
  /* 
   * Free searchParams input arguement structure from InitializeRoutine
   * Free elements first 
   */

  if (dirpath)
    XLALDestroyCHARVector(dirpath);
  if (cachefile)
    XLALDestroyCHARVector(cachefile);
  if (params->channelName)
    XLALFree(params->channelName);
  if (params->auxlabel)
    XLALFree(params->auxlabel);
  if (params->numSlaves)
    XLALFree(params->numSlaves);
  if (params->injectMapCache)
    XLALFree(params->injectMapCache);
  if (params->injectSingleMap)
    XLALFree(params->injectSingleMap);
  if (params)
    XLALFree(params);   
  if (injectParams->injectTxtFile)
    XLALDestroyCHARVector(injectParams->injectTxtFile);
  if (injectParams)
    XLALFree(injectParams);
  /* 
   * Done freeing for search params struct memory
   */
  /*LALCheckMemoryLeaks();*/
  return 0;
}
/* End Main Code Body */



/* ********************************************************************** */
/*
 * These are the SemiPrivate functions for LALapps tracksearch
 * Useful and permanent to local lalapps code
 *
 */

/*
 * Function that prepared data for analysis
 * Can do whitening and calibration on the segments
 * TINA -- Setup the removal of COHERENT LINES here!
 */
void LALappsTrackSearchPrepareData( LALStatus        *status,
				    REAL4TimeSeries  *dataSet,
				    REAL4TimeSeries  *injectSet,
				    TSSegmentVector  *dataSegments,
				    TSSearchParams    params)
/* Add option NULL or with data called REAL4TimeSeries injectSet */
{
  if (params.verbosity >= printFiles)
    {
      print_lalUnit(dataSet->sampleUnits,"Pre_DataCondEntireInputDataSet_Units.diag");
      print_real4tseries(dataSet,"Pre_DataCondEntireInputDataSet.diag");
    }
  /*
   * Calibrate entire input data set
   */
  if (params.calibrate)
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"Performing calibration on input data\n");
	  fflush(stdout);
	}

      if (params.verbosity >= printFiles)
	{
	  print_lalUnit(dataSet->sampleUnits,"Pre_DataCalibrationTimeDomain_Units.diag");
	  print_real4tseries(dataSet,"Pre_DataCalibrationTimeDomain.diag");
	}
      LALappsTrackSearchCalibrate(status,dataSet,params);
      if (params.verbosity >= printFiles)
	{
	  print_lalUnit(dataSet->sampleUnits,"Post_DataCalibrationTimeDomain_Units.diag");
	  print_real4tseries(dataSet,"Post_DataCalibrationTimeDomain.diag");
	}
    }
  else
    {
      if (params.verbosity >= verbose)
	fprintf(stdout,"Calibration not requested, assuming no calibration needed\n");
    }
  /* 
   *End calibration conditional 
   */

  if ((params.highPass > 0)||(params.lowPass > 0))
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"Bandpassing requested.\n");
	  fflush(stdout);
	}
      if (params.verbosity >= printFiles)
	print_real4tseries(dataSet,"Pre_ButterworthFiltered_TimeDomain.diag");
      LALappsTrackSearchBandPassing(status,dataSet,params);
      if (params.verbosity >= printFiles)
	print_real4tseries(dataSet,"Post_ButterworthFiltered_TimeDomain.diag");
    }
  else if (params.verbosity >= verbose)
    {
      fprintf(stdout,"No bandpassing requensted.\n");
      fflush(stdout);
    }

  /*
   * Perform injections when inject structure has data.
   */
  if (injectSet != NULL)
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"Making requested injections.\n");
	  fprintf(stderr,"If frame was REAL8 data converted to REAL4 data.\n");
	  fprintf(stderr,"Potential problem with possibly prefactored input data.\n");
	  fflush(stdout);
	}
      if (params.verbosity >= printFiles)
	print_real4tseries(dataSet,"Pre_SoftwareInjectDataSet.diag");
      LALappsTrackSearchPerformInjection(status,dataSet,injectSet,params);
      if (params.verbosity >= printFiles)
	print_real4tseries(dataSet,"Post_SoftwareInjectDataSet.diag");
    }

  /*
   * NOTE:
   * Tina, due to the comment in LSD about removing lines in data
   * stretches seconds to minutes we will now attempt line removal
   * on the individual segments.  For ASQ data and other high fs
   * data these stretches will be minutes and not tens of minutes in 
   * duration... This should work (I HOPE).
   */

  if (params.numLinesToRemove > 0)
    {
      if (1)
	/* 0 Do individual segment line removal */
	/* 1 Do entire data set line removal */
	{
	  if (params.verbosity >= verbose)
	    {
	      fprintf(stdout,
		      "Preparing to remove %i lines from entire input data at one time.\n",params.numLinesToRemove);
	      fflush(stdout);
	    }
	  if (params.verbosity >= printFiles)
	    print_real4tseries(dataSet,"Pre_LineRemoval_TimeDomain.diag");
	  LALappsTracksearchRemoveHarmonics(status,dataSet,params);
	  if (params.verbosity >= printFiles)
	    print_real4tseries(dataSet,"Post_LineRemoval_TimeDomain.diag");
	  /*
	   * Split incoming data into Segment Vector
	   * Adjust the segmenter to make orignal segments with
	   * SegBufferPoints on each side of the segments.
	   */
	  /*DO I CONVERT TO XLALish type?*/
	  LAL_CALL(LALTrackSearchDataSegmenter(status,
					       dataSet,
					       dataSegments,
					       params),
		   status);

	}
      else
	{
	  /*
	   * Split incoming data into Segment Vector
	   * Adjust the segmenter to make orignal segments with
	   * SegBufferPoints on each side of the segments.
	   */
	  /*DO I CONVERT TO XLALish type?*/
	  LAL_CALL(LALTrackSearchDataSegmenter(status,
					       dataSet,
					       dataSegments,
					       params),
		   status);

	  if (params.verbosity >= verbose)
	    fprintf(stdout,
		    "Preparing to remove %i lines from individual data segments.\n",
		    params.numLinesToRemove);
	  LALappsTracksearchRemoveHarmonicsFromSegments(status,
							dataSet,
							dataSegments,
							params);
	}
    }
  else 
    {
	  /* REMEMBER: If we don't remove line we still need to
	   * split the data!!!
	   * Split incoming data into Segment Vector
	   * Adjust the segmenter to make orignal segments with
	   * SegBufferPoints on each side of the segments.
	   */
	  LAL_CALL(LALTrackSearchDataSegmenter(status,
					       dataSet,
					       dataSegments,
					       params),
		   status);
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"Lines and Harmonics will NOT be removed.\n");
	  fflush(stdout);
	}
    }

  /*
   * If we are to whiten first let us calculate the 
   * average PSD using the non segment data structure
   */
  if (params.whiten !=0 )
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"Preparing to whiten data segments.\n");
	  fflush(stdout);
	}
      LALappsTrackSearchWhitenSegments(status,dataSet,dataSegments,params);
    }
  else if (params.verbosity >= verbose)
    {
      fprintf(stdout,"No data whitening will be done.\n");
      fflush(stdout);
    }

  if (params.verbosity >= verbose)
    {
      fprintf(stdout,"Done preparing input data.\n");
      fflush(stdout);
    }
  return;
}
  /* End the LALappsTrackSearchPrepareData routine */

/*
 * Setup params structure by parsing the command line
 */
void LALappsTrackSearchInitialize(
				  LALStatus          *status,
				  int                 argc,
				  char               *argv[], 
				  TSSearchParams     *params,
				  TSappsInjectParams *injectParams,
				  CHARVector        **dPtrCachefile,
				  CHARVector        **dPtrDirPath
				  )
{
  /*Local Variables*/
  LIGOTimeGPS     tempGPS;
  INT4 len;
  REAL4 optWidth=0;
  const char lineDelimiters[]=",";
  char *token;
  CHARVector *lineTokens=NULL;

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
      {"printPGM",            no_argument,        0,    'A'},
      {"color_map",           required_argument,  0,    'B'},
      {"inject_map_cache",    required_argument,  0,    'C'},
      {"inject_map",          required_argument,  0,    'D'},
      {"sample_rate",         required_argument,  0,    'E'},
      {"smooth_average_spectrum", required_argument,    0,    'F'},
      {"filter_high_pass",    required_argument,  0,    'G'},
      {"filter_low_pass",     required_argument,  0,    'H'},
      {"verbosity",           required_argument,  0,    'I'},
      {"inject_offset",       required_argument,  0,    'J'},
      {"inject_count",        required_argument,  0,    'K'},
      {"inject_space",        required_argument,  0,    'L'},
      {"inject_file",         required_argument,  0,    'M'},
      {"inject_scale",        required_argument,  0,    'N'},
      {"channel_type",        required_argument,  0,    'O'},
      {"bin_buffer",          required_argument,  0,    'P'},
      {"zcontrast",           required_argument,  0,    'Q'},
      {"zlength",             required_argument,  0,    'R'},
      {"remove_line",         required_argument,  0,    'S'},
      {"max_harmonics",       required_argument,  0,    'T'},
      {"snr",                 required_argument,  0,    'U'},
      {"heterodyne_frequency",required_argument,  0,    'V'},
      {"heterodyne_sample_rate",required_argument, 0,   'W'},
      {0,                     0,                  0,      0}
    };
  
  int              C;
  /*
   * Set default values for optional arguments 
   * These options will be used if omitted from the command line 
   * Required arguements will be initialized as well to simply debugging
   */
  params->tSeriesAnalysis=1;/*Default Yes this is tseries to analyze*/
  params->LineWidth = 0;
  params->TransformType = Spectrogram;
  params->NumSeg = 1;
  params->SegLengthPoints = 0; /* Set for debuging this with overlap # */
  params->SegBufferPoints =0; /* Set to zero to test software */
  params->colsToClip=0;/*Number of cols to clip off ends of TFR */
  params->overlapFlag = 0; /* Change this to an INT4 for overlapping */
  params->TimeLengthPoints = 0;
  params->discardTLP = 0;
  params->window = Hann; /*Welch*/
  params->whiten = 0;
  params->avgSpecMethod = -1; /*Will trigger error*/
  params->avgSpecWindow = Rectangular; /* Default no window */
  params->smoothAvgPSD = 0;/*0 means no smoothing*/
  params->highPass=-1;/*High pass filter freq*/
  params->lowPass=-1;/*Low pass filter freq*/
  params->FreqBins = 0;
  params->TimeBins = 0;
  params->windowsize = 0;/*256*/ /*30*/
  params->StartThresh =-1; /* We choose invalid values to check incoming args*/
  params->LinePThresh =-1; /* We choose invalid values to check incoming args*/
  params->MinPower = 0;
  params->MinLength = 0;
  params->MinSNR = 0;
  params->GPSstart.gpsSeconds = 0;
  params->GPSstart.gpsNanoSeconds = 0;
  /* 
   * In future...
   * Must include error checking for these values 
   */
  params->TimeLengthPoints = 0;
  params->SamplingRate = 1;
  params->HeterodyneFrequency=0;
  params->HeterodyneSamplingRate=0;
  params->SamplingRateOriginal = params->SamplingRate;
  params->makenoise = -1;
  params->calChannelType=SimDataChannel;/*See lsd*/
  params->channelName=NULL; /* Leave NULL to avoid frame read problems */
  params->channelNameType=LAL_ADC_CHAN; /*Default for data analysis */
  params->dataDirPath=NULL;
  params->singleDataCache=NULL;
  params->detectorPSDCache=NULL;
  params->channelNamePSD=NULL;
  params->auxlabel=NULL;;
  params->joinCurves=FALSE; /*default no curve joining */
  params->dataSegVec=NULL;
  params->numSlaves=0;
  params->calFrameCache=NULL;
  params->printPGM=0; /*Don't output PGMs of map*/
  params->pgmColorMapFile=NULL;
  params->verbosity=quiet;/*printFiles;*/
  params->calibrate=0;/*No calibration default */
  params->calCatalog=NULL;
  strcpy(params->calibrateIFO,"L1");
  params->injectMapCache=NULL;
  params->injectSingleMap=NULL;
  params->autoThresh=-1;/*Default no automatic lambda runs*/
  params->relaThresh=-1;/*Default no automatic lambda runs*/
  params->autoLambda=0;/*False assume no automatic*/
  /*params.listLinesToRemove*/
  params->numLinesToRemove=0;/*Num of coherent lines to remove */
  params->maxHarmonics=1;/*Num of harmonics to try less than fnyq */
  /* Inject params defaults */
  injectParams->startTimeOffset=0;
  injectParams->numOfInjects=0;
  injectParams->injectSpace=0;
  injectParams->scaleFactor=0;
  injectParams->injectTxtFile=NULL;
  injectParams->sampleRate=params->SamplingRate;
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
			   "a:b:c:d:e:f:h:i:j:k:l:m:o:p:q:r:s:t:u:v:w:x:y:z:A:B:C:D:E:F:G:H",
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
	    /* Allow intepreting float values */
	    REAL8 startTime = atof(optarg);
	    /*XLALFloatToGPS(&tempGPS,startTime);*/
	    LAL_CALL(LALFloatToGPS(status,
				   &tempGPS,
				   &startTime),
		     status);
	    params->GPSstart.gpsSeconds = tempGPS.gpsSeconds;
	    params->GPSstart.gpsNanoSeconds = tempGPS.gpsNanoSeconds;
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
	    name = (CHAR*) XLALMalloc(strlen(optarg)+1);
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
	  }
	  break;
	  
	case 'e':
	  /* Setting Start Threshold parameter */
	  {
	    params->StartThresh = atof(optarg);
	  }
	  break;
	  
	case 'f':
	  /* Setting member point threshold paramter */
	  {
	    params->LinePThresh = atof(optarg);
	  }
	  break;
	  
	case 'g':
	  /* Setting Thresholding on line length for search */
	  {
	    params->MinLength = atoi(optarg);
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
	    params->channelName = (CHAR*) XLALMalloc(strlen(optarg)+1);
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
	    fprintf(stderr,"Option not available:Can not specify number of maps\n");
	    fprintf(stderr,TRACKSEARCHC_MSGEMEM);
	    exit(TRACKSEARCHC_EMEM);
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
	    if (!strcmp(optarg,"Rectangular"))
	      params->window=Rectangular;
	    else if (!strcmp(optarg,"Hann"))
	      params->window=Hann;
	    else if (!strcmp(optarg,"Welch"))
	      params->window=Welch;
	    else if (!strcmp(optarg,"Bartlett"))
	      params->window=Bartlett;
	    else if (!strcmp(optarg,"Parzen"))
	      params->window=Parzen;
	    else if (!strcmp(optarg,"Papoulis"))
	      params->window=Papoulis;
	    else if (!strcmp(optarg,"Hamming"))
	      params->window=Hamming;
	    else if (!strcmp(optarg,"Kaiser"))
	      params->window=Kaiser;
	    else if (!strcmp(optarg,"Creighton"))
	      params->window=Creighton;
	    else
	      {
		fprintf(stderr,"Invalid window option: using Rectangular window");
		params->window=Rectangular;
	      };
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
	  }
	  break;

	case 'p':
	  /* Choose transform type */
	  {
	    if (strcmp(optarg,"Spectrogram")==0)
	      {
		params->TransformType=Spectrogram;
	      }
	    else if (strcmp(optarg,"RSpectrogram")==0)
	      {
		params->TransformType=RSpectrogram;
	      }
	    else if (strcmp(optarg,"WignerVille")==0)
	      {
		params->TransformType=WignerVille;
	      }
	    else if (strcmp(optarg,"PSWignerVille")==0)
	      {
		params->TransformType=PSWignerVille;
	      }
	    else
	      {
		fprintf(stdout,"\n Invald TF Transform selected defaulting to Spectrogram!\n");
		fflush(stdout);
		params->TransformType=Spectrogram;
	      };
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
	    params->auxlabel = (CHAR*) XLALMalloc(strlen(optarg)+1);
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
	    /*Setup structure for path or cachefile name temp*/
	    *(dPtrDirPath)=XLALCreateCHARVector(maxFilenameLength);
	    len= strlen(optarg) +1;
	    strncpy((*dPtrDirPath)->data,optarg,len);
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
	    *dPtrCachefile=XLALCreateCHARVector(maxFilenameLength);
	    len = strlen(optarg) +1;
	    memcpy((*dPtrCachefile)->data,optarg,len);
	  }
	  break;

	case 'w':
	  {
	    params->joinCurves=TRUE;
	  }
	  break;

	case 'x':
	  {
	    len=strlen(optarg)+1;
	    params->detectorPSDCache=(CHAR *)
	      XLALCalloc(len,sizeof(CHAR));
	    memcpy(params->detectorPSDCache,optarg,len);
	  }
	  break;
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
	  }
	  break;
	case 'z':
	  { /* Setting up aux label if present */
	    /* Insert string copy code here */
	    len = strlen(optarg) +1;
	    params->calFrameCache = (CHAR *) XLALCalloc(len,sizeof(CHAR));
	    memcpy(params->calFrameCache,optarg,len);
	    params->calibrate=1;
	  }
	  break;
	case 'A':
	  {
	    params->printPGM=1;/*Enable PGM printing*/
	  }
	  break;
	case 'B':
	  {
	    len = strlen(optarg) +1;
	    params->pgmColorMapFile = (CHAR *) XLALCalloc(len,sizeof(CHAR));
	    memcpy(params->pgmColorMapFile,optarg,len);
	  }
	case 'C':
	  {
	    params->injectMapCache=
	      (CHAR*) XLALCalloc(strlen(optarg)+1,sizeof(CHAR));
	    memcpy(params->injectMapCache,optarg,strlen(optarg)+1);
	    /*
	     * Turn off tseries flag
	     */
	    params->tSeriesAnalysis=0;
	  }
	  break;

	case 'D':
	  {
	    params->injectSingleMap=
	      (CHAR*) XLALCalloc(strlen(optarg)+1,sizeof(CHAR));
	    memcpy(params->injectSingleMap,optarg,strlen(optarg)+1);
	    /*
	     * Turn off tseries flag
	     */
	    params->tSeriesAnalysis=0;
	  }
	  break;

	case 'E':
	  {
	    params->SamplingRate=atof(optarg);
	    params->SamplingRateOriginal=atof(optarg);
	    injectParams->sampleRate=params->SamplingRate;
	  }
	  break;

	case 'F':
	  {
	    params->smoothAvgPSD =(UINT4) atoi(optarg); 
	    /* 
	     * >0 Means do running median with this block size
	     */
	  }
	  break;

	case 'G':
	  {
	    params->highPass=atof(optarg);
	  }
	  break;

	case 'H':
	  {
	    params->lowPass=atof(optarg);
	  }
	  break;

	case 'I':
	  {
	    if(!strcmp(optarg,"quiet"))
	      params->verbosity=quiet;
	    else if(!strcmp(optarg,"verbose"))
	      params->verbosity=verbose;
	    else if(!strcmp(optarg,"printFiles"))
	      params->verbosity=printFiles;
	    else if(!strcmp(optarg,"all"))
	      params->verbosity=all;
	    else
	      {
		fprintf(stderr,"Invalid option:verbosity: assuming quiet\n");
		params->verbosity=quiet;
	      };
	  }
	  break;
	
	case 'J':
	  {
	    injectParams->startTimeOffset=atof(optarg);
	  }
	  break;

	case 'K':
	  {
	    injectParams->numOfInjects=atoi(optarg);
	  }
	  break;

	case 'L':
	  {
	    injectParams->injectSpace=atof(optarg);
	  }
	  break;

	case 'M':
	  {
	    injectParams->injectTxtFile=XLALCreateCHARVector(strlen(optarg)+1);
	    memcpy(injectParams->injectTxtFile->data,optarg,strlen(optarg)+1);
	  }
	  break;

	case 'N':
	  {
	    injectParams->scaleFactor=atof(optarg);
	  }
	  break;

	case 'O':
	  {
	    if(!strcmp(optarg,"LAL_ADC_CHAN"))
	      params->channelNameType=LAL_ADC_CHAN;
	    else if(!strcmp(optarg,"LAL_SIM_CHAN"))
	      params->channelNameType=LAL_SIM_CHAN;
	    else if(!strcmp(optarg,"LAL_PROC_CHAN"))
	      params->channelNameType=LAL_PROC_CHAN;
	    else
	      {
		fprintf(stderr,"Invalid channel type specified assuming LAL_ADC_CHAN\n");
		params->channelNameType=LAL_ADC_CHAN;
	      };
	  }
	  break;

	case 'P':
	  {
	    params->colsToClip=atoi(optarg);
	  }
	  break;

	case 'Q':
	  {
	    params->autoThresh=atof(optarg);
	    params->autoLambda=1;
	  }
	  break;

	case 'R':
	  {
	    params->relaThresh=atof(optarg);
	  }
	  break;

	case 'S':
	  {
	    lineTokens=XLALCreateCHARVector(maxFilenameLength);
	    len = strlen(optarg) +1;
	    strncpy(lineTokens->data,optarg,len);
	    token=strtok(lineTokens->data,lineDelimiters);
	    while (token!=NULL)
	      {
		params->numLinesToRemove=params->numLinesToRemove+1;
		params->listLinesToRemove[params->numLinesToRemove-1]=atof(token);
		token=strtok(NULL,lineDelimiters);
	      }
	    if (lineTokens)
	      XLALDestroyCHARVector(lineTokens);
	  }
	  break;

	case 'T':
	  {
	    params->maxHarmonics=atoi(optarg);
	  }
	  break;

	case 'U':
	  {
	    params->MinSNR=atof(optarg);
	  }
	  break;

	case 'V':
	  {
	    params->HeterodyneFrequency=atof(optarg);
	  }
	  break;

	case 'W':
	  {
	    params->HeterodyneSamplingRate=atof(optarg);
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
   * moving error catching code out of the case statement in some cases
   */
  /* 
   * Tabulate the number of segments to create for specified overlap
   * Check for valid overlap number which < segment length number 
   */
  /* 
   *If doing a MAP analysis skip these checks these options should
   *not be invoked
   */
  if ((params->HeterodyneFrequency > 0) || (params->HeterodyneSamplingRate > 0))
    {
      fprintf(stdout,"WARNING: All measure such as segment length are assumed relative to the effective sampling rate of the heterodyned data.\n");
      fflush(stdout);
    }
  if ( params->tSeriesAnalysis == 1)
    {
      fprintf(stderr,"Checking tSeries parameters.\n");
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
	  params->discardTLP=(params->TimeLengthPoints)%(params->SegLengthPoints);
	  params->TimeLengthPoints=params->TimeLengthPoints-params->discardTLP;
	}
      else
	{
	  /* Determine the number of maps */
	  params->NumSeg=floor((params->TimeLengthPoints-params->overlapFlag)
			       /(params->SegLengthPoints - params->overlapFlag));
	  params->discardTLP=(params->TimeLengthPoints-params->overlapFlag)%(params->SegLengthPoints - params->overlapFlag);
	  /* 
	   * Reset points to process by N-discardTLP, this gives us
	   * uniform segments to make maps with
	   */
	  params->TimeLengthPoints=params->TimeLengthPoints-params->discardTLP;
	}
      if (params->verbosity > quiet)
	{
	  fprintf(stdout,
		  "Need to trim %i subsegment lengths.\n",
		  params->NumSeg);
	  fprintf(stdout,
		  "Total data points requested :%i\n",
		  params->TimeLengthPoints+params->discardTLP);
	  fprintf(stdout,
		  "Subsegment length in data points requested :%i\n",
		  params->SegLengthPoints);
	  fprintf(stdout,
		  "Segment overlap in points :%i\n",
		  params->overlapFlag);
	  fprintf(stdout,
		  "Points to discard from total data requested : %i\n",
		  params->discardTLP);
	  fprintf(stdout,
		  "Total data points that will be processed :%i\n",
		      params->TimeLengthPoints);
	  fflush(stdout);
	}
      /* 
       * Determine the number of additional points required to negate
       * edge effects in first and last TFR
       */
      params->SegBufferPoints=((params->TimeLengthPoints/params->TimeBins)*params->colsToClip);
      params->SegBufferPoints=((params->SegLengthPoints/params->TimeBins)*params->colsToClip);
      if ((params->verbosity >= verbose))
	{
	  fprintf(stdout,"As a result of bin buffering requests an additional %i points are needed to negate TFR edge effects, %i points before and then after the data set to analyze.\n",2*params->SegBufferPoints,params->SegBufferPoints);
	  fflush(stdout);
	}
      if (params->whiten > 2 )
	{
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  exit(TRACKSEARCHC_EARGS);
	}
      if (params->SamplingRate < 0 )
	{
	  fprintf(stderr,"Requested Sampling rate less than zero!\n");
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  fflush(stderr);
	  exit(TRACKSEARCHC_EARGS);
	}
      if ((params->HeterodyneFrequency < 0 ) || (2*params->HeterodyneFrequency > params->SamplingRate))
	{
	  fprintf(stderr,"Heterodyne frequency invalid!\n");
	  fprintf(stderr,"Heterodyne frequency %9.2f Hz must lie between 0 and %9.3f.\n",params->HeterodyneFrequency,0.5*params->SamplingRate);
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  fflush(stderr);
	  exit(TRACKSEARCHC_EARGS);
	}
      if (2*params->HeterodyneSamplingRate > params->SamplingRate)
	{
	  fprintf(stderr,"WARNING: New Heterodyne sampling rate should never exceed half of the frame data sampling rate of %9.3f.\n",params->SamplingRate);
	  fprintf(stderr,"The Heterodyned time series will have some undesired artifacts for frequencies above %9.3f and the frequency axis of the time-frequency plot will be wrong above this frequency.\n",(params->SamplingRate)/2); 
	  fflush(stderr);
	}
    }  /*
   * Do following checks
   */
  /* If not auto lambda choosen then */
  if (params->autoLambda == 0)
    {
      if (params->StartThresh < params->LinePThresh)
	{
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  exit(TRACKSEARCHC_EARGS);
	}
      if (params->StartThresh < params->LinePThresh)
	{
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  exit(TRACKSEARCHC_EARGS);
	}
    }
  else
    { /*Assuming we need to autocalibrate here*/
      if ((params->relaThresh < 0) || (params->relaThresh > 1))
	{
	  params->autoLambda=0;
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  exit(TRACKSEARCHC_EARGS);
	}
    }

  if (params->MinLength < 3)
    {
      fprintf(stderr,"Minimum length threshold invalid!\n");
      fprintf(stderr,TRACKSEARCHC_MSGEARGS);
      exit(TRACKSEARCHC_EARGS);
    }
  if (/*Convolution kernel sigma must \geq 1*/
      params->LineWidth/(2*sqrt(3)) < 1 
      )
    {
      fprintf(stderr,"Sigma value too small to be useful try increasing line width to at least %i\n",
	      (INT4) ceil(sqrt(3)*2.0));
      fprintf(stderr,TRACKSEARCHC_MSGEARGS);
      exit(TRACKSEARCHC_EARGS);
    }
      
  if (/* Minimum line width to use: sigma \geq (line_width/(2*\sqrt(3))) */
      ((params->LineWidth) < ceil((2.0/MASK_SIZE)*(2.0*sqrt(3))))
      )
    {
      optWidth=0;
      optWidth=ceil((2.0/MASK_SIZE)*(2.0*sqrt(3)));
      fprintf(stderr,"Line width inappropriate try mimimum of: %i\n",
	      (INT4) optWidth);
      fprintf(stderr,TRACKSEARCHC_MSGEARGS);
      exit(TRACKSEARCHC_EARGS);
    }
  return;
}
/* End initialization subroutine for search */

/* 
 * This is the frame reading routine taken from power.c
 * it has been modified very slightly
 */
void LALappsGetFrameData(LALStatus*          status,
			 TSSearchParams*     params,
			 REAL4TimeSeries*    DataIn,
			 CHARVector*         dirname,
			 CHAR*               cachefile
			 )
{
  FrStream             *stream = NULL;
  FrCache              *frameCache = NULL;
  FrChanIn              channelIn;
  LALTYPECODE           dataTypeCode=0;
  LIGOTimeGPS           bufferedDataStartGPS;
  REAL8                 bufferedDataStart=0;
  LIGOTimeGPS           bufferedDataStopGPS;
  REAL8                 bufferedDataStop=0;
  REAL8                 bufferedDataTimeInterval=0;
  UINT4                 errcode=0;
  UINT4                 index=0;
  /* Set all variables from params structure here */
  channelIn.name = params->channelName;
  channelIn.type = params->channelNameType;
  if(dirname)
    {
      /* Open frame stream */
      stream=XLALFrOpen(dirname->data,"*.gwf");
    }
  else if (cachefile)
    {
      /* Open frame cache */
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Importing cache file: %s\n",cachefile);
	  fflush(stdout);
	}
      lal_errhandler = LAL_ERR_EXIT;
      frameCache=XLALFrImportCache(cachefile);
      if (frameCache == NULL)
	{
	  fprintf(stderr,"Error importing frame cache file!\n");
	  fprintf(stderr,"%s\n",TRACKSEARCHC_MSGEDATA);
	  fflush(stderr);
	  exit(TRACKSEARCHC_EDATA);
	}
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Creating data stream.\n");
	  fflush(stdout);
	}
      stream=XLALFrCacheOpen(frameCache);	  
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Clearing imported cache file.\n");
	  fflush(stdout);
	}
      XLALFrDestroyCache(frameCache);
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Data stream ready.\n");
	  fflush(stdout);
	}
    }
  else 
    {
      fprintf(stderr,"Cache file or path to frame files not specified!\n");
      fprintf(stderr,"%s\n",TRACKSEARCHC_MSGEVAL);
      exit(TRACKSEARCHC_EVAL);
    }
  lal_errhandler = LAL_ERR_EXIT;
  /* Set verbosity of stream so user sees frame read problems! */
  XLALFrSetMode(stream,LAL_FR_VERBOSE_MODE);
  /*DataIn->epoch SHOULD and MUST equal params->startGPS - (params.SegBufferPoints/params.SamplingRate)*/
  memcpy(&bufferedDataStartGPS,&(DataIn->epoch),sizeof(LIGOTimeGPS));
  LAL_CALL(LALGPStoFloat(status,&bufferedDataStart,&bufferedDataStartGPS),
	   status);
  /* 
   * Seek to end of requested data makes sure that all stream is complete!
   */
  bufferedDataStop=bufferedDataStart+(DataIn->data->length * DataIn->deltaT);
  LAL_CALL(LALFloatToGPS(status,
			 &bufferedDataStopGPS,
			 &bufferedDataStop),
	   status);
  bufferedDataTimeInterval=bufferedDataStop-bufferedDataStart;
  if (params->verbosity >= verbose)
    {
      fprintf(stderr,"Checking frame stream spans requested data interval, including the appropriate data buffering!\n");
      fprintf(stderr,"Start           : %f\n",bufferedDataStart);
      fprintf(stderr,"Stop            : %f\n",bufferedDataStop);
      fprintf(stderr,"Interval length : %f\n",bufferedDataTimeInterval);
    }

  errcode=XLALFrSeek(stream,&bufferedDataStopGPS);
  if (errcode!=0)
    {
      fprintf(stderr,"Error seeking stream to %f\n",bufferedDataStart);
      fflush(stderr);
      exit(errcode);
    }
  errcode=XLALFrSeek(stream,&bufferedDataStartGPS);
  if (errcode!=0)
    {
      fprintf(stderr,"Error seeking stream to %f\n",bufferedDataStop);
      fflush(stderr);
      exit(errcode);
    }
  /*
   * Determine the variable type of data in the frame file.
   */
  dataTypeCode=XLALFrGetTimeSeriesType(channelIn.name, stream);
  if (params->verbosity >= verbose)
    {
      fprintf(stdout,"Checking data stream variable type for :%s \n",channelIn.name);
      fflush(stdout);
    }

  switch (dataTypeCode)
    {
    case LAL_D_TYPE_CODE:
      {
	if (params->verbosity >= verbose)
	  {
	    fprintf(stdout,"Converting frame stream data from REAL8.\n");
	    fflush(stdout);
	  }
	errcode=XLALFrGetREAL8FrameConvertToREAL4TimeSeries(DataIn,stream);
	if (errcode!=0)
	  {
	    fprintf(stderr,"Conversion of stream data FAILED!\n");
	    fflush(stderr);
	  }
      }
      break;

    case LAL_S_TYPE_CODE:
      {
	if (params->verbosity >= verbose)
	  {
	    fprintf(stdout,"Converting frame stream data from REAL4.\n");
	    fflush(stdout);
	  }
	errcode=XLALFrGetREAL4FrameConvertToREAL4TimeSeries(DataIn,stream);
	if (errcode!=0)
	  {
	    fprintf(stderr,"Conversion of stream data FAILED!\n");
	    fflush(stderr);
	  }
      }
      break;

    case LAL_I2_TYPE_CODE:
      {
	if (params->verbosity >= verbose)
	  {
	    fprintf(stdout,"Converting frame stream data from INT2.\n");
	    fflush(stdout);
	  }
	errcode=XLALFrGetINT2FrameConvertToREAL4TimeSeries(DataIn,stream);
	if (errcode!=0)
	  {
	    fprintf(stderr,"Conversion of stream data FAILED!\n");
	    fflush(stderr);
	  }
      }
      break;

    case LAL_I4_TYPE_CODE:
      {
	if (params->verbosity >= verbose)
	  {
	    fprintf(stdout,"Converting frame stream data from INT4.\n");
	    fflush(stdout);
	  }
	errcode=XLALFrGetINT4FrameConvertToREAL4TimeSeries(DataIn,stream);
	if (errcode!=0)
	  {
	    fprintf(stderr,"Conversion of stream data FAILED!\n");
	    fflush(stderr);
	  }
      }
      break;
	  
    case LAL_I8_TYPE_CODE:
      {
	if (params->verbosity >= verbose)
	  {
	    fprintf(stdout,"Converting frame stream data from INT8.\n");
	    fflush(stdout);
	  }
	errcode=XLALFrGetINT8FrameConvertToREAL4TimeSeries(DataIn,stream);
	if (errcode!=0)
	  {
	    fprintf(stderr,"Conversion of stream data FAILED!\n");
	    fflush(stderr);
	  }
      }
      break;

    default:
      {
	fprintf(stderr,"Data type code found can't be loaded.\n");
	fprintf(stderr,"Data type code value %i\n.",dataTypeCode);
	fprintf(stderr,"%s\n",TRACKSEARCHC_MSGEREAD);
	fflush(stderr);
	exit(TRACKSEARCHC_EREAD);
      }
    }
  if (stream)
    {
      /*Close the frame stream if found open*/
      errcode=XLALFrClose(stream);
      if (errcode != 0)
	{
	  fprintf(stderr,"Problem trying to close the frame stream!\n");
	  exit(errcode);
	}
    }
  /*
   * Heterodyne the signal if requested by user.
   */
  if (params->HeterodyneFrequency > 0)
    {
      if (params->verbosity >= verbose)
	{
	  fprintf(stdout,"Heterodying with reference signal %9.3f.\n",params->HeterodyneFrequency);
	  fprintf(stdout,"Preparing to heterodyne the loaded time series.\n");
	  fflush(stdout);
	}
      /*H-dyne*/
      errcode=LALappsQuickHeterodyneTimeSeries(DataIn,params->HeterodyneFrequency);
      /*
       * Resample the data to the relative rate specified.
       */
      if ((1/DataIn->deltaT) > params->HeterodyneSamplingRate)
	{
	  if (params->verbosity >= verbose)
	    {
	      fprintf(stdout,"A resampling of frame data is required.\n");
	      fprintf(stdout,"Current sample rate           :%9.3f\n",1/DataIn->deltaT);
	      fprintf(stdout,"Heterodyned down to frequency :%9.3f\n",params->HeterodyneFrequency);
	      fprintf(stdout,"New effective sample rate     :%9.3f\n",params->HeterodyneSamplingRate);
	      fprintf(stdout,"Will resample to requested rate.\n");
	      fflush(stdout);
	    }
	  errcode=XLALResampleREAL4TimeSeries(DataIn,(1/params->HeterodyneSamplingRate));
	  params->SamplingRateOriginal=params->SamplingRate;
	  params->SamplingRate=params->HeterodyneSamplingRate;

	  if (DataIn->data->length != params->TimeLengthPoints)
	    {
	      fprintf(stdout,"Warning some REAL4TimeSeries resizeing required.\n");
	      fprintf(stdout,"Vector has %i points but should actually have %i.\n",
		      DataIn->data->length,
		      params->TimeLengthPoints);
	      fprintf(stdout,"Zero padding where appropriate.\n");
	      DataIn=XLALResizeREAL4TimeSeries(DataIn,0,params->TimeLengthPoints);
	      for (index = DataIn->data->length+1; 
		   index < params->TimeLengthPoints;
		   index++)
		DataIn->data->data[index]=0;
	      fprintf(stdout,"Resizing completed.\n");
	    }
	}
      /*
       * Dump the heterodyned and resampled data to disk if printFiles
       */
      if (params->verbosity >= printFiles)
	{
	  print_real4tseries(DataIn,"HeterodynedAndResampleTimeSeriesFrameData.diag");
	  print_lalUnit(DataIn->sampleUnits,"HeterodynedAndResampleTimeSeriesFrameData_Units.diag");
	}
    }
}
  /* End frame reading code */


/*
 * Routine to allow use of acsii input rather than frames
 */
void LALappsGetAsciiData(LALStatus*          status,
			 TSSearchParams*     params,
			 REAL4TimeSeries*    DataIn,
			 CHARVector*         dirname
			 )
{
  FILE                 *fp=NULL;
  INT4                  i;

  /* File opening via an absolute path */
  fp = fopen(dirname->data,"r");
  if (!fp)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEREAD);
      exit(TRACKSEARCHC_EREAD);
    };
  for(i=0;i<(INT4)params->TimeLengthPoints;i++)
    {
      fscanf(fp,"%f\n",&(DataIn->data->data[i]));
    }
  fclose(fp);
  if (DataIn->data->length != params->TimeLengthPoints)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEREAD);
      exit(TRACKSEARCHC_EREAD);
    }
  return;
}
/* End of Ascii reader function */


/*
 * Local routine to actually carry out search
 */
void LALappsDoTrackSearch(
			  LALStatus                     *status,
			  TimeFreqRep                   *tfmap,
			  TrackSearchParams              tsInputs,
			  TrackSearchMapMarkingParams    tsMarkers,
			  TSSearchParams                 params)
{
  TrackSearchOut         outputCurves;/*Native LAL format for results*/
  TrackSearchOut         outputCurvesThreshold; /*Curve data thresholded */
  CHARVector            *outputFilename=NULL;
  CHARVector            *outputFilenameMask=NULL;
  CHARVector            *outputCandidateFilename=NULL;
  int                    errCode=0;
  /*************************************************************/
  /* 
   * The LALTracksearch seems to map the 
   * width to tCol and the height to fRow
   */
  outputFilename=XLALCreateCHARVector(maxFilenameLength);
  outputFilenameMask=XLALCreateCHARVector(maxFilenameLength);
  /*
   * Setup file mask for output filenames
   */
  sprintf(outputFilenameMask->data,"CandidateList_Start_%i_%i_Stop_%i_%i",
	  tsMarkers.mapStartGPS.gpsSeconds,
	  tsMarkers.mapStartGPS.gpsNanoSeconds,
	  tsMarkers.mapStopGPS.gpsSeconds,
	  tsMarkers.mapStopGPS.gpsNanoSeconds);
  outputCurves.curves = NULL;
  outputCurves.numberOfCurves=0;
  outputCurves.minPowerCut=params.MinPower;
  outputCurves.minLengthCut=params.MinLength;
  outputCurves.minSNRCut=params.MinSNR;
  outputCurves.startThreshCut=tsInputs.high;
  outputCurves.linePThreshCut=tsInputs.low;

  /*
   * The power thresholding is not done in the LALroutine
   * This information is forwarded to this function for a post processing 
   * Should be default be givin minimum curve length parameter 
   * We want to apply thresholds in a seperate routine 
   */
  /*
   * DO THE AUTO ADJUSTMENTS!!!
   */
  if (params.autoLambda)
    {
      /* Do the calculate of Lh given parameters */
      lal_errhandler = LAL_ERR_RTRN;
      /*      errCode = LAL_CALL( LALTracksearchFindLambdaMean(status,*tfmap,&params),
	      status);*/
      errCode = LAL_CALL( LALTracksearchFindLambdaMedian(status,*tfmap,&params),
			  status);

      if ( errCode != 0 )
	{
	  fprintf(stderr,"Error calling automagic Lambda selection routine.\n");
	  fprintf(stderr,"%s\n",status->statusDescription);
	  fflush(stderr);
	}
      tsInputs.high=params.StartThresh;
      tsInputs.low=params.LinePThresh;
      /*Reset to show the auto selected values.*/
      outputCurves.startThreshCut=tsInputs.high;
      outputCurves.linePThreshCut=tsInputs.low;

    }
  /* Perform the analysis on the data seg given.*/
  tsInputs.allocFlag = 1;
  lal_errhandler = LAL_ERR_RTRN;
  errCode = LAL_CALL( LALSignalTrackSearch(status,&outputCurves,tfmap,&tsInputs),
	    status);
  if ( errCode != 0 )
    {
      fprintf(stderr,"Error calling LALSignalTrackSearch!\n");
      fprintf(stderr,"Error Code: %s\n",status->statusDescription);
      fprintf(stderr,"Function  : %s\n",status->function);
      fprintf(stderr,"File      : %s\n",status->file);
      fprintf(stderr,"Line      : %i\n",status->line);
      fprintf(stderr,"Aux information.\n");
      fprintf(stderr,"Start  point count: %i\n",outputCurves.store.numLStartPoints);
      fprintf(stderr,"Member point count: %i\n",outputCurves.store.numLPoints);
      fprintf(stderr,"Curves returned   : %i\n",outputCurves.numberOfCurves);
      fflush(stderr);
    }
  
  /* 
   * Call tracksearch again to free any temporary ram in 
   * variable outputCurves which is no longer required
   */
  tsInputs.allocFlag = 2;
  LAL_CALL(  LALSignalTrackSearch(status,&outputCurves,tfmap,&tsInputs),
	     status);
  /*
   * Setup for call to function to do map marking
   * We mark maps is convert Tbin and Fbin to 
   * Hz and GPSseconds
   */
  /*Height -> time bins */
  /*Width -> freq bins */
  /* Need to configure mapMarkerParams correctly!*/
  lal_errhandler = LAL_ERR_RTRN;
  errCode = LAL_CALL(  LALTrackSearchInsertMarkers(status,
					 &outputCurves,
					 &tsMarkers),
	     status);
  if ( errCode != 0 )
    {
      fprintf(stderr,"Error calling insert markers routine.\n");
      fprintf(stderr,"%s\n",status->statusDescription);
      fflush(stderr);
    }
  /* 
   *Call the connect curve routine if argument is specified 
   */
  if (params.joinCurves)
    {
      if (params.verbosity > quiet)
	{
	  fprintf(stdout,"Connecting found tracks.\n");
	  fflush(stdout);
	}
      fprintf(stdout,"Subroutines for segment connection NOT TESTED\n");
      fprintf(stdout,"as of Fri-Jun-19-2009:200906191018\n");
      lal_errhandler = LAL_ERR_RTRN;
      errCode=LAL_CALL( LALTrackSearchConnectSigma(status,
						   &outputCurves,
						   *tfmap,
						   tsInputs),
			status);
      if (errCode !=0)
	{
	  fprintf(stderr,"Call to connect trigger subroutine failed!.\n");
	  fprintf(stderr,"%s\n",status->statusDescription);
	  fflush(stderr);
	}
    };
  /*
   * Write Pre-Threshold Results to Disk
   */
  if (params.verbosity >= printFiles)
    {
      sprintf(outputFilename->data,"Pre-%s",outputFilenameMask->data);
      LALappsWriteCurveList(status,
			    outputFilename->data,
			    outputCurves,
			    NULL);
    }
  /*
   * Apply the user requested thresholds 
   */
  outputCurvesThreshold.curves = NULL;
  lal_errhandler = LAL_ERR_RTRN;
  errCode = LAL_CALL( LALTrackSearchApplyThreshold(status,
					 &outputCurves,
					 &outputCurvesThreshold,
					 params),
	    status);
  if ( errCode != 0 )
    {
      fprintf(stderr,"Error calling trigger threshold selection routine.\n");
      fprintf(stderr,"%s\n",status->statusDescription);
      fflush(stderr);
    }
  /*
   * Record user request thresholds into output data structure
   * as a simple reference
   */
  outputCurvesThreshold.minPowerCut=params.MinPower;
  outputCurvesThreshold.minLengthCut=params.MinLength;
  outputCurvesThreshold.minSNRCut=params.MinSNR;
  outputCurvesThreshold.startThreshCut=tsInputs.high;
  outputCurvesThreshold.linePThreshCut=tsInputs.low;
  /* 
   * Dump out list of surviving candidates
   */
  LALappsDetermineFilename(status,
			   tsMarkers,
			   &outputCandidateFilename,
			   ".candidates");
  LALappsWriteSearchResults(status,
			    outputCandidateFilename->data,
			    outputCurvesThreshold);

/*   LALAPPSWRITESEARCHRESULTS(status, */
/* 			    outputCandidateFilename->data, */
/* 			    outputCurves); */

  if (params.verbosity >= printFiles)
    LALappsWriteCurveList(status,
			  outputFilenameMask->data,
			  outputCurvesThreshold,
			  &params);
  /*
   * General Memory Cleanup
   */
  LALappsDestroyCurveDataSection(status,
				 &(outputCurves.curves),
				 outputCurves.numberOfCurves);
  LALappsDestroyCurveDataSection(status,
				 &(outputCurvesThreshold.curves),
				 outputCurvesThreshold.numberOfCurves);
  if (outputFilename)
    XLALDestroyCHARVector(outputFilename);
  if (outputFilenameMask)
    XLALDestroyCHARVector(outputFilenameMask);
  if (outputCandidateFilename)
    XLALDestroyCHARVector(outputCandidateFilename);

  /*************************************************************/
}
/* 
 * End Do_tracksearch function
 */

void
LALappsDoTSeriesSearch(LALStatus         *status,
		       REAL4TimeSeries   *signalSeries,
		       TSSearchParams     params,
		       INT4               callNum)
{
  CreateTimeFreqIn       tfInputs;/*Input params for map making*/
  INT4                   j;
  LALWindowParams        windowParams;/*Needed to generate windowing funcs*/
  REAL4Window           *tempWindow = NULL;
  TimeFreqParam         *autoparams = NULL;/*SelfGenerated values for TFmap*/
  TimeFreqRep           *tfmap = NULL;/*TF map of dataset*/
  TrackSearchMapMarkingParams  mapMarkerParams;/*Struct of params for marking*/
  TrackSearchParams      inputs;/*Native LAL format for tracksearch module*/
  UINT4                  injectionRun=0;/*1=yes*/
  CHARVector            *binaryFilename=NULL;
  TSAMap                *mapBuilder=NULL;/* point to write vars for saving*/
  TSAMap                *tmpTSA=NULL;/*point to correct vars for cropping*/
  REAL8                 signalStop=0;
  REAL8                 signalStart=0;
  REAL8                 cropDeltaT=0;
  INT4                  tmpSegDataPoints=0;
  INT4                  tmpMapTimeBins=0;
  INT4                  errCode=0;
  /*
   * Error checking section
   */
  LALappsTSassert((callNum >= 0),
		  TRACKSEARCHC_EVAL,
		  TRACKSEARCHC_MSGEVAL);
  if ((signalSeries == NULL) && ((params.injectMapCache !=NULL) ||
				 (params.injectSingleMap !=NULL)))
    injectionRun=1;
  /*
   * If we are creating our maps for analysis via tseries data
   */
  if (injectionRun == 0)
    { 
      tfInputs.type = params.TransformType;
      tfInputs.fRow = 2*params.FreqBins;
      /*
       * TimeBin information needs to have buffer TimeBins added here
       * this should be transparent to most users unless they invoke
       * a --bin_buffer flag.  Then the assumed 2 Tbins becomes
       * something else.
       */
      tfInputs.tCol = params.TimeBins+(2*params.colsToClip);
      tfInputs.wlengthF = params.windowsize;
      tfInputs.wlengthT = params.windowsize;

      errCode=LAL_CALL(LALCreateTimeFreqParam(status,&autoparams,&tfInputs),
		       status);
      if (errCode != 0)
	{
	  fprintf(stderr,"Error calling LALCreateTimeFreqParam\n");
	  fprintf(stderr,"Error Code: %s\n",status->statusDescription);
	  fprintf(stderr,"Function  : %s\n",status->function);
	  fprintf(stderr,"File      : %s\n",status->file);
	  fprintf(stderr,"Line      : %i\n",status->line);
	}
      errCode=LAL_CALL(LALCreateTimeFreqRep(status,&tfmap,&tfInputs),
		       status);
      if (errCode != 0)
	{
	  fprintf(stderr,"Error calling LALCreateTimeFreqRep\n");
	  fprintf(stderr,"Error Code: %s\n",status->statusDescription);
	  fprintf(stderr,"Function  : %s\n",status->function);
	  fprintf(stderr,"File      : %s\n",status->file);
	  fprintf(stderr,"Line      : %i\n",status->line);
	}
      /*
       * There is an issue with the overlapping of fft windows used to
       * construct the TFR.  Example:
       * SegLength = 32
       * Tbins = 32 
       * then tfmap->timeInstant[i]=i
       * is Ok but consider the case of we now want a TFR with only
       * 8 time bins.  Then the time indices of TFR variable are no
       * longer right.  We must reindex timeInstant into
       * 2,6,10,14,18,22,26,30
       * This will then make the correct TFR
       * SegPoints = T, TimeBins=B
       * Indices  ::: floor(T/2B+(i*(T/B)))
       * i goes 0,1,2,...B 
       * Bug found Date 26Oct06 Thur
       */
      /* Tue-Nov-06-2007:200711061138 */
      /* Add in timeInstant calculations with buffered dat */
      tmpSegDataPoints=params.SegLengthPoints+(2*params.SegBufferPoints);
      tmpMapTimeBins=(params.TimeBins+(2*params.colsToClip));
      for (j=0;j<tfmap->tCol;j++)
	{
	  tfmap->timeInstant[j]=floor(
				      (tmpSegDataPoints/(2*tmpMapTimeBins))
				      +
				      ((j*(tmpSegDataPoints/tmpMapTimeBins)))
				      );
	}
      windowParams.length = params.windowsize;
      windowParams.type = params.window;
      /*Create selected type of window*/
      switch (params.window) {
	case Rectangular:
		tempWindow = XLALCreateRectangularREAL4Window(params.windowsize);
		break;

	case Hann:
		tempWindow = XLALCreateHannREAL4Window(params.windowsize);
		break;

	case Welch:
		tempWindow = XLALCreateWelchREAL4Window(params.windowsize);
		break;

	case Bartlett:
		tempWindow = XLALCreateBartlettREAL4Window(params.windowsize);
		break;

	case Parzen:
		tempWindow = XLALCreateParzenREAL4Window(params.windowsize);
		break;

	case Papoulis:
		tempWindow = XLALCreatePapoulisREAL4Window(params.windowsize);
		break;

	case Hamming:
		tempWindow = XLALCreateHammingREAL4Window(params.windowsize);
		break;

	case Kaiser:
		tempWindow = XLALCreateKaiserREAL4Window(params.windowsize, 1000);
		fprintf(stderr,"For Kaiser beta hard wired to 1,000\n");
		fflush(stderr);
		break;

	case Creighton:
		tempWindow = XLALCreateCreightonREAL4Window(params.windowsize, 1000);
		fprintf(stderr,"For Creighton beta hard wired to 1,000\n");
		fflush(stderr);
		break;

	case Tukey:
		tempWindow = XLALCreateTukeyREAL4Window(params.windowsize, 1000);
		fprintf(stderr,"For Tukey beta hard wired to 1,000\n");
		break;

	default:
	  fprintf(stderr,"Window specified not allowed.\n");
	  fflush(stderr);
	  exit(TRACKSEARCHC_EVAL);
	}

      /* 
       * Case statment that look to do the various TF Reps 
       */
      if (params.verbosity > quiet)
	{
	  fprintf(stdout,"Creating TFR: ");
	  fflush(stdout);
	}

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
	    if (params.verbosity > quiet)
	      {
		fprintf(stdout,"Spectrogram\n");
		fflush(stdout);
	      }
	    /* Required from deprication of LALWindow function */
	    memcpy(autoparams->windowT->data,
		   tempWindow->data->data,
		   (windowParams.length * sizeof(REAL4)));
	    errCode=LAL_CALL( LALTfrSp(status,signalSeries->data,tfmap,autoparams),
			      status);
	    if ( errCode != 0 )
	      {
		fprintf(stderr,"Error calling LALTfrSp!\n");
		fprintf(stderr,"Error Code: %s\n",status->statusDescription);
		fprintf(stderr,"Function  : %s\n",status->function);
		fprintf(stderr,"File      : %s\n",status->file);
		fprintf(stderr,"Line      : %i\n",status->line);
	      }
	  }
	  break;
	case WignerVille:
	  {

	    LAL_CALL( LALTfrWv(status,signalSeries->data,tfmap,autoparams),
		      status);
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
	    LAL_CALL( LALTfrPswv(status,signalSeries->data,tfmap,autoparams),
		      status);
	  }
	  break;
	case RSpectrogram:
	  {
	    /* Required from deprication of LALWindow function */
	    memcpy(autoparams->windowT->data,
		   tempWindow->data->data,
		   (windowParams.length * sizeof(REAL4)));
	    LAL_CALL( LALTfrRsp(status,signalSeries->data,tfmap,autoparams),
		      status);
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
      XLALDestroyREAL4Window(tempWindow);
      LAL_CALL( LALDestroyTimeFreqParam(status,&autoparams),
		status);
    }
  /*
   * End preparing map from time series data
   */
  /*
   * Setup map variable via reading in the params with the 
   * map file name
   */
  if (injectionRun==1)
    {
      /*
       * Read in file with a map
       */
    }
  /* 
   * Determine the 'true' cropped map starting epoch.  This should in
   * pratice match the epoch the user originally intended without the
   * implicit time bin addition and subsequent cropping
   */
  /*Translate buffer points to GPS delta T */
  /*Add this time to the time of the buffered data segment */
  /*Report this as the true start time of cropped map */
  cropDeltaT=signalSeries->deltaT*params.SegBufferPoints;

  mapMarkerParams.mapStartGPS=signalSeries->epoch;
  signalStart=XLALGPSGetREAL8(&signalSeries->epoch);
  /*
   * Fix the signalStop time stamp to be without the buffer points.  It
   * should be the stop time of the clipped TFR.
   */
  signalStop=(signalSeries->deltaT*(signalSeries->data->length))+signalStart;
  XLALGPSSetREAL8(&(mapMarkerParams.mapStopGPS),signalStop);
  mapMarkerParams.mapTimeBins=tfmap->tCol;
  mapMarkerParams.mapFreqBins=((tfmap->fRow/2)+1);
  /*This is the map time resolution*/
  mapMarkerParams.deltaT=(signalStop-signalStart)/mapMarkerParams.mapTimeBins;
  mapMarkerParams.dataDeltaT=signalSeries->deltaT;
  /* 
   * Properly CROP TFR due to buffered segments used to create the
   * TFR.  MAKE SURE to account for proper TFR dims and time/freq
   * information.
   */
  tmpTSA=XLALMalloc(sizeof(TSAMap));
  tmpTSA->imageCreateParams=tfInputs;
  tmpTSA->clippedWith=0;
  tmpTSA->imageBorders=mapMarkerParams;
  tmpTSA->imageRep=tfmap;
  if ((params.verbosity > quiet) && (params.colsToClip > 0))
    {
      fprintf(stdout,"Cropping TFR ");
      fflush(stdout);
    }
  LALappsTSACropMap(status,&tmpTSA,params.colsToClip);
    
  /* 
   *Copy information from cropping procedure to relevant structures.
   */
  memcpy(&tfInputs,&(tmpTSA->imageCreateParams),sizeof(CreateTimeFreqIn));
  memcpy(&mapMarkerParams,&(tmpTSA->imageBorders),sizeof(TrackSearchMapMarkingParams));
  if (params.verbosity >= printFiles)
    LALappsTSAWritePGM(status,tmpTSA,NULL); 
  tfmap=tmpTSA->imageRep;
  XLALFree(tmpTSA);

  /*
   * Fill in LALSignalTrackSearch params structure via the use of
   * the TSSearch huge struct elements.  These options have been
   * adjusted properly in the Crop subroutine call.
   */
  inputs.sigma=(params.LineWidth)/(2*sqrt(3));
  inputs.high=params.StartThresh;
  inputs.low=params.LinePThresh;
  inputs.width=((tfmap->fRow/2)+1);
  inputs.height=tfmap->tCol;

  /* If requested to the auto lambda determination */
  if (params.verbosity > quiet)
    {
      fprintf(stdout,"Created TFR. \n");
      fflush(stdout);
    }
  /*
   * Call subroutine to run the search
   */
  LALappsDoTrackSearch(status,
		       tfmap,
		       inputs,
		       mapMarkerParams,
		       params);

  if (params.verbosity > quiet)
    {
      fprintf(stdout,"Analyzed TFR. \n");
      fflush(stdout);
    }
  /*
   * Assemble the TSAmap to write to disk
   */
  if (1==1)
    {
      binaryFilename=XLALCreateCHARVector(maxFilenameLength);
      /*
       * Set filename
       */
      /* tsaMap_gpsSeconds_gpsNanoseconds.map*/
      /*
       * Assemble a tsaMap structure 
       */
      mapBuilder=XLALMalloc(sizeof(TSAMap));
      mapBuilder->imageCreateParams=tfInputs;
      mapBuilder->clippedWith=0;
      mapBuilder->imageBorders=mapMarkerParams;
      mapBuilder->imageRep=tfmap;
      sprintf(binaryFilename->data,
	      "tsaMap_Start_%i_%i_Stop_%i_%i.map",
	      mapMarkerParams.mapStartGPS.gpsSeconds,
	      mapMarkerParams.mapStartGPS.gpsNanoSeconds,
	      mapMarkerParams.mapStopGPS.gpsSeconds,
	      mapMarkerParams.mapStopGPS.gpsNanoSeconds);
      LALappsTSAWriteMapFile(status,
			     mapBuilder,
			     NULL);
      if (binaryFilename)
	XLALDestroyCHARVector(binaryFilename);

      XLALFree(mapBuilder);
    }
  /* 
   * Destroy the memory holding the actual TF map
   */
  LAL_CALL( LALDestroyTimeFreqRep(status,&tfmap),
	    status);
}
/*
 * End LALappsDoTSeriesSearch
 */

void
LALappsDoTimeSeriesAnalysis(LALStatus          *status,
			    TSSearchParams      params,
			    TSappsInjectParams  injectParams,
			    CHAR               *cachefile,
			    CHARVector         *dirpath)
{
  TSSegmentVector  *SegVec = NULL;/*Holds array of data sets*/
  TSCreateParams    SegVecParams;/*Hold info for memory allocation*/
  REAL4TimeSeries  *dataset = NULL;/*Structure holding entire dataset to analyze*/
  REAL4TimeSeries  *injectSet = NULL;/*Struct for holding inject data*/
  UINT4             i;  /* Counter for data breaking */
  INT4              j;
  UINT4             productioncode = 1; /* Variable for ascii or frame */
  UINT4             newPointLength=0;
  LIGOTimeGPS       edgeOffsetGPS;
  REAL8             originalFloatTime=0;
  REAL8             revisedSamplingRate=0;/*Account for Fhet*/
  REAL8             newFloatTime=0;
  CHARVector       *dataLabel=NULL;
  /* Set to zero to use ascii files */
  /*
   * Check for nonNULL directory path Ptr to frame files
   */
  LALappsTSassert(((dirpath!=NULL)||(cachefile!=NULL)),
		  TRACKSEARCHC_ENULL,
		  TRACKSEARCHC_MSGENULL);
  LALappsTSassert((status!=NULL),
		  TRACKSEARCHC_ENULL,
		  TRACKSEARCHC_MSGENULL);
  /* 
   * This is the data reading section of code
   * Sources: Frame file or Single Ascii file (1C)
   * All pathing is relative for this code
   * 0 - use frames
   * 1 - use Ascii file
   */
  /*
   * Allocate space for dataset time series
   */
  /* create and initialize the time series vector */
  /* 
   *Adjust value of edgeOffsetGPS to be new start point 
   * including col buffer data
   */
    originalFloatTime=XLALGPSGetREAL8(&(params.GPSstart));
    newFloatTime=originalFloatTime-(params.SegBufferPoints/params.SamplingRate);
    LAL_CALL(LALFloatToGPS(status,
			   &edgeOffsetGPS,
			   &newFloatTime),
	     status);
    /*
     * To handle special case where data is hetordyned and resampled
     * we determine the interval of time that the data would have spanned 
     * at params->HeterodyneSamplingRate then we determine how many points
     * are needed to properly read in this amount of data at the a sampling
     * rate of params->SamplingRate
     */
    newPointLength=params.TimeLengthPoints+(2*params.SegBufferPoints);
    if ((params.HeterodyneFrequency > 0) && 
	(params.HeterodyneSamplingRate > 0))
      {
	newPointLength=floor(
			     (params.TimeLengthPoints+
			      (2*params.SegBufferPoints))
			     /params.HeterodyneSamplingRate*
			     params.SamplingRate);
      }
    dataset=XLALCreateREAL4TimeSeries(params.channelName,
				      &edgeOffsetGPS,
				      0,
				      1/params.SamplingRate,
				      &lalADCCountUnit,
				      newPointLength);
  if (productioncode == 1) /* Use production frame file inputs */
    {
      if (params.makenoise < 1 )
	{  
	  if (params.verbosity >= verbose)
	    {
	      fprintf(stdout,"Reading frame files.\n");
	      fflush(stdout);
	    }
	  LALappsGetFrameData(status,&params,dataset,dirpath,cachefile);
	  if (params.verbosity >= verbose)
	    {
	      fprintf(stdout,"Reading frame complete.\n");
	      fflush(stdout);
	    }
	}
      else
	{
	  /* 
	   * Fake Data Generation routing generates UnitGaussian 
	   * random values 
	   *  the output array is desire length in a 
	   * time series structure 
	   */
	  fakeDataGeneration(status,dataset,
			     (params.NumSeg*params.SegLengthPoints),
			     params.makenoise);
	};
    }
  else
    {
      LALappsGetAsciiData(status,&params,dataset,dirpath);
    }
  /*
   * If injections were requested load them up!  We return a NULL time
   * series variable if no injection is to be done!
   */
  if (injectParams.injectTxtFile != NULL)
    {
      /* Wed-Oct-24-2007:200710241529 
       * Adjust starting point of injection to avoid
       * cropping first injection from TFR
       */
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"Preparing the injection data!\n");
	  fprintf(stdout,"Inject Offset    %f\n",injectParams.startTimeOffset);
	  fprintf(stdout,"Inject Space     %f\n",injectParams.injectSpace);
	  fprintf(stdout,"Inject Count     %i\n",injectParams.numOfInjects);
	  fprintf(stdout,"Inject Scale     %f\n",injectParams.scaleFactor);
	  fprintf(stdout,"Inject Sampling  %f\n",injectParams.sampleRate);
	  fflush(stdout);
	}

      LALappsCreateInjectableData(status,
				  &injectSet,
				  injectParams);
    }
  if ((params.verbosity >= printFiles) && (injectSet != NULL))
    print_real4tseries(injectSet,"CreatedSoftwareInjectableData.diag");
  /* 
   * Prepare the data for the call to Tracksearch 
   */
  SegVecParams.dataSegmentPoints = params.SegLengthPoints;
  SegVecParams.numberDataSegments = params.NumSeg;
  SegVecParams.SegBufferPoints = params.SegBufferPoints;
  /*
   * Allocate structure to hold all segment data
   */
  LAL_CALL(LALCreateTSDataSegmentVector(status,&SegVec,&SegVecParams),
	   status);

  /*
   * Wed-Oct-24-2007:200710241539 
   * Edit function to fill buffered segments
   */
  LALappsTrackSearchPrepareData(status,dataset,injectSet,SegVec,params);
  /*
   * Remove the dataset variable to make room in RAM for TFRs
   */
  if (params.verbosity > quiet)
    {
      fprintf(stdout,"Freeing RAM associated with original time series data.\n");
      fflush(stdout);
    }

  if (dataset)
    {
      XLALDestroyREAL4TimeSeries(dataset);
    }
  if (injectSet)
    {
      XLALDestroyREAL4TimeSeries(injectSet);
    }

  j=0;
  if (params.verbosity >= verbose)
    {
      fprintf(stdout,"Analyzing a total of %i subsegments of data.\n",params.NumSeg);
      fflush(stdout);
    }
  for(i = 0;i < params.NumSeg;i++)
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"\n");
	  fprintf(stdout,"Analyzing Segment %i of %i :",i+1,params.NumSeg);
	  fflush(stdout);
	}
      /*
       * Call to prepare tSeries search
       * Rewrite functions inside to CROP maps and adjust time marks
       * according to the uncropped parts of TFR. Ignore segment
       * buffers....
       */
      if (params.verbosity >= printFiles)
	{
	  dataLabel=XLALCreateCHARVector(maxFilenameLength);
	  sprintf(dataLabel->data,"SegmentAnalyzed_%i.diag",i);
	  print_real4tseries(SegVec->dataSeg[j],dataLabel->data);
	  sprintf(dataLabel->data,"SegmentAnalyzed_%i_Units.diag",i);
	  print_lalUnit(SegVec->dataSeg[j]->sampleUnits,dataLabel->data);
	  XLALDestroyCHARVector(dataLabel);
	}
      LALappsDoTSeriesSearch(status,SegVec->dataSeg[j],params,j);
      j++;
    };
  fprintf(stdout,"\n");
  fflush(stdout);
  /* Free some of the memory used to do the analysis */
  if (params.dataSegVec)
    {
      LAL_CALL(LALDestroyTSDataSegmentVector(status,(params.dataSegVec)),
	       status);
    }
  if (SegVec)
    {
      LAL_CALL(LALDestroyTSDataSegmentVector(status,SegVec),
	       status);
    }
}
/*
 * End LALappsDoTimeSeriesAnalysis
 */

void
LALappsDoTSAMapAnalysis(LALStatus        *status,
			TSSearchParams    params)
{
  UINT4                  i=0;
  TSAMap                *currentMap=NULL;
  TSAcache              *mapCache=NULL;
  CHARVector            *mapFilenameVec=NULL;
  /*
   * Error Checking
   */
  LALappsTSassert((params.injectMapCache != NULL),
		  TRACKSEARCHC_ENULL,
		  TRACKSEARCHC_MSGENULL);
  /*
   * Only configuration options we need are
   * minLength
   * minPower
   * lineWidth
   */
  /* 
   * Open a map cache and sequentially loop through the entries 
   * performing the analysis
   */
  mapFilenameVec=XLALCreateCHARVector(maxFilenameLength);

  strcpy(mapFilenameVec->data,params.injectMapCache);

  LALappsTSALoadCacheFile(status,
			  mapFilenameVec,
			  &mapCache);

  XLALDestroyCHARVector(mapFilenameVec);

  if (mapCache->numMapFilenames < 1)
    {
      fprintf(stderr,"Cache file has no entries!\n");
      LALappsTSassert(0,
		      TRACKSEARCHC_EARGS,
		      TRACKSEARCHC_MSGEARGS);
    }
  for (i=0;i<mapCache->numMapFilenames;i++)
    {
      /*
       * Perform the search for each individual map
       * Create the map to process
       */
      LALappsTSAReadMapFile(status,
			    &currentMap,
			    mapCache->filename[i]);
      /*
       * Execute complete search and write results
       */
      LALappsDoTSAMapSearch(status,
			    currentMap,
			    &params,
			    0);
      /*
       * Free RAM
       */
      LALappsTSADestroyMap(status,
			   &currentMap);
    }
  /*
   * Free map cache structure
   */
  LALappsTSADestroyCache(status,
			 &mapCache);
}
/*
 * End LALappsDoTSAMapAnalysis
 */

void
LALappsDoTSAMapSearch(LALStatus          *status,
		      TSAMap             *tfmap,
		      TSSearchParams     *params,
		      INT4                callNum)
{
  TrackSearchParams      inputs;/*Native LAL format for tracksearch
				  module*/
  LALappsTSassert((callNum >= 0),
		  TRACKSEARCHC_EVAL,
		  TRACKSEARCHC_MSGEVAL);
  /*
   * Setup search parameter structure
   */
  inputs.sigma=(params->LineWidth)/(2*sqrt(3));
  inputs.high=params->StartThresh;
  inputs.low=params->LinePThresh;
  inputs.width=((tfmap->imageRep->fRow/2)+1);
  inputs.height=tfmap->imageRep->tCol;
  LALappsDoTrackSearch(status,
		       tfmap->imageRep,
		       inputs,
		       tfmap->imageBorders,
		       *params);
}
/*
 * End LALappsDoTSAMapSearch
 */

void
LALappsWriteCurveList(LALStatus            *status,
		      CHAR                 *filename,
		      TrackSearchOut        outCurve,
		      TSSearchParams       *params)
{

  CHARVector      *totalName=NULL;
  CHARVector      *breveName=NULL;
  CHARVector      *configName=NULL;
  /*
   * Error checks
   */
  LALappsTSassert((status != NULL),
		  TRACKSEARCHC_ENULL,
		  TRACKSEARCHC_MSGENULL);
  LALappsTSassert((filename !=NULL),
		  TRACKSEARCHC_ENULL,
		  TRACKSEARCHC_MSGENULL);


  /* Output Breve file */
  breveName=XLALCreateCHARVector(maxFilenameLength);
  sprintf(breveName->data,"%s.breve",filename);
  LALappsWriteBreveResults(status,
			   breveName->data,
			   outCurve);
  if (breveName)
    XLALDestroyCHARVector(breveName);

  /* Output Total file */
  totalName=XLALCreateCHARVector(maxFilenameLength);

  sprintf(totalName->data,"%s.full",filename);
  LALappsWriteSearchResults(status,
			    totalName->data,
			    outCurve);
  if (totalName)
    XLALDestroyCHARVector(totalName);

  /* If possible output configuration information */
  if (params!=NULL)
    {
      configName=XLALCreateCHARVector(maxFilenameLength);
      sprintf(configName->data,"%s.config",filename);
      LALappsWriteSearchConfig(status->statusPtr,
			       configName->data,
			       *params);
    }     
  if (configName)
    XLALDestroyCHARVector(configName);

  return;
}
/*
 * End LALappsWriteCurveList
 */

void
LALappsWriteSearchConfig(LALStatus          *status,
			 const CHAR*         myFilename,
			 TSSearchParams      myParams)
{
  FILE            *configFile=NULL;
  INT4             TB=0;
  INT4             FB=0;
  
  configFile=fopen(myFilename,"w");
  fprintf(configFile,"tSeriesAnalysis\t: %i\n",myParams.tSeriesAnalysis);
  fprintf(configFile,"searchMaster\t: %i\n",myParams.searchMaster);
  fprintf(configFile,"haveData\t: %i\n",myParams.haveData);
  fprintf(configFile,"numSlaves\t: NO USED\n");
  fprintf(configFile,"GPSstart\t: %i,%i\n",
	  myParams.GPSstart.gpsSeconds,
	  myParams.GPSstart.gpsNanoSeconds);
  fprintf(configFile,"TimeLengthPoints\t: %i\n",
	  myParams.TimeLengthPoints);
  fprintf(configFile,"discardTLP\t: %i\n",myParams.discardTLP);
  fprintf(configFile,"SegLengthPoints\t: %i\n",
	  myParams.SegLengthPoints);
  fprintf(configFile,"NumSeg\t: %i\n",myParams.NumSeg);
  fprintf(configFile,"SamplingRate\t: %8.3f\n",myParams.SamplingRate);
  fprintf(configFile,"OriginalSamplingRate\t: %8.3f\n",
	  myParams.SamplingRateOriginal);
  fprintf(configFile,"Tlength\t: %i,%i\n",
	  (INT4)myParams.Tlength.gpsSeconds,
	  (INT4)myParams.Tlength.gpsNanoSeconds);
  fprintf(configFile,"TransformType\t: %i\n",
	  (myParams.TransformType));
  fprintf(configFile,"LineWidth\t: %i\n",myParams.LineWidth);
  fprintf(configFile,"MinLength\t: %i\n",myParams.MinLength);
  fprintf(configFile,"MinPower\t: %f\n",myParams.MinPower);
  fprintf(configFile,"overlapFlag\t: %i\n",myParams.overlapFlag);
  fprintf(configFile,"whiten\t: %i\n",myParams.whiten);
  fprintf(configFile,"avgSpecMethod\t: %i\n",
	  (myParams.avgSpecMethod));
  fprintf(configFile,"avgSpecWindow\t: %i\n",
	  (myParams.avgSpecWindow));
  fprintf(configFile,"multiResolution\t: %i\n",
	  myParams.multiResolution);
  FB=myParams.FreqBins;
  TB=myParams.TimeBins;
  fprintf(configFile,"FreqBins\t: %i\n",FB);
  fprintf(configFile,"TimeBins\t: %i\n",TB);
  fprintf(configFile,"windowsize\t: %i\n",myParams.windowsize);
  fprintf(configFile,"window\t: %i\n",myParams.window);
  fprintf(configFile,"numEvents\t: %i\n",myParams.numEvents);
  if (myParams.channelName == NULL)
    fprintf(configFile,"channelName\t: NULL\n");
  else
    fprintf(configFile,"channelName\t: %s\n",myParams.channelName);
  if (myParams.dataDirPath == NULL)
    fprintf(configFile,"channelName\t: NULL\n");
  else
    fprintf(configFile,"dataDirPath\t: %s\n",myParams.dataDirPath);
  if (myParams.singleDataCache == NULL)
    fprintf(configFile,"singleDataCache\t: NULL\n");
  else
    fprintf(configFile,"singleDataCache\t: %s\n",
	    myParams.singleDataCache);
  if (myParams.detectorPSDCache == NULL)
    fprintf(configFile,"detectorPSDCache\t: NULL\n");
  else
    fprintf(configFile,"detectorPSDCache\t: %s\n",
	    myParams.detectorPSDCache);
  if (myParams.channelNamePSD == NULL)
    fprintf(configFile,"channelNamePSD\t: NULL\n");
  else
    fprintf(configFile,"channelNamePSD\t: %s\n",
	    myParams.channelNamePSD);

  fprintf(configFile,"calChannelType\t: %i\n",
	  myParams.calChannelType);
  if (myParams.calFrameCache == NULL)
    fprintf(configFile,"calFrameCache\t: \n");
  else
    fprintf(configFile,"calFrameCache\t: %s\n",
	    myParams.calFrameCache);
  fprintf(configFile,"calibrate\t: %i\n",
	  myParams.calibrate);
  fprintf(configFile,"calibrateIFO\t: %s\n",
	  myParams.calibrateIFO);
  if (myParams.calCatalog == NULL)
    fprintf(configFile,"calCatalog\t: NULL\n");
  else
    fprintf(configFile,"calCatalog\t: %s\n",
	    myParams.calCatalog);
  fprintf(configFile,"dataSegVect\t: MEMORY SPACE\n");
  fprintf(configFile,"currentSeg\t: %i\n",
	  myParams.currentSeg);
  fprintf(configFile,"makenoise\t: %i\n",
	  myParams.makenoise);
  if (myParams.auxlabel == NULL)
    fprintf(configFile,"auxlabel\t: NULL\n");
  else
    fprintf(configFile,"auxlabel\t: %s\n",
	    myParams.auxlabel);
  fprintf(configFile,"joinCurves\t: %i\n",
	  myParams.joinCurves);
  fprintf(configFile,"verbosity\t: %i\n",
	  myParams.verbosity);
  fprintf(configFile,"printPGM\t: %i\n",
	  myParams.printPGM);
  fprintf(configFile,"pgmColorMapFile\t: %s\n",
	  myParams.pgmColorMapFile);
  fprintf(configFile,"injectMapCache\t: %s\n",
	  myParams.injectMapCache);
  fprintf(configFile,"injectSingleMap\t: %s\n",
	  myParams.injectSingleMap);
  fprintf(configFile,"High Pass Freq\t: %f\n",
	  myParams.highPass);
  fprintf(configFile,"Low Pass Freq\t: %f\n",
	  myParams.lowPass);
  fclose(configFile);
  return;
  
}
/*
 * End LALappsWriteSearchConfig
 */

void
LALappsWriteSearchResults(LALStatus      *status,
			  const CHAR*     myFilename,
			  TrackSearchOut  outCurve)
{
  FILE            *totalFile=NULL;
  INT4             i=0;
  INT4             j=0;

  totalFile=fopen(myFilename,"w");
  fprintf(totalFile,"# Total Curves,Lh,Ll: %i,%e,%e\n",outCurve.numberOfCurves,
	  outCurve.startThreshCut,
	  outCurve.linePThreshCut);
  fprintf(totalFile,"# Legend: Col,Row;gpsSec,gpsNanoSec,Freq,depth\n");
  /*Form of solution FreqIndex,TimeIndex,GPSSec,GPSNano,Power*/
  for (i = 0;i < outCurve.numberOfCurves;i++)
    {
      fprintf(totalFile,"Curve number,length,power:%i,%i,%6.18f,%6.1f\n",
	      i,
	      outCurve.curves[i].n,
	      outCurve.curves[i].totalPower,
	      outCurve.curves[i].snrEstimate);
      for (j = 0;j < outCurve.curves[i].n;j++)
	{ /*Long info*/
	  fprintf(totalFile,"%i,%i;%i,%i,%f,%6.18f",
		  outCurve.curves[i].col[j],
		  outCurve.curves[i].row[j],
		  outCurve.curves[i].gpsStamp[j].gpsSeconds,
		  outCurve.curves[i].gpsStamp[j].gpsNanoSeconds,
		  outCurve.curves[i].fBinHz[j],
		  outCurve.curves[i].depth[j]);
	  if (j+1 < outCurve.curves[i].n)
	    fprintf(totalFile,":");
	}
      fprintf(totalFile,"\n");
    }
  /*
   * Close files and clear mem
   */
  fclose(totalFile);
}
/*
 * End LALappsWriteSearchResults
 */

void
LALappsWriteBreveResults(LALStatus      *status,
			 const CHAR*     myFilename,
			 TrackSearchOut  outCurve)
{
  FILE       *breveFile=NULL;
  INT4        i=0;
  REAL8            startStamp=0;
  REAL8            stopStamp=0;

  breveFile=fopen(myFilename,"w");
  /*Short info*/
  fprintf(breveFile,"%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	  "CurveNum",
	  "Power",
	  "Length",
	  "StartF",
	  "FinishF",
	  "StartT",
	  "FinishT",
	  "FreqStart",
	  "FreqStop",
	  "GPSstart",
	  "GPSstop");
  for (i = 0;i < outCurve.numberOfCurves;i++)
    {
      startStamp=XLALGPSGetREAL8(&outCurve.curves[i].gpsStamp[0]);
      stopStamp=XLALGPSGetREAL8(&outCurve.curves[i].gpsStamp[outCurve.curves[i].n-1]);
      fprintf(breveFile,
	      "%12i %12e %12i %12i %12i %12i %12i %12.3f %12.3f %12.3f %12.3f\n",
	      i,
	      outCurve.curves[i].totalPower,
	      outCurve.curves[i].n,
	      outCurve.curves[i].col[0],
	      outCurve.curves[i].col[outCurve.curves[i].n - 1],
	      outCurve.curves[i].row[0],
	      outCurve.curves[i].row[outCurve.curves[i].n - 1],
	      outCurve.curves[i].fBinHz[0],
	      outCurve.curves[i].fBinHz[outCurve.curves[i].n -1],
	      startStamp,
	      stopStamp
	      );

    }
  fclose(breveFile);
}
/*
 * End LALappsWriteBreveResults
 */

/*
 * Local routine to allocate memory for the curve structure
 */
void LALappsCreateCurveDataSection(LALStatus    *status,
				   Curve        **curvein)
{
  INT4    counter;
  Curve  *curve;
 
  INITSTATUS (status, "fixCurveStrucAllocation", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  *curvein = curve = XLALMalloc(MAX_NUMBER_OF_CURVES * sizeof(Curve));
  for (counter = 0; counter < MAX_NUMBER_OF_CURVES;counter++)
    {
      curve[counter].n = 0;
      curve[counter].junction = 0;
      curve[counter].totalPower = 0;
      curve[counter].row = NULL;
      curve[counter].col = NULL;
      curve[counter].depth = NULL;
      curve[counter].row = XLALMalloc(MAX_CURVE_LENGTH * sizeof(INT4));
      curve[counter].col = XLALMalloc(MAX_CURVE_LENGTH * sizeof(INT4));
      curve[counter].depth = XLALMalloc(MAX_CURVE_LENGTH * sizeof(REAL4));
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
      /* Need to gracefully exit here by unallocating later */
    };
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/* End of curve structure allocation function */


/*
 * Local routine to clean up memory once used to hold
 * collection of curve candidates
 */
void LALappsDestroyCurveDataSection(LALStatus    *status,
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
      XLALFree(curve[counter].fBinHz);
      XLALFree(curve[counter].row);
      XLALFree(curve[counter].col);
      XLALFree(curve[counter].depth);
      XLALFree(curve[counter].gpsStamp);
      /* Need to gracefully exit here by unallocating later */
      /* Just in case there is an error in deallocation */
    };
  if (numCurves > 0)
    {
      XLALFree(curve);
    }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/* End function to deallocate the curve structures */


void LALappsCreateInjectableData(LALStatus           *status,
				 REAL4TimeSeries    **injectSet,
				 TSappsInjectParams   params)
{
  UINT4  i=0;
  UINT4  j=0;
  UINT4  k=0;
  UINT4 timePoints=0;
  UINT4 offSetPoints=0;
  UINT4 lineCount=0;
  UINT4 newLineCount=0;
  UINT4 pointLength=0;

  REAL4Sequence   *domain=NULL;
  REAL4Sequence   *range=NULL;
  REAL4Sequence   *waveDomain=NULL;
  REAL4Sequence   *waveRange=NULL;
  REAL4TimeSeries *ptrInjectSet=NULL;
  REAL8            fileDuration=0;
  FILE            *fp=NULL;
  CHAR             c;
  const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;
  /*
   * Try and load the file
   */
  fp = fopen(params.injectTxtFile->data,"r");
  if (!fp)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEFILE);
      exit(TRACKSEARCHC_EFILE);
    }
  /*
   * Count the lines in the file
   */
  while ((c = fgetc(fp)) != EOF)
    {
      if (c == '\n')
	lineCount++;
    }

  if (!(lineCount > 0))
    {
      fprintf(stderr,TRACKSEARCHC_MSGEREAD);
      exit(TRACKSEARCHC_EREAD);
    }
 
  fclose(fp);
  /*
   * Allocate RAM to load up values
   */
  domain=XLALCreateREAL4Vector(lineCount);
  range=XLALCreateREAL4Vector(lineCount);

  /*
   * Expecting 2C data tstep and amp
   */
  /* Reopen file to read in contents! */
  fp = fopen(params.injectTxtFile->data,"r");
  if (!fp)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEFILE);
      exit(TRACKSEARCHC_EFILE);
    }
  for (i=0;i<lineCount;i++)
    {
      if ((fscanf(fp,"%f %f \n",&domain->data[i],&range->data[i])) == EOF)
	{
	  fprintf(stderr,TRACKSEARCHC_MSGEREAD);
	  exit(TRACKSEARCHC_EREAD);
	}
      range->data[i]=(range->data[i]*params.scaleFactor);
    }
  fclose(fp);
  fileDuration=domain->data[domain->length-1]-domain->data[0];
  newLineCount=params.sampleRate*fileDuration;
  offSetPoints=params.sampleRate*params.startTimeOffset;
  timePoints=params.sampleRate*params.injectSpace;
  waveDomain=XLALCreateREAL4Vector(newLineCount);
  waveRange=XLALCreateREAL4Vector(newLineCount);

  for (i=0;i<waveDomain->length;i++)
    {
      waveDomain->data[i]=i*(1/params.sampleRate);
      waveRange->data[i]=1.0;
    }
  /*
   * Call interpolating routine to create resample waveform
   */

  LAL_CALL(LALSVectorPolynomialInterpolation(status,
					     waveDomain,
					     waveRange,
					     domain,
					     range),
	   status);

  pointLength=(params.numOfInjects*
	       (waveDomain->length+timePoints))+offSetPoints;
  LALappsTSassert((pointLength>0),
		  TRACKSEARCHC_EARGS,
		  TRACKSEARCHC_MSGEARGS);
  *injectSet=XLALCreateREAL4TimeSeries("Injections",
				      &gps_zero,
				      0,
				      1/params.sampleRate,
				      &lalADCCountUnit,
				      pointLength);

  /*
   * Copy the data into a new variable to be returned
   */
  i=0;
  ptrInjectSet = *injectSet;
  if (ptrInjectSet != NULL)
    {
      /*Inject Offset silence*/
      for (j=0;j<offSetPoints;j++)
	{
	  (*injectSet)->data->data[i] = 0;
	  i++;
	}
      /*Loop over the requested injections and silence*/
      for (k=0;k<params.numOfInjects;k++)
	{
	  for (j=0;j<waveDomain->length;j++)
	    {
	      (*injectSet)->data->data[i] = waveRange->data[j];
	      i++;
	    }
	  for (j=0;j<timePoints;j++)
	    {
	      (*injectSet)->data->data[i] = 0;
	      i++;
	    }
	}
    }
  if (domain)
    XLALDestroyREAL4Vector(domain);

  if (range)
    XLALDestroyREAL4Vector(range);

  if (waveRange)
    XLALDestroyREAL4Vector(waveRange);

  if (waveDomain)
    XLALDestroyREAL4Vector(waveDomain);

  return;
}
/*
 * End of function LALappsCreateInjectableData to create injectable
 * data structure for calibration etc
 */

/* 
 * Function to perform calibration of input data if necessary
 */
void LALappsTrackSearchCalibrate( LALStatus          *status,
				  REAL4TimeSeries    *dataSet,
				  TSSearchParams      params)
{
  FrCache                  *calcache = NULL;
  COMPLEX8FrequencySeries  *response=NULL;
  REAL4FFTPlan             *dataSetPlan=NULL;
  COMPLEX8FrequencySeries  *dataSetFFT=NULL;
  const LALUnit     strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  UINT4                     segmentPoints=0;
  CalibrationUpdateParams   calfacts;
  UINT4                     i=0;
  UINT4                     j=0;

  /* 
   * Create response function based on starting epoch from 
   * first data segment.  We assume it is constant over 
   * all the other data segments available.
   */
  segmentPoints=params.SegLengthPoints;
  /* 
   * Implicity set Epoch, Unit and df for Response extract
   */
  response=XLALCreateCOMPLEX8FrequencySeries("tmpResponse",
					     &(dataSet->epoch),
					     0,
					     1/(dataSet->deltaT*dataSet->data->length),
					     &strainPerCount,
					     segmentPoints/2+1);

  LAL_CALL(LALFrCacheImport(status,
			    &calcache, 
			    params.calFrameCache),
	   status);
  /*
   * Extract response function from calibration cache
   */
  /* Configure calfacts */
  memset(&calfacts, 0, sizeof(calfacts));
  calfacts.ifo = (CHAR*) LALMalloc(3*sizeof(CHAR));
  strcpy(calfacts.ifo,params.calibrateIFO);


  LAL_CALL( LALExtractFrameResponse(status,
				    response, 
				    calcache, 
				    &calfacts),
	    status);
  if (params.verbosity >= printFiles)
    {
      print_lalUnit(response->sampleUnits,"extractedResponseFunction_Units.diag");
      print_complex8fseries(response,"extractedResponseFunction.diag");
    }
  /*
   * Destroy IFO text pointer
   */
  XLALFree(calfacts.ifo);
  /*
   * Destroy frame cache
   */
  if (calcache)
    {
      LAL_CALL( LALDestroyFrCache(status, &calcache),status);
    }
  /* Done getting response function */
  /*
   * Allocate frequency series for Set FFT
   */
  dataSetFFT=XLALCreateCOMPLEX8FrequencySeries("Entire data FFTed",
					       &(dataSet->epoch),
					       0,
					       (dataSet->deltaT*dataSet->data->length),
					       &(dataSet->sampleUnits),
					       dataSet->data->length/2+1);
  /*
   * FFT input unsegmented data set
   */
  dataSetPlan=XLALCreateForwardREAL4FFTPlan(dataSet->data->length,0);

  LAL_CALL( LALForwardREAL4FFT(status,
			       dataSetFFT->data,
			       dataSet->data,
			       dataSetPlan),
	    status);
  if (dataSetPlan)
    XLALDestroyREAL4FFTPlan(dataSetPlan);


  if (params.verbosity >= printFiles)
    {
      print_complex8fseries(dataSetFFT,"Pre_CalibrateDataSetFFT.diag");
      print_lalUnit(dataSetFFT->sampleUnits,"Pre_CalibrateDataSetFFT_Units.diag");
    }
  LAL_CALL(LALTrackSearchCalibrateCOMPLEX8FrequencySeries(status,
							  dataSetFFT,
							  response),
	   status);
  if (params.verbosity >= printFiles)
    {
      print_complex8fseries(dataSetFFT,"Post_CalibrateDataSetFFT.diag");
      print_lalUnit(dataSetFFT->sampleUnits,"Post_CalibrateDataSetFFT_Units.diag");
    }
  /*
   * Inverse FFT the entire dataSet
   */
  fprintf(stderr,"Hack on Response function...Don't forget\n");
  fprintf(stderr,"Will cause memory leak for some reason?\n");
  dataSetFFT->data->data[dataSetFFT->data->length-1].im=0;
  dataSetFFT->data->data[dataSetFFT->data->length].im=0;
      
  dataSetPlan=XLALCreateReverseREAL4FFTPlan(dataSet->data->length,0);

  LAL_CALL(LALReverseREAL4FFT(status,
			      dataSet->data,
			      dataSetFFT->data,
			      dataSetPlan),
	   status);
  /* 
   * Migrate units fields from F basis to T basis
   */
  dataSet->sampleUnits=dataSetFFT->sampleUnits;
  /*
   * Apply missing factor 1/n
   */
  for (j=0;j<dataSet->data->length;j++)
    dataSet->data->data[i]=dataSet->data->data[i]/dataSet->data->length;

  if (dataSetPlan)
    LAL_CALL(LALDestroyREAL4FFTPlan(status,&dataSetPlan),status);
  /*
   * Free the FFT memory space
   */
  if (dataSetFFT)
    XLALDestroyCOMPLEX8FrequencySeries(dataSetFFT);

  if (response)
    XLALDestroyCOMPLEX8FrequencySeries(response);

  return;
}
/*
 * End function to perform calibration
 */

/*
 * Butterworth Band passing
 */
void LALappsTrackSearchBandPassing( LALStatus           *status,
				    REAL4TimeSeries     *dataSet,
				    TSSearchParams       params)
{
  PassBandParamStruc       bandPassParams;
  int                      errcode=0;
  if (params.lowPass > 0)
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"You requested a low pass filter of the data at %f Hz\n",params.lowPass);
	  fflush(stdout);
	}

      bandPassParams.name=NULL;
      bandPassParams.nMax=20;
      /* F < f1 kept, F > f2 kept */
      bandPassParams.f1=params.lowPass;
      bandPassParams.f2=0;
      bandPassParams.a1=0.9;
      bandPassParams.a2=0;
      /*
       * Band pass is achieved by low pass first then high pass!
       * Call the low pass filter function.
       */
      errcode=XLALButterworthREAL4TimeSeries(dataSet,&bandPassParams);
      if (errcode !=0)
	{
	  fprintf(stderr,"Error trying to Butterworth filter data (lowPass).\n");
	  exit(errcode);
	}
    }
  /*
   * Call the high pass filter function.
   */
  if (params.highPass > 0)
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"You requested a high pass filter of the data at %f Hz\n",params.highPass);
	  fflush(stdout);
	}
      bandPassParams.name=NULL;
      bandPassParams.nMax=20;
      /* F < f1 kept, F > f2 kept */
      bandPassParams.f1=0;
      bandPassParams.f2=params.highPass;
      bandPassParams.a1=0;
      bandPassParams.a2=0.9;
      errcode=XLALButterworthREAL4TimeSeries(dataSet,&bandPassParams);
      if (errcode !=0)
	{
	  fprintf(stderr,"Error trying to Butterworth filter data (lowPass).\n");
	  exit(errcode);
	}
    }
  return;
}
/*
 * End Butterworth Band passing
 */
/* 
 * Variation that removes harmonics from individual data segments
 */ 
void LALappsTracksearchRemoveHarmonicsFromSegments(LALStatus       *status,
						   REAL4TimeSeries *dataSet,
						   TSSegmentVector *dataSegments,
						   TSSearchParams   params)
{
  UINT4             j=0;
  REAL4TimeSeries  *tmpSegmentPtr=NULL;
  CHARVector       *dataLabel=NULL;

  for (j=0;j<dataSegments->length;j++)
    {
      if (params.verbosity > quiet)
	{
	  fprintf(stdout,"Removing lines from segment %i.\n",j+1);
	  fflush(stdout);
	}
      tmpSegmentPtr=(dataSegments->dataSeg[j]);
      dataLabel=XLALCreateCHARVector(128);
	sprintf(dataLabel->data,"Pre_LineRemovalTimeDomainDataSeg_%i.diag",j);
      if (params.verbosity >= printFiles)
	print_real4tseries(tmpSegmentPtr,dataLabel->data);
      sprintf(dataLabel->data,"Pre_LineRemovalTimeDomainDataSeg_%i_Units.diag",j);
      if (params.verbosity >= printFiles)
	print_lalUnit(tmpSegmentPtr->sampleUnits,dataLabel->data);
      LALappsTracksearchRemoveHarmonics(status,tmpSegmentPtr,params);
      sprintf(dataLabel->data,"Post_LineRemovalTimeDomainDataSeg_%i.diag",j);
      if (params.verbosity >= printFiles)
	print_real4tseries(tmpSegmentPtr,dataLabel->data);
      sprintf(dataLabel->data,"Post_LineRemovalTimeDomainDataSeg_%i_Units.diag",j);
      if (params.verbosity >= printFiles)
	print_lalUnit(tmpSegmentPtr->sampleUnits,dataLabel->data);
      XLALDestroyCHARVector(dataLabel);
      tmpSegmentPtr=NULL;
    }
}

/*
 * End Harmonics line removal variant
 */
/*
 * Removing harmonic lines from data
 */
void LALappsTracksearchRemoveHarmonics( LALStatus             *statusX,
				        REAL4TimeSeries       *dataSet,
				        TSSearchParams         params)
{
  COMPLEX8Vector           *referenceSignal=NULL;
  CHARVector               *dataLabel=NULL;
  UINT4                     tmpHarmonicCount=0;
  REAL4                     tmpLineFreq=0.0;
  COMPLEX8FrequencySeries  *signalFFT_harmonic=NULL;
  REAL4FFTPlan             *forwardPlan=NULL;
  INT4Vector               *harmonicIndex=NULL;
  INT4Vector               *harmonicIndexCompliment=NULL;
  REAL4TVectorCLR          *inputTVectorCLR=NULL;
  REAL4FVectorCLR          *inputFVectorCLR=NULL;
  REAL4Vector              *inputPSDVector=NULL;
  REAL4Vector              *cleanData=NULL;
  UINT4                     i=0;
  UINT4                     j=0;
  UINT4                     l=0;
  UINT4                     nanCount=0;
  const LIGOTimeGPS         gps_zero = LIGOTIMEGPSZERO;
  UINT4                     planLength=0;
  int                       code=0;
  LALStatus                 status=blank_status;


  /*lal_errhandler=LAL_ERR_RTRN;*/
  if (params.numLinesToRemove > 0)
    {
      /*
       * Perform the removal of COHERENT LINES specified by user to entire
       * data set.  We want to remove up to n harmonics which have been
       * specified at the command line.  
       */
      if (params.verbosity > quiet)
	{
	  fprintf(stdout,"Removing requested coherent lines, if possible.\n");
	  fflush(stdout);
	}
      for (l=0;l<params.numLinesToRemove;l++)
	{
	  /*Setup the proper number of harmonics that should be done!*/
	  /* Create reference signal IFF Line to remove is less than Fnyq*/
	  tmpHarmonicCount=params.maxHarmonics;
	  tmpLineFreq=params.listLinesToRemove[l];
	  if (params.listLinesToRemove[l] < (params.SamplingRate/2))
	    {
	      if ((params.SamplingRate/(2*tmpLineFreq)) <  tmpHarmonicCount)
		{
		  tmpHarmonicCount=floor(params.SamplingRate/(2*tmpLineFreq));
		  if ((tmpHarmonicCount*tmpLineFreq) >= (params.SamplingRate/2))
		    tmpHarmonicCount=tmpHarmonicCount-1;
		}
	      if (params.verbosity > quiet)
		{
		  fprintf(stdout,"Reference %f, Max Harmonics to remove %i ",tmpLineFreq,tmpHarmonicCount);
		  fflush(stdout);
		}
	      /* Setup both T and F domain input data structures */
	      inputTVectorCLR=(REAL4TVectorCLR*)XLALMalloc(sizeof(REAL4TVectorCLR));
	      inputFVectorCLR=(REAL4FVectorCLR*)XLALMalloc(sizeof(REAL4FVectorCLR));
  
	      /* Take data to F domain */
	      planLength=dataSet->data->length;
	      referenceSignal=XLALCreateCOMPLEX8Vector(planLength);
	      signalFFT_harmonic=XLALCreateCOMPLEX8FrequencySeries("tmpSegPSD",
								   &gps_zero,
								   0,
								   1/(dataSet->data->length*dataSet->deltaT),
								   &(dataSet->sampleUnits),
								   planLength/2+1);
	      inputPSDVector=XLALCreateREAL4Vector(planLength/2+1);
	      forwardPlan=XLALCreateForwardREAL4FFTPlan(planLength,0);
	      code=LAL_CALL( LALForwardREAL4FFT(&status,
					   signalFFT_harmonic->data,
					   dataSet->data,
					   forwardPlan),
			&status);
	      if (params.verbosity >= printFiles)
		{
		  dataLabel=XLALCreateCHARVector(maxFilenameLength);
		  sprintf(dataLabel->data,"LineRemoval_Fdomain_R_%f.diag",tmpLineFreq);
		  print_complex8fseries(signalFFT_harmonic,dataLabel->data);
		  sprintf(dataLabel->data,"LineRemoval_Fdomain_R_%f_Units.diag",tmpLineFreq);
		  print_lalUnit(signalFFT_harmonic->sampleUnits,dataLabel->data);
		  XLALDestroyCHARVector(dataLabel);
		}

	      if (code != 0)
		{
		  fprintf(stdout,"FAILED!\n");
		  fflush(stdout);
		  fprintf(stderr,"Error: LALForwardReal4FFT\n");
		  exit(TRACKSEARCHC_EMISC);
		}
	      /* May need to use XLALREAL[4,8]PowerSpectrum */
	      code=XLALREAL4PowerSpectrum(inputPSDVector,
					  dataSet->data,
					  forwardPlan);
	      if (code != 0)
		{
		  fprintf(stdout,"FAILED!\n");
		  fflush(stdout);
		  fprintf(stderr,"Error: XLALREAL4PowerSpectrum\n");
		  exit(TRACKSEARCHC_EMISC);
		}

	      /*Assign CLR data structures*/
	      /* The  CLR Time Vector  */
	      inputTVectorCLR->length = dataSet->data->length;
	      inputTVectorCLR->data = dataSet->data->data;
	      inputTVectorCLR->deltaT = dataSet->deltaT; /* inverse of the sampling frequency */
	      inputTVectorCLR->fLine = 0.0;
	      /* The  CLR Frequency Vector  */
	      inputFVectorCLR->length = inputPSDVector->length;
	      inputFVectorCLR->data = inputPSDVector->data; 
	      inputFVectorCLR->deltaF = 1.0/(dataSet->data->length*dataSet->deltaT);
	      inputFVectorCLR->fLine =  inputTVectorCLR->fLine;

	      inputTVectorCLR->fLine = tmpLineFreq;
	      inputFVectorCLR->fLine = tmpLineFreq;

	      for (i=0;i<referenceSignal->length;i++)
		{
		  referenceSignal->data[i].re=0;
		  referenceSignal->data[i].im=0;
		}
	      harmonicIndex=XLALCreateINT4Vector(tmpHarmonicCount);
	      if  (harmonicIndex == NULL)
		{
		  fprintf(stderr,"Memory allocation problem in line removal subroutine.\n");
		}
	      harmonicIndexCompliment=XLALCreateINT4Vector(3*tmpHarmonicCount);
	      if  (harmonicIndexCompliment == NULL)
		{
		  fprintf(stderr,"Memory allocation problem in line removal subroutine.\n");
		}
	      for (i=0;i<harmonicIndex->length;i++)
		harmonicIndex->data[i]=i+1;
	      /*Unknown code failure for following function call!!!!*/
	      /*Introduced simple change to code to fix this contacted */
	      /*author via email to notify of changes.Thu-Jun-26-2008:200806261734 */
	      code=LAL_CALL( LALHarmonicFinder(&status,
					       harmonicIndexCompliment,
					       inputFVectorCLR,
					       harmonicIndex),
			&status);
	      if (code != 0)
		{
		  fprintf(stdout,"FAILED!\n");
		  fflush(stdout);
		  fprintf(stderr,"Error: LALHarmonicFinder\n");
		  exit(TRACKSEARCHC_EMISC);
		}
	      code=LAL_CALL( LALRefInterference(&status,
					   referenceSignal,
					   signalFFT_harmonic->data,
					   harmonicIndexCompliment),
			&status);

	      if (code != 0)
		{
		  fprintf(stdout,"FAILED!\n");
		  fflush(stdout);
		  fprintf(stderr,"Error: LALRefInterference\n");
		  exit(TRACKSEARCHC_EMISC);
		}
	      if (harmonicIndexCompliment)
		  XLALDestroyINT4Vector(harmonicIndexCompliment);
	      if (harmonicIndex)
		XLALDestroyINT4Vector(harmonicIndex);

	      /* Clean up input time series with global reference signal */
	      cleanData=XLALCreateREAL4Vector(dataSet->data->length);

	      code=LAL_CALL( LALCleanAll(&status,
					 cleanData,
					 referenceSignal,
					 inputTVectorCLR),
			     &status);

	      if (code != 0)
		{
		  fprintf(stdout,"FAILED!\n");
		  fflush(stdout);
		  fprintf(stderr,"Error: LALCleanAll\n");
		  exit(TRACKSEARCHC_EMISC);
		}

	      /*Copy the clean data back into the variable dataSet */
	      nanCount=0;
	      /* Check that the cleaned data has no NaN values.  If it
	       * does print a warning to the user that this frequency was
	       * not removed.
	       */
	      for (j=0;j<cleanData->length;j++)
		{
		  if (cleanData->data[j] != cleanData->data[j])
		    nanCount=nanCount+1;
		}
	      if (nanCount > 0)
		{
		  fprintf(stderr,"Error with cleaned data, NaN encountered.\n");
		  fprintf(stderr,"For cleaning frequency %f, %i elements of %i elements in array contain NaN value!\n",tmpLineFreq,nanCount,cleanData->length);
		}
	      else
		{
		  for (j=0;j<dataSet->data->length;j++)
		    {
		      dataSet->data->data[j]=cleanData->data[j];
		    }
		}

	      if (cleanData)
		XLALDestroyREAL4Vector(cleanData);

	      if (signalFFT_harmonic)
		XLALDestroyCOMPLEX8FrequencySeries(signalFFT_harmonic);

	      if (inputPSDVector)
		XLALDestroyREAL4Vector(inputPSDVector);

	      if (forwardPlan)
		XLALDestroyREAL4FFTPlan(forwardPlan);

	      if (referenceSignal)
		XLALDestroyCOMPLEX8Vector(referenceSignal);

	      if (inputTVectorCLR)
		XLALFree(inputTVectorCLR);

	      if (inputFVectorCLR)
		XLALFree(inputFVectorCLR);
	      if (params.verbosity > quiet)
		{
		  if (nanCount == 0)
		    {
		      fprintf(stdout,"Removed!\n");
		      fflush(stdout);
		    }
		  else
		    {
		      fprintf(stdout,"Warning, Line not found!\n");
		      fflush(stdout);
		    }
		}
	    }
	  else
	    {
	      fprintf(stderr,"Ignoring line to remove can not be present in data! F_nyq:%f F_line: %f\n",(params.SamplingRate/2),tmpLineFreq);
	    }
	  /* END OF NEW LOOP TO REMOVE ONE LINE AT A TIME */
	  /*
	   * END COHERENT LINE removal section
	   */
	}
    }
  /*lal_errhandler=LAL_ERR_DFLT;*/
  return;
}
/*
 * End line removal
 */
/*
 * Perform software injections if data is available
 */
void LALappsTrackSearchPerformInjection(LALStatus        *status,
					REAL4TimeSeries  *dataSet,
					REAL4TimeSeries  *injectSet,
					TSSearchParams     params)
{
  /*
   * NOTE TO SELF: Tina add code to compare any factor in the dataSet
   * units structure and factor the injected data appropriately to
   * match the input data (dataSet).  This is bug affected by the
   * REAL8->REAL4 conversions in frame reading code only.
   */
  UINT4                     i=0;
  UINT4                     j=0;

  for (i=0,j=0;
       ((i<injectSet->data->length)&&((j+params.SegBufferPoints)<dataSet->data->length));
       i++,j++)
    if (j<params.SegBufferPoints)
      dataSet->data->data[j]=dataSet->data->data[i];
    else
      dataSet->data->data[j]=
	injectSet->data->data[i-params.SegBufferPoints]+
	dataSet->data->data[i];
  return;
}
/*
 * End software injections
 */
/*
 * Whiten data segments
 */
void LALappsTrackSearchWhitenSegments( LALStatus        *status,
				       REAL4TimeSeries  *dataSet,
				       TSSegmentVector  *dataSegments,
				       TSSearchParams     params)
{
  UINT4                    planLength=0;
  REAL4TimeSeries          *tmpSignalPtr=NULL;
  LALUnit                   originalFrequecyUnits;
  AverageSpectrumParams     avgPSDParams;
  CHARVector               *dataLabel=NULL;
  COMPLEX8FrequencySeries  *signalFFT=NULL;
  REAL4FFTPlan             *reversePlan=NULL;
  REAL4FFTPlan             *forwardPlan=NULL;
  REAL4FFTPlan             *averagePSDPlan=NULL;
  REAL4FrequencySeries     *averagePSD=NULL;
  REAL4Window              *windowPSD=NULL;
  UINT4                     i=0;
  UINT4                     j=0;
  UINT4                     stride=0;
  INT4                      segCount=0;
  INT4                      originalDataLength=0;
  UINT4                     segmentLength=0;
  int                       errcode=0;
  const LIGOTimeGPS         gps_zero = LIGOTIMEGPSZERO;

  if (params.verbosity >= verbose )
    {
      fprintf(stdout,"Estimating PSD for whitening\n");
      fflush(stdout);
    }
  averagePSD=XLALCreateREAL4FrequencySeries("averagePSD",
					    &gps_zero,
					    0,
					    1/((params.SegLengthPoints+(2*params.SegBufferPoints))*dataSet->deltaT),
					    &lalDimensionlessUnit,
					    (params.SegLengthPoints+(2*params.SegBufferPoints))/2+1);
  /*
   * The given units above need to be correct to truly reflect the
   * units stored in the frame file
   */
  switch (params.avgSpecWindow) 
    {
    case Rectangular:
      windowPSD = XLALCreateRectangularREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)));
      break;

    case Hann:
      windowPSD = XLALCreateHannREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)));
      break;

    case Welch:
      windowPSD = XLALCreateWelchREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)));
      break;

    case Bartlett:
      windowPSD = XLALCreateBartlettREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)));
      break;

    case Parzen:
      windowPSD = XLALCreateParzenREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)));
      break;

    case Papoulis:
      windowPSD = XLALCreatePapoulisREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)));
      break;

    case Hamming:
      windowPSD = XLALCreateHammingREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)));
      break;

    case Kaiser:
      windowPSD = XLALCreateKaiserREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)), 1000);
      fprintf(stderr,"For Kaiser beta hard wired to 1,000\n");
      break;

    case Creighton:
      windowPSD = XLALCreateCreightonREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)), 1000);
      fprintf(stderr,"For Creighton beta hard wired to 1,000\n");
      break;

    case Tukey:
      windowPSD = XLALCreateTukeyREAL4Window((params.SegLengthPoints+(2*params.SegBufferPoints)), 1000);
      fprintf(stderr,"For Tukey beta hard wired to 1,000\n");
      break;

    default:
      fprintf(stderr,"Window specified not allowed.\n");
      exit(TRACKSEARCHC_EVAL);
    }

  /* If we only do one map we need to something special here*/
  if (params.NumSeg < 2)
    avgPSDParams.overlap=1;
  else /* use same as segment overlap request*/
    {
      avgPSDParams.overlap=params.overlapFlag;
      /*
       * Shift around the overlap for PSD estimation to accomodate
       * the extra data points of 2*params.SegBufferPoints.  This
       * not guarantee that we will have the same number of
       * segments in our PSD estimate as we will process.
       */
      if (params.SegBufferPoints > 0)
	avgPSDParams.overlap=2*params.SegBufferPoints;
    }
      
  avgPSDParams.window=windowPSD;
  avgPSDParams.method=params.avgSpecMethod;
  /* Set PSD plan length and plan*/
  averagePSDPlan=XLALCreateForwardREAL4FFTPlan(params.SegLengthPoints+(2*params.SegBufferPoints),
					       0);

  avgPSDParams.plan=averagePSDPlan;
 
  /*
   * Stride is fixed by 
   * params.NumSeg = 1 + (dataSet->data->length - params.SegLengthPoints)/stride
   * We have not direct control over stride length (correlation amount) of
   * each individual PSD estimate to the mean PSD requested.  Advantage
   * we will make the least correlated PSD estimate possible.  So 
   * numSegs for PSD estimate <= numSegs to analyze
   */
/*   stride=params.SegLengthPoints; */
/*   segCount=1+(dataSet->data->length - params.SegLengthPoints)/stride; */

  switch (params.avgSpecMethod)
    {
    case useMean:
      segmentLength=params.SegLengthPoints+(2*params.SegBufferPoints);
      if ((params.NumSeg%2 !=0))
	{
	  /* Add one to odd number segCount to make it even for this method*/
	  segCount=params.NumSeg+1;
	}
      else
	{
	  segCount=params.NumSeg;
	}
      stride=floor((dataSet->data->length - segmentLength)
		   /((segCount)-1));
      
      /* 
       * If stride is LESS THAN segmentLength/2 up the segCount by 2
       * by adjusting the value of our stride variable
       */
      if (stride < segmentLength/2)
	{
	  segCount=segCount-2;
	  stride=floor((dataSet->data->length - segmentLength)
                       /((segCount)-1));
	}
      /* Workaround "HACK" for proper spanning of data */
      /* See AverageSpectrum.c line 845 */
      if ((dataSet->data->length-((segCount-1)*stride+segmentLength))>=1)
	{
	  fprintf(stdout,"Notice: We can not estimate the PSD using XLALREAL4AverageSpectrumMedianMean because the segment length, bin buffer, and segment overlaps options combine to form a data interval that isn't spanned correctly.  See documentation for the above XLAL function.  We are off by %i points to span the data properly.\n",
		  (dataSet->data->length-((segCount-1)*stride+segmentLength)));
	  fflush(stdout);
	}

      if (/*Check and use a tolerance for mismatches to use medianmean estimator*/
	  (dataSet->data->length-((segCount-1)*stride+segmentLength)
	   <
	   floor(segmentLength*0.005))
	&&
	  (dataSet->data->length-((segCount-1)*stride+segmentLength)
	   >
	   0)
	  )
	{
	  fprintf(stderr,"WARNING! To estimate PSD for whitening we are ignoring the last data %i points in the input data used to make the PSD estimate.\n",
		  (dataSet->data->length-((segCount-1)*stride+segmentLength)));
	  originalDataLength=dataSet->data->length;	
	  dataSet->data->length=originalDataLength-(dataSet->data->length-((segCount-1)*stride+segmentLength));
	  errcode=XLALREAL4AverageSpectrumMedianMean(averagePSD,
						     dataSet,
						     segmentLength,
						     stride,
						     windowPSD,
						     averagePSDPlan);
	  dataSet->data->length=originalDataLength;
	}
      else
	{
	  errcode=XLALREAL4AverageSpectrumMedianMean(averagePSD,
						     dataSet,
						     segmentLength,
						     stride,
						     windowPSD,
						     averagePSDPlan);
	}
	
      if (errcode !=0)
	{
	  fprintf(stderr,"Problem calculating average PSD using mean method.\n");
	  fprintf(stderr,"SegCount Caculated :%i\n Segments :%i\n Dataset :%i\n PSD :%i\n Seglen :%i\n Seglen with buffer bin points :%i\n Stride Caculated:%i\n",
		  segCount,
		  params.NumSeg,
		  dataSet->data->length,
		  averagePSD->data->length,
		  params.SegLengthPoints,
		  params.SegLengthPoints+(2*params.SegBufferPoints),
		  stride);
	  exit(errcode);
	}
      break;
    
    case useMedian:
      /* This case not yet debugged!!! */
      /*Fri-Jun-13-2008:200806131546 */
      fprintf(stderr,"The option useMedian, not yet debugged. Behavior unknown!\n");
      segmentLength=params.SegLengthPoints+(2*params.SegBufferPoints);

      errcode=XLALREAL4AverageSpectrumMedian(averagePSD,
					     dataSet,
					     params.SegLengthPoints,
					     stride,
					     windowPSD,
					     averagePSDPlan);
      if (errcode !=0)
	{
	  fprintf(stderr,"Problem calculating average PSD using median method.\n");
	  fprintf(stderr," Dataset :%i\n PSD :%i\n Seglen :%i\n Seglen + Buff :%i\n Stride :%i\n",
		  dataSet->data->length,
		  averagePSD->data->length,
		  params.SegLengthPoints,
		  params.SegLengthPoints+(2*params.SegBufferPoints),
		  stride);
	  exit(errcode);
	}
      break;
 
    case useUnity:
      fprintf(stderr,"This method has been removed!\n");
      exit(TRACKSEARCHC_EVAL);
      break;

    default:
	fprintf(stderr,"Method requested to estimate PSD is not allowed.\n");
      exit(TRACKSEARCHC_EVAL);
    }

  if (averagePSDPlan)
    XLALDestroyREAL4FFTPlan(averagePSDPlan);

  if (params.smoothAvgPSD > 2)
    {
      /* Error check that median block size less than =
	 length of power psd
      */
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"We will do the running median smoothing of the average PSD with block size of %i\n",
		  params.smoothAvgPSD);
	  fflush(stdout);
	}
      if (params.verbosity >= printFiles)
	{
	  print_real4fseries(averagePSD,"Pre_SmoothingAveragePSD.diag");
	  print_lalUnit(averagePSD->sampleUnits,
			"Pre_SmoothingAveragePSD_Units.diag");
	}
      LALappsSmoothWithRunningMedian(averagePSD->data,(UINT4) params.smoothAvgPSD);
    }
  if (params.verbosity == printFiles)
    {
      print_real4tseries(dataSet,"timeDomainAllDataSet.diag");
      print_real4fseries(averagePSD,"Post_SmoothingAveragePSD.diag");
      print_lalUnit(averagePSD->sampleUnits,"Post_SmoothingAveragePSD_Units.diag");
    }
  /*
   * Setup global FFT plans 
   */
  planLength=params.SegLengthPoints+(2*params.SegBufferPoints);
  
  forwardPlan=XLALCreateForwardREAL4FFTPlan(planLength,
					    0);

  reversePlan=XLALCreateReverseREAL4FFTPlan(planLength,
					    0);

  /*
   * Allocate Frequency Struture to use as memory space for 
   * subsequent ffts of segments
   */
  signalFFT=XLALCreateCOMPLEX8FrequencySeries("tmpSegPSD",
					      &gps_zero,
					      0,
					      1/(dataSegments->dataSeg[0]->deltaT * dataSegments->dataSeg[0]->data->length),
					      &(dataSegments->dataSeg[0]->sampleUnits),
					      planLength/2+1);
  /*
   * Grab original units for manipulations later
   */
  originalFrequecyUnits=signalFFT->sampleUnits;

  /* Sanity check those units above should be counts or strain */
  /* 
   * Actually whiten each data segment
   */
  if (params.verbosity >= printFiles)
    {
      tmpSignalPtr=NULL;
      for (i=0;i<dataSegments->length;i++)
	{
	  tmpSignalPtr=(dataSegments->dataSeg[i]);
	  dataLabel=XLALCreateCHARVector(128);
	  sprintf(dataLabel->data,"Pre_WhitenTimeDomainDataSeg_%i.diag",i);
	  print_real4tseries(tmpSignalPtr,dataLabel->data);
	  sprintf(dataLabel->data,"Pre_WhitenTimeDomainDataSeg_%i_Units.diag",i);
	  print_lalUnit(tmpSignalPtr->sampleUnits,dataLabel->data);
	  XLALDestroyCHARVector(dataLabel);
	}
      tmpSignalPtr=NULL;

      for (i=0;i<dataSegments->length;i++)
	{
	  tmpSignalPtr = (dataSegments->dataSeg[i]);
	  /*
	   * FFT segment
	   */
	  errcode = XLALREAL4ForwardFFT(signalFFT->data,
					tmpSignalPtr->data,
					forwardPlan);
	  if (errcode !=0)
	    {
	      fprintf(stderr,"Error trying to fft input!\n");
	      exit(errcode);
	    }
	  /*
	   * Whiten
	   */
	  if (params.verbosity >= printFiles)
	    {
	      dataLabel=XLALCreateCHARVector(128);
	      sprintf(dataLabel->data,"Pre_whitenSignalFFT_%i.diag",i);
	      print_complex8fseries(signalFFT,dataLabel->data);
	      XLALDestroyCHARVector(dataLabel);
	    }
	  if (params.verbosity >= verbose)
	    {
	      fprintf(stdout,"Whitening data segment: %i of %i\n",i+1,dataSegments->length);
	      fflush(stdout);
	    }
	  LAL_CALL(LALTrackSearchWhitenCOMPLEX8FrequencySeries(status,
							       signalFFT,
							       averagePSD,
							       params.whiten),
		   status);
	  if (params.verbosity >= printFiles)
	    {
	      dataLabel=XLALCreateCHARVector(128);
	      sprintf(dataLabel->data,"Post_whitenSignalFFT_%i.diag",i);
	      print_complex8fseries(signalFFT,dataLabel->data);
	      XLALDestroyCHARVector(dataLabel);
	    }
		
	    
	  /*
	   * Reverse FFT
	   */
	  errcode=XLALREAL4ReverseFFT(tmpSignalPtr->data,
				      signalFFT->data,
				      reversePlan);
	  if (errcode != 0)
	    {
	      fprintf(stderr,"Error trying to inverse fft input!\n");
	      exit(errcode);
	    }

	  /* 
	   * Normalize the IFFT by 1/n factor
	   * See lsd-5 p259 10.1 for explaination
	   */
	  for (j=0;j<tmpSignalPtr->data->length;j++)
	    tmpSignalPtr->data->data[j]= 
	      tmpSignalPtr->data->data[j]/tmpSignalPtr->data->length;
	}
      /*
       * Reset the tmp frequency series units
       */
      signalFFT->sampleUnits=originalFrequecyUnits;
      tmpSignalPtr=NULL;
      
      /*
       * Free temporary Frequency Series
       */
      if (signalFFT)
	XLALDestroyCOMPLEX8FrequencySeries(signalFFT);

      if (windowPSD)
	XLALDestroyREAL4Window(windowPSD);

      if (averagePSD)
	XLALDestroyREAL4FrequencySeries(averagePSD);

      if (forwardPlan)
	XLALDestroyREAL4FFTPlan(forwardPlan);

      if (reversePlan)
	XLALDestroyREAL4FFTPlan(reversePlan);

    }
  if (params.verbosity >= printFiles)
    {
      tmpSignalPtr=NULL;
      for (i=0;i<dataSegments->length;i++)
	{
	  tmpSignalPtr=(dataSegments->dataSeg[i]);
	  dataLabel=XLALCreateCHARVector(128);
	  sprintf(dataLabel->data,"Post_WhitenTimeDomainDataSeg_%i.diag",i);
	  print_real4tseries(tmpSignalPtr,dataLabel->data);
	  sprintf(dataLabel->data,"Post_WhitenTimeDomainDataSeg_%i_Units.diag",i);
	  print_lalUnit(tmpSignalPtr->sampleUnits,dataLabel->data);
	  XLALDestroyCHARVector(dataLabel);
	}
    }
  return;
}
/*
 * End Whiten data segments
 */

int LALappsQuickHeterodyneTimeSeries(REAL4TimeSeries *data,
				     REAL8            fHet)
{
  REAL4 a=0;
  REAL4 b=0;
  REAL4 c=0;
  REAL4 d=0;
  REAL4 t=0;
  COMPLEX8 z;
  UINT4 index=0;

  if ((1/data->deltaT) < 1/fHet)
    {
      fprintf(stderr,"Error the heterodyne frequency %9.2f Hz exceeds the input data sampling rate %9.2f.\n",fHet,(1/data->deltaT));
      return TRACKSEARCHC_ESUB;
    }

  for (index=0;index<data->data->length;index++)
    {
      t=index*data->deltaT;
      a=data->data->data[index];
      b=0;
      c=cos(2*LAL_PI*fHet*t);
      d=sin(2*LAL_PI*fHet*t);
      z.re=((a*c)-(b*d));
      z.im=((b*d)+(b*c));
      /*
       *Keeping only REAL(Data*LocalOscillator)
       */
      data->data->data[index]=z.re;
    }
  return 0;
}

int LALappsSmoothWithRunningMedian(REAL4Vector *data,
				   UINT4   blocksize)
{
  UINT4                halfBlock=0;
  REAL4Vector         *smoothedDataVector=NULL;
  REAL4Vector         *tmpExtendedDataVector=NULL;
  REAL8                medianBiasFactor=0;
  LALRunningMedianPar  smoothingMedianParams;
  UINT4                i=0;
  LALStatus            status;
  UINT4                errcode=0;
  /*
   * Error check that median block size < data->length
   *
   * Perform running median on data with blocksize but
   * extend the data vector by floor(blocksize/2) to 
   * wrap data for smoother median estimates
   */
  memset(&status, 0, sizeof(status));
  halfBlock=(UINT4) floor(blocksize/2);
  smoothedDataVector=XLALCreateREAL4Vector(data->length+halfBlock);
  tmpExtendedDataVector=XLALCreateREAL4Vector(smoothedDataVector->length+blocksize-1);
  /*
   * Use tmp vector to make a buffer of smoothed input data vector
   *
   * Insert time --> The point shift equals nPoint \def block/2
   * points wraped on input data vector
   */
  for (i=0;i<halfBlock;i++)
    tmpExtendedDataVector->data[i]=data->data[halfBlock-i-1];
  for (i=0;i<data->length;i++)
    tmpExtendedDataVector->data[i+halfBlock]=data->data[i];
  for (i=0;i<blocksize-1;i++)
    tmpExtendedDataVector->data[i+data->length+halfBlock]=
      data->data[data->length-i-1];
  /*
   * Setup running median parameters
   */
  smoothingMedianParams.blocksize=blocksize;
  errcode=LAL_CALL(LALSRunningMedian(&status,
			     smoothedDataVector,
			     tmpExtendedDataVector,
			     smoothingMedianParams),
	   &status);
  if (errcode!=0)
    {
      fprintf(stderr,"Error calling LALSRunningMedian\n");
      fprintf(stderr,"%s\n",status.statusDescription);
      fflush(stderr);
      return TRACKSEARCHC_ESUB;
    }
  errcode=LAL_CALL(LALRngMedBias(&status,
			 &medianBiasFactor,
			 blocksize),
	   &status);
  if (errcode!=0)
    {
      fprintf(stderr,"Error calling LALRngMedBias\n");
      fprintf(stderr,"%s\n",status.statusDescription);
      fflush(stderr);
      return TRACKSEARCHC_ESUB;
    }

  /*
   * Fix possible bias of the running median
   */
  for (i=0;i<smoothedDataVector->length;i++)
    smoothedDataVector->data[i]=
      medianBiasFactor*smoothedDataVector->data[i];
  /*
   * Remove the shifted data points from the tmpExtendedDataVector
   */
  for (i=0;i<data->length;i++)
    data->data[i]=smoothedDataVector->data[i];
  if (smoothedDataVector)
    XLALDestroyREAL4Vector(smoothedDataVector);

  if (tmpExtendedDataVector)
    XLALDestroyREAL4Vector(tmpExtendedDataVector);

  return 0;
}

/*
 * End of Semi Private functions
 */


/****************************************************************************/

/* ********************************************************************** */

/*
 * These are local scope subroutines mean for use 
 * by tracksearch only if it may be globally useful then we will place
 * it in the LAL libraries.
 */


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
  fprintf(fp,"Data Sampling Rate:            %f\n",params.SamplingRate);
  fprintf(fp,"TF Transform Type used:        %d\n",params.TransformType);
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
  int                errcode = 0;

  INITSTATUS (status, "makefakenoise", TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  seed = 0;
  fd=XLALCreateREAL4Vector(pointnum);

  RP=XLALCreateRandomParams(seed);
  errcode=XLALNormalDeviates(fd,RP);
  if (errcode != 0)
    {
      fprintf(stderr,"Error generating random data for injection.\n");
      exit(errcode);
    }
  XLALDestroyRandomParams(RP);
  fakeData->data = NULL;
  fakeData->data=XLALCreateREAL4Vector(pointnum);
  for (k = 0;k < ((INT4) fd->length);k++)
    {
      fakeData->data->data[k] = fd->data[k];
    }
  XLALDestroyREAL4Vector(fd);
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
 
/*
 * End local scratch functions
 * These functions may eventually be integrated more into
 * lalapps tracksearch code
 */

/****************************************************************************/
