/*
 * Author: Torres C. (Univ of TX at Brownsville)
 */

#include "tracksearch.h"

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



/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */

/* TEMPORARY */
/* Non Compliant code taken from EPSearch.c */
static void print_real4tseries(const REAL4TimeSeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALSPrintTimeSeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  LALStatus   status=blank_status;
  REAL8   timeT;
  size_t i;
  LAL_CALL(LALGPStoFloat(&status,
			 &timeT,
			 &(fseries->epoch)),
	   &status);
  if(fp) 
    {
      for(i = 0; i < fseries->data->length; i++)
	fprintf(fp, "%f\t%g\n", (i * fseries->deltaT)+timeT, fseries->data->data[i]);
      fclose(fp);
    }
#endif
}

static void print_real4fseries(const REAL4FrequencySeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALCPrintFrequencySeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  size_t i;

  if(fp) {
    for(i = 0; i < fseries->data->length; i++)
	fprintf(fp, "%f\t%g\n", (i * fseries->deltaF), fseries->data->data[i]);
    fclose(fp);
  }
#endif
}

static void print_complex8fseries(const COMPLEX8FrequencySeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALCPrintFrequencySeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  size_t i;

  if(fp) {
    for(i = 0; i < fseries->data->length; i++)
      fprintf(fp, "%f\t%g\n", i * fseries->deltaF, sqrt(fseries->data->data[i].re * fseries->data->data[i].re + fseries->data->data[i].im * fseries->data->data[i].im));
    fclose(fp);
  }
#endif
}

static void print_complex8_RandC_fseries(const COMPLEX8FrequencySeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALCPrintFrequencySeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  size_t i;

  if(fp) {
    for(i = 0; i < fseries->data->length; i++)
      fprintf(fp, "%f\t%g\t%g\n", i * fseries->deltaF, 
	      fseries->data->data[i].re,
	      fseries->data->data[i].im);
    fclose(fp);
  }

#endif
}

static void print_lalUnit(LALUnit unit,const char *file)
{
  FILE *fp = fopen(file,"w");
  CHARVector *unitString=NULL;
  LALStatus  status=blank_status;

  LAL_CALL(LALCHARCreateVector(&status,&unitString,maxFilenameLength),&status);
  LAL_CALL(LALUnitAsString(&status,unitString,&unit),&status);
  if (fp)
    fprintf(fp,"%s\n",unitString->data);
  fclose(fp);
  LAL_CALL(LALCHARDestroyVector(&status,&unitString),&status);
}
/*
 * End diagnostic code
 */

/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */



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
  LALStatus         status;/* Error containing structure */
  TSSearchParams   *params;/*Large struct for all params*/
  TSappsInjectParams *injectParams;/*Struct for performing injects*/
  CHARVector       *cachefile=NULL;/* name of file with frame cache info */
  CHAR            *cachefileDataField=NULL;/*field of cachefile*/
  CHARVector       *dirpath=NULL;

  /*
   *Sleep for Attaching DDD 
   */
  unsigned int doze = 8;
  pid_t myPID;
  myPID = getpid( );
  fprintf( stdout, "pid %d sleeping for %d seconds\n", myPID, doze );
  fflush( stdout );
  sleep( doze );
  fprintf( stdout, "pid %d awake\n", myPID );
  fflush( stdout );

  /* End Global Variable Declarations */

  /* SET LAL DEBUG STUFF */
  set_debug_level("NONE");
  memset(&status, 0, sizeof(status));
  lal_errhandler = LAL_ERR_EXIT;
  lal_errhandler = LAL_ERR_RTRN;
  set_debug_level("MEMDBG");
  /*  set_debug_level("ERROR | WARNING | MEMDBG");*/
  /*  set_debug_level("ERROR | WARNING | MEMDBG");*/
  /*  set_debug_level("ERROR | WARNING ");*/
  /*    set_debug_level("ALLDBG");*/

  /*
   * Initialize status structure 
   */
  
  params=LALMalloc(sizeof(TSSearchParams));
  injectParams=LALMalloc(sizeof(TSappsInjectParams));

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
      LALappsDoTimeSeriesAnalysis(&status,
				  *params,
				  *injectParams,
				  cachefileDataField,
				  dirpath);
    }
  else
    {
      /*
       * Assume we have preformed maps to use and do a map
       * analysis instead
       */
      LALappsDoTSAMapAnalysis(&status,*params);
    }
  /* 
   * Free searchParams input arguement structure from InitializeRoutine
   * Free elements first 
   */

  if (dirpath)
    LAL_CALL(LALCHARDestroyVector(&status,&dirpath),&status);
  if (cachefile)
    LAL_CALL(LALCHARDestroyVector(&status,&cachefile),&status);
  if (params->channelName)
    LALFree(params->channelName);
  if (params->auxlabel)
    LALFree(params->auxlabel);
  if (params->numSlaves)
    LALFree(params->numSlaves);
  if (params->injectMapCache)
    LALFree(params->injectMapCache);
  if (params->injectSingleMap)
    LALFree(params->injectSingleMap);
  if (params)
    LALFree(params);   
  if (injectParams->injectTxtFile)
    LAL_CALL(LALCHARDestroyVector(&status,&(injectParams->injectTxtFile)),
	     &status);
  if (injectParams)
    LALFree(injectParams);
  /* 
   * Done freeing for search params struct memory
   */
  LALCheckMemoryLeaks();
  return 0;
}
/* End Main Code Body */


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* ********************************************************************** */
/*
 * These are the SemiPrivate functions for LALapps tracksearch
 * Useful and permanent to local lalapps code
 *
 */

/*
 * Function that prepared data for analysis
 * Can do whitening and calibration on the segments
 */
void LALappsTrackSearchPrepareData( LALStatus*        status,
				    REAL4TimeSeries  *dataSet,
				    REAL4TimeSeries  *injectSet,
				    TSSegmentVector  *dataSegments,
				    TSSearchParams    params)
     /* Add option NULL or with data called REAL4TimeSeries injectSet */
{
  AverageSpectrumParams     avgPSDParams;
  CHARVector               *dataLabel=NULL;
  COMPLEX8FrequencySeries  *response=NULL;
  COMPLEX8FrequencySeries  *signalFFT=NULL;
  CalibrationUpdateParams   calfacts;
  FrCache                  *calcache = NULL;
  LALUnit                   originalFrequecyUnits;
  LALWindowParams           windowParamsPSD;
  REAL4FFTPlan             *forwardPlan=NULL;
  REAL4FFTPlan             *reversePlan=NULL;
  REAL4FFTPlan             *averagePSDPlan=NULL;
  REAL4FFTPlan             *dataSetPlan=NULL;
  COMPLEX8FrequencySeries  *dataSetFFT=NULL;
  REAL4FrequencySeries     *averagePSD=NULL;
  REAL4Vector              *smoothedAveragePSD=NULL;
  REAL4Vector              *tmpExtendedAveragePSD=NULL;
  REAL8                     smoothingAveragePSDBias=0;
  LALRunningMedianPar       smoothingPSDparams;
  REAL4TimeSeries          *tmpSignalPtr=NULL;
  REAL4Window              *windowPSD=NULL;
  UINT4                     i=0;
  UINT4                     j=0;
  UINT4                     segmentPoints=0;
  UINT4                    planLength=0;
  PassBandParamStruc       bandPassParams;
  const LALUnit     strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;


  if (params.verbosity == printFiles)
    {
      print_lalUnit(dataSet->sampleUnits,"EntireInputDataSet_Units.diag");
      print_real4tseries(dataSet,"EntireInputDataSet.diag");
    }

  /*
   * Calibrate entire input data set
   */
  if (params.calibrate)
    {
      if (params.verbosity >= verbose)
	fprintf(stdout,"Performing calibration on input data\n");
      /* 
       * Create response function based on starting epoch from 
       * first data segment.  We assume it is constant over 
       * all the other data segments available.
       */
      segmentPoints=params.SegLengthPoints;
      /* 
       * Implicity set Epoch, Unit and df for Response extract
       */
      LAL_CALL( LALCreateCOMPLEX8FrequencySeries(status,
						 &response,
						 "tmpResponse",
						 dataSet->epoch,
						 0,
						 1/(dataSet->deltaT*dataSet->data->length),
						 strainPerCount,
						 segmentPoints/2+1),
		status);

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
      if (params.verbosity == printFiles)
	{
	  print_lalUnit(response->sampleUnits,"responseFunction_Units.diag");
	  print_complex8fseries(response,"responseFunction.diag");
	}
      /*
       * Destroy IFO text pointer
       */
      LALFree(calfacts.ifo);
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
      LAL_CALL(
	       LALCreateCOMPLEX8FrequencySeries(status,
						&dataSetFFT,
						"Entire data FFTed",
						dataSet->epoch,
						0,
						(dataSet->deltaT*dataSet->data->length),
						dataSet->sampleUnits,
						dataSet->data->length/2+1),
	       status);
      /*
       * FFT input unsegmented data set
       */
      LAL_CALL( LALCreateForwardREAL4FFTPlan(status,
					     &dataSetPlan,
					     dataSet->data->length,
					     0),
		status);
      LAL_CALL( LALForwardREAL4FFT(status,
				   dataSetFFT->data,
				   dataSet->data,
				   dataSetPlan),
		status);
      if (dataSetPlan)
	LAL_CALL(LALDestroyREAL4FFTPlan(status,&dataSetPlan),status);

      if (params.verbosity == printFiles)
	{
	  print_complex8fseries(dataSetFFT,"dataSetFFT.diag");
	  print_lalUnit(dataSetFFT->sampleUnits,"dataSetFFT_Units.diag");
	}
      /*
       * Calibrate
       */
      LAL_CALL(LALTrackSearchCalibrateCOMPLEX8FrequencySeries(status,
							      dataSetFFT,
							      response),
	       status);

      if (params.verbosity == printFiles)
	{
	  print_complex8fseries(dataSetFFT,"dataSetFFTCalibrate.diag");
	  print_lalUnit(dataSetFFT->sampleUnits,"dataSetFFTCalibrate_Units.diag");
	}
      /*
       * Inverse FFT the entire dataSet
       */
      fprintf(stderr,"Hack on Response function...Don't forget\n");
      fprintf(stderr,"Will cause memory leak for some reason?\n");
      dataSetFFT->data->data[dataSetFFT->data->length-1].im=0;
      dataSetFFT->data->data[dataSetFFT->data->length].im=0;
      
      LAL_CALL(LALCreateReverseREAL4FFTPlan(status,
					    &dataSetPlan,
					    dataSet->data->length,
					    0),
	       status);
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
      
      if (params.verbosity == printFiles)
	{
	  print_real4tseries(dataSet,"dataSetCalibrate.diag");
	  print_lalUnit(dataSet->sampleUnits,"dataSetCalibrate_Units.diag");
	}
      /*
       * Free the FFT memory space
       */
	if (dataSetFFT)
	LAL_CALL(LALDestroyCOMPLEX8FrequencySeries(status,dataSetFFT),
	status);
      if (response)
	LAL_CALL(LALDestroyCOMPLEX8FrequencySeries(status,response),
		 status);
    }
  else
    {
      if (params.verbosity >= verbose)
	fprintf(stdout,"Calibration not requested, assuming no calibration needed\n");
    }
  /* 
   *End calibration conditional 
   */

  /* STRIP OUT LATER */
  /* Set all data points to 1 and see what happens */
  /*  for (j=0;j<dataSet->data->length;j++)*/
  /*    dataSet->data->data[j]=1;*/

  /*
   * Perform low pass filtering on the input data stream if requested
   */
  if ((params.highPass > 0) || (params.lowPass > 0))
    {
      if (params.verbosity >= verbose)
	{
	  fprintf(stdout,"You requested a high pass filter of the data at %f Hz\n",params.highPass);
	  fprintf(stdout,"You requested a low pass filter of the data at %f Hz\n",params.lowPass);
	}
      bandPassParams.name=NULL;
      bandPassParams.nMax=10;
      /* F < f1 kept, F > f2 kept */
      bandPassParams.f1=params.lowPass;
      bandPassParams.f2=params.highPass;
      bandPassParams.a1=0.9;
      bandPassParams.a2=0.9;
      /*
       * Call the low pass filter function.
       */
      LAL_CALL(LALButterworthREAL4TimeSeries(status,
					     dataSet,
					     &bandPassParams),
	       status);
      if (params.verbosity >= printFiles)
	print_real4tseries(dataSet,"ButterworthFiltered.diag");
    }
  /*
   * End the high pass filter
   */
  /*
   * Add in injections if injection variable is not NULL copy up to 
   * end of dataSet series
   */
  if (injectSet != NULL)
    {
      if (injectSet->data->length >= dataSet->data->length)
	for (j=0;j<dataSet->data->length;j++)
	  dataSet->data->data[j]=
	    dataSet->data->data[j]+injectSet->data->data[j];
      else
	for (j=0;j<injectSet->data->length;j++)
	  dataSet->data->data[j]=
	    dataSet->data->data[j]+injectSet->data->data[j];
    }
  /*
   * Split incoming data into Segment Vector
   */
  LAL_CALL(LALTrackSearchDataSegmenter(status,
				       dataSet,
				       dataSegments,
				       params),
	   status);
  /*
   * If injection requested inject a single waveform into each segment
   * scaled to fit the segment sample rate and use specified scale factor
   * We don't bandpass the injected waveform!!
   */
  /*
   * If we are to whiten first let us calculate the 
   * average PSD using the non segment data structure
   */

  if (params.whiten !=0 )
    {
      if (params.verbosity >= verbose )
	fprintf(stdout,"Estimating PSD for whitening\n");
      LAL_CALL(
	       LALCreateREAL4FrequencySeries(status,
					     &averagePSD,
					     "averagePSD",
					     gps_zero,
					     0,
					     1/((params.SegLengthPoints)*dataSet->deltaT),
					     lalDimensionlessUnit,
					     params.SegLengthPoints/2+1),
	       status);
      /*
       * The given units above need to be correct to truly reflect the
       * units stored in the frame file
       */
      windowParamsPSD.length=params.SegLengthPoints;
      windowParamsPSD.type=params.avgSpecWindow;
      LAL_CALL( LALCreateREAL4Window(status,
				     &(windowPSD),
				     &windowParamsPSD),
		status);
      /* If we only do one map we need to something special here*/
      if (params.NumSeg < 2)
	avgPSDParams.overlap=1;
      else /* use same as segment overlap request*/
	avgPSDParams.overlap=params.overlapFlag;

      avgPSDParams.window=windowPSD;
      avgPSDParams.method=params.avgSpecMethod;
      /* Set PSD plan length and plan*/
      LAL_CALL( LALCreateForwardREAL4FFTPlan(status,
					     &averagePSDPlan,
					     params.SegLengthPoints,
					     0),
		status);
      avgPSDParams.plan=averagePSDPlan;

      LAL_CALL(  LALREAL4AverageSpectrum(status,
					 averagePSD,
					 dataSet,
					 &avgPSDParams),
		status);
      if (params.smoothAvgPSD > 2)
	{
	  /* Error check that median block size less than =
	     length of power psd
	  */
	  if (params.verbosity >= verbose)
	    {
	      fprintf(stdout,"We will do the running median smoothing of the average PSD with block size of %i\n",
		      params.smoothAvgPSD);
	    }
	  if (params.verbosity >= printFiles)
	    {
	      print_real4fseries(averagePSD,"PreSmoothing-averagePSD.diag");
	      print_lalUnit(averagePSD->sampleUnits,
			    "PreSmoothing-averagePSD_Units.diag");
	    }
	  
	  /*
	   * Perform running median on the average PSD using
	   * blocks of size n
	   */
	  LAL_CALL(LALSCreateVector(status,
				    &smoothedAveragePSD,
				    averagePSD->data->length),
		   status);
	  LAL_CALL(LALSCreateVector(status,
				    &tmpExtendedAveragePSD,
				    averagePSD->data->length+params.smoothAvgPSD-1),
		   status);
	  /*
	   * Build a buffered PSD rep so that median estimate
	   * returned will have m points equal to the n points of the 
	   * original average PSD
	   */
	  for (i=0;i<averagePSD->data->length;i++)
	    tmpExtendedAveragePSD->data[i]=averagePSD->data->data[i];
	  for (i=averagePSD->data->length-1;
	       i<tmpExtendedAveragePSD->length;
	       i++)
	    tmpExtendedAveragePSD->data[i]=averagePSD->data->data[averagePSD->data->length-1];
	  /*
	   * Setup running median parameters
	   */
	  smoothingPSDparams.blocksize=params.smoothAvgPSD;
	  LAL_CALL(LALSRunningMedian(status,
				     smoothedAveragePSD,
				     tmpExtendedAveragePSD,
				     smoothingPSDparams),
		   status);
	  LALPrintVector(tmpExtendedAveragePSD);
	  LALPrintVector(smoothedAveragePSD);
	  LAL_CALL(LALRngMedBias(status,
				 &smoothingAveragePSDBias,
				 smoothingPSDparams.blocksize),
		   status);
	  /*
	   * Fix possible bias of the running median
	   */
	  for (i=0;i<smoothedAveragePSD->length;i++)
	    smoothedAveragePSD->data[i]=
	      smoothingAveragePSDBias*smoothedAveragePSD->data[i];
	  /*
	   * Copy newly smoothed PSD to frequency series
	   */
	  for (i=0;i<averagePSD->data->length;i++)
	   averagePSD->data->data[i]=smoothedAveragePSD->data[i];

	  if (smoothedAveragePSD)
	    LAL_CALL(LALSDestroyVector(status,
				       &smoothedAveragePSD),
		     status);

	  if (tmpExtendedAveragePSD)
	    LAL_CALL(LALSDestroyVector(status,
				       &tmpExtendedAveragePSD),
		     status);
	  
	  
	}
      if (params.verbosity == printFiles)
	{
	  print_real4tseries(dataSet,"dataSet.diag");
	  print_real4fseries(averagePSD,"averagePSD.diag");
	  print_lalUnit(averagePSD->sampleUnits,"averagePSD_Units.diag");
	}
      /*
       * Setup global FFT plans 
       */
      planLength=params.SegLengthPoints;
      LAL_CALL(  LALCreateForwardREAL4FFTPlan(status,
					      &forwardPlan,
					      planLength,
					      0),
		 status);
      LAL_CALL( LALCreateReverseREAL4FFTPlan(status,
					     &reversePlan,
					     planLength,
					     0),
		status);
      /*
       * Allocate Frequency Struture to use as memory space for 
       * subsequent ffts of segments
       */
      LAL_CALL( 
	       LALCreateCOMPLEX8FrequencySeries(
						status,
						&signalFFT,
						"tmpSegPSD",
						gps_zero,
						0,
						1/(dataSegments->dataSeg[0]->deltaT * dataSegments->dataSeg[0]->data->length),
						dataSegments->dataSeg[0]->sampleUnits,
						planLength/2+1),
	       status);
      /*
       * Grab original units for manipulations later
       */
      originalFrequecyUnits=signalFFT->sampleUnits;

      /* Sanity check those units above should be counts or strain */
      /* 
       * Actually whiten each data segment
       */
      for (i=0;i<dataSegments->length;i++)
	{
	  tmpSignalPtr = (dataSegments->dataSeg[i]);
	  if (params.verbosity == printFiles)
	    {
	      print_lalUnit(tmpSignalPtr->sampleUnits,"InputTimeSeries_Units.diag");
	    }
	  /*
	   * FFT segment
	   */
	  LAL_CALL( LALForwardREAL4FFT(status,
				       signalFFT->data,
				       tmpSignalPtr->data,
				       forwardPlan),
		    status);
	  if (params.verbosity == printFiles)
	    {
	      print_complex8fseries(signalFFT,"signalFFT.diag");
	      print_lalUnit(signalFFT->sampleUnits,"signalFFT_Units.diag");
	    }
	  /*
	   * Whiten
	   */
	  if (params.whiten != 0)
	    {
	      LAL_CALL(LALTrackSearchWhitenCOMPLEX8FrequencySeries(status,
								   signalFFT,
								   averagePSD,
								   params.whiten),
		       status);
	      if (params.verbosity == printFiles)
		{
		  print_complex8fseries(signalFFT,"signalFFTWhiten.diag");
		}
	    }
	  /*
	   * Reverse FFT
	   */
	  LAL_CALL( LALReverseREAL4FFT(status,
				       tmpSignalPtr->data,
				       signalFFT->data,
				       reversePlan),
		    status);
	  /* 
	   * Normalize the IFFT by 1/n factor
	   * See lsd-5 p259 10.1 for explaination
	   */
	  for (j=0;j<tmpSignalPtr->data->length;j++)
	    tmpSignalPtr->data->data[j]= 
	      tmpSignalPtr->data->data[j]/tmpSignalPtr->data->length;
	  if (params.verbosity == printFiles)
	    {
	      LAL_CALL(LALCHARCreateVector(status,&dataLabel,128),
		       status);
	      sprintf(dataLabel->data,"signalFFT_%i.diag",i);
	      print_complex8fseries(signalFFT,dataLabel->data);
	      LAL_CALL(LALCHARDestroyVector(status,&dataLabel),
		       status);
	    }
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
	{
	  LAL_CALL( LALDestroyCOMPLEX8FrequencySeries(status,signalFFT),
		    status);
	}
      if (windowPSD)
	LAL_CALL( LALDestroyREAL4Window(status,
					&windowPSD),
		  status);
      if (averagePSDPlan)
	LAL_CALL( LALDestroyREAL4FFTPlan(status,&averagePSDPlan),
		  status);
      if (averagePSD)
	LAL_CALL(LALDestroyREAL4FrequencySeries(status,averagePSD),
		 status);
      if (forwardPlan)
	{
	  LAL_CALL(
		   LALDestroyREAL4FFTPlan(status,&forwardPlan),
		   status);
	}
      if (reversePlan)
	{
	  LAL_CALL( LALDestroyREAL4FFTPlan(status,&reversePlan),
		    status);
	}
    }
  if (params.verbosity == printFiles)
    {
      tmpSignalPtr=NULL;
      for (i=0;i<dataSegments->length;i++)
	{
	  tmpSignalPtr=(dataSegments->dataSeg[i]);
	  LAL_CALL(LALCHARCreateVector(status,&dataLabel,128),
		   status);
	  sprintf(dataLabel->data,"dataSeg_%i.diag",i);
	  print_real4tseries(tmpSignalPtr,dataLabel->data);
	  sprintf(dataLabel->data,"dataSeg_%i_Units.diag",i);
	  print_lalUnit(tmpSignalPtr->sampleUnits,dataLabel->data);
	  LAL_CALL(LALCHARDestroyVector(status,&dataLabel),
		   status);
	}
    }
  return;
}
/* End the PrepareData routine */

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
      {0,                     0,                  0,      0}
    };
  
  int              C;

  /*  INITSTATUS (status, "getparams", TRACKSEARCHC);*/
  /*  ATTATCHSTATUSPTR (status);*/

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
  params->overlapFlag = 0; /* Change this to an INT4 for overlapping */
  params->TimeLengthPoints = 0;
  params->discardTLP = 0;
  params->window = Hann; /*Welch*/
  params->whiten = 0;
  params->avgSpecMethod = -1; /*Will trigger error*/
  params->avgSpecWindow = Rectangular; /* Default no window */
  params->smoothAvgPSD = 0;/*0 means no smoothing*/
  params->highPass=0;/*High pass filter freq*/
  params->lowPass=0;/*Low pass filter freq*/
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
  params->SamplingRateOriginal = params->SamplingRate;
  params->makenoise = -1;
  params->calChannelType=SimDataChannel;/*See lsd*/
  params->channelName=NULL; /* Leave NULL to avoid frame read problems */
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
	    LAL_CALL(
		     LALFloatToGPS(status,
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
	    params->window = atoi(optarg);
	    /* Error check FIX THIS LATER*/
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
	    /*Setup structure for path or cachefile name temp*/
	    LAL_CALL(LALCHARCreateVector(status,dPtrDirPath,maxFilenameLength),
		     status);

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
	    LAL_CALL(LALCHARCreateVector(status,dPtrCachefile,maxFilenameLength),
		     status);
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
	      LALCalloc(len,sizeof(CHAR));
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
	    params->calFrameCache = (CHAR *) LALCalloc(len,sizeof(CHAR));
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
	    params->pgmColorMapFile = (CHAR *) LALCalloc(len,sizeof(CHAR));
	    memcpy(params->pgmColorMapFile,optarg,len);
	  }
	case 'C':
	  {
	    params->injectMapCache=
	      (CHAR*) LALCalloc(strlen(optarg)+1,sizeof(CHAR));
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
	      (CHAR*) LALCalloc(strlen(optarg)+1,sizeof(CHAR));
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
	    LAL_CALL(LALCHARCreateVector(status,
					 &(injectParams->injectTxtFile),
					 strlen(optarg)+1),
		     status);
	    memcpy(injectParams->injectTxtFile->data,optarg,strlen(optarg)+1);
	  }
	  break;

	case 'N':
	  {
	    injectParams->scaleFactor=atof(optarg);
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
	}
      else
	{
	  /* Determine the number of maps */
	  params->NumSeg=floor((params->TimeLengthPoints-params->overlapFlag)
			       /(params->SegLengthPoints - params->overlapFlag));
	  /* Determine number of points to throw away (not process) */
	  /*params->discardTLP=floor(params->NumSeg*
	   *		       (params->SegLengthPoints-params->overlapFlag)
	   *		       /(params->TimeLengthPoints-params->SegLengthPoints));
	   */
	  params->discardTLP=(params->TimeLengthPoints-params->overlapFlag)%(params->SegLengthPoints - params->overlapFlag);
	  /* 
	   * Reset points to process by N-discardTLP, this gives us
	   * uniform segments to make maps with
	   */
	  params->TimeLengthPoints=params->TimeLengthPoints-params->discardTLP;
	  if ((params->verbosity >= verbose))
	    fprintf(stdout,"We need to throw away %i trailing data points for %i, evenly matched data segments\n",params->discardTLP,params->NumSeg);
	  
	};
      if (params->whiten > 2 )
	{
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  exit(TRACKSEARCHC_EARGS);
	}
    }
  /*
   * Do following checks
   */
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
  if (params->MinLength < 2)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEARGS);
      exit(TRACKSEARCHC_EARGS);
    }
  if ((params->LineWidth) < 1)
    {
      fprintf(stderr,TRACKSEARCHC_MSGEARGS);
      exit(TRACKSEARCHC_EARGS);
    }

  /*  DETATCHSTATUSPTR(status);*/
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
  REAL4TimeSeries      *tmpData=NULL;
  UINT4                 loadPoints=0;
  UINT4                 i=0;
  ResampleTSParams      resampleParams;

  /* Set all variables from params structure here */
  channelIn.name = params->channelName;
  memcpy(DataIn->name,
	 params->channelName,
	 (strlen(params->channelName)*sizeof(CHAR)));
 
  memset( DataIn->data->data, 0, DataIn->data->length*sizeof(REAL4) );
  /* Need to use NULL DataIN->name so we will skip this line of code
     Nullifying the DataIn->name pointer by hand hope it works */

  /* only try to load frame if name is specified */
  if (dirname || cachefile)
    {
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
      /*
       * Load the metadata to check the frame sampling rate
       */
      LAL_CALL( LALFrGetREAL4TimeSeriesMetadata( status, 
						 DataIn, 
						 &channelIn, 
						 stream), 
		status);
      /*
       * Determine the sampling rate and assuming the params.sampling rate
       * how many of these points do we load to get the same time duration
       * of data points
       */
      params->SamplingRateOriginal=1/(DataIn->deltaT);
      loadPoints=params->SamplingRateOriginal*(params->TimeLengthPoints/params->SamplingRate);
      LAL_CALL(
	       LALCreateREAL4TimeSeries(status,
					&tmpData,
					"Higher Sampling Tmp Data",
					params->GPSstart,
					0,
					1/params->SamplingRateOriginal,
					lalADCCountUnit,
					loadPoints),
	       status);
      /* get the data */
      LAL_CALL( LALFrSeek(status, &(tmpData->epoch), stream), status);
      LAL_CALL( LALFrGetREAL4TimeSeries( status, 
					 tmpData, 
					 &channelIn, 
					 stream), 
		status);
      if (params->verbosity > 0)
	print_real4tseries(tmpData,"OriginalInputTimeSeries.diag");
      /*
       * Prepare for the resample if needed or just copy the data so send
       * back
       */
      if (params->SamplingRate < params->SamplingRateOriginal)
	{
	  fprintf(stdout,"We will resample the data from %6.3f to %6.3f.\n",
		  1/tmpData->deltaT,
		  params->SamplingRate);
	  resampleParams.deltaT=1/params->SamplingRate;
	  resampleParams.filterType=defaultButterworth;
	  /*resampleParams->filterParams*/
	  LAL_CALL(
		   LALResampleREAL4TimeSeries(status,
					      tmpData,
					      &resampleParams),
		   status);
	  if (params->verbosity > 0)
	    print_real4tseries(tmpData,"ResampledlInputTimeSeries.diag");
	  /*
	   * Copy only the valid data and fill the returnable metadata
	   */
	  DataIn->deltaT=(1/params->SamplingRate);
	  LALappsTSassert((tmpData->data->length >=
			   DataIn->data->length),
			  TRACKSEARCHC_EDATA,
			  TRACKSEARCHC_MSGEDATA);
	  for (i=0;i<DataIn->data->length;i++)
	    DataIn->data->data[i]=tmpData->data->data[i];
	}
      else
	{
	  fprintf(stdout,"Resampling to %f, not possible we have %f sampling rate in the frame file.\n",
		  params->SamplingRate,
		  params->SamplingRateOriginal);
	  /*
	   * Straight copy the data over to returnable structure
	   */
	  params->SamplingRate=params->SamplingRateOriginal;
	  DataIn->deltaT=1/params->SamplingRateOriginal;
	  for (i=0;i<DataIn->data->length;i++)
	    DataIn->data->data[i]=tmpData->data->data[i];
	}
      if (params->verbosity > 0)
    	print_real4tseries(tmpData,"ActualReDoneTimeSeries.diag");
      /*
       * Release the memory for the temporary time series
       */
      if (tmpData)
	LAL_CALL(LALDestroyREAL4TimeSeries(status,tmpData),
		 status);
      /* close the frame stream */
      LAL_CALL( LALFrClose( status, &stream ), status);
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
 /*************************************************************/
  /* 
   * The LALTracksearch seems to map the 
   * width to tCol and the height to fRow
   */
  LAL_CALL(LALCHARCreateVector(status,&outputFilename,maxFilenameLength),status);
  LAL_CALL(LALCHARCreateVector(status,&outputFilenameMask,maxFilenameLength),status);
  /*
   * Setup file mask for output filenames
   */
   sprintf(outputFilenameMask->data,"CandidateList_Start_%i_%i_Stop_%i_%i",
	   tsMarkers.mapStartGPS.gpsSeconds,
	   tsMarkers.mapStartGPS.gpsNanoSeconds,
	   tsMarkers.mapStopGPS.gpsSeconds,
	   tsMarkers.mapStopGPS.gpsNanoSeconds);
  /* Flag to prepare structure inside LALTrackSearch */
  tsInputs.allocFlag = 1; 
  outputCurves.curves = NULL;
  /*
   * The power thresholding is not done in the LALroutine
   * This information is forwarded to this function for a post processing 
   * Should be default be givin minimum curve length parameter 
   * We want to apply thresholds in a seperate routine 
   */
  /* Catch the error code here and abort */
  LAL_CALL( LALSignalTrackSearch(status,&outputCurves,tfmap,&tsInputs),
	    status);

  /* 
   * Call tracksearch again to free any temporary ram in 
   * variable outputCurves which is no longer required
   */
  tsInputs.allocFlag = 2;
    LAL_CALL(  LALSignalTrackSearch(status,&outputCurves,tfmap,&tsInputs),
	     status);
  /*
   * Write PGM to disk if requested
   */
  if (params.verbosity == printFiles)
    DumpTFImage(tfmap->map,"DumpMap0",tsInputs.height,tsInputs.width,1);
  /*
   * Setup for call to function to do map marking
   * We mark maps is convert Tbin and Fbin to 
   * Hz and GPSseconds
   */
  /*Height -> time bins */
  /*Width -> freq bins */
  /* Need to configure mapMarkerParams correctly!*/
  LAL_CALL(  LALTrackSearchInsertMarkers(status,
					 &outputCurves,
					 &tsMarkers),
	     status);
  /* 
   *Call the connect curve routine if argument is specified 
   */
  if (params.joinCurves)
    {
      LAL_CALL( LALTrackSearchConnectSigma(status,
					   &outputCurves,
					   *tfmap,
					   tsInputs),
		status);
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
  outputCurvesThreshold.numberOfCurves = 0;
  LAL_CALL( LALTrackSearchApplyThreshold(status,
					 &outputCurves,
					 &outputCurvesThreshold,
					 params),
	    status);

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
    LAL_CALL(LALCHARDestroyVector(status,&outputFilename),
	     status);
  if (outputFilenameMask)
    LAL_CALL(LALCHARDestroyVector(status,&outputFilenameMask),
	     status);
  if (outputCandidateFilename)
    LAL_CALL(LALCHARDestroyVector(status,&outputCandidateFilename),
	     status);

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
  REAL8                 signalStop=0;
  REAL8                 signalStart=0;
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
      tfInputs.tCol = params.TimeBins;
      /*  tfInputs.wlengthT = params.FreqBins+1;*/
      /*  tfInputs.wlengthF = params.FreqBins+1;*/
      tfInputs.wlengthF = params.windowsize;
      tfInputs.wlengthT = params.windowsize;

      LAL_CALL(LALCreateTimeFreqParam(status,&autoparams,&tfInputs),
	       status);

      LAL_CALL(LALCreateTimeFreqRep(status,&tfmap,&tfInputs),
	       status);

      for (j=0;j<tfmap->tCol;j++)
	{
	  tfmap->timeInstant[j]=j;
	}
      windowParams.length = params.windowsize;
      windowParams.type = params.window;

      LAL_CALL(LALCreateREAL4Window(status,&(tempWindow),&(windowParams)),
	       status);

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
	    LAL_CALL( LALTfrSp(status,signalSeries->data,tfmap,autoparams),
		      status);
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
      LAL_CALL( LALDestroyREAL4Window(status,&tempWindow),
		status);
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
   * Fill in LALSignalTrackSearch params structure via the use of
   * the TSSearch huge struct elements
   */
  inputs.sigma=params.LineWidth;
  inputs.high=params.StartThresh;
  inputs.low=params.LinePThresh;
  inputs.width=((tfmap->fRow/2)+1);
  inputs.height=tfmap->tCol;
  /*
   * Fill the map marking parameter structure
   */
  mapMarkerParams.deltaT=signalSeries->deltaT;
  mapMarkerParams.mapTimeBins=inputs.height;
  mapMarkerParams.mapFreqBins=inputs.width;
  mapMarkerParams.mapStartGPS.gpsSeconds=signalSeries->epoch.gpsSeconds;
  mapMarkerParams.mapStartGPS.gpsNanoSeconds=
    signalSeries->epoch.gpsNanoSeconds;
  LAL_CALL(
	   LALGPStoFloat(status,
			 &signalStart,
			 &mapMarkerParams.mapStartGPS),
	   status);
  signalStop=(signalSeries->deltaT*signalSeries->data->length)+signalStart;
  LAL_CALL(
	   LALFloatToGPS(status,
      			 &(mapMarkerParams.mapStopGPS),
			 &signalStop),
			 
	   status);

  /*
   * Call subroutine to run the search
   */
  LALappsDoTrackSearch(status,
		       tfmap,
		       inputs,
		       mapMarkerParams,
		       params);
  /*
   * Assemble the TSAmap to write to disk
   */
 if (1==1)
    {
      LAL_CALL(LALCHARCreateVector(status,&binaryFilename,maxFilenameLength),
	       status);
      /*
       * Set filename
       */
      /* tsaMap_gpsSeconds_gpsNanoseconds.map*/
      /*
       * Assemble a tsaMap structure 
       */
      mapBuilder=LALMalloc(sizeof(TSAMap));
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
      /*
       *LALappsTSAWriteMapFile(status,
       *		     mapBuilder,
       *		     binaryFilename);
       */
      LALappsTSAWriteMapFile(status,
			     mapBuilder,
			     NULL);
      LAL_CALL(LALCHARDestroyVector(status,&binaryFilename),
	       status);
      LALFree(mapBuilder);
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
  LAL_CALL(
	   LALCreateREAL4TimeSeries(status,
				    &dataset,
				    "Uninitialized",
				    params.GPSstart,
				    0,
				    1/params.SamplingRate,
				    lalADCCountUnit,
				    params.TimeLengthPoints),
	   status);
  if (productioncode == 1) /* Use production frame file inputs */
    {
      if (params.makenoise < 1 )
	{  
	  if (params.verbosity >= verbose)
	    fprintf(stdout,"Reading frame files.\n");
	  LALappsGetFrameData(status,&params,dataset,dirpath,cachefile);
	  if (params.verbosity >= verbose)
	    fprintf(stdout,"Reading frame complete.\n");

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
      if (params.verbosity >= verbose)
	printf("Preparing the injection data!\n");
      LALappsCreateInjectableData(status,
				  &injectSet,
				  injectParams);
    }
  if ((params.verbosity >= printFiles) && (injectSet != NULL))
    print_real4tseries(injectSet,"PreparedInjectData.diag");
  /* 
   * Prepare the data for the call to Tracksearch 
   */
  SegVecParams.dataSegmentPoints = params.SegLengthPoints;
  SegVecParams.numberDataSegments = params.NumSeg;
  LAL_CALL(LALCreateTSDataSegmentVector(status,&SegVec,&SegVecParams),
	   status);

  if (params.verbosity >= printFiles)
    print_real4tseries(dataset,"InputReal4TimeSeriesPrePrepare.diag");
  LALappsTrackSearchPrepareData(status,dataset,injectSet,SegVec,params);
  if (params.verbosity >= printFiles)
    print_real4tseries(dataset,"InputReal4TimeSeriesPostPrepare.diag");

  j=0;
  for(i = 0;i < params.NumSeg;i++)
    {
      printf(".");
      /*
       * Call to prepare tSeries search
       */
      LALappsDoTSeriesSearch(status,SegVec->dataSeg[j],params,j);
      j++;
    };
  printf("\n");
  /* Free some of the memory used to do the analysis */
  if (params.dataSegVec)
    {
      LALDestroyTSDataSegmentVector(status,(params.dataSegVec));
    }
  if (SegVec)
    LAL_CALL(LALDestroyTSDataSegmentVector(status,SegVec),
	     status);
  if (dataset)
    LAL_CALL(LALDestroyREAL4TimeSeries(status,dataset),
	     status);
  if (injectSet)
    LAL_CALL(LALDestroyREAL4TimeSeries(status,injectSet),
	     status);
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
  LAL_CALL(LALCHARCreateVector(status,&mapFilenameVec,maxFilenameLength),
	   status);
  strcpy(mapFilenameVec->data,params.injectMapCache);
  LALappsTSALoadCacheFile(status,
			  mapFilenameVec,
			  &mapCache);
  LAL_CALL(LALCHARDestroyVector(status,&mapFilenameVec),
	   status);
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
  inputs.sigma=params->LineWidth;
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
  LAL_CALL(LALCHARCreateVector(status,&breveName,maxFilenameLength),status);
  sprintf(breveName->data,"%s.breve",filename);
  LALappsWriteBreveResults(status,
			   breveName->data,
			   outCurve);
  if (breveName)
    LAL_CALL(LALCHARDestroyVector(status,&breveName),status);

  /* Output Total file */
  LAL_CALL(LALCHARCreateVector(status,&totalName,maxFilenameLength),status);
  sprintf(totalName->data,"%s.full",filename);
  LALappsWriteSearchResults(status,
			    totalName->data,
			    outCurve);
  if (totalName)
    LAL_CALL(LALCHARDestroyVector(status,&totalName),status);
  /* If possible output configuration information */
  if (params!=NULL)
    {
      LAL_CALL(LALCHARCreateVector(status,&configName,maxFilenameLength),
	       status);
      sprintf(configName->data,"%s.config",filename);
      LALappsWriteSearchConfig(status->statusPtr,
			       configName->data,
			       *params);
    }     
  if (configName)
    LAL_CALL(LALCHARDestroyVector(status,&configName),status);

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
  fprintf(totalFile,"# Total Curves: %i\n",outCurve.numberOfCurves);
  fprintf(totalFile,"# Legend: Col,Row;gpsSec,gpsNanoSec,Freq,depth\n");
  /*Form of solution FreqIndex,TimeIndex,GPSSec,GPSNano,Power*/
  for (i = 0;i < outCurve.numberOfCurves;i++)
    {
      fprintf(totalFile,"Curve number,length,power:%i,%i,%f\n",
	      i,
	      outCurve.curves[i].n,
	      outCurve.curves[i].totalPower);
    for (j = 0;j < outCurve.curves[i].n;j++)
      { /*Long info*/
	fprintf(totalFile,"%i,%i;%i,%i,%f,%f",
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
      LAL_CALL(
	       LALGPStoFloat(status,
			     &startStamp,
			     &(outCurve.curves[i].gpsStamp[0]))
	       ,status);
      LAL_CALL(
	       LALGPStoFloat(status,
			     &stopStamp,
			     &(outCurve.curves[i].gpsStamp[outCurve.curves[i].n -1]))
	       ,status);
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
    if (c == '\n')
      lineCount++;
  /*printf("Lines found in input wave file: %i\n",lineCount);*/
  rewind(fp);
  /*
   * Allocate RAM to load up values
   */
  LAL_CALL(LALCreateVector(status,&domain,lineCount),status);
  LAL_CALL(LALCreateVector(status,&range,lineCount),status);
  /*
   * Expecting 2C data tstep and amp
   */
  for (i=0;i<lineCount;i++)
    {
      fscanf(fp,"%f %f\n",&domain->data[i],&range->data[i]);
      range->data[i]=(range->data[i]*params.scaleFactor);
    }
  fileDuration=domain->data[domain->length-1]-domain->data[0];
  newLineCount=params.sampleRate*fileDuration;
  offSetPoints=params.sampleRate*params.startTimeOffset;
  timePoints=params.sampleRate*params.injectSpace;
  LAL_CALL(LALCreateVector(status,&waveDomain,newLineCount),status);
  LAL_CALL(LALCreateVector(status,&waveRange,newLineCount),status);

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
  LAL_CALL(
	   LALCreateREAL4TimeSeries(status,
				    injectSet,
				    "Injections",
				    gps_zero,
				    0,
				    1/params.sampleRate,
				    lalADCCountUnit,
				    pointLength),
	   status);
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
    LAL_CALL(LALDestroyVector(status,&domain),status);
  if (range)
    LAL_CALL(LALDestroyVector(status,&range),status);
  if (waveRange)
    LAL_CALL(LALDestroyVector(status,&waveRange),status);
  if (waveDomain)
    LAL_CALL(LALDestroyVector(status,&waveDomain),status);
return;
}
/*
 * End of function LALappsCreateInjectableData to create injectable
 * data structure for calibration etc
 */

/*
 * End of Semi Private functions
 */


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
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
 * End local scratch functions
 * These functions may eventually be integrated more into
 * lalapps tracksearch code
 */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/




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
