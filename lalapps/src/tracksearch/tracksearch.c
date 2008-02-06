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
 * Author: Torres Cristina (Univ of TX at Brownsville)
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
  set_debug_level("NONE");
  memset(&status, 0, sizeof(status));
  /*lal_errhandler = LAL_ERR_EXIT;*/
  lal_errhandler = LAL_ERR_ABRT;
  /*lal_errhandler = LAL_ERR_RTRN;*/
  set_debug_level("MEMDBG");
  set_debug_level("ERROR");
  set_debug_level("3");
  /*  set_debug_level("ERROR | WARNING | MEMDBG");*/
  /*  set_debug_level("ERROR | WARNING | MEMDBG");*/
  /*  set_debug_level("ERROR | WARNING ");*/
  /*  set_debug_level("ALLDBG"); */

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
 * TINA -- Setup the removal of COHERENT LINES here!
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
  UINT4                     k=0;
  UINT4                     segmentPoints=0;
  UINT4                    planLength=0;
  PassBandParamStruc       bandPassParams;
  const LALUnit     strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;
  INT4Vector              *harmonicIndex=NULL;
  INT4Vector              *harmonicIndexCompliment=NULL;
  COMPLEX8Vector          *referenceSignal=NULL;
  COMPLEX8Vector          *tmpReferenceSignal=NULL;
  REAL4TVectorCLR         *inputTVectorCLR=NULL;
  REAL4FVectorCLR         *inputFVectorCLR=NULL;
  REAL4Vector             *inputPSDVector=NULL;
  REAL4Vector             *cleanData=NULL;
  UINT4                    tmpHarmonicCount=0;
  REAL4                    tmpLineFreq=0.0;
  
  /* 31May07
   * Note to self: Tina please double check that data in units
   * structure is properly migrated though function call...
   */
  if (params.verbosity == printFiles)
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
	  print_lalUnit(response->sampleUnits,"extractedResponseFunction_Units.diag");
	  print_complex8fseries(response,"extractedResponseFunction.diag");
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
	  print_complex8fseries(dataSetFFT,"Pre_CalibrateDataSetFFT.diag");
	  print_lalUnit(dataSetFFT->sampleUnits,"Pre_CalibrateDataSetFFT_Units.diag");
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
	  print_real4tseries(dataSet,"Post_CalibrateTimeDomainDataSet.diag");
	  print_lalUnit(dataSet->sampleUnits,"Post_CalibrateTimeDomainDataSe__Units.diag");
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

  /*****************************************
   * Perform low and or high  pass filtering on 
   * the input data stream if requested
   */
  if ((params.verbosity >= printFiles) && ((params.highPass > 0)||(params.lowPass)))
    print_real4tseries(dataSet,"Pre_ButterworthFiltered_AllDataset.diag");

  if (params.lowPass > 0)
    {
      if (params.verbosity >= verbose)
	fprintf(stdout,"You requested a low pass filter of the data at %f Hz\n",params.lowPass);

      bandPassParams.name=NULL;
      bandPassParams.nMax=10;
      /* F < f1 kept, F > f2 kept */
      bandPassParams.f1=params.lowPass;
      bandPassParams.f2=0;
      bandPassParams.a1=0.9;
      bandPassParams.a2=0;
      /*
       * Band pass is achieved by low pass first then high pass!
       * Call the low pass filter function.
       */
      LAL_CALL(LALButterworthREAL4TimeSeries(status,
					     dataSet,
					     &bandPassParams),
	       status);
    }
  /*
   * Call the high pass filter function.
   */
  if (params.highPass > 0)
    {
      if (params.verbosity >= verbose)
	fprintf(stdout,"You requested a high pass filter of the data at %f Hz\n",params.highPass);

      bandPassParams.name=NULL;
      bandPassParams.nMax=10;
      /* F < f1 kept, F > f2 kept */
      bandPassParams.f1=0;
      bandPassParams.f2=params.highPass;
      bandPassParams.a1=0;
      bandPassParams.a2=0.9;
      LAL_CALL(LALButterworthREAL4TimeSeries(status,
					     dataSet,
					     &bandPassParams),
	       status);
    }
  if ((params.verbosity >= printFiles) && ((params.highPass > 0)||(params.lowPass)))
    print_real4tseries(dataSet,"Post_ButterworthFiltered_AllDataset.diag");
  /*
   * End the Butterworth filtering
   *****************************************
   */
  if (params.numLinesToRemove > 0)
    {
      /*
       * Perform the removal of COHERENT LINES specified by user to entire
       * data set.  We want to remove up to n harmonics which have been
       * specified at the command line.  
       */
      if (params.verbosity > quiet)
	{
	  fprintf(stdout,"Removing requested coherent lines.\n");
	  for (i=0;i<params.numLinesToRemove;i++)
	    fprintf(stdout," %f, ",params.listLinesToRemove[i]);
	  fprintf(stdout,"\n");
	  fprintf(stdout,"Up to %i harmonics for %i lines will be removed.\n",params.maxHarmonics,params.numLinesToRemove);
	}
      /* Loop over all lines */
      /* Setup tmp reference signal */
      /* Add tmp reference to global Reference */
      /* Use global reference signal to clean input data */
      /* Take input data to Fourier Domain */
  
      /* Setup both T and F domain input data structures */
      inputTVectorCLR=(REAL4TVectorCLR*)LALMalloc(sizeof(REAL4TVectorCLR));
      inputFVectorCLR=(REAL4FVectorCLR*)LALMalloc(sizeof(REAL4FVectorCLR));
  
      /* Take data to F domain */
      planLength=dataSet->data->length;
      LAL_CALL(LALCCreateVector(status,
				&tmpReferenceSignal,
				planLength),
	       status);

      LAL_CALL(LALCCreateVector(status,
				&referenceSignal,
				planLength),
	       status);

      LAL_CALL(
	       LALCreateCOMPLEX8FrequencySeries(
						status,
						&signalFFT,
						"tmpSegPSD",
						gps_zero,
						0,
						1/(dataSet->data->length*dataSet->deltaT),
						dataSet->sampleUnits,
						planLength/2+1),
	       status);


      LAL_CALL(  LALSCreateVector(status,
				  &inputPSDVector,
				  planLength/2+1),
		 status);

      LAL_CALL(  LALCreateForwardREAL4FFTPlan(status,
					      &forwardPlan,
					      planLength,
					      0),
		 status);

  
      LAL_CALL( LALForwardREAL4FFT(status,
				   signalFFT->data,
				   dataSet->data,
				   forwardPlan),
		status);


      LAL_CALL( LALRealPowerSpectrum(status,
				     inputPSDVector,
				     dataSet->data,
				     forwardPlan),
		status);
      /*Assign CLR data structures*/
      /* The  CLR Time Vector  */
      inputTVectorCLR->length = dataSet->data->length;
      inputTVectorCLR->data = dataSet->data->data;
      inputTVectorCLR->deltaT = dataSet->deltaT; /* inverse of the sampling frequency */
      inputTVectorCLR->fLine = 0.0;
      /* The  CLR Frequency Vector  */
      inputFVectorCLR->length = inputPSDVector->length;
      inputFVectorCLR->data = inputPSDVector->data; 
      inputFVectorCLR->deltaF = signalFFT->deltaF;
      inputFVectorCLR->fLine =  inputTVectorCLR->fLine;

      for (i=0;i<referenceSignal->length;i++)
	{
	  referenceSignal->data[i].re=0;
	  referenceSignal->data[i].im=0;
	  tmpReferenceSignal->data[i].re=0;
	  tmpReferenceSignal->data[i].im=0;
	}
      for (j=0;j<params.numLinesToRemove;j++)
	{
	  if (params.verbosity > quiet)
	    {
	      fprintf(stdout,"Creating reference signal :%f\n",params.listLinesToRemove[j]);
	    }
	  /* Create reference signal */
	  tmpHarmonicCount=params.maxHarmonics;
	  tmpLineFreq=params.listLinesToRemove[j];
	  if ((params.SamplingRate/tmpLineFreq) < tmpHarmonicCount)
	    tmpHarmonicCount=floor(params.SamplingRate/tmpLineFreq);
      
	  LAL_CALL(LALI4CreateVector(status,
				     &harmonicIndex,
				     tmpHarmonicCount),
		   status);
	  LAL_CALL(LALI4CreateVector(status,
				     &harmonicIndexCompliment,
				     3*tmpHarmonicCount),
		   status);
	  for (i=0;i<tmpHarmonicCount;i++)
	    harmonicIndex->data[i]=i+1;

	  inputTVectorCLR->fLine = tmpLineFreq;
	  inputFVectorCLR->fLine = tmpLineFreq;

	  LAL_CALL( LALHarmonicFinder(status,
				      harmonicIndexCompliment,
				      inputFVectorCLR,
				      harmonicIndex),
		    status);

	  LAL_CALL( LALRefInterference(status,
				       tmpReferenceSignal,
				       signalFFT->data,
				       harmonicIndexCompliment),
		    status);



	  if (harmonicIndexCompliment)
	    {
	      LAL_CALL(LALI4DestroyVector(status,
					  &harmonicIndexCompliment),
		       status);
	    }

	  if (harmonicIndex)
	    {
	      LAL_CALL(LALI4DestroyVector(status,
					  &harmonicIndex),
		       status);
	    }
	  /* Copy this temp reference into global reference signal variable??*/
	  for (k=0;k < referenceSignal->length;k++)
	    {
	      referenceSignal->data[k].re=
		referenceSignal->data[k].re +
		tmpReferenceSignal->data[k].re;
	      referenceSignal->data[k].im=
		referenceSignal->data[k].im +
		tmpReferenceSignal->data[k].im;
	    }
	}
      /* Clean up input time series with global reference signal */
      LAL_CALL(  LALSCreateVector(status,
				  &cleanData,
				  dataSet->data->length),
		 status);
  
      LAL_CALL( LALCleanAll(status,
			    cleanData,
			    referenceSignal,
			    inputTVectorCLR),
		status);

      /*Copy the clean data back into the variable dataSet */
      for (j=0;j<dataSet->data->length;j++)
	dataSet->data->data[j]=cleanData->data[j];
  

      if (cleanData)
	{
	  LAL_CALL(
		   LALDestroyVector(status,
				    &cleanData),
		   status);
	}

      if (signalFFT)
	{
	  LAL_CALL( LALDestroyCOMPLEX8FrequencySeries(status,
						      signalFFT),
		    status);
	}

      if (inputPSDVector)
	{
	  LAL_CALL(
		   LALDestroyVector(status,
				    &inputPSDVector),
		   status);
	}
      if (forwardPlan)
	{
	  LAL_CALL(
		   LALDestroyREAL4FFTPlan(status,
					  &forwardPlan),
		   status);
	}
      if (tmpReferenceSignal)
	{
	  LAL_CALL(LALCDestroyVector(status,
				     &tmpReferenceSignal),
		   status);
	}

      if (referenceSignal)
	{
	  LAL_CALL(LALCDestroyVector(status,
				     &referenceSignal),
		   status);
	}

      if (inputTVectorCLR)
	LALFree(inputTVectorCLR);

      if (inputFVectorCLR)
	LALFree(inputFVectorCLR);

      /*
       * END COHERENT LINE removal section
       */
    }
  /*
   * Perform injections when inject structure has data.
   */
  if (injectSet != NULL)
    {
      if (params.verbosity > quiet)
	{
	  printf("Making requested injections.\n");
	  fprintf(stderr,"If frame was REAL8 data converted to REAL4 data.\n");
	  fprintf(stderr,"Potential problem with possibly prefactored input data.\n");
	}
      if (params.verbosity >= printFiles)
	print_real4tseries(dataSet,"Pre_SoftwareInjectDataSet.diag");
      /*
       * NOTE TO SELF: Tina add code to compare any factor in the dataSet
       * units structure and factor the injected data appropriately to
       * match the input data (dataSet).  This is bug affected by the
       * REAL8->REAL4 conversions in frame reading code only.
       */
      for (i=0,j=0;
	   ((i<injectSet->data->length)&&((j+params.SegBufferPoints)<dataSet->data->length));
	   i++,j++)
	if (j<params.SegBufferPoints)
	  dataSet->data->data[j]=dataSet->data->data[i];
	else
	  dataSet->data->data[j]=
	    injectSet->data->data[i-params.SegBufferPoints]+
	    dataSet->data->data[i];

      if (params.verbosity >= printFiles)
	print_real4tseries(dataSet,"Post_SoftwareInjectDataSet.diag");
    }
  /*
   * Split incoming data into Segment Vector
   * Adjust the segmenter to make orignal segments with
   * SegBufferPoints on each side of the segments.
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
					     1/((params.SegLengthPoints+(2*params.SegBufferPoints))*dataSet->deltaT),
					     lalDimensionlessUnit,
					     (params.SegLengthPoints+(2*params.SegBufferPoints))/2+1),
	       status);
      /*
       * The given units above need to be correct to truly reflect the
       * units stored in the frame file
       */
      windowParamsPSD.length=(params.SegLengthPoints+(2*params.SegBufferPoints));
      windowParamsPSD.type=params.avgSpecWindow;
      LAL_CALL( LALCreateREAL4Window(status,
				     &(windowPSD),
				     &windowParamsPSD),
		status);
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
      LAL_CALL( LALCreateForwardREAL4FFTPlan(status,
					     &averagePSDPlan,
					     params.SegLengthPoints+(2*params.SegBufferPoints),
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
	      print_real4fseries(averagePSD,"Pre_SmoothingAveragePSD.diag");
	      print_lalUnit(averagePSD->sampleUnits,
			    "Pre_SmoothingAveragePSD_Units.diag");
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
	  print_real4tseries(dataSet,"timeDomainAllDataSet.diag");
	  print_real4fseries(averagePSD,"Post_SmoothingAveragePSD.diag");
	  print_lalUnit(averagePSD->sampleUnits,"Post_SmoothingAveragePSD_Units.diag");
	}
      /*
       * Setup global FFT plans 
       */
      planLength=params.SegLengthPoints+(2*params.SegBufferPoints);
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
      if (params.verbosity >= printFiles)
	{
	  tmpSignalPtr=NULL;
	  for (i=0;i<dataSegments->length;i++)
	    {
	      tmpSignalPtr=(dataSegments->dataSeg[i]);
	      LAL_CALL(LALCHARCreateVector(status,&dataLabel,128),
		       status);
	      sprintf(dataLabel->data,"Pre_WhitenTimeDomainDataSeg_%i.diag",i);
	      print_real4tseries(tmpSignalPtr,dataLabel->data);
	      sprintf(dataLabel->data,"Pre_WhitenTimeDomainDataSeg_%i_Units.diag",i);
	      print_lalUnit(tmpSignalPtr->sampleUnits,dataLabel->data);
	      LAL_CALL(LALCHARDestroyVector(status,&dataLabel),
		       status);
	    }
	  tmpSignalPtr=NULL;

	  for (i=0;i<dataSegments->length;i++)
	    {
	      tmpSignalPtr = (dataSegments->dataSeg[i]);
	      /*
	       * FFT segment
	       */
	      LAL_CALL( LALForwardREAL4FFT(status,
					   signalFFT->data,
					   tmpSignalPtr->data,
					   forwardPlan),
			status);
	      /*
	       * Whiten
	       */
	      if (params.verbosity >= printFiles)
		{
		  LAL_CALL(LALCHARCreateVector(status,&dataLabel,128),
			   status);
		  sprintf(dataLabel->data,"Pre_whitenSignalFFT_%i.diag",i);
		  print_complex8fseries(signalFFT,dataLabel->data);
		  LAL_CALL(LALCHARDestroyVector(status,&dataLabel),
			   status);
		}
	      if (params.verbosity >= verbose)
		fprintf(stdout,"Whitening data segment: %i of %i\n",i,dataSegments->length);
	      LAL_CALL(LALTrackSearchWhitenCOMPLEX8FrequencySeries(status,
								   signalFFT,
								   averagePSD,
								   params.whiten),
		       status);
	      if (params.verbosity >= printFiles)
		{
		  LAL_CALL(LALCHARCreateVector(status,&dataLabel,128),
			   status);
		  sprintf(dataLabel->data,"Post_whitenSignalFFT_%i.diag",i);
		  print_complex8fseries(signalFFT,dataLabel->data);
		  LAL_CALL(LALCHARDestroyVector(status,&dataLabel),
			   status);
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
      if (params.verbosity >= printFiles)
	{
	  tmpSignalPtr=NULL;
	  for (i=0;i<dataSegments->length;i++)
	    {
	      tmpSignalPtr=(dataSegments->dataSeg[i]);
	      LAL_CALL(LALCHARCreateVector(status,&dataLabel,128),
		       status);
	      sprintf(dataLabel->data,"Post_WhitenTimeDomainDataSeg_%i.diag",i);
	      print_real4tseries(tmpSignalPtr,dataLabel->data);
	      sprintf(dataLabel->data,"Post_WhitenTimeDomainDataSeg_%i_Units.diag",i);
	      print_lalUnit(tmpSignalPtr->sampleUnits,dataLabel->data);
	      LAL_CALL(LALCHARDestroyVector(status,&dataLabel),
		       status);
	    }
	}
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
	    params->listLinesToRemove[params->numLinesToRemove]=
	      atof(optarg);
	    params->numLinesToRemove=params->numLinesToRemove+1;
	  }
	  break;

	case 'T':
	  {
	    params->maxHarmonics=atoi(optarg);
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
	  params->discardTLP=(params->TimeLengthPoints-params->overlapFlag)%(params->SegLengthPoints - params->overlapFlag);
	  /* 
	   * Reset points to process by N-discardTLP, this gives us
	   * uniform segments to make maps with
	   */
	  params->TimeLengthPoints=params->TimeLengthPoints-params->discardTLP;
	  if ((params->verbosity >= verbose))
	    fprintf(stdout,"We need to throw away %i trailing data points for %i, evenly matched data segments\n",params->discardTLP,params->NumSeg);
	  
	};
      /* 
       * Determine the number of additional points required to negate
       * edge effects in first and last TFR
       */
      params->SegBufferPoints=((params->TimeLengthPoints/params->TimeBins)*params->colsToClip);
      params->SegBufferPoints=((params->SegLengthPoints/params->TimeBins)*params->colsToClip);
      if ((params->verbosity >= verbose))
	fprintf(stdout,"As a result of bin buffering requests an additional %i points are needed to negate TFR edge effects, %i points before and then after the data set to analyze.\n",2*params->SegBufferPoints,params->SegBufferPoints);

      if (params->whiten > 2 )
	{
	  fprintf(stderr,TRACKSEARCHC_MSGEARGS);
	  exit(TRACKSEARCHC_EARGS);
	}
    }
  /*
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

  if (params->MinLength < 2)
    {
      fprintf(stderr,"Minimum length threshold invalid!\n");
      fprintf(stderr,TRACKSEARCHC_MSGEARGS);
      exit(TRACKSEARCHC_EARGS);
    }
  if (
      ((params->LineWidth) < ceil((2.0/MASK_SIZE)*(2.0*sqrt(3))))
    ||
      (params->LineWidth/(2.0*sqrt(3)) < 1)
      )
    {
      optWidth=0;
      if ((2.0*sqrt(3)) < ceil((2.0/MASK_SIZE)))
	{
	  optWidth=ceil((2.0/MASK_SIZE)*(2.0*sqrt(3)));
	}
      else
	{
	  optWidth=ceil(2.0*sqrt(3));
	}
      fprintf(stderr,"Line width inappropriate try mimimum of: %i\n",
	      (INT4) optWidth);
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
  REAL8TimeSeries      *convertibleREAL8Data=NULL;
  REAL8TimeSeries      *tmpREAL8Data=NULL;
  INT2TimeSeries       *convertibleINT2Data=NULL;
  INT2TimeSeries       *tmpINT2Data=NULL;
  INT4TimeSeries       *convertibleINT4Data=NULL;
  INT4TimeSeries       *tmpINT4Data=NULL;
  PassBandParamStruc    bandPassParams;
  UINT4                 loadPoints=0;
  UINT4                 i=0;
  LALTYPECODE           dataTypeCode=0;
  ResampleTSParams      resampleParams;
  INT4                  errCode=0;
  LIGOTimeGPS           bufferedDataStartGPS;
  REAL8                 bufferedDataStart=0;
  LIGOTimeGPS           bufferedDataStopGPS;
  REAL8                 bufferedDataStop=0;
  REAL8                 bufferedDataTimeInterval=0;

  /* Set all variables from params structure here */
  channelIn.name = params->channelName;
  channelIn.type = params->channelNameType;
  /* only try to load frame if name is specified */
  if (dirname || cachefile)
    {
      if(dirname)
	{
	  /* Open frame stream */
	  LAL_CALL( LALFrOpen( status, &stream, dirname->data, "*.gwf" ), status);
	}
      else if (cachefile)
	{
	  /* Open frame cache */
	  lal_errhandler = LAL_ERR_EXIT;
	  LAL_CALL( LALFrCacheImport( status, &frameCache, cachefile ), status);
	  LAL_CALL( LALFrCacheOpen( status, &stream, frameCache ), status);
	  LAL_CALL( LALDestroyFrCache( status, &frameCache ), status );
	}
      lal_errhandler = LAL_ERR_EXIT;
      /* Set verbosity of stream so user sees frame read problems! */
      LAL_CALL( LALFrSetMode(status,LAL_FR_VERBOSE_MODE,stream),status);
      /*DataIn->epoch SHOULD and MUST equal params->startGPS - (params.SegBufferPoints/params.SamplingRate)*/
      memcpy(&bufferedDataStartGPS,&(DataIn->epoch),sizeof(LIGOTimeGPS));
      LAL_CALL(
	       LALGPStoFloat(status,
			     &bufferedDataStart,
			     &bufferedDataStartGPS),
	       status);
      /* 
       * Seek to end of requested data makes sure that all stream is complete!
       */
      bufferedDataStop=bufferedDataStart+(DataIn->data->length * DataIn->deltaT);
      LAL_CALL(
	       LALFloatToGPS(status,
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

      LAL_CALL( LALFrSeek(status, &(bufferedDataStopGPS),stream),status);

      LAL_CALL( LALFrSeek(status, &(bufferedDataStartGPS), stream), status);
      /*
       * Determine the variable type of data in the frame file.
       */
      LAL_CALL( LALFrGetTimeSeriesType(status, &dataTypeCode, &channelIn, stream),status);
      

      if (dataTypeCode == LAL_S_TYPE_CODE)
	{
	  /* Proceed as usual reading in REAL4 data type information */
	  /* Load the metadata to check the frame sampling rate */
	  if (params->verbosity >= verbose)
	    fprintf(stderr,"NO conversion of frame REAL4 data needed.\n");

	  lal_errhandler = LAL_ERR_EXIT;
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
	  /*
	   * Add the bin_buffer points to the loadpoints variable.
	   * This will load the correct number of points.
	   */
	  loadPoints=loadPoints+(params->SamplingRateOriginal*2*(params->SegBufferPoints/params->SamplingRate));

	  LAL_CALL(
		   LALCreateREAL4TimeSeries(status,
					    &tmpData,
					    "Higher Sampling Tmp Data",
					    bufferedDataStartGPS,
					    0,
					    1/params->SamplingRateOriginal,
					    lalADCCountUnit,
					    loadPoints),
		   status);
	  /* get the data */
	  LAL_CALL( LALFrSeek(status, &(tmpData->epoch), stream), status);
	  lal_errhandler = LAL_ERR_RTRN;
	  errCode=LAL_CALL( LALFrGetREAL4TimeSeries( status, 
						     tmpData, 
						     &channelIn, 
						     stream), 
			    status);
	  if (errCode != 0)
	    {
	      fprintf(stderr,"The span of REAL4TimeSeries data is not available in frame stream!\n");
	      exit(errCode);
	    }
	} /*End of reading type REAL4TimeSeries*/
      else if (dataTypeCode == LAL_D_TYPE_CODE)
	{
	  /*
	   * Proceed with a REAL8 Version of above and 
	   * copy to REAL4 Struct in the end
	   */
	  /* Create temporary space use input REAL4TimeSeries traits!*/
	  if (params->verbosity >= verbose)
	    fprintf(stderr,"Must convert frame REAL8 data to REAL4 data for analysis!\n");
	  LAL_CALL(
		   LALCreateREAL8TimeSeries(status,
					    &tmpREAL8Data,
					    "tmp space",
					    bufferedDataStartGPS,
					    0,
					    DataIn->deltaT,
					    lalADCCountUnit,
					    DataIn->data->length),
		   status);

	  /* Load the metadata to check the frame sampling rate */
	  lal_errhandler = LAL_ERR_EXIT;
	  LAL_CALL( LALFrGetREAL8TimeSeriesMetadata( status, 
						     tmpREAL8Data, 
						     &channelIn, 
						     stream), 
		    status);
	  /*
	   * Determine the sampling rate and assuming the params.sampling rate
	   * how many of these points do we load to get the same time duration
	   * of data points
	   */
 	  params->SamplingRateOriginal=1/(tmpREAL8Data->deltaT);
 	  loadPoints=params->SamplingRateOriginal*(params->TimeLengthPoints/params->SamplingRate);
	  loadPoints=loadPoints+(params->SamplingRateOriginal*2*(params->SegBufferPoints/params->SamplingRate));
	  LAL_CALL(
		   LALCreateREAL8TimeSeries(status,
					    &convertibleREAL8Data,
					    "Higher Sampling Tmp Data",
					    bufferedDataStartGPS,
					    0,
					    1/params->SamplingRateOriginal,
					    lalADCCountUnit,
					    loadPoints),
		   status);
	  /* get the data */
	  LAL_CALL( LALFrSeek(status, &(tmpREAL8Data->epoch), stream), status);
	  lal_errhandler = LAL_ERR_RTRN;
	  errCode=LAL_CALL( LALFrGetREAL8TimeSeries( status, 
						     convertibleREAL8Data, 
						     &channelIn, 
						     stream), 
			    status);
	  if (errCode != 0)
	    {
	      fprintf(stderr,"The span of REAL8TimeSeries data is not available in framestream!\n");
	      exit(errCode);
	    }
	  /*
	   * REAL8-->REAL4 Workaround (tmp solution factor into units struct!)
	   * eventually we will rework this technique to use REAL8
	   * effectively!!!! In order to recover higher freq
	   * componenet in time series we will high pass the REAL8 data
	   * first before factoring it.
	   */
	  /******************************************/
	  /*
	   * Perform low and or high  pass filtering on 
	   * the input data stream if requested
	   */
	  if (params->verbosity >= verbose)
	    fprintf(stdout,"Frame Reader needs to high pass input for casting!\n");

	  if (params->lowPass > 0)
	    {
	      if (params->verbosity >= verbose)
		fprintf(stdout,"FRAME READER: You requested a low pass filter of the data at %f Hz\n",params->lowPass);
	      bandPassParams.name=NULL;
	      bandPassParams.nMax=10;
	      /* F < f1 kept, F > f2 kept */
	      bandPassParams.f1=params->lowPass;
	      bandPassParams.f2=0;
	      bandPassParams.a1=0.9;
	      bandPassParams.a2=0;
	      /*
	       * Band pass is achieved by low pass first then high pass!
	       * Call the low pass filter function.
	       */
	      LAL_CALL(LALButterworthREAL8TimeSeries(status,
						     convertibleREAL8Data, 
						     &bandPassParams),
		       status);
	    }
	  if (params->highPass > 0)
	    {
	      if (params->verbosity >= verbose)
		fprintf(stdout,"FRAME READER: You requested a high pass filter of the data at %f Hz\n",params->highPass);
	      bandPassParams.name=NULL;
	      bandPassParams.nMax=10;
	      /* F < f1 kept, F > f2 kept */
	      bandPassParams.f1=0;
	      bandPassParams.f2=params->highPass;
	      bandPassParams.a1=0;
	      bandPassParams.a2=0.9;
	      LAL_CALL(LALButterworthREAL8TimeSeries(status,
						     convertibleREAL8Data, 
						     &bandPassParams),
		       status);
	      /*
	       * End Butterworth filtering
	       */
	    }
	  if (params->highPass < 30)
	    {
	      fprintf(stdout,"WARNING! Input data should be high pass filtered!\n");
	      fprintf(stdout,"Use knee frequency of 30Hz minimum!\n");
	      fprintf(stderr,"PSD estimate and results are questionable!\n");
	    };
	  if ((params->verbosity >= printFiles) && ((params->highPass > 0)||(params->lowPass)))
	    {
	      print_real8tseries(convertibleREAL8Data,"FRAME_READER_1_ButterworthFiltered.diag");
	      print_lalUnit(convertibleREAL8Data->sampleUnits,"FRAME_READER_1_ButterworthFiltered_Units.diag");
	    };
	  /*
	   * End the high pass filter
	   *****************************************
	   */
	  /*
	   * Factor the data by 10^20! See Xavi email.
	   */
	  if (params->verbosity >= printFiles)
	    {
	      print_real8tseries(convertibleREAL8Data,"FRAME_READER_2_PreFactoredREAL8.diag");
	      print_lalUnit(convertibleREAL8Data->sampleUnits,"FRAME_READER_2_PreFactoredREAL8_Units.diag");
	    }
	  if (params->verbosity >= verbose)
	    {	    
	      fprintf(stderr,"We are factoring the input data for casting.\n");
	      fprintf(stderr,"Factoring out : 10e20\n");
	      fprintf(stderr,"Factor place in data units structure.\n");
	    }
	  /* Factoring out 10^-20 and to unitsstructure! */
	  for (i=0;i<convertibleREAL8Data->data->length;i++)
	    convertibleREAL8Data->data->data[i]=
	      convertibleREAL8Data->data->data[i]*pow(10,20);
	  convertibleREAL8Data->sampleUnits.powerOfTen=-20;
	  if (params->verbosity >= printFiles)
	    {
	      print_real8tseries(convertibleREAL8Data,"FRAME_READER_3_FactoredREAL8.diag");
	      print_lalUnit(convertibleREAL8Data->sampleUnits,"FRAME_READER_3_FactoredREAL8_Units.diag");
	      LALappsPSD_Check(convertibleREAL8Data);

	    }

	  /* Prepare to copy/cast REAL8 data into REAL4 structure */
	  /* Copy REAL8 data into new REAL4 Structure */
	  LALappsCreateR4FromR8TimeSeries(status,&tmpData,convertibleREAL8Data);
	  if (tmpREAL8Data)
	    LAL_CALL(LALDestroyREAL8TimeSeries(status,tmpREAL8Data),
		     status); 
	  if (convertibleREAL8Data)
	    LAL_CALL(LALDestroyREAL8TimeSeries(status,convertibleREAL8Data),
		     status); 
	}/*End of reading type REAL8TimeSeries*/
      else if  (dataTypeCode == LAL_I2_TYPE_CODE)
	{
	  /*Proceed with a INT2 Version of above and copy to REAL4 Struct in the end*/
	  /* Create temporary space use input REAL4TimeSeries traits!*/
	  if (params->verbosity >= verbose)
	    fprintf(stderr,"Must convert frame INT2 data to REAL4 data for analysis!\n");
	  LAL_CALL( 
		   LALCreateINT2TimeSeries(status,
					   &tmpINT2Data,
					   "tmp space",
					   bufferedDataStartGPS,
					   0,
					   DataIn->deltaT,
					   lalADCCountUnit,
					   DataIn->data->length),
		   status);

	  /* Load the metadata to check the frame sampling rate */
	  lal_errhandler = LAL_ERR_EXIT;
	  LAL_CALL( LALFrGetINT2TimeSeriesMetadata( status, 
						    tmpINT2Data, 
						    &channelIn, 
						    stream), 
		    status);
	  /*
	   * Determine the sampling rate and assuming the params.sampling rate
	   * how many of these points do we load to get the same time duration
	   * of data points
	   */
	  params->SamplingRateOriginal=1/(tmpINT2Data->deltaT);
	  loadPoints=params->SamplingRateOriginal*(params->TimeLengthPoints/params->SamplingRate);
	  loadPoints=loadPoints+(params->SamplingRateOriginal*2*(params->SegBufferPoints/params->SamplingRate));
	  LAL_CALL(
		   LALCreateINT2TimeSeries(status,
					   &convertibleINT2Data,
					   "Tmp Data",
					   bufferedDataStartGPS,
					   0,
					   1/params->SamplingRateOriginal,
					   lalADCCountUnit,
					   loadPoints),
		   status);
	  /* get the data */
	  LAL_CALL( LALFrSeek(status, &(tmpINT2Data->epoch), stream), status);
	  lal_errhandler = LAL_ERR_RTRN;
	  errCode=LAL_CALL( LALFrGetINT2TimeSeries( status, 
						    convertibleINT2Data, 
						    &channelIn, 
						    stream), 
			    status);
	  if (errCode != 0)
	    {
	      fprintf(stderr,"The span of INT2TimeSeries data is not available in framestream!\n");
	      exit(errCode);
	    }
	  /* Create high sampled REAL4 variable */
	  LAL_CALL(
		   LALCreateREAL4TimeSeries(status,
					    &tmpData,
					    "Higher Sampling REAL4 Data",
					    bufferedDataStartGPS,
					    0,
					    1/params->SamplingRateOriginal,
					    lalADCCountUnit,
					    loadPoints),
		   status);
	  /* Prepare to copy/cast REAL8 data into REAL4 structure */
	  for (i=0;i<tmpData->data->length;i++)
	    tmpData->data->data[i] = (REAL4) convertibleINT2Data->data->data[i];
	  tmpData->sampleUnits=convertibleINT2Data->sampleUnits;
	  if (convertibleINT2Data)
	    LAL_CALL(LALDestroyINT2TimeSeries(status,convertibleINT2Data),
		     status); 
	  if (tmpINT2Data)
	    LAL_CALL(LALDestroyINT2TimeSeries(status,tmpINT2Data),
		     status); 
	}
      else if  (dataTypeCode == LAL_I4_TYPE_CODE)
	{
	  /*Proceed with a INT2 Version of above and copy to REAL4 Struct in the end*/
	  /* Create temporary space use input REAL4TimeSeries traits!*/
	  if (params->verbosity >= verbose)
	    fprintf(stderr,"Must convert frame INT4 data to REAL4 data for analysis!\n");
	  LAL_CALL( 
		   LALCreateINT4TimeSeries(status,
					   &tmpINT4Data,
					   "REAL8 tmp space",
					   bufferedDataStartGPS,
					   0,
					   DataIn->deltaT,
					   lalADCCountUnit,
					   DataIn->data->length),
		   status);

	  /* Load the metadata to check the frame sampling rate */
	  lal_errhandler = LAL_ERR_EXIT;
	  LAL_CALL( LALFrGetINT4TimeSeriesMetadata( status, 
						    tmpINT4Data, 
						    &channelIn, 
						    stream), 
		    status);
	  /*
	   * Determine the sampling rate and assuming the params.sampling rate
	   * how many of these points do we load to get the same time duration
	   * of data points
	   */
	  params->SamplingRateOriginal=1/(tmpINT4Data->deltaT);
	  loadPoints=params->SamplingRateOriginal*(params->TimeLengthPoints/params->SamplingRate);
	  loadPoints=loadPoints+(params->SamplingRateOriginal*2*(params->SegBufferPoints/params->SamplingRate));
	  LAL_CALL(
		   LALCreateINT4TimeSeries(status,
					   &convertibleINT4Data,
					   "Tmp Data",
					   bufferedDataStartGPS,
					   0,
					   1/params->SamplingRateOriginal,
					   lalADCCountUnit,
					   loadPoints),
		   status);
	  /* get the data */
	  LAL_CALL( LALFrSeek(status, &(tmpINT4Data->epoch), stream), status);
	  lal_errhandler = LAL_ERR_RTRN;
	  errCode=LAL_CALL( LALFrGetINT4TimeSeries( status, 
						    convertibleINT4Data, 
						    &channelIn, 
						    stream), 
			    status);
	  if (errCode != 0)
	    {
	      fprintf(stderr,"The span of INT2TimeSeries data is not available in framestream!\n");
	      exit(errCode);
	    }
	  /* Create high sampled REAL4 variable */
	  LAL_CALL(
		   LALCreateREAL4TimeSeries(status,
					    &tmpData,
					    "Higher Sampling REAL4 Data",
					    bufferedDataStartGPS,
					    0,
					    1/params->SamplingRateOriginal,
					    lalADCCountUnit,
					    loadPoints),
		   status);
	  /* Prepare to copy/cast REAL8 data into REAL4 structure */
	  for (i=0;i<tmpData->data->length;i++)
	    tmpData->data->data[i] = (REAL4) convertibleINT4Data->data->data[i];
	  tmpData->sampleUnits=convertibleINT4Data->sampleUnits;
	  if (convertibleINT4Data)
	    LAL_CALL(LALDestroyINT4TimeSeries(status,convertibleINT4Data),
		     status); 
	  if (tmpINT4Data)
	    LAL_CALL(LALDestroyINT4TimeSeries(status,tmpINT4Data),
		     status); 
	}
      else
	{
	  fprintf(stderr,"This type of frame data(%o) can not be used.\n",dataTypeCode);
	  fprintf(stderr,"Expecting INT2, INT4, REAL4 or REAL8\n");
	  fprintf(stderr,"Contact tracksearch authors for further help.\n");
	  LALappsTSassert(1==0,
			  TRACKSEARCHC_EDATA,
			  TRACKSEARCHC_MSGEDATA);
	}
      /* End of error for invalid data types in frame file.*/
      /* close the frame stream */
      LAL_CALL( LALFrClose( status, &stream ), status);

      if (params->verbosity >= printFiles)
	{
	  print_real4tseries(tmpData,"RawOriginalInputTimeSeries.diag");
	  print_lalUnit(tmpData->sampleUnits,"RawOriginalInputTimeSeries_Units.diag");
	}
      /*
       * Prepare for the resample if needed or just copy the data so send
       * back
       */
      if (params->SamplingRate < params->SamplingRateOriginal)
	{
	  if (params->verbosity >= verbose)
	    fprintf(stderr,"We will resample the data from %6.3f to %6.3f.\n",
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
	  if (params->verbosity >= printFiles)
	    {
	      print_real4tseries(tmpData,"ResampledOriginalInputTimeSeries.diag");
	      print_lalUnit(tmpData->sampleUnits,"ResampledlOriginalInputTimeSeries_Units.diag");
	    }
	  /*
	   * Copy only the valid data and fill the returnable metadata
	   */
	  DataIn->deltaT=(1/params->SamplingRate);
	  DataIn->sampleUnits=tmpData->sampleUnits;
	  /*
	   * Store new 'epoch' information into DataIn due to possible
	   * shift in epoch information due to resampling
	   */
	  memcpy(&(DataIn->epoch),&(tmpData->epoch),sizeof(LIGOTimeGPS));
	  DataIn->epoch=bufferedDataStartGPS;
	  LALappsTSassert((tmpData->data->length >=
			   DataIn->data->length),
			  TRACKSEARCHC_EDATA,
			  TRACKSEARCHC_MSGEDATA);
	  for (i=0;i<DataIn->data->length;i++)
	    DataIn->data->data[i]=tmpData->data->data[i];
	}
      else
	{
	  if (params->verbosity >= verbose)
	    fprintf(stdout,"Resampling to %f, not possible we have %f sampling rate in the frame file.\n",
		    params->SamplingRate,
		    params->SamplingRateOriginal);
	  /*
	   * Straight copy the data over to returnable structure
	   */
	  params->SamplingRate=params->SamplingRateOriginal;
	  DataIn->deltaT=1/params->SamplingRateOriginal;
	  DataIn->sampleUnits=tmpData->sampleUnits;
	  memcpy(&(DataIn->epoch),&(tmpData->epoch),sizeof(LIGOTimeGPS));
	  for (i=0;i<DataIn->data->length;i++)
	    DataIn->data->data[i]=tmpData->data->data[i];
	}
      if (params->verbosity >= printFiles)
	{
	  print_real4tseries(DataIn,"CopyForUse_ResampledTimeSeries.diag");
	  print_lalUnit(DataIn->sampleUnits,"CopyForUse_ResampledTimeSeries_Units.diag");
	}
      /*
       * Release the memory for the temporary time series
       */
      if (tmpData)
	LAL_CALL(LALDestroyREAL4TimeSeries(status,tmpData),
		 status);
    }
}
/* End frame reading code */

/* MISSING CLOSE } */


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
  /*
   * if the tssearchparams 'params' variable contains the directive
   * to apply auto determined lambda thresholds, calculate them for
   * the map and use them to run the analysis
   */
  /*
   * DO THE AUTO ADJUSTMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */
  if (params.autoLambda)
    {
      /* Do the calculate of Lh given parameters */
      /* Just pending inclussion of gaussian integral */
      LAL_CALL( LALTracksearchFindLambda(status,*tfmap,&params),
		status);
      tsInputs.high=params.StartThresh;
      tsInputs.low=params.LinePThresh;
    }
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
  if (params.verbosity >= printFiles)
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
  LAL_CALL( LALTrackSearchApplyThreshold(status,
					 &outputCurves,
					 &outputCurvesThreshold,
					 params),
	    status);

  /*
   * Record user request thresholds into output data structure
   * as a simple reference
   */
  outputCurvesThreshold.minPowerCut=params.MinPower;
  outputCurvesThreshold.minLengthCut=params.MinLength;
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
  TSAMap                *tmpTSA=NULL;/*point to correct vars for cropping*/
  REAL8                 signalStop=0;
  REAL8                 signalStart=0;
  REAL8                 cropDeltaT=0;
  INT4                  tmpSegDataPoints=0;
  INT4                  tmpMapTimeBins=0;
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

      LAL_CALL(LALCreateTimeFreqParam(status,&autoparams,&tfInputs),
	       status);

      LAL_CALL(LALCreateTimeFreqRep(status,&tfmap,&tfInputs),
	       status);
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
	    lal_errhandler = LAL_ERR_EXIT;
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
   * Determine the 'true' cropped map starting epoch.  This should in
   * pratice match the epoch the user originally intended without the
   * implicit time bin addition and subsequent cropping
   */
  /*Translate buffer points to GPS delta T */
  /*Add this time to the time of the buffered data segment */
  /*Report this as the true start time of cropped map */
  cropDeltaT=signalSeries->deltaT*params.SegBufferPoints;

  mapMarkerParams.mapStartGPS=signalSeries->epoch;
  LAL_CALL(
	   LALGPStoFloat(status,
			 &signalStart,
			 &signalSeries->epoch),
	   status);
  /*
   * Fix the signalStop time stamp to be without the buffer points.  It
   * should be the stop time of the clipped TFR.
   */
  signalStop=(signalSeries->deltaT*(signalSeries->data->length))+signalStart;
  LAL_CALL(
	   LALFloatToGPS(status,
			 &(mapMarkerParams.mapStopGPS),
			 &signalStop),
			 
	   status);
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
  tmpTSA=LALMalloc(sizeof(TSAMap));
  tmpTSA->imageCreateParams=tfInputs;
  tmpTSA->clippedWith=0;
  tmpTSA->imageBorders=mapMarkerParams;
  tmpTSA->imageRep=tfmap;
/*   /\* */
/*    * Dump out PGM of map before we crop it! */
/*    *\/ */
/*   if (params.verbosity >= verbose) */
/*     LALappsTSAWritePGM(status,tmpTSA,NULL);  */
  LALappsTSACropMap(status,&tmpTSA,params.colsToClip);
    
  /* 
   *Copy information from cropping procedure to relevant structures.
   */
  memcpy(&tfInputs,&(tmpTSA->imageCreateParams),sizeof(CreateTimeFreqIn));
  memcpy(&mapMarkerParams,&(tmpTSA->imageBorders),sizeof(TrackSearchMapMarkingParams));
  if (params.verbosity >= printFiles)
    LALappsTSAWritePGM(status,tmpTSA,NULL); 
  tfmap=tmpTSA->imageRep;
  LALFree(tmpTSA);

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
/*   if (params.verbosity >= verbose) */
/*     { */
/*       fprintf(stderr,"Map related search parameters\n"); */
/*       fprintf(stderr,"Sigma\t\t%f\n",inputs.sigma); */
/*       fprintf(stderr,"High\t\t%f\n",inputs.high); */
/*       fprintf(stderr,"Low\t\t%f\n",inputs.low); */
/*       fprintf(stderr,"Width\t\t%i\t\t%i\n",inputs.width,(tfmap->fRow/2)+1); */
/*       fprintf(stderr,"Height\t\t%i\t\t%i\n",inputs.height,tfmap->tCol); */
/*     } */

  /* If requested to the auto lambda determination */

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
  LIGOTimeGPS       edgeOffsetGPS;
  REAL8             originalFloatTime=0;
  REAL8             newFloatTime=0;
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
    LAL_CALL(
	     LALGPStoFloat(status,
			   &originalFloatTime,
			   &(params.GPSstart)),
	     status);
    newFloatTime=originalFloatTime-(params.SegBufferPoints/params.SamplingRate);
    LAL_CALL(
	     LALFloatToGPS(status,
			   &edgeOffsetGPS,
			   &newFloatTime),
	     status);
  LAL_CALL(
	   LALCreateREAL4TimeSeries(status,
				    &dataset,
				    "Uninitialized",
				    edgeOffsetGPS,
				    0,
				    1/params.SamplingRate,
				    lalADCCountUnit,
				    (params.TimeLengthPoints+(2*params.SegBufferPoints))),
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
      /* Wed-Oct-24-2007:200710241529 
       * Adjust starting point of injection to avoid
       * cropping first injection from TFR
       */
      if (params.verbosity >= verbose)
	{
	  printf("Preparing the injection data!\n");
	  printf("Inject Offset    %f\n",injectParams.startTimeOffset);
	  printf("Inject Space     %f\n",injectParams.injectSpace);
	  printf("Inject Count     %i\n",injectParams.numOfInjects);
	  printf("Inject Scale     %f\n",injectParams.scaleFactor);
	  printf("Inject Sampling  %f\n",injectParams.sampleRate);
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
  if (dataset)
    {
      LAL_CALL(LALDestroyREAL4TimeSeries(status,dataset),
	       status);
      dataset=NULL;
    }
  if (injectSet)
    {
      LAL_CALL(LALDestroyREAL4TimeSeries(status,injectSet),
	       status);
      injectSet=NULL;
    }

  j=0;
  for(i = 0;i < params.NumSeg;i++)
    {
      printf(".");
      /*
       * Call to prepare tSeries search
       * Rewrite functions inside to CROP maps and adjust time marks
       * according to the uncropped parts of TFR. Ignore segment
       * buffers....
       */
      LALappsDoTSeriesSearch(status,SegVec->dataSeg[j],params,j);
      j++;
    };
  printf("\n");
  /* Free some of the memory used to do the analysis */
  if (params.dataSegVec)
    {
      LALDestroyTSDataSegmentVector(status,(params.dataSegVec));
      params.dataSegVec=NULL;
    }
  if (SegVec)
    {
      LAL_CALL(LALDestroyTSDataSegmentVector(status,SegVec),
	       status);
      SegVec=NULL;
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
  fprintf(totalFile,"# Total Curves,Lh,Ll: %i,%f,%f\n",outCurve.numberOfCurves,outCurve.startThreshCut,outCurve.linePThreshCut);
  fprintf(totalFile,"# Legend: Col,Row;gpsSec,gpsNanoSec,Freq,depth\n");
  /*Form of solution FreqIndex,TimeIndex,GPSSec,GPSNano,Power*/
  for (i = 0;i < outCurve.numberOfCurves;i++)
    {
      fprintf(totalFile,"Curve number,length,power:%i,%i,%6.18f\n",
	      i,
	      outCurve.curves[i].n,
	      outCurve.curves[i].totalPower);
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
  LAL_CALL(LALCreateVector(status,&domain,lineCount),status);
  LAL_CALL(LALCreateVector(status,&range,lineCount),status);
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
  LALappsTSassert((pointLength>0),
		  TRACKSEARCHC_EARGS,
		  TRACKSEARCHC_MSGEARGS);
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
 
