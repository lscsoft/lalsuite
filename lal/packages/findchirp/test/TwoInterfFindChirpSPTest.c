/*----------------------------------------------------------------------- 
 * 
 * File Name: TwoInterfFindChirpSPTest.c
 *
 * Author: Bose, S., Noel, J. S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h>
#include <lal/TwoInterfFindChirp.h>

static REAL4 cartesianInnerProduct(REAL4 x[3], REAL4 y[3])
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

NRCSID(MAIN,"$Id$");

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

static void 
Usage (const char *program, int exitflag);

static void 
ParseOptions (int argc, char *argv[]);

void ReadResponse(LALStatus *status, COMPLEX8FrequencySeries *resp, char* respFile, BOOLEAN invert);

typedef enum{
  impulse, gaussian, file
    }
InputDataType;

extern char    *optarg;
extern int      optind;

int lalDebugLevel = 1;

static InputDataType    inputDataType        = gaussian; /* other options: impulse, file*/
static BOOLEAN          twoInterfRhosqout    = 0;
static BOOLEAN          rhosqout             = 0;
static BOOLEAN          verbose              = 0; 

static BOOLEAN          invertResponse       = 0;

static INT4             numPoints            = 32768;
static INT4             numSegments          = 1;
static INT4             numTmplts            = 1;
static INT4             numChisqBins         = 8;
static INT4             srate                = 8192;
static REAL4            sigmasq              = 64.0;
static REAL4            sigmasq1             = 64.0;
static REAL4            sigmasq2             = 64.0;
static INT4             invSpecTrunc         = 0;
static REAL4            fLow                 = 150.0;
static REAL4            rhosqThresh[2]       = {0.0, 0.0};
static REAL4            chisqThresh[2]       = {0.010, 0.010};
static REAL4            twoInterfRhosqThresh = 50.0;
static REAL4            mass                 = 1.4;
static REAL4            dynRange1            = 1.0; /*other options: 69.0; 1.0*/
static REAL4            dynRange2            = 1.0; /*other options: 69.0; 1.0*/
static UINT4            site0                = 0;
static UINT4            site1                = 1;
static UINT4            startSeconds         = 10000;


int
main (int argc, char *argv[])
{
  static LALStatus              status;

  INT4                          i;
  INT4                          j;
  INT4                          k;
  INT4                          l;
  INT4                          n;

  INT4                          seed[2] = {1,123};

  INT4                          flag = 0;

  FILE                         *fpData[2];
  FILE                         *fpRhosq[2];
  FILE                         *fpTwoInterfRhosq;
  
  REAL4                         sigma;
  REAL4                         Sfk, Sfk1, Sfk2;
  REAL4                         respRe; 
  REAL4                         respIm;
  REAL4                         deltaT;
  REAL4                         deltaF;   
  REAL4                         dynRange[2] = {dynRange1,dynRange2}; 
  REAL4Vector                  *noiseVec[2] = {NULL,NULL};

  RandomParams                 *randParams[2] = {NULL,NULL};
  LALDetectorPair               detectors;
  DataSegmentVector            *dataSegVec = NULL;
  InspiralTemplate             *tmplt[2] = {NULL,NULL};

  FindChirpSPDataParams        *dataParams   = NULL;
  FindChirpSPTmpltParams       *tmpltParams[2]  = {NULL,NULL};
  FindChirpFilterParams        *filterParams = NULL;
  FindChirpFilterInput         *filterInput  = NULL;
  FindChirpInitParams          *initParams[2]={NULL,NULL}; 
  
  TwoInterfInspiralEvent                *twoInterfInspiralEvent = NULL;
  TwoInterfFindChirpSPDataParamsVector  *twoInterfDataParamsVec = NULL;
  TwoInterfFindChirpInitParams          *twoInterfInitParams = NULL;
  TwoInterfFindChirpFilterParams        *twoInterfFilterParams = NULL;
  TwoInterfDataSegmentVector            *twoInterfDataSegVec = NULL;
  TwoInterfFindChirpSegmentVector       *twoInterfFcSegVec = NULL;
  TwoInterfFindChirpFilterInputVector   *twoInterfFilterInputVec = NULL;

  RealFFTPlan                  *powerFFTPlan = NULL;

  REAL4                         s[3];
  REAL4                         distance;
 
  /*
   *
   * parse options, allocate memory, init params and set values
   *
   */
  
  ParseOptions (argc, argv);
  
  /*
   *
   * Write network and individual detector rho-squares 
   *
   */
  
  /* override numSegments if outputting rhosq */
  if ( rhosqout )
    {
      numSegments = 1; 
      numTmplts = 1;
      fpRhosq[0] = fopen ("rhosq1.dat", "w");
      fpRhosq[1] = fopen ("rhosq2.dat", "w");
    }

  if ( twoInterfRhosqout )
    {
      numSegments = 1; 
      numTmplts = 1;
      fpTwoInterfRhosq= fopen ("twoInterfRhosq.dat", "w");
    }
  
  twoInterfInitParams = (TwoInterfFindChirpInitParams *) 
    LALMalloc (sizeof(TwoInterfFindChirpInitParams));
  
  twoInterfInitParams->numDetectors            = 2; /*hardwired to 2*/
  twoInterfInitParams->numSegments             = numSegments;
  twoInterfInitParams->numPoints               = numPoints;
  twoInterfInitParams->numChisqBins            = numChisqBins;
  twoInterfInitParams->createRhosqVec          = rhosqout;
  twoInterfInitParams->createTwoInterfRhosqVec = twoInterfRhosqout;


  if (inputDataType == gaussian) {
    for ( n = 0; n < 2; ++n )
      {
	LALCreateVector (&status, &noiseVec[n], numPoints);
	TestStatus (&status, "0", 1);
	ClearStatus (&status);
	
	LALCreateRandomParams (&status, &randParams[n], seed[n]);
	TestStatus (&status, "0", 1);
	ClearStatus (&status);
      }
  }  

  /*
   *
   * allocate memory for segments for DETECTOR 1
   *
   */
  
  LALCreateTwoInterfDataSegmentVector (&status, &twoInterfDataSegVec, 
				       twoInterfInitParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  LALCreateTwoInterfFindChirpSegmentVector (&status, &twoInterfFcSegVec, 
					    twoInterfInitParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  LALCreateTwoInterfFindChirpInputVector (&status, &twoInterfFilterInputVec, 
					  twoInterfInitParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  /*allocate memory to structures needed to create stationary phase templates*/
  tmplt[0] = twoInterfFilterInputVec->filterInput[0].tmplt = 
    (InspiralTemplate *) LALMalloc (sizeof(InspiralTemplate));
  memset( tmplt[0], 0, sizeof(InspiralTemplate) );
  
  tmplt[1] = twoInterfFilterInputVec->filterInput[1].tmplt = 
    (InspiralTemplate *) LALMalloc (sizeof(InspiralTemplate));
  memset( tmplt[1], 0, sizeof(InspiralTemplate) );

  /*
   *
   * initialize functions and create parameter structures for DETECTOR 1
   *
   */
  
  
  LALTwoInterfFindChirpFilterInit (&status, &twoInterfFilterParams, 
				   twoInterfInitParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  LALTwoInterfFindChirpSPDataInit (&status, &twoInterfDataParamsVec, 
				   twoInterfInitParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  if (inputDataType == file)
    /* create fft plan for generating power spectra */
    LALCreateForwardRealFFTPlan(&status, &powerFFTPlan, numPoints, 0);

  for ( n = 0 ; n < 2 ; ++n)
    {
      initParams[n] = (FindChirpInitParams *) 
	LALMalloc (sizeof(FindChirpInitParams));
      initParams[n]->numSegments       = twoInterfInitParams->numSegments;
      initParams[n]->numPoints         = twoInterfInitParams->numPoints;
      initParams[n]->numChisqBins      = twoInterfInitParams->numChisqBins;
      initParams[n]->createRhosqVec    = rhosqout;
      
      LALFindChirpSPTemplateInit (&status, &tmpltParams[n], initParams[n]);
      TestStatus (&status, "0", 1);
      ClearStatus (&status);

      /* moved from 'abort' section */
      LALFree (initParams[n]);
      initParams[n] = NULL;
    }

  /* moved from 'abort' section */
  LALFree (twoInterfInitParams);
  twoInterfInitParams = NULL;

  /*
   *
   * fill segments with data, spectrum and response
   *
   */
  

  /* analytical stuff */
  sigma  = sqrt( sigmasq );
  deltaT = 1.0 / (float) srate;
  deltaF = 1.0 / (numPoints * deltaT);
  Sfk    = 2.0 * sigmasq * deltaT;
  Sfk1    = 2.0 * sigmasq1* deltaT;
  Sfk2    = 2.0 * sigmasq2* deltaT;
  respRe = 1.0;
  respIm = 0.0;

  /* print input parameters */
  fprintf( stdout, "\n              srate = %d\n          numPoints = %d\n",
	   srate, numPoints );
  fprintf( stdout, "        numSegments = %d\n\n", numSegments );
  fprintf( stdout, "            sigma^2 = %5.2f\n", sigmasq );
  fprintf( stdout, "                Sfk = %10.8f\n", Sfk );
  fprintf( stdout, "                Sfk1 = %10.8f\n", Sfk1 );
  fprintf( stdout, "                Sfk2 = %10.8f\n", Sfk2 );
  fprintf( stdout, "             deltaT = %e\n             deltaF = %e\n\n", 
	   deltaT, deltaF);
  fprintf( stdout, "       invSpecTrunc = %d\n", invSpecTrunc );
  fprintf( stdout, "               fLow = %5.3f\n\n", fLow );
  fprintf( stdout, "     rhosqThreshold = %5.3f\n\n", rhosqThresh[0]);
  fprintf( stdout, "     rhosqThreshold = %5.3f\n\n", rhosqThresh[1]);
  fprintf( stdout, "    chisqThreshold1 = %5.3f\n\n", chisqThresh[0]);
  fprintf( stdout, "    chisqThreshold2 = %5.3f\n\n", chisqThresh[1]);

  fprintf( stdout, "     twoInterfRhosqThresh = %5.3f\n", twoInterfRhosqThresh);
  fprintf( stdout, "               mass = %5.2f\n\n", mass );

  /*
   *
   * fill segments for DETECTOR 1
   *
   */
  
  dataSegVec = twoInterfDataSegVec->data;
  
  for ( i = 0; i < dataSegVec[0].length; ++i )
    {
      REAL4 time;
      dataSegVec[0].data[i].chan->deltaT = (REAL8) deltaT;
      dataSegVec[0].data[i].spec->deltaF = (REAL8) deltaF;
      
      time = numPoints * deltaT * (REAL4) i + startSeconds;
      dataSegVec[0].data[i].chan->epoch.gpsSeconds     = (INT4) floor( time );
      time = (time - floor( time )) * 1.0E9;
      dataSegVec[0].data[i].chan->epoch.gpsNanoSeconds = (INT4) floor( time );
      
      dataSegVec[0].data[i].resp->epoch = dataSegVec[0].data[i].chan->epoch;

      /* impulse */
      if ( inputDataType == impulse )
	{
	  memset( dataSegVec[0].data[i].chan->data->data, 0, 
		  numPoints * sizeof(REAL4) );
	  
	  dataSegVec[0].data[i].chan->data->data[0] = 1.0;
	  
	  /* spectrum and response */
	  for ( k = 0; k < numPoints/2 + 1; ++k )
	    {
	      dataSegVec[0].data[i].spec->data->data[k]    = 1.0;
	      dataSegVec[0].data[i].resp->data->data[k].re = 1.0;
	      dataSegVec[0].data[i].resp->data->data[k].im = 0.0;
	    }
	  
	}
      /* gaussian noise */
      else if ( inputDataType == gaussian )
	{
	  LALNormalDeviates (&status, noiseVec[0], randParams[0]);
	  TestStatus (&status, "0", 1);
	  ClearStatus (&status);
	  
	  /* ifodmro */
	  for ( j = 0; j < numPoints; ++j )
	    {
	      REAL4 noise;
	      REAL4 ifodmro;
	      
	      noise = floor( 0.5 + sigma * noiseVec[0]->data[j] );
	      
	      if ( noise > -2048 && noise < 2047 ) 
		ifodmro= noise;
	      else if ( noise < -2048 ) 
		ifodmro = -2048.0;
	      else 
		ifodmro=2047.0;
	      
	      dataSegVec[0].data[i].chan->data->data[j] = ifodmro;
	    }
	  
	  /* spectrum and response */
	  for ( k = 0; k < numPoints/2 + 1; ++k )
	    {
	      dataSegVec[0].data[i].spec->data->data[k]    = Sfk;
	      dataSegVec[0].data[i].resp->data->data[k].re = 1.0;
	      dataSegVec[0].data[i].resp->data->data[k].im = 0.0;
	    }
	}
      else  if ( inputDataType == file )


	/* this should eventually read in start epoch */
	{
	  /* open the input files */
	  if ( !(fpData[0] = fopen( "data1.dat", "r" )) )
	    {
	      fprintf( stdout, "unable to open the data file from detector 1 for reading\n" );
	      fflush( stdout );
	      goto abort;
	    }
	  
	  /* read in ifodmro data */ 
	  for ( j = 0; j < numPoints; ++j )
	    {
	      if (( (flag = fscanf( fpData[0], "%f\n", 
				    &(dataSegVec[0].data[i].chan->data->data[j]) )) != 1 || flag == EOF ) 
		  && j < numPoints ) 
		{
		  fprintf( stdout, "error reading input data\n" );
		  fflush( stdout );
		  fclose( fpData[0] );
		  goto abort;
		}
	    }
	  
	  ReadResponse(&status, dataSegVec[0].data[i].resp, "response1.dat", invertResponse);

	  /* close the file */
	  fclose( fpData[0] );
	
	  LALRealPowerSpectrum(&status, dataSegVec[0].data[i].spec->data,
			       dataSegVec[0].data[i].chan->data, powerFFTPlan);
          /* CHANGE BACK TO ORIGINAL NORMALIZATION -- JC */
          {
            REAL4Vector *myvector = dataSegVec[0].data[i].spec->data;
            UINT4 mybin;
            for ( mybin = 1; mybin < myvector->length - 1; ++mybin )
              myvector->data[mybin] *= 0.5;
          }


	
	  /* scale the power spectrum to the findchirp definition */
	  {
	    UINT4 j;
	    for(j=0; j<dataSegVec[1].data[i].spec->data->length; j++)

#ifdef CALCULATE_PSD
	      dataSegVec[0].data[i].spec->data->data[j] *= 4.0 * deltaT / numPoints;
#endif

#ifndef CALCULATE_PSD
	    dataSegVec[0].data[i].spec->data->data[j] = Sfk1;
#endif
	  }

	}
      
    }

  /*
   *
   * fill segments for DETECTORS 2
   *
   */
  
  for ( i = 0; i < dataSegVec[1].length; ++i )
    {
      REAL4 time;

      dataSegVec[1].data[i].chan->deltaT = (REAL8) deltaT;
      dataSegVec[1].data[i].spec->deltaF = (REAL8) deltaF;
      
      time = numPoints * deltaT * (REAL4) i + startSeconds;
      dataSegVec[1].data[i].chan->epoch.gpsSeconds     = (INT4) floor( time );
      time = (time - floor( time )) * 1.0E9;
      dataSegVec[1].data[i].chan->epoch.gpsNanoSeconds = (INT4) floor( time );
      
      dataSegVec[1].data[i].resp->epoch = dataSegVec[1].data[i].chan->epoch;

      /* impulse */
      if ( inputDataType == impulse )
	{
	  memset( dataSegVec[1].data[i].chan->data->data, 0, 
		  numPoints * sizeof(REAL4) );
	  
	  dataSegVec[1].data[i].chan->data->data[0] = 1.0;
	  
	  /* spectrum and response */
	  for ( k = 0; k < numPoints/2 + 1; ++k )
	    {
	      dataSegVec[1].data[i].spec->data->data[k]    = 1.0;
	      dataSegVec[1].data[i].resp->data->data[k].re = 1.0;
	      dataSegVec[1].data[i].resp->data->data[k].im = 0.0;
	    }
	  
	}
      /* gaussian noise */
      else if ( inputDataType == gaussian )
	{
	  LALNormalDeviates (&status, noiseVec[1], randParams[1]);
	  TestStatus (&status, "0", 1);
	  ClearStatus (&status);
	  
	  /* ifodmro */
	  for ( j = 0; j < numPoints; ++j )
	    {
	      REAL4 noise;
	      REAL4 ifodmro;
	      
	      noise = floor( 0.5 + sigma * noiseVec[1]->data[j] );
	      
	      if ( noise > -2048 && noise < 2047 ) 
		ifodmro= noise;
	      else if ( noise < -2048 ) 
		ifodmro = -2048.0;
	      else 
		ifodmro=2047.0;
	      
	      dataSegVec[1].data[i].chan->data->data[j] = ifodmro;
	    }
	  
	  /* spectrum and response */
	  for ( k = 0; k < numPoints/2 + 1; ++k )
	    {
	      dataSegVec[1].data[i].spec->data->data[k]    = Sfk;
	      dataSegVec[1].data[i].resp->data->data[k].re = 1.0;
	      dataSegVec[1].data[i].resp->data->data[k].im = 0.0;
	    }
	}
      else if ( inputDataType == file )
	
	/* this should eventually read in start epoch somehow */
	{
	  /* open the input files */
	  if ( !(fpData[1] = fopen( "data2.dat", "r" )) )
	    {
	      fprintf( stdout, "unable to open the data file from detector 1 for reading\n" );
	      fflush( stdout );
	      goto abort;
	    }

	  /* read in ifodmro data */ 
	  for ( j = 0; j < numPoints; ++j )
	    {
	      if (( (flag = fscanf( fpData[1], "%f\n", 
				    &(dataSegVec[1].data[i].chan->data->data[j]) )) != 1 || flag == EOF ) 
		  && j < numPoints ) 
		{
		  fprintf( stdout, "error reading input data\n" );
		  fflush( stdout );
		  fclose( fpData[1] );
		  goto abort;
		}
	    }

	  ReadResponse(&status, dataSegVec[1].data[i].resp, "response2.dat", invertResponse);
	  
	  /* close the files */
	  fclose( fpData[1] );

	  LALRealPowerSpectrum(&status, dataSegVec[1].data[i].spec->data,
			       dataSegVec[1].data[i].chan->data, powerFFTPlan);
          /* CHANGE BACK TO ORIGINAL NORMALIZATION -- JC */
          {
            REAL4Vector *myvector = dataSegVec[1].data[i].spec->data;
            UINT4 mybin;
            for ( mybin = 1; mybin < myvector->length - 1; ++mybin )
              myvector->data[mybin] *= 0.5;
          }



	  /* scale the power spectrum to corresponds with the findchirp package definition */
	 
	  {
	    UINT4 j;
	    for(j=0; j<dataSegVec[1].data[i].spec->data->length; j++)
#ifdef CALCULATE_PSD
	      dataSegVec[1].data[i].spec->data->data[j] *= 4.0 * deltaT / numPoints;
#endif

#ifndef CALCULATE_PSD
	    dataSegVec[1].data[i].spec->data->data[j] = Sfk2;
#endif 
	  }
	 
	}
    }

  if (inputDataType == file)
    LALDestroyRealFFTPlan(&status, &powerFFTPlan);
  
  /*
   *
   * condition data from DETECTORS for stationary phase chirps
   *
   */
  
  dataParams = twoInterfDataParamsVec->data;
  
  for ( n = 0 ; n < 2 ; ++n )
    {
      /* dynamic range */
      dynRange[n] = pow( 2.0, dynRange[n] );

      /* set parameters */
      dataParams[n].deltaT       = deltaT;
      dataParams[n].fLow         = fLow;
      dataParams[n].dynRange     = dynRange[n];
      dataParams[n].invSpecTrunc = invSpecTrunc;
    }/* ends loop over detectors */
  
  LALTwoInterfFindChirpSPData (&status, twoInterfFcSegVec, twoInterfDataSegVec, twoInterfDataParamsVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  /* moved from 'abort' section */
  LALTwoInterfFindChirpSPDataFinalize (&status, &twoInterfDataParamsVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALDestroyTwoInterfDataSegmentVector (&status, &twoInterfDataSegVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  /* template generation parameters */

  for ( n = 0 ; n < 2 ; ++n )
    {
      tmpltParams[n]->deltaT   = deltaT;
      tmpltParams[n]->dynRange = dynRange[n];
      tmpltParams[n]->fLow     = fLow;
    }/* ends loop over detectors */
  

  /* network filter parameters */
  
  detectors.detectorOne = lalCachedDetectors[site0];
  detectors.detectorTwo = lalCachedDetectors[site1];

  /* calculate separation vector between sites */
  for ( j=0; j<3; j++) {
    s[j] = (REAL4) ( detectors.detectorOne.location[j] -
		     detectors.detectorTwo.location[j]);
  }
  /* calculate distance between sites (in meters) */
  distance = sqrt( cartesianInnerProduct(s,s) );

  twoInterfFilterParams->twoInterfRhosqThresh = twoInterfRhosqThresh;
  twoInterfFilterParams->detectors            = &detectors;

  filterParams = twoInterfFilterParams->paramsVec->filterParams;
  
  /*  PrintAllocs(0, "before Chisq veto init" ); */

  for ( n = 0; n < 2; ++n )
    {
      /* filter parameters for each detector */
      filterParams[n].deltaT               = deltaT;
      filterParams[n].rhosqThresh          = rhosqThresh[n];
      filterParams[n].chisqThresh          = chisqThresh[n];
      filterParams[n].computeNegFreq       = 0;
      filterParams[n].maximiseOverChirp    = 1; 

      LALTwoInterfFindChirpChisqVetoInit (&status, filterParams[n].chisqParams, 
				 numChisqBins, numPoints);
      TestStatus (&status, "0", 1); 
      ClearStatus (&status);
    }/* ends loop over detectors */

  /* PrintAllocs(0, "after Chisq veto init" ); */

  /*
   *
   * loop over templates
   *
   */

  filterInput = twoInterfFilterInputVec->filterInput;

  for ( l = 0; l < numTmplts; ++l, mass +=0.001 )
    {
      /* loop over detectors */
      for ( n = 0; n < 2; ++n )
	{
	  /*
	   *
	   * generate stationary phase template
	   *
	   */
	  
	  /* generate dummy template parameters */
	  REAL4 m1 = mass;
	  REAL4 m2 = mass;
	  
	  tmplt[n]->mass1     = m1;
	  tmplt[n]->mass2     = m2;
	  tmplt[n]->totalMass = m1 + m2;
	  tmplt[n]->mu        = m1 * m2 / tmplt[n]->totalMass;
	  tmplt[n]->eta       = tmplt[n]->mu / tmplt[n]->totalMass;
	  
	  LALFindChirpSPTemplate (&status, filterInput[n].fcTmplt, tmplt[n], tmpltParams[n]);
	  TestStatus (&status, "0", 1);
	  ClearStatus (&status);
	  
	  LALFindChirpSPTemplateFinalize (&status, &tmpltParams[n]);
	  TestStatus (&status, "0", 1);
	  ClearStatus (&status);

	  /*
	   *
	   * loop over segments
	   *
	   */
	  
	  for ( i = 0; i < numSegments; ++i )
	    {
	      twoInterfFilterInputVec->filterInput[n].segment = twoInterfFcSegVec->data[n].data + i;
	      
	    }/* end loop over number of segments */
	
	}/* end loop over detectors */
      
      twoInterfInspiralEvent = NULL;

      LALTwoInterfFindChirpFilterSegment (&status, &twoInterfInspiralEvent, twoInterfFilterInputVec, twoInterfFilterParams); 
      TestStatus (&status, "0", 1);
      
      if ( twoInterfInspiralEvent )
	{
	  fprintf( stdout, "\nEvents found in segment!\n\n" );
	  while (twoInterfInspiralEvent  )
	    {
	      TwoInterfInspiralEvent *thisEvent = twoInterfInspiralEvent;
	      INT4  delayPts = (INT4) thisEvent->timeIndex2 - (INT4)thisEvent->timeIndex;
	      REAL4 delay = delayPts * deltaT;

	      twoInterfInspiralEvent = thisEvent->next;
	      
	      fprintf( stdout, "event id                              = %d\n\n", thisEvent->twoInterfId+1 );
	      
	      fprintf( stdout, "delay from d1 to d2                   = %.2f ms\n\n", delay * 1000);

	      if (site0 != site1) {
		fprintf( stdout, "Half cone angle                       = %.2f\n", thisEvent->twoInterfAngle );
		fprintf( stdout, "Half cone central axis                = %f radians RA\n", thisEvent->twoInterfAxisRa);
		fprintf( stdout, "                                        %f radians dec\n\n", thisEvent->twoInterfAxisDec);
              } 
	      else
		fprintf( stdout, "Cannot determined half cone angle, because detectors 1 and 2 are coincident\n\n");


	      fprintf( stdout, "network snrsq                         = %.2f\n", thisEvent->snrsq );

	      fprintf( stdout, "detector 1 snrsq                      = %.2f\n", thisEvent->eventIn1->snrsq );
	      fprintf( stdout, "detector 2 snrsq                      = %.2f\n\n", thisEvent->eventIn2->snrsq );
	      
	      fprintf( stdout, "source distance estimated from det. 1 = %.4f\n", thisEvent->eventIn1->effDist * 1000); 
	      fprintf( stdout, "source distance estimated from det. 2 = %.4f\n\n", thisEvent->eventIn2->effDist * 1000); 
	      
	      fprintf( stdout, "chisq in detector 1                   = %.4f\n", thisEvent->chisq1);
	      fprintf( stdout, "chisq in detector 2                   = %.4f\n\n", thisEvent->chisq2);
		
	      fprintf( stdout, "event GPS time in DETECTOR 1          = %d.%d\n", 
		       thisEvent->time.gpsSeconds, thisEvent->time.gpsNanoSeconds );
	      fprintf( stdout, "'network' timeIndex                   = %d\n", thisEvent->timeIndex );
	      
	      fprintf( stdout, "DETECTOR 2 timeIndex                  = %d\n\n", thisEvent->eventIn2->timeIndex );
	      
	      fprintf( stdout, "distance between detectors (in meters)= %e\n\n", distance );
	      
	      fprintf( stdout, "m1                                    = %f\n", thisEvent->tmplt.mass1 );
	      fprintf( stdout, "m2                                    = %f\n\n", thisEvent->tmplt.mass2 );

	      fflush( stdout );
	      
	      LALFree( thisEvent->eventIn1 );
	      LALFree( thisEvent->eventIn2 );
	      LALFree( thisEvent );
	    }
	  
	}
      else
	{
	  fprintf( stdout, "\nNo events found in segment for mass = %.4f solar mass template!\n", mass );
	}
    }/* end loop over templates */
  
  
  
  /* write rhosq of detectors 1 and 2 to file */
  if ( rhosqout )
    {
      for ( j = 0; j < filterParams[0].rhosqVec->data->length; ++j )
	{
	  fprintf( fpRhosq[0], "%d\t%e\n", j, filterParams[0].rhosqVec->data->data[j] );
	}
      
      for ( j = 0; j < filterParams[1].rhosqVec->data->length; ++j )
	{
	  fprintf( fpRhosq[1], "%d\t%e\n", j, filterParams[1].rhosqVec->data->data[j] );
	}
    }
  
  /* write network rhosq to file */
  if ( twoInterfRhosqout )
    {
      for ( j = 0; j < twoInterfFilterParams->twoInterfRhosqVec->data->length; ++j )
	{
	  fprintf( fpTwoInterfRhosq, "%d\t%e\n", j, twoInterfFilterParams->twoInterfRhosqVec->data->data[j] );
	}
    }
  
  /*
   *
   * finalize functions and destroy parameter structures
   *
   */
  
 abort:
  
  for (n= 0; n<2; ++n)
    {
      LALDestroyRandomParams (&status, &randParams[n]);
      TestStatus (&status, "0", 1);
      ClearStatus (&status);
      
      LALDestroyVector (&status, &noiseVec[n]);
      TestStatus (&status, "0", 1);
      ClearStatus (&status);
      
      LALTwoInterfFindChirpChisqVetoFinalize (&status, filterParams[n].chisqParams,
					      numChisqBins);
      TestStatus (&status, "0", 1);
      ClearStatus (&status);
    }
  
  LALTwoInterfFindChirpFilterFinalize (&status, &twoInterfFilterParams);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);

  LALDestroyTwoInterfFindChirpSegmentVector (&status, &twoInterfFcSegVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  LALDestroyTwoInterfFindChirpInputVector (&status, &twoInterfFilterInputVec);
  TestStatus (&status, "0", 1);
  ClearStatus (&status);
  
  for( n = 0 ; n < 2 ; ++n )
    {
      LALFree (tmplt[n]);
      tmplt[n] = NULL;
    }
  
  LALCheckMemoryLeaks ();
  
  return 0;
}


/*
 * TestStatus ()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus (
    LALStatus  *status, 
    const char *ignored, 
    int         exitcode
           )
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    REPORTSTATUS (status);
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}


/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus (
    LALStatus   *status
            )
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}


static void
Usage (
    const char *program, 
    int         exitcode
      )
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options (defaults shown in square brackets):\n");
  fprintf (stderr, "  General:\n");
  fprintf (stderr, "    -h                print this message\n");
  fprintf (stderr, "    -V                verbose\n");
  fprintf (stderr, "    -d debuglevel     LAL status debug level [1]\n");
  fprintf (stderr, "    -o                write single-detector rhosq output files [0]\n");
  fprintf (stderr, "    -u                write network rhosq (maximized over time delay) output files\n");
  fprintf (stderr, "  Data Creation:\n");
  fprintf (stderr, "    -I(i|g|f)         input data is (i)mpulse, (g)aussian, (f)file [g]\n");
  fprintf (stderr, "    -C                response files given in counts / strain\n");
  fprintf (stderr, "    -n numPoints      number of points in a segment [1048576]\n");
  fprintf (stderr, "    -s numSegments    number of data segments [1]\n");
  fprintf (stderr, "    -r srate          sampling rate of the data [8192]\n");
  fprintf (stderr, "    -v sigmasq        rms of the gaussian noise [8.0]\n");
  fprintf (stderr, "    -k sigmasq1       rms of the gaussian noise in detector 1[8.0]\n");
  fprintf (stderr, "    -K sigmasq2       rms of the gaussian noise in detector 2[8.0]\n");
  fprintf (stderr, "    -m mass           mass of each binary component in solar masses [1.4]\n");
  fprintf (stderr, "  Data Conditioning:\n");
  fprintf (stderr, "    -i invSpecTrunc   number of points to truncate S^{-1}(f) [0]\n");
  fprintf (stderr, "    -f fLow           low frequency cutoff for S^{-1}(f) [150.0]\n");
  fprintf (stderr, "    -y dynRange1       log_2( dynamicRange1 ) [1.0]\n");
  fprintf (stderr, "    -z dynRange2       log_2( dynamicRange2 ) [1.0]\n");
  fprintf (stderr, "  Filtering:\n");
  fprintf (stderr, "    -N numTmplts      number templates to filter against [1]\n");
  fprintf (stderr, "    -b numChisqBins   number of bins for chi squared test [8]\n");
  fprintf (stderr, "    -t rhosqThresh    signal to noise squared threshold [50.0]\n");
  fprintf (stderr, "    -T twoInterfRhosqThresh    signal to noise squared threshold for network [50.0]\n");
  fprintf (stderr, "    -x chisqThresh1   chi squared threshold for detector 1 [0.01]\n");
  fprintf (stderr, "    -X chisqThresh2   chi squared threshold for detector 2 [0.01]\n");
  fprintf (stderr, "    -p site0          first detector's site number [0]\n");
  fprintf (stderr, "    -q site1          second detector's site number [1]\n");
  exit (exitcode);
}


static void
ParseOptions (
    int         argc, 
    char       *argv[]
             )
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "Vhd:""n:""s:""y:""z:""r:""v:""k:""K:""i:""f:""R:""t:""T:""b:""m:""N:""x:""X:""oI:""uI:""p:""q:""C:");
    if (c == -1)
    {
      break;
    }
    
    switch (c)
      {
      case 'd': /* set debuglevel */
        lalDebugLevel = atoi (optarg);
        break;
      case 'V': /* set verbosity */
        verbose = 1;
        break;
      case 'o': /* sets flag to write rhosq */
        rhosqout = 1;
        break;
      case 'u': /* sets flag to write twoInterfRhosq */
        twoInterfRhosqout = 1;
        break;
      case 'n': /* set number of points in a segment */
        numPoints = atoi (optarg);
        break;
      case 's': /* set number of segments */
        numSegments = atoi (optarg);
        break;
      case 'r': /* set sampling rate */
        srate = atoi (optarg);
        break;
      case 'i': /* set invSpecTrunc */
        invSpecTrunc = atoi (optarg);
        break;
      case 'N': /* set numChisqBins */
        numTmplts = atoi (optarg);
        break;
      case 'b': /* set numChisqBins */
        numChisqBins = atoi (optarg);
        break;
      case 'f': /* set fLow */
        fLow = (REAL4) atof (optarg);
        break;
      case 'y': /* set dynRange1 */
        dynRange1 = (REAL4) atof (optarg);
        break;
      case 'z': /* set dynRange2 */
        dynRange2 = (REAL4) atof (optarg);
        break;
      case 'm': /* mass */
        mass = (REAL4) atof (optarg);
        break;
      case 't': /* set rhosq threshold */
        rhosqThresh[0] = rhosqThresh[1] = (REAL4) atof (optarg);
        break;
      case 'T': /* set twoInterfRhosq threshold */
        twoInterfRhosqThresh = (REAL4) atof (optarg);
        break;
      case 'x': /* set chisq threshold */
        chisqThresh[0] = (REAL4) atof (optarg);
        break;
      case 'X': /* set chisq threshold */
        chisqThresh[1] = (REAL4) atof (optarg);
        break;
      case 'v': /* set variance */
        sigmasq = (float) atof (optarg) * (float) atof(optarg);
        break;
      case 'k': /* set noise variance in detector 2*/
        sigmasq1 = (float) atof (optarg) * (float) atof(optarg);
        break;
      case 'K': /* set noise variance in detector 2*/
        sigmasq2 = (float) atof (optarg) * (float) atof(optarg);
        break;
      case 'I': /* input data */
        switch (*optarg)
	  {
          case 'i':
            inputDataType = impulse;
            break;
          case 'g':
            inputDataType = gaussian;
            break;
          case 'f':
            inputDataType = file;
            break;
          default:
            inputDataType = gaussian;
	  }
        break;

      case 'C': /* invert response function? 
                   uses opposite convention as injection code 
                   (hence the inversion after input here)     */
	invertResponse = atoi(optarg);
	break;

      case 'h':
        Usage (argv[0], 0);
        break;
      
      case 'p':  /*site value of detector 1 */
	site0 = (UINT4) atof (optarg);
	break;
      case 'q':  /*site value of detector 2 */
	site1 = (UINT4) atof (optarg);
	break;

      default:
        Usage (argv[0], 1);
      }
  }
  
  if (optind < argc)
    {
      Usage (argv[0], 1);
    }
  
  return;
}


/*
 *
 * graphREAL4 ()
 *
 * function to graph an array of REAL4's for debugging
 *
 */
void
graphREAL4 (
    REAL4      *array, 
    INT4        n,
    INT4        spacing
           ) 
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; i++ )
    fprintf( fp, "%d\t%e\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

  return;
}


/*
 *
 * graphINT2 ()
 *
 * function to graph an array of INT2's for debugging
 *
 */
void 
graphINT2 (
    INT2       *array, 
    INT4        n,
    INT4        spacing
          ) 
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; i++ )
    fprintf( fp, "%d\t%d\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

  return;
}


/*
 *
 * graphINT4 ()
 *
 * function to graph an array of INT4's for debugging
 *
 */
void 
graphINT4 (
    INT4       *array, 
    INT4        n,
    INT4        spacing
          ) 
{
  FILE *fp;
  INT4 i;

  /* open a file for writing */
  if ( !(fp = fopen( "temp.graph", "w" )) )
  {
    printf( "couldn't open file\n" );
  }

  /* print data into the file */
  for ( i = 0; i < n; i++ )
    fprintf( fp, "%d\t%d\n", i, array[i * spacing] );

  /* close the file */
  fclose( fp );

  /* start up graphing program with data in the file */
  /* system( "xmgr temp.graph 1>/dev/null 2>&1 &" ); */

  return;
}


/* read a response file of format 'frequency  real response  imaginary response' */
/* considers any line that it can't convert to floating point numbers to be a comment */

/* sets f0 to the first frequency point, deltaF to the difference of the first
 * and second frequency, and reads the COMPLEX8 values of the response from the file */

/* The 'invert' flag will invert the response if set. This is to make up for the
 * inconsistent definition of the response/transfer function */

/* does not initialize the name, units, nor epoch of the response */ 

void ReadResponse(LALStatus *status, COMPLEX8FrequencySeries *resp, char* respFile, BOOLEAN invert) {
  UINT4 fLength, i, nScanned;
  REAL4 frequency;
  FILE  *respPtr;
  CHAR  garbage;

  fLength = resp->data->length;
  respPtr = fopen(respFile, "r");

  if ( ! respPtr ) {
    fprintf(stderr, "Error reading response function from %s\n", respFile);
    exit(0);
  }
  
  for(i=0; i<fLength; i++) {
    /* read data into frequency series */
    nScanned = fscanf(respPtr, "%g %g %g\n", &frequency, &resp->data->data[i].re, 
	   &resp->data->data[i].im);

    /* read start frequency */
    if (i==0) 
      resp->f0 = frequency;

    /* determine frequency sampling */
    else if (i==1)
      resp->deltaF = frequency - resp->f0;

    /* if conversion failed, this is a comment line, so skip it */
    if (nScanned == 0) {
      i--;
      garbage = 'a';
      while( garbage != EOF && garbage != '\n')
	fscanf(respPtr, "%c", &garbage);
    }

    /* if EOF or other error, exit */
    else if ( nScanned != 3 ) {
      fprintf(stderr, "Error reading %s\n", respFile);
      exit(0);
    }
  }

  /* if response file was specified incorrectly, invert it */
  if ( invert ) {
    REAL4 modsq;

    for(i=0; i<fLength; i++) { 
      if ( resp->data->data[i].re == 0.0 && resp->data->data[i].im == 0.0 ) {
	fprintf(stderr, "Error: response function is zero at index %d\n", i);
	exit(0);
      }

      modsq = resp->data->data[i].re * resp->data->data[i].re +
	resp->data->data[i].im * resp->data->data[i].im;

      resp->data->data[i].re /= modsq; 
      resp->data->data[i].im /= 0.0 - modsq; 
    }
  }

  fclose(respPtr);
}
