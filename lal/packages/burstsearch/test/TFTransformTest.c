/*----------------------------------------------------------------------- 
 * 
 * File Name: TFTransformTest.c
 * 
 * Author: Eanna Flanagan
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * main()
 *
 * SYNOPSIS 
 * 
 * DESCRIPTION 
 * Test suite for functions in TFTransform.c
 * 
 * DIAGNOSTICS
 * Writes PASS or FAIL to stdout as tests are passed or failed.
 *
 * CALLS
 *
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "LALStdlib.h"
#include "SeqFactories.h"
#include "VectorOps.h"
#include "TFTransform.h"
#include "PrintVector.h"


#define _CODES(x) #x
#define CODES(x) _CODES(x)


NRCSID (MAIN, "$Id$");


extern char *optarg;
extern int   optind;

INT4 debuglevel = 1;   /* set to 2 to get full status information for tests */
INT4 verbose    = 1;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (Status *status, const char *expectedCodes, int exitCode);

static REAL4 ff(REAL4 w);   /* simple function used to construct a waveform */


int
main (int argc, char *argv[])
{
  const INT4 ntot   = 1000;   /* total number of points in time domain */
  const REAL8 alpha = 0.27;   /* ln(nt)/ln(ntot) */
  const REAL8 beta  = 0.2;    /* nf_actual / nf_total */

  static Status                 status; 
  TFPlaneParams                 params;
  VerticalTFTransformIn         transformparams;
  HorizontalTFTransformIn       transformparams1;
  REAL4TimeSeries               tseries;
  COMPLEX8FrequencySeries       fseries;
  COMPLEX8TimeFrequencyPlane    *tfp=NULL;
  RealDFTParams                 *dftparams1=NULL;


  INT4                          i;
  INT4                          tseglength;
  INT4                          fseglength;
  INT4                          nt;
  INT4                          nf;
  INT4                          nforig;

  /*
   *
   * Parse the command line options
   *
   */

  ParseOptions (argc, argv);



  /* compute parameters */

  nt = (INT4)(exp( alpha * log( (REAL8)(ntot))));
  nforig = ntot / (2 * nt);
  nf = (INT4)(beta * (REAL8)(nforig));

  if(verbose)
    {
      printf("Total number of data points :    %d\n",ntot);
      printf("Number of time bins         :    %d\n",nt);
      printf("Number of frequency bins    :    %d\n",nf);
      printf("Original # of freq bins     :    %d\n",nforig);
    }


    

  /* 
   *  
   *  Set up input time series
   *
   */

  tseries.epoch.gpsSeconds=0;
  tseries.epoch.gpsNanoSeconds=0;
  tseries.deltaT = 0.001;  /* 1 kHz sampling */
  tseries.f0 = 0.0;
  tseries.name = NULL;
  tseries.sampleUnits=NULL;
  tseries.data=NULL;


  SCreateVector (&status, &(tseries.data), ntot);
  TestStatus (&status, CODES(0), 1);

  for(i=0; i< tseries.data->length; i++)
    {
      tseries.data->data[i] = ff( (REAL4)(i)/ (REAL4)(ntot));
    };

  /*
   *
   *   Set up time-frequency plane structure
   *
   */

  /* setup parameters structure for creating TF plane */
  params.timeBins = nt;
  params.freqBins = nf;
  params.deltaT = tseries.deltaT * 2.0 * (REAL8)(nforig);
  params.flow = 0.0;

  /* Create TF plane  */
  CreateTFPlane( &status, &tfp, &params);
  TestStatus (&status, CODES(0), 1);


  /*  
   *
   *
   *  Test of vertical TFTransform
   *
   *
   */
  
  {
    LALWindowParams winParams;
    winParams.type=Rectangular;

    tseglength = 2 * (INT4)( (params.deltaT) / (2.0*(tseries.deltaT)) );
    winParams.length=tseglength;

    /* setup input structure for computing TF transform */
    transformparams.startT=0;
    transformparams.dftParams=NULL;
    CreateRealDFTParams( &status, &(transformparams.dftParams), &winParams, 1); 
    TestStatus (&status, CODES(0), 1);


    /* Compute TF transform */
    if(verbose)
      {
	printf("Computing vertical time-frequency plane\n");
      }
    TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
    TestStatus (&status, CODES(0), 1);


    /* Destroy stuff */
    DestroyRealDFTParams( &status, &(transformparams.dftParams));
    TestStatus (&status, CODES(0), 1);
  }





  /*  
   *
   *
   *  Test of horizontal TFTransform
   *
   *
   */

  {
    LALWindowParams winParams;
    LALWindowParams winParams1;
    winParams.type=Rectangular;
    winParams.length=ntot;
    
    fseries.data=NULL;
    CCreateVector( &status, &(fseries.data), ntot/2+1);
    TestStatus (&status, CODES(0), 1);
      
    CreateRealDFTParams( &status, &dftparams1, &winParams, 1);
    TestStatus (&status, CODES(0), 1);
    
    if(verbose)
      {
	printf("Computing FFT of time domain data\n");
      }
    ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
    TestStatus (&status, CODES(0), 1);

    fseglength = (INT4)( 0.5+1/(tfp->params->deltaT * fseries.deltaF));

    /* setup input structure for computing TF transform */
    transformparams1.startT=0;  /* not used for horizontal transforms */
    transformparams1.dftParams=NULL;

    winParams1.type=Rectangular;
    winParams1.length=fseglength;

    CreateComplexDFTParams( &status, &(transformparams1.dftParams),
                            &winParams1,-1); 
    TestStatus (&status, CODES(0), 1);


    
    if(verbose)
      {
	printf("Computing horizontal time-frequency plane\n");
      }
    /* Compute TF transform */
    FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
    TestStatus (&status, CODES(0), 1);


    /* Destroy stuff */
    
    DestroyComplexDFTParams( &status, &(transformparams1.dftParams));
    TestStatus (&status, CODES(0), 1);
    
    DestroyRealDFTParams( &status, &dftparams1);
    TestStatus (&status, CODES(0), 1);

    CDestroyVector( &status, &(fseries.data));
    TestStatus (&status, CODES(0), 1);
  }


  
  DestroyTFPlane( &status, &tfp);
  TestStatus (&status, CODES(0), 1);

  SDestroyVector (&status, &(tseries.data) );
  TestStatus (&status, CODES(0), 1);





  /*************************************************************************
   *                                                                       * 
   *                                                                       *
   *  Now check to make sure that correct error codes are generated.       *
   *                                                                       *
   *                                                                       * 
   *************************************************************************/


  if (verbose || debuglevel)
  {
    printf ("\n===== Check Errors =====\n");
  }

  /* 
   *
   *  Test functions CreateRealDFTParams() and DestroyRealDFTParams()
   *
   */
  {
    LALWindowParams winParams;

    if (verbose)
      {
	printf ("\n--- Testing CreateRealDFTParams() and DestroyRealDFTParams() \n\n");
      }

    winParams.type=Rectangular;
    winParams.length=tseglength;

    CreateRealDFTParams( &status, NULL, &winParams, 1); 
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    winParams.length=0;
    CreateRealDFTParams( &status, &(transformparams.dftParams), &winParams, 1);  
    TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
    winParams.length=tseglength;
    
    transformparams.dftParams=NULL;
    DestroyRealDFTParams( &status, &(transformparams.dftParams)); 
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    /* next few tests require a valid DFTParams */
    
    CreateRealDFTParams( &status, &(transformparams.dftParams), &winParams, 1);
    TestStatus (&status, CODES(0), 1);

    CreateRealDFTParams( &status, &(transformparams.dftParams), &winParams, 1);
    TestStatus (&status, CODES(TFTRANSFORM_EALLOCP), 1);

    DestroyRealDFTParams( &status, NULL);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    {
      RealFFTPlan *p;
      p = transformparams.dftParams->plan;
      transformparams.dftParams->plan = NULL;
      DestroyRealDFTParams( &status, &(transformparams.dftParams));
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams.dftParams->plan = p;
    }

    {
      REAL4Vector *p;
      p = transformparams.dftParams->window;
      transformparams.dftParams->window = NULL;
      DestroyRealDFTParams( &status, &(transformparams.dftParams));
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams.dftParams->window = p;
    }

    DestroyRealDFTParams( &status, &(transformparams.dftParams));
    TestStatus (&status, CODES(0), 1);
  }


  /* 
   *
   *  Test functions CreateComplexDFTParams() and DestroyComplexDFTParams()
   *
   */
  {
    LALWindowParams winParams;

    if (verbose)
      {
	printf ("\n--- Testing CreateComplexDFTParams() and DestroyComplexDFTParams() \n\n");
      }

    winParams.type=Rectangular;
    winParams.length=tseglength;

    CreateComplexDFTParams( &status, NULL, &winParams, 1); 
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    winParams.length=0;
    CreateComplexDFTParams( &status, &(transformparams1.dftParams), &winParams, 1);  
    TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
    winParams.length=tseglength;
    
    transformparams.dftParams=NULL;
    DestroyComplexDFTParams( &status, &(transformparams1.dftParams)); 
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    /* next few tests require a valid DFTParams */
    
    CreateComplexDFTParams( &status, &(transformparams1.dftParams), &winParams, 1);
    TestStatus (&status, CODES(0), 1);

    CreateComplexDFTParams( &status, &(transformparams1.dftParams), &winParams, 1);
    TestStatus (&status, CODES(TFTRANSFORM_EALLOCP), 1);


    DestroyComplexDFTParams( &status, NULL);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);


    {
      ComplexFFTPlan *p;
      p = transformparams1.dftParams->plan;
      transformparams1.dftParams->plan = NULL;
      DestroyComplexDFTParams( &status, &(transformparams1.dftParams));
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams1.dftParams->plan = p;
    }

    {
      REAL4Vector *p;
      p = transformparams1.dftParams->window;
      transformparams1.dftParams->window = NULL;
      DestroyComplexDFTParams( &status, &(transformparams1.dftParams));
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams1.dftParams->window = p;
    }

    DestroyComplexDFTParams( &status, &(transformparams1.dftParams));
    TestStatus (&status, CODES(0), 1);

  }




  /* 
   *
   *  Test function ComputeFrequencySeries() 
   *
   */
  {
    LALWindowParams winParams;
    winParams.type=Rectangular;
    winParams.length=ntot;

    if (verbose)
      {
	printf("\n--- Testing ComputeFrequencySeries()\n\n");
      }

    CCreateVector( &status, &(fseries.data), ntot/2+1);
    TestStatus (&status, CODES(0), 1);
      
    SCreateVector (&status, &(tseries.data), ntot);
    TestStatus (&status, CODES(0), 1);

    CreateRealDFTParams( &status, &dftparams1, &winParams, 1);
    TestStatus (&status, CODES(0), 1);
    
    ComputeFrequencySeries ( &status, &fseries, &tseries, NULL);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    ComputeFrequencySeries ( &status, &fseries, NULL, dftparams1);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    ComputeFrequencySeries ( &status, NULL, &tseries, dftparams1);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    {
      COMPLEX8Vector *p;
      p = fseries.data;
      fseries.data=NULL;
      ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      fseries.data=p;
    }

    {
      REAL4Vector *p;
      p = tseries.data;
      tseries.data=NULL;
      ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tseries.data=p;

      p = dftparams1->window;
      dftparams1->window=NULL;
      ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      dftparams1->window=p;
    }

    {
      RealFFTPlan *p;
      p = dftparams1->plan;
      dftparams1->plan=NULL;
      ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      dftparams1->plan=p;
    }

    {
      INT4 n;
      n = dftparams1->plan->size;
      dftparams1->plan->size=0;
      ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      dftparams1->plan->size=n;
    }

    tseries.data->length--;
    ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    tseries.data->length++;
    
    fseries.data->length--;
    ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    fseries.data->length++;

    dftparams1->window->length--;
    ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    dftparams1->window->length++;

    {
      REAL8 p;
      p = tseries.deltaT;
      tseries.deltaT=0.0;
      ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      tseries.deltaT=p;

      p = dftparams1->sumofsquares;
      dftparams1->sumofsquares=0;
      ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      dftparams1->sumofsquares=p;
    }

    DestroyVector (&status, &(tseries.data));
    TestStatus (&status, CODES(0), 1);

    DestroyRealDFTParams( &status, &dftparams1);
    TestStatus (&status, CODES(0), 1);

    CDestroyVector( &status, &(fseries.data));
    TestStatus (&status, CODES(0), 1);
  }



  /* 
   *
   *  Test functions CreateTFPlane() and DestroyTFPlane()
   *
   */
  {
    if (verbose)
      {
	printf("\n--- Testing CreateTFPlane() and DestroyTFPlane()\n\n");
      }

    CreateTFPlane( &status, &tfp, NULL);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    CreateTFPlane( &status, NULL, &params);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    DestroyTFPlane( &status, NULL);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    DestroyTFPlane( &status, &tfp);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    {
      INT4 n;
      n = params.timeBins;
      params.timeBins=0;
      CreateTFPlane( &status, &tfp, &params);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      params.timeBins=n;

      n = params.freqBins;
      params.freqBins=0;
      CreateTFPlane( &status, &tfp, &params);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      params.freqBins=n;
    }

    {
      REAL8 p;
      p = params.deltaT;
      params.deltaT = 0.0;
      CreateTFPlane( &status, &tfp, &params);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      params.deltaT=p;
    }

      
  
    /* next few tests require a valid TF plane */
  
    CreateTFPlane( &status, &tfp, &params);
    TestStatus (&status, CODES(0), 1);

    CreateTFPlane( &status, &tfp, &params);
    TestStatus (&status, CODES(TFTRANSFORM_EALLOCP), 1);
    
    {
      COMPLEX8 *p;
      p = tfp->data;
      tfp->data = NULL;
      DestroyTFPlane( &status, &tfp);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tfp->data = p;
    }

    {
      TFPlaneParams *p;
      p = tfp->params;
      tfp->params=NULL;
      DestroyTFPlane( &status, &tfp);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tfp->params = p;
    }
    
    DestroyTFPlane( &status, &tfp);
    TestStatus (&status, CODES(0), 1);
  }






  /* 
   *
   *  Test function TimeSeriesToTFPlane() 
   *
   */
  {
    LALWindowParams winParams;
    winParams.type=Rectangular;
    winParams.length=tseglength;

    if (verbose)
      {
	printf("\n--- Testing TimeSeriesToTFPlane()\n\n");
      }

    CreateRealDFTParams( &status, &(transformparams.dftParams), &winParams, 1);    TestStatus (&status, CODES(0), 1);

    SCreateVector (&status, &(tseries.data), ntot);
    TestStatus (&status, CODES(0), 1);

    for(i=0; i< tseries.data->length; i++)
      {
	tseries.data->data[i] = ff( (REAL4)(i)/ (REAL4)(ntot));
      };

    CreateTFPlane( &status, &tfp, &params);
    TestStatus (&status, CODES(0), 1);


    /* Now start checking errors */

    TimeSeriesToTFPlane( &status, tfp, &tseries, NULL);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
    TimeSeriesToTFPlane( &status, tfp, NULL, &transformparams);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
    TimeSeriesToTFPlane( &status, NULL, &tseries, &transformparams);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    {
      REAL4Vector *p;
      p = tseries.data;
      tseries.data=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tseries.data=p;

      p = transformparams.dftParams->window;
      transformparams.dftParams->window=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams.dftParams->window=p;

    }

    {
      REAL4 *p;
      p = tseries.data->data;
      tseries.data->data=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tseries.data->data=p;

      p = transformparams.dftParams->window->data;
      transformparams.dftParams->window->data=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams.dftParams->window->data=p;
    }

    {
      COMPLEX8 *p;
      p = tfp->data;
      tfp->data=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tfp->data=p;
    }

    {
      RealDFTParams *p;
      p = transformparams.dftParams;
      transformparams.dftParams=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams.dftParams=p;
    }

    {
      RealFFTPlan *p;
      p = transformparams.dftParams->plan;
      transformparams.dftParams->plan=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams.dftParams->plan=p;
    }

    {
      void *p;
      p = transformparams.dftParams->plan->plan;
      transformparams.dftParams->plan->plan=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams.dftParams->plan->plan=p;
    }

    {
      TFPlaneParams *p;
      p = tfp->params;
      tfp->params=NULL;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tfp->params=p;
    }

    {
      INT4 p;
      p = tfp->params->timeBins;
      tfp->params->timeBins=0;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      tfp->params->timeBins=p;

      p = tfp->params->freqBins;
      tfp->params->freqBins=0;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);

      tfp->params->freqBins=10000000;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
      tfp->params->freqBins=p;

      p = transformparams.startT;
      transformparams.startT = 1000000;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
      transformparams.startT = p;
    }

    {
      REAL4 p;
      p = transformparams.dftParams->sumofsquares;
      transformparams.dftParams->sumofsquares=0;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      transformparams.dftParams->sumofsquares=p;
    }

    {
      REAL8 p;
      p = tseries.deltaT;
      tseries.deltaT=0.0;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      tseries.deltaT=p;

      p = tfp->params->deltaT;
      tfp->params->deltaT=0.0;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      tfp->params->deltaT = 0.1 * tseries.deltaT;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
      tfp->params->deltaT=p;

      p = tseries.f0;
      tseries.f0 = -1.0;
      TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      tseries.f0=p;
    }

    transformparams.dftParams->plan->size--;
    TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    transformparams.dftParams->plan->size++;

    transformparams.dftParams->window->length--;
    TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    transformparams.dftParams->window->length++;

    transformparams.dftParams->plan->sign=-1;
    TimeSeriesToTFPlane( &status, tfp, &tseries, &transformparams);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    transformparams.dftParams->plan->sign=1;

    
    /* clean up */

    DestroyTFPlane( &status, &tfp);
    TestStatus (&status, CODES(0), 1);

    SDestroyVector (&status, &(tseries.data) );
    TestStatus (&status, CODES(0), 1);

    DestroyRealDFTParams( &status, &(transformparams.dftParams));
    TestStatus (&status, CODES(0), 1);
  }






  /* 
   *
   *  Test function FreqSeriesToTFPlane() 
   *
   */
  {
    LALWindowParams winParams;
    LALWindowParams winParams1; 

    winParams.type=Rectangular;
    winParams.length=ntot;
    winParams1.type=Rectangular;
    winParams1.length=fseglength;

    if (verbose)
      {
	printf("\n--- Testing FreqSeriesToTFPlane()\n\n");
      }


    SCreateVector (&status, &(tseries.data), ntot);
    TestStatus (&status, CODES(0), 1);

    for(i=0; i< tseries.data->length; i++)
      {
	tseries.data->data[i] = ff( (REAL4)(i)/ (REAL4)(ntot));
      };

    CreateTFPlane( &status, &tfp, &params);
    TestStatus (&status, CODES(0), 1);

    CCreateVector( &status, &(fseries.data), ntot/2+1);
    TestStatus (&status, CODES(0), 1);
      

    CreateRealDFTParams( &status, &dftparams1, &winParams,1);
    TestStatus (&status, CODES(0), 1);
    
    ComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
    TestStatus (&status, CODES(0), 1);

    CreateComplexDFTParams( &status, &(transformparams1.dftParams), 
                            &winParams1, -1); 
    TestStatus (&status, CODES(0), 1);

    FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
    TestStatus (&status, CODES(0), 1);



    /* Now start checking errors */

    FreqSeriesToTFPlane( &status, tfp, &fseries, NULL);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
    FreqSeriesToTFPlane( &status, tfp, NULL, &transformparams1);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
    FreqSeriesToTFPlane( &status, NULL, &fseries, &transformparams1);
    TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);

    {
      COMPLEX8Vector *p;
      p = fseries.data;
      fseries.data=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      fseries.data=p;
    }
    
    {
      REAL4Vector *p;
      p = transformparams1.dftParams->window;
      transformparams1.dftParams->window=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams1.dftParams->window=p;
    }

    {
      COMPLEX8 *p;
      p = fseries.data->data;
      fseries.data->data=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      fseries.data->data=p;

      p = tfp->data;
      tfp->data=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tfp->data=p;
    }

    {
      REAL4 *p;
      p = transformparams1.dftParams->window->data;
      transformparams1.dftParams->window->data=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams1.dftParams->window->data=p;
    }

    {
      ComplexDFTParams *p;
      p = transformparams1.dftParams;
      transformparams1.dftParams=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams1.dftParams=p;
    }


    {
      ComplexFFTPlan *p;
      p = transformparams1.dftParams->plan;
      transformparams1.dftParams->plan=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams1.dftParams->plan=p;
    }

    {
      void *p;
      p = transformparams1.dftParams->plan->plan;
      transformparams1.dftParams->plan->plan=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      transformparams1.dftParams->plan->plan=p;
    }

    {
      TFPlaneParams *p;
      p = tfp->params;
      tfp->params=NULL;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_ENULLP), 1);
      tfp->params=p;
    }

    {
      INT4 p;
      p = tfp->params->timeBins;
      tfp->params->timeBins=0;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      tfp->params->timeBins=p;

      p = tfp->params->freqBins;
      tfp->params->freqBins=0;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);

      tfp->params->freqBins=10000000;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
      tfp->params->freqBins=p;
    }

    {
      REAL4 p;
      p = transformparams1.dftParams->sumofsquares;
      transformparams1.dftParams->sumofsquares=0;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      transformparams1.dftParams->sumofsquares=p;
    }

    {
      REAL8 p;
      p = fseries.deltaF;
      fseries.deltaF=0.0;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      fseries.deltaF=p;

      p = tfp->params->deltaT;
      tfp->params->deltaT=0.0;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      tfp->params->deltaT = 20.0 / fseries.deltaF;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
      tfp->params->deltaT=p;

      p = fseries.f0;
      fseries.f0 = -1.0;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EPOSARG), 1);
      fseries.f0 = 10000000.0;
      FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
      TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
      fseries.f0=p;
    }

    transformparams1.dftParams->plan->size--;
    FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    transformparams1.dftParams->plan->size++;

    transformparams1.dftParams->window->length--;
    FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    transformparams1.dftParams->window->length++;

    transformparams1.dftParams->plan->sign=1;
    FreqSeriesToTFPlane( &status, tfp, &fseries, &transformparams1);
    TestStatus (&status, CODES(TFTRANSFORM_EINCOMP), 1);
    transformparams1.dftParams->plan->sign=-1;




    /* Now clean up */

    DestroyComplexDFTParams( &status, &(transformparams1.dftParams));
    TestStatus (&status, CODES(0), 1);
    
    DestroyRealDFTParams( &status, &dftparams1);
    TestStatus (&status, CODES(0), 1);

    CDestroyVector( &status, &(fseries.data));
    TestStatus (&status, CODES(0), 1);

    DestroyTFPlane( &status, &tfp);
    TestStatus (&status, CODES(0), 1);

    SDestroyVector (&status, &(tseries.data) );
    TestStatus (&status, CODES(0), 1);

  }


  LALCheckMemoryLeaks ();

  if(verbose)  printf("PASS: all tests\n");
  
  return 0;
}


static REAL4 ff(REAL4 w)
{
  /* simple waveform function used for testing */
  REAL4 t,s,f;
  REAL4 sigma = 0.4;
  t = 2.0* (w - 0.5);
  f = 70.0 + 30.0*t;
  s = sin(f*t)*exp(-t*t/( 2.0 * sigma * sigma));
  return(s);
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
TestStatus (Status *status, const char *ignored, int exitcode)
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    /*REPORTSTATUS (status);*/
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
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h         print this message\n");
  fprintf (stderr, "  -q         quiet: run silently\n");
  fprintf (stderr, "  -v         verbose: print extra information\n");
  fprintf (stderr, "  -d level   set debuglevel to level\n");
  exit (exitcode);
}


/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debug level */
        debuglevel = atoi (optarg);
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
        break;

      case 'h':
        Usage (argv[0], 0);
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





