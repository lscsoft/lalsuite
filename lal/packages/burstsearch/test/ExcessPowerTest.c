/*----------------------------------------------------------------------- 
 * 
 * File Name: ExcessPowerTest.c
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
 * Test suite for functions in ExcessPower.c
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

#include <lal/ExcessPower.h>
#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/Thresholds.h>
#include <lal/VectorOps.h>


#define _CODES(x) #x
#define CODES(x) _CODES(x)


NRCSID (MAIN, "$Id$");


extern char *optarg;
extern int   optind;

INT4 lalDebugLevel = 2;   /* set to 2 to get full status information for tests */
INT4 verbose    = 0;


static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);


static REAL4 ff(REAL4 w);   /* simple function used to construct a waveform */

int
main (int argc, char *argv[])
{
  static LALStatus             status; 
  
  const INT4                ntot=512;
  const REAL8               factor=1.1;        /* controls SNR of signal */

  TFTiling                  *tfTiling=NULL;
  COMPLEX8FrequencySeries   fseries;
  REAL4TimeSeries           tseries;
  ComputeExcessPowerIn      input;
  INT4                      ntotT;
  INT4                      ntotF;
  INT4                      i;
  INT4                      numEvents;
  REAL8                     lambda;
  REAL8                     snr2;

  


  /*
   *
   * Parse the command line options
   *
   */


  ParseOptions (argc, argv);


  /* 
   *  
   *  Set up input time series
   *
   */

  {
    ntotT = 2*ntot +2;
   
    tseries.epoch.gpsSeconds=0;
    tseries.epoch.gpsNanoSeconds=0;
    tseries.deltaT = 0.001;  /* 1 kHz sampling */
    tseries.f0 = 0.0;
    /*
     * OMITTED
     *
    tseries.name = NULL;
    tseries.sampleUnits=NULL;
     */
    tseries.data=NULL;

    LALSCreateVector (&status, &(tseries.data), ntotT);
    TestStatus (&status, CODES(0), 1);
    
    for(i=0; i< (INT4)tseries.data->length; i++)
      {
	tseries.data->data[i] = factor * ff( (REAL4)(i)/ (REAL4)(ntotT));
      };

  }



  /*
   *
   *  Take DFT of time series to give frequency series
   *
   */

#define PATRICK 0 
#ifdef PATRICK
  {
    RealDFTParams                   *dftparams1=NULL;
    
    ntotF = ntotT/2+1;
    fseries.data=NULL;
    LALCCreateVector( &status, &(fseries.data), ntotF);
    TestStatus (&status, CODES(0), 1);
      
    {
      LALWindowParams winParams;
      winParams.type=Rectangular;
      winParams.length=ntotT;
      LALCreateRealDFTParams( &status, &dftparams1, &winParams, 1);
      TestStatus (&status, CODES(0), 1);
    }

    if(verbose)
      {
	printf("Computing FFT of time domain data\n");
      }
    LALComputeFrequencySeries ( &status, &fseries, &tseries, dftparams1);
    TestStatus (&status, CODES(0), 1);

    LALDestroyRealDFTParams (&status, &dftparams1);
    TestStatus (&status, CODES(0), 1);
        

  }
#endif


  /*
   *
   *   Compute matched filtering SNR
   *
   */

  {
    
    snr2=0.0;
    for(i=0;i<(INT4)fseries.data->length;i++)
      {
	COMPLEX8 z=fseries.data->data[i];
	snr2 += z.re*z.re + z.im*z.im;
      }
  }
    


  /*
   *
   *  Add white noise to frequency series
   *
   */
  {

    LALAddWhiteNoise (&status, fseries.data, 1.0);
    TestStatus (&status, CODES(0), 1);

    /*  
     *  1st and last point must be purely real to be able to FFT back
     *  to time domain, so ...
     */
    fseries.data->data[0].im=0.0;
    fseries.data->data[ntotF-1].im=0.0;
    
  }



  /*
   *
   *  DFT back to time domain to get modified data
   * 
   */
  {
    RealFFTPlan            *pinv = NULL;

    LALCreateReverseRealFFTPlan (&status, &pinv, ntotT, 0);
    TestStatus (&status, CODES(0), 1);

    LALReverseRealFFT (&status, tseries.data, fseries.data, pinv);
    TestStatus (&status, CODES(0), 1);

    LALDestroyRealFFTPlan     (&status, &pinv);
    TestStatus (&status, CODES(0), 1);

    /* correct normalization to match original time domain data */
    for(i=0; i<(INT4)tseries.data->length; i++)
      {
	tseries.data->data[i] /= sqrt((REAL8)(ntotT));
      }



  }    


  /*
   *
   *  set up TFTiling
   *
   */

  {
    CreateTFTilingIn          params;
    
    params.overlapFactor = 3;
    params.minFreqBins = 1;
    params.minTimeBins = 1;
    params.flow = 0.0;
    params.deltaF = fseries.deltaF;
    params.length = ntot;
    params.maxTileBand = 64.0;
  
    LALCreateTFTiling (&status, &tfTiling, &params);
    TestStatus (&status, CODES(0), 1);
  }





  /*
   *
   *  Now do the computations
   *
   */

  input.numSigmaMin=2.0;
  input.alphaDefault=0.5;

  LALComputeTFPlanes (&status, tfTiling, &fseries);
  TestStatus (&status, CODES(0), 1);

  LALComputeExcessPower (&status, tfTiling, &input);
  TestStatus (&status, CODES(0), 1);

  LALComputeLikelihood  (&status, &lambda, tfTiling);
  TestStatus (&status, CODES(0), 1);

  LALSortTFTiling (&status, tfTiling);
  TestStatus (&status, CODES(0), 1);

  {
	TFTile *tile;
	numEvents= 0;
	for(tile = tfTiling->firstTile; tile && (tile->alpha <= 0.000000001); tile = tile->nextTile)
		numEvents++;
  }


  if(verbose)
    {
    REAL8 sigma2;

    /* 
     *  if lambda<1 computation of effective number of sigma fails
     *  so trap this
     */

    if(lambda<=1.0) 
      {
	sigma2 = 0.0;
      }
    else
      {
	Chi2ThresholdIn locinput;

	locinput.dof = 1.0;
	locinput.falseAlarm = 1.0 / lambda;
	
	LALChi2Threshold (&status, &sigma2, &locinput);
	TestStatus (&status, CODES(0), 1);
      }

    LALPrintTFTileList (&status, stdout, tfTiling, 4);
    TestStatus (&status, CODES(0), 1);
    printf("\n");
    printf("  ***    Matched filtering SNR           : %f\n",sqrt(snr2));
    printf("  ***    Number of tiles above threshold : %d out of %d\n",
           numEvents,tfTiling->numTiles);
    printf("  ***    Average likelihood              : %f\n",lambda);
    printf("  ***    Effective number of sigma       : %f\n\n",sqrt(sigma2));
    
    }

  /*
   *
   *  Clean up memory used
   *
   */


  LALCDestroyVector (&status, &(fseries.data) );
  TestStatus (&status, CODES(0), 1);

  LALDestroyTFTiling (&status, &tfTiling);
  TestStatus (&status, CODES(0), 1);

  LALDestroyVector (&status, &(tseries.data) );
  TestStatus (&status, CODES(0), 1);






  /*************************************************************************
   *                                                                       * 
   *                                                                       *
   *  Now check to make sure that correct error codes are generated.       *
   *                                                                       *
   *                                                                       * 
   *************************************************************************/

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
  if (verbose || lalDebugLevel)
  {
    printf ("\n===== Check Errors =====\n");
  }

  /* 
   *
   *  Test function LALAddWhiteNoise()
   *
   */
  {
    COMPLEX8Vector *v=NULL;
    INT4 ii;
    COMPLEX8 *p=NULL;

    if (verbose)
      {
	printf ("\n--- Testing LALAddWhiteNoise() \n\n");
      }
    
    LALCCreateVector( &status, &v, 10);
    TestStatus (&status, CODES(0), 1);
    for(ii=0; ii<10; ii++)
      {
	v->data[ii].re=0.0;
	v->data[ii].im=0.0;
      }

    LALAddWhiteNoise( &status, NULL, 1.0);
    TestStatus( &status, CODES(LAL_NULL_ERR), 1);

    p = v->data;
    v->data=NULL;
    LALAddWhiteNoise( &status, v, 1.0);
    TestStatus( &status, CODES(LAL_NULL_ERR), 1);
    v->data=p;

    ii=v->length;
    v->length=0;
    LALAddWhiteNoise( &status, v, 1.0);
    TestStatus( &status, CODES(LAL_RANGE_ERR), 1);
    v->length=ii;

    LALCDestroyVector (&status, &v);
    TestStatus (&status, CODES(0), 1);
  }



  /* 
   *
   *  Test functions LALCreateTFTiling() and LALDestroyTFTiling()
   *
   */
  {
    CreateTFTilingIn          params;

    if (verbose)
      {
	printf ("\n--- Testing LALCreateTFTiling() and LALDestroyTFTiling() \n\n");
      }
    
    params.overlapFactor = 3;
    params.minFreqBins = 1;
    params.minTimeBins = 1;
    params.flow = 0.0;
    params.deltaF = 1.0;
    params.length = 16;
    params.maxTileBand = 64.0;
  
    LALCreateTFTiling (&status, &tfTiling, NULL);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    LALCreateTFTiling (&status, NULL, &params);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    LALDestroyTFTiling (&status, &tfTiling);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);


    {
      INT4 p;
      
      p = params.overlapFactor;
      params.overlapFactor=0;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.overlapFactor=p;

      p = params.length;
      params.length=0;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.length=17;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_BADPARM_ERR), 1);
      params.length=p;

      p = params.minFreqBins;
      params.minFreqBins=0;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.minFreqBins=100;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.minFreqBins=p;

      p = params.minTimeBins;
      params.minTimeBins=0;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.minTimeBins=100;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.minTimeBins=p;
    }
      

    {
      REAL8 p;
      
      p=params.deltaF;
      params.deltaF=0.0;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.deltaF=p;
      
      p=params.flow;
      params.flow=-1.0;
      LALCreateTFTiling (&status, &tfTiling, &params);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      params.flow=p;
    }

    LALCreateTFTiling (&status, &tfTiling, &params);
    TestStatus (&status, CODES(0), 1);
 
    LALCreateTFTiling (&status, &tfTiling, &params);
    TestStatus (&status, CODES(LAL_NNULL_ERR), 1);

    LALDestroyTFTiling (&status, NULL);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    {
      COMPLEX8TimeFrequencyPlane **p;
      p = tfTiling->tfp;
      tfTiling->tfp=NULL;
      LALDestroyTFTiling (&status, &tfTiling);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->tfp=p;
    }

    {
      TFTile *p;
      p = tfTiling->firstTile;
      tfTiling->firstTile=NULL;
      LALDestroyTFTiling (&status, &tfTiling);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->firstTile=p;
    }

    {
      INT4 p;
      p = tfTiling->numPlanes;
      tfTiling->numPlanes=0;
      LALDestroyTFTiling (&status, &tfTiling);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      tfTiling->numPlanes=p;
    }

    LALDestroyTFTiling (&status, &tfTiling);
    TestStatus (&status, CODES(0), 1);
  }




  /* 
   *
   *  Test functions LALComputeTFPlanes()
   *
   */
  {
    CreateTFTilingIn          params;
    COMPLEX8FrequencySeries   locfseries;

    if (verbose)
      {
	printf ("\n--- Testing LALComputeTFPlanes() \n\n");
      }
    
    params.overlapFactor = 3;
    params.minFreqBins = 1;
    params.minTimeBins = 1;
    params.flow = 0.0;
    params.deltaF = 1.0;
    params.length = 16;
    params.maxTileBand = 64.0;
  
    LALCreateTFTiling (&status, &tfTiling, &params);
    TestStatus (&status, CODES(0), 1);

    locfseries.epoch.gpsSeconds=0;
    locfseries.epoch.gpsNanoSeconds=0;
    locfseries.deltaF = 1.0;
    locfseries.f0 = 0.0;
    /*
     * OMITTED
     *
    locfseries.name = NULL;
    locfseries.sampleUnits=NULL;
     */
    locfseries.data=NULL;

    LALCCreateVector( &status, &(locfseries.data), 1000);
    TestStatus (&status, CODES(0), 1);
    

    /* now start checking errors */

    LALComputeTFPlanes( &status, tfTiling, NULL);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    LALComputeTFPlanes( &status, NULL, &locfseries);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    {
      COMPLEX8Vector *p;
      p = locfseries.data;
      locfseries.data=NULL;
      LALComputeTFPlanes( &status, tfTiling, &locfseries);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      locfseries.data=p;
    }

    {
      COMPLEX8 *p;
      p = locfseries.data->data;
      locfseries.data->data=NULL;
      LALComputeTFPlanes( &status, tfTiling, &locfseries);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      locfseries.data->data=p;
    }

    {
      TFTile *p;
      p=tfTiling->firstTile;
      tfTiling->firstTile=NULL;
      LALComputeTFPlanes( &status, tfTiling, &locfseries);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->firstTile=p;
    }

    {
      INT4 p;
      p = tfTiling->numPlanes;
      tfTiling->numPlanes=0;
      LALComputeTFPlanes( &status, tfTiling, &locfseries);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      tfTiling->numPlanes=p;
    }

    {
      COMPLEX8TimeFrequencyPlane **p;
      p = tfTiling->tfp;
      tfTiling->tfp=NULL;
      LALComputeTFPlanes( &status, tfTiling, &locfseries);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->tfp=p;
    }

    {
      ComplexDFTParams **p;
      p = tfTiling->dftParams;
      tfTiling->dftParams=NULL;
      LALComputeTFPlanes( &status, tfTiling, &locfseries);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->dftParams=p;
    }

    {
      INT4 ii;
      COMPLEX8TimeFrequencyPlane **thisPlane;
      ComplexDFTParams      **thisdftParams;
      COMPLEX8TimeFrequencyPlane *p;
      ComplexDFTParams *p1;
            
      for(ii=0;ii<tfTiling->numPlanes; ii++)
	{
	  thisPlane = tfTiling->tfp+ii;
	  p=*thisPlane;
	  *thisPlane=NULL;
	  LALComputeTFPlanes( &status, tfTiling, &locfseries);
	  TestStatus (&status, CODES(LAL_NULL_ERR), 1);
	  *thisPlane=p;

	  thisdftParams = tfTiling->dftParams+ii;
	  p1=*thisdftParams;
	  *thisdftParams=NULL;
	  LALComputeTFPlanes( &status, tfTiling, &locfseries);
	  TestStatus (&status, CODES(LAL_NULL_ERR), 1);
	  *thisdftParams=p1;
	}
    }

    /* clean up */

    LALCDestroyVector( &status, &(locfseries.data));
    TestStatus (&status, CODES(0), 1);

    LALDestroyTFTiling (&status, &tfTiling);
    TestStatus (&status, CODES(0), 1);
  }







  /* 
   *
   *  Test function LALComputeExcessPower()
   *
   */
  {
    CreateTFTilingIn          params;
    COMPLEX8FrequencySeries   locfseries;

    if (verbose)
      {
	printf ("\n--- Testing LALComputeExcessPower() \n\n");
      }
    
    params.overlapFactor = 3;
    params.minFreqBins = 1;
    params.minTimeBins = 1;
    params.flow = 0.0;
    params.deltaF = 1.0;
    params.length = 16;
    params.maxTileBand = 64.0;
  
    LALCreateTFTiling (&status, &tfTiling, &params);
    TestStatus (&status, CODES(0), 1);

    locfseries.epoch.gpsSeconds=0;
    locfseries.epoch.gpsNanoSeconds=0;
    locfseries.deltaF = 1.0;
    locfseries.f0 = 0.0;
    /*
     * OMITTED
     *
    locfseries.name = NULL;
    locfseries.sampleUnits=NULL;
     */
    locfseries.data=NULL;

    LALCCreateVector( &status, &(locfseries.data), 1000);
    TestStatus (&status, CODES(0), 1);
    
    LALComputeExcessPower( &status, tfTiling, &input);	
    TestStatus (&status, CODES(EXCESSPOWERH_EORDER), 1);

    LALComputeTFPlanes( &status, tfTiling, &locfseries);
    TestStatus (&status, CODES(0), 1);



    /* now start checking errors */

    LALComputeExcessPower( &status, NULL, &input);	
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    LALComputeExcessPower( &status, tfTiling, NULL);	
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);


    {
      COMPLEX8TimeFrequencyPlane **p;
      p = tfTiling->tfp;
      tfTiling->tfp=NULL;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->tfp=p;
    }

    {
      ComplexDFTParams **p;
      p = tfTiling->dftParams;
      tfTiling->dftParams=NULL;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->dftParams=p;
    }

    {
      TFTile *p;
      p=tfTiling->firstTile;
      tfTiling->firstTile=NULL;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->firstTile=p;
    }

    {
      INT4 p;
      p = tfTiling->numPlanes;
      tfTiling->numPlanes=0;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      tfTiling->numPlanes=p;
    }

    {
      INT4 ii;
      COMPLEX8TimeFrequencyPlane **thisPlane;
      COMPLEX8TimeFrequencyPlane *p;
      COMPLEX8 *p2;
      TFPlaneParams *p3;
            
      for(ii=0;ii<tfTiling->numPlanes; ii++)
	{
	  thisPlane = tfTiling->tfp+ii;
	  p=*thisPlane;
	  *thisPlane=NULL;
	  LALComputeExcessPower( &status, tfTiling, &input);
	  TestStatus (&status, CODES(LAL_NULL_ERR), 1);
	  *thisPlane=p;
	  
	  p2=(*thisPlane)->data;
	  (*thisPlane)->data=NULL;
	  LALComputeExcessPower( &status, tfTiling, &input);
	  TestStatus (&status, CODES(LAL_NULL_ERR), 1);
	  (*thisPlane)->data=p2;

	  p3=(*thisPlane)->params;
	  (*thisPlane)->params=NULL;
	  LALComputeExcessPower( &status, tfTiling, &input);
	  TestStatus (&status, CODES(LAL_NULL_ERR), 1);
	  (*thisPlane)->params=p3;
	}
    }

    {
      INT4 p;
      TFTile *t=tfTiling->firstTile;
      COMPLEX8TimeFrequencyPlane *tfPlane=*(tfTiling->tfp);

      p = t->whichPlane;
      t->whichPlane=-1;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      t->whichPlane = tfTiling->numPlanes;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      t->whichPlane=p;
      
      p = tfPlane->params->timeBins;
      tfPlane->params->timeBins=0;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      tfPlane->params->timeBins=p;

      p = tfPlane->params->freqBins;
      tfPlane->params->freqBins=0;
      LALComputeExcessPower( &status, tfTiling, &input);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      tfPlane->params->freqBins=p;
    }
    
    /* clean up */

    LALCDestroyVector( &status, &(locfseries.data));
    TestStatus (&status, CODES(0), 1);

    LALDestroyTFTiling (&status, &tfTiling);
    TestStatus (&status, CODES(0), 1);
  }







  /* 
   *
   *  Test function LALComputeLikelihood()
   *
   */
  {
    CreateTFTilingIn          params;
    COMPLEX8FrequencySeries   locfseries;
    REAL8                     loclambda;

    if (verbose)
      {
	printf ("\n--- Testing LALComputeLikelihood() \n\n");
      }
    
    params.overlapFactor = 3;
    params.minFreqBins = 1;
    params.minTimeBins = 1;
    params.flow = 0.0;
    params.deltaF = 1.0;
    params.length = 16;
    params.maxTileBand = 64.0;
  
    LALCreateTFTiling (&status, &tfTiling, &params);
    TestStatus (&status, CODES(0), 1);

    locfseries.epoch.gpsSeconds=0;
    locfseries.epoch.gpsNanoSeconds=0;
    locfseries.deltaF = 1.0;
    locfseries.f0 = 0.0;
    /*
     * OMITTED
     *
    locfseries.name = NULL;
    locfseries.sampleUnits=NULL;
     */
    locfseries.data=NULL;

    LALCCreateVector( &status, &(locfseries.data), 1000);
    TestStatus (&status, CODES(0), 1);
    memset( locfseries.data->data, 0,
        locfseries.data->length * sizeof( *locfseries.data->data ) );
    
    LALComputeTFPlanes( &status, tfTiling, &locfseries);
    TestStatus (&status, CODES(0), 1);

    LALComputeLikelihood( &status, &loclambda, tfTiling);	
    TestStatus (&status, CODES(EXCESSPOWERH_EORDER), 1);

    LALComputeExcessPower( &status, tfTiling, &input);
    TestStatus (&status, CODES(0), 1);
      

    /* now start checking errors */

    LALComputeLikelihood( &status, NULL, tfTiling);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    LALComputeLikelihood( &status, &loclambda, NULL);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    {
      TFTile *p;
      p=tfTiling->firstTile;
      tfTiling->firstTile=NULL;
      LALComputeLikelihood( &status, &loclambda, tfTiling);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->firstTile=p;
    }

    /* clean up */

    LALCDestroyVector( &status, &(locfseries.data));
    TestStatus (&status, CODES(0), 1);

    LALDestroyTFTiling (&status, &tfTiling);
    TestStatus (&status, CODES(0), 1);
  }



  /* 
   *
   *  Test functions LALSortTFTiling() and LALPrintTFTileList()
   *
   */
  {
    CreateTFTilingIn          params;
    COMPLEX8FrequencySeries   locfseries;
    INT4                      ii;

    if (verbose)
      {
	printf ("\n--- Testing LALSortTFTiling() \n\n");
      }
    
    params.overlapFactor = 3;
    params.minFreqBins = 1;
    params.minTimeBins = 1;
    params.flow = 0.0;
    params.deltaF = 1.0;
    params.length = 16;
    params.maxTileBand = 64.0;
  
    LALCreateTFTiling (&status, &tfTiling, &params);
    TestStatus (&status, CODES(0), 1);

    locfseries.epoch.gpsSeconds=0;
    locfseries.epoch.gpsNanoSeconds=0;
    locfseries.deltaF = 1.0;
    locfseries.f0 = 0.0;
    /*
     * OMITTED
     *
    locfseries.name = NULL;
    locfseries.sampleUnits=NULL;
     */
    locfseries.data=NULL;

    LALCCreateVector( &status, &(locfseries.data), 1000);
    TestStatus (&status, CODES(0), 1);
    
    for(ii=0; ii<1000; ii++)
      {
	locfseries.data->data[ii].re = 1.1;
	locfseries.data->data[ii].im = 1.1;
      }


    LALComputeTFPlanes( &status, tfTiling, &locfseries);
    TestStatus (&status, CODES(0), 1);

    LALSortTFTiling( &status, tfTiling);
    TestStatus (&status, CODES(EXCESSPOWERH_EORDER), 1);

    LALPrintTFTileList( &status, stdout, tfTiling, 2);
    TestStatus (&status, CODES(EXCESSPOWERH_EORDER), 1);

    LALComputeExcessPower( &status, tfTiling, &input);
    TestStatus (&status, CODES(0), 1);

    

    /* now start checking errors */

    LALSortTFTiling( &status, NULL);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    {
      COMPLEX8TimeFrequencyPlane **p;
      p = tfTiling->tfp;
      tfTiling->tfp=NULL;
      LALSortTFTiling( &status, tfTiling);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->tfp=p;
    }

    {
      ComplexDFTParams **p;
      p = tfTiling->dftParams;
      tfTiling->dftParams=NULL;
      LALSortTFTiling( &status, tfTiling);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->dftParams=p;
    }

    {
      TFTile *p;
      p=tfTiling->firstTile;
      tfTiling->firstTile=NULL;
      LALSortTFTiling( &status, tfTiling);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->firstTile=p;
    }

    {
      INT4 p;
      p = tfTiling->numPlanes;
      tfTiling->numPlanes=0;
      LALSortTFTiling( &status, tfTiling);
      TestStatus (&status, CODES(LAL_RANGE_ERR), 1);
      tfTiling->numPlanes=p;
    }






    if (verbose)
      {
	printf ("\n--- Testing PrintTFTiling() \n\n");
      }

    LALPrintTFTileList( &status, NULL, tfTiling, 2);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);

    LALPrintTFTileList( &status, stdout, NULL, 2);
    TestStatus (&status, CODES(LAL_NULL_ERR), 1);
    
    {
      TFTile *p;
      p = tfTiling->firstTile;
      tfTiling->firstTile=NULL;
      LALPrintTFTileList( &status, stdout, tfTiling, 2);
      TestStatus (&status, CODES(LAL_NULL_ERR), 1);
      tfTiling->firstTile=p;
    }

    /* clean up */

    LALCDestroyVector( &status, &(locfseries.data));
    TestStatus (&status, CODES(0), 1);

    LALDestroyTFTiling (&status, &tfTiling);
    TestStatus (&status, CODES(0), 1);
  }


  }  /* end if laldebuglevel */
#endif



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
TestStatus (LALStatus *status, const char *ignored, int exitcode)
{
  char  str[64];
  char *tok;

  /*  if (verbose)
  {
    REPORTSTATUS (status);
    }*/

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
  fprintf (stderr, "  -d level   set lalDebugLevel to level\n");
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
        lalDebugLevel = atoi (optarg);
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











