#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include "FindChirpTDTest.h"

#define TEST_STATUS( ps ) \
  ( ( ps )->statusCode && ( exit( 1 ), \
    fprintf( stderr, "LAL Routine failed near line %d\n", __LINE__ ), 1 ) )

LALStatus status;
int lalDebugLevel = 1;

int main( void )
{
  const UINT4 numSegments  = 1;
  const UINT4 numPoints    = 262144;
  const UINT4 numChisqBins = 8;
  const UINT4 invSpecTrunc = 0;
  const REAL4 mass1        = 1.4;    /* solar masses */
  const REAL4 mass2        = 1.4;    /* solar masses */
  const REAL4 srate        = 16384;  /* Hz */
  const REAL4 fmin         = 100;    /* Hz */
  const REAL4 fmax         = 1000;   /* Hz */
  const REAL4 dynRange     = 1;

  FindChirpInitParams initParams; /* need to populate this by hand */

  /* these are created by Start() */
  FindChirpFilterInput   *filterInput   = NULL;
  FindChirpFilterParams  *filterParams  = NULL;
  FindChirpSegmentVector *fcSegVec      = NULL;
  DataSegmentVector      *dataSegVec    = NULL;

  /* these are required for SPFilter(); they are created by SPInit() */
  FindChirpTmpltParams *spTmpltParams = NULL;
  FindChirpDataParams  *spDataParams  = NULL;

  /* this is required for TDFilter(); it is created by TDInit() */
  FindChirpDataParams  *tdDataParams  = NULL;


  /* set initialization parameters */
  initParams.numSegments    = numSegments;
  initParams.numPoints      = numPoints;
  initParams.numChisqBins   = numChisqBins;
  initParams.approximant    = TaylorF2;
  initParams.createRhosqVec = 1;


  /* create generic objects needed by both SP and TD filters */
  Start( &dataSegVec, &filterInput, &filterParams, &fcSegVec, &initParams );


  /* set filter parameters, e.g., thresholds for events */
  filterParams->deltaT         = 1 / srate;
  filterParams->rhosqThresh    = 1e-6;
  filterParams->chisqThresh    = 1e+6;


  /* create some fake data */
  MakeData( dataSegVec, mass1, mass2, srate, fmin, fmax );


  /*
   * initialize the SP- and TD-specific parameters
   */

  SPInit( &spTmpltParams, &spDataParams, &initParams, srate, fmin, dynRange,
     invSpecTrunc );

  TDInit( &tdDataParams, &initParams, srate, fmin, dynRange, invSpecTrunc );


  /*
   * do the SP- and TD-filtering... could loop over templates and data here
   * without re-initializing
   */

  SPFilter( dataSegVec, mass1, mass2, filterInput, filterParams, fcSegVec,
      spTmpltParams, spDataParams );

  TDFilter( dataSegVec, mass1, mass2, fmax, filterInput, filterParams, fcSegVec,
      tdDataParams );


  /* print out the correct effective distance for the injected template */
  fprintf( stdout, "\nExpected:\n\td_Mpc = %e\n",
      ( mass1 + mass2 ) * LAL_MRSUN_SI / ( 1e6 * LAL_PC_SI ) );


  /* clean up memory and exit */

  SPFini( &spTmpltParams, &spDataParams );
  TDFini( &tdDataParams );
  Stop( &dataSegVec, &filterInput, &filterParams, &fcSegVec, numChisqBins );

  LALCheckMemoryLeaks();
  return 0;
}



/*
 *
 * Start(), Stop()
 *
 * Start() creates and initializes various structures needed by FindChirp.
 * Stop() destroys the allocated memory.
 *
 */

int Start(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    FindChirpInitParams     *initParams
    )
{
  LALCreateFindChirpSegmentVector( &status, fcSegVec, initParams );
  TEST_STATUS( &status );

  LALCreateFindChirpInput( &status, filterInput, initParams );
  TEST_STATUS( &status );

  LALFindChirpFilterInit( &status, filterParams, initParams );
  TEST_STATUS( &status );

  LALFindChirpChisqVetoInit( &status, (*filterParams)->chisqParams,
      initParams->numChisqBins, initParams->numPoints );
  TEST_STATUS( &status );

  LALCreateDataSegmentVector( &status, dataSegVec, initParams );
  TEST_STATUS( &status );

  return 0;
}

int Stop(
    DataSegmentVector      **dataSegVec,
    FindChirpFilterInput   **filterInput,
    FindChirpFilterParams  **filterParams,
    FindChirpSegmentVector **fcSegVec,
    UINT4 numChisqBins
    )
{
  LALDestroyFindChirpSegmentVector( &status, fcSegVec );
  TEST_STATUS( &status );

  LALDestroyFindChirpInput( &status, filterInput );
  TEST_STATUS( &status );

  LALFindChirpChisqVetoFinalize( &status, (*filterParams)->chisqParams,
     numChisqBins );
  TEST_STATUS( &status );

  LALFindChirpFilterFinalize( &status, filterParams );
  TEST_STATUS( &status );

  LALDestroyDataSegmentVector( &status, dataSegVec );
  TEST_STATUS( &status );

  return 0;
}


/*
 *
 * SPInit(), SPFini()
 *
 * SPInit() creates and initializes various structures needed for SP filtering.
 * SPFini() destroys the allocated memory.
 *
 */

int SPInit(
    FindChirpTmpltParams **spTmpltParams,
    FindChirpDataParams  **spDataParams,
    FindChirpInitParams     *initParams,
    REAL4 srate,
    REAL4 fmin,
    REAL4 dynRange,
    UINT4 trunc
    )
{
  LALFindChirpTemplateInit( &status, spTmpltParams, initParams );
  TEST_STATUS( &status );

  (*spTmpltParams)->deltaT   = 1 / srate;
  (*spTmpltParams)->fLow     = fmin;
  (*spTmpltParams)->dynRange = dynRange;

  LALFindChirpDataInit( &status, spDataParams, initParams );
  TEST_STATUS( &status );

  (*spDataParams)->fLow         = fmin;
  (*spDataParams)->dynRange     = dynRange;
  (*spDataParams)->invSpecTrunc = trunc;

  return 0;
}

int SPFini(
    FindChirpTmpltParams **spTmpltParams,
    FindChirpDataParams  **spDataParams
    )
{
  LALFindChirpTemplateFinalize( &status, spTmpltParams );
  TEST_STATUS( &status );

  LALFindChirpDataFinalize( &status, spDataParams );
  TEST_STATUS( &status );

  return 0;
}



/*
 *
 * TDInit(), TDFini()
 *
 * TDInit() creates and initializes various structures needed for TD filtering.
 * TDFini() destroys the allocated memory.
 *
 */

int TDInit(
    FindChirpDataParams **tdDataParams,
    FindChirpInitParams    *initParams,
    REAL4 srate,
    REAL4 fmin,
    REAL4 dynRange,
    UINT4 trunc
    )
{
  UINT4 i;
  LALFindChirpDataInit( &status, tdDataParams, initParams );
  TEST_STATUS( &status );

  for ( i = 0; i < (*tdDataParams)->ampVec->length; ++i )
    (*tdDataParams)->ampVec->data[i] = 1;

  (*tdDataParams)->fLow         = fmin;
  (*tdDataParams)->dynRange     = dynRange;
  (*tdDataParams)->invSpecTrunc = trunc;

  return 0;
}

int TDFini(
    FindChirpDataParams **tdDataParams
    )
{
  LALFindChirpDataFinalize( &status, tdDataParams );
  TEST_STATUS( &status );
  return 0;
}



/*
 *
 * SPFilter()
 *
 * Filter the data segment vector through a specified template.  Use a SP
 * template to do the filtering.
 *
 */

int SPFilter(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    FindChirpFilterInput *filterInput,
    FindChirpFilterParams *filterParams,
    FindChirpSegmentVector *fcSegVec,
    FindChirpTmpltParams *tmpltParams,
    FindChirpDataParams *dataParams
    )
{
  InspiralTemplate tmplt;
  UINT4 segment, i;

  tmplt.mass1     = mass1;
  tmplt.mass2     = mass2;
  tmplt.totalMass = tmplt.mass1 + tmplt.mass2;
  tmplt.mu        = tmplt.mass1 * tmplt.mass2 / tmplt.totalMass;
  tmplt.eta       = tmplt.mu / tmplt.totalMass;

  LALFindChirpSPData( &status, fcSegVec, dataSegVec, dataParams );
  TEST_STATUS( &status );

  LALFindChirpSPTemplate( &status, filterInput->fcTmplt, &tmplt, tmpltParams );
  TEST_STATUS( &status );

  for ( segment = 0; segment < fcSegVec->length; ++segment )
  {
    FILE *fp;
    SnglInspiralTable *event = NULL;

    filterInput->segment = fcSegVec->data + segment;

    LALFindChirpFilterSegment( &status, &event, filterInput, filterParams );
    TEST_STATUS( &status );

    fp = fopen( "sp_rhosq.out", "w" );
    for ( i = 0; i < filterParams->rhosqVec->data->length; ++i )
      fprintf( fp, "%e\n", filterParams->rhosqVec->data->data[i] );
    fclose( fp );

    while ( event )
    {
      SnglInspiralTable *thisEvent = event;
      event = thisEvent->next;
      fprintf( stdout, "\nSP Observed Event:\n" );
      fprintf( stdout, "\tsnr = %e\n", thisEvent->snr );
      fprintf( stdout, "\tchisq = %e\n", thisEvent->chisq );
      fprintf( stdout, "\tsigma_sq = %e\n", thisEvent->sigmasq );
      fprintf( stdout, "\td_Mpc = %e\n", thisEvent->eff_distance );
      LALFree( thisEvent );
    }
  }

  return 0;
}



/*
 *
 * TDFilter()
 *
 * Filter the data segment vector through a specified template.  Use a TD
 * template to do the filtering.
 *
 * NOTE: Chi-Squared Stuff is wrong!!!  Need to re-compute bins!!!
 *
 */

int TDFilter(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    REAL4 fmax,
    FindChirpFilterInput *filterInput,
    FindChirpFilterParams *filterParams,
    FindChirpSegmentVector *fcSegVec,
    FindChirpDataParams *dataParams
    )
{
  REAL4Vector    *signal = NULL;
  COMPLEX8Vector *stilde = NULL;
  InspiralTemplate tmplt;
  UINT4 segment;
  UINT4 numPoints;
  REAL4 segNormSum;
  UINT4 n;
  UINT4 k;
  UINT4 i;

  numPoints = dataSegVec->data->chan->data->length;

  LALFindChirpSPData( &status, fcSegVec, dataSegVec, dataParams );
  TEST_STATUS( &status );

  tmplt.mass1           = mass1;
  tmplt.mass2           = mass2;
  tmplt.totalMass       = tmplt.mass1 + tmplt.mass2;
  tmplt.mu              = tmplt.mass1 * tmplt.mass2 / tmplt.totalMass;
  tmplt.eta             = tmplt.mu / tmplt.totalMass;
  tmplt.ieta            = 1;
  tmplt.massChoice      = m1Andm2;
  tmplt.startTime       = 0;
  tmplt.startPhase      = 0;
  tmplt.tSampling       = 1 / fcSegVec->data->deltaT;
  tmplt.fLower          = fcSegVec->data->fLow;
  tmplt.fCutoff         = fmax;
  tmplt.signalAmplitude = 1;
  tmplt.nStartPad       = 0;
  tmplt.nEndPad         = 0;
  tmplt.order           = twoPN;
  tmplt.approximant     = TaylorT2;
  tmplt.massChoice      = m1Andm2;
  tmplt.OmegaS          = 0;
  tmplt.Theta           = 0;

  LALInspiralParameterCalc( &status, &tmplt );
  TEST_STATUS( &status );

  LALInspiralWaveLength( &status, &n, tmplt );
  TEST_STATUS( &status );
  if ( n > numPoints )
  {
    fprintf( stderr, "Chirp is too long!\n" );
    exit( 1 );
  }

  LALSCreateVector( &status, &signal, numPoints );
  TEST_STATUS( &status );

  LALCCreateVector( &status, &stilde, numPoints / 2 + 1 );
  TEST_STATUS( &status );

  LALInspiralWave( &status, signal, &tmplt );
  TEST_STATUS( &status );

  /* shift chirp to end of vector */
  n = numPoints;
  while ( signal->data[--n] == 0 )
    ;
  ++n;
  memmove( signal->data + numPoints - n, signal->data,
      n * sizeof( *signal->data ) );
  memset( signal->data, 0, ( numPoints - n ) * sizeof( *signal->data ) );

  /* fft chirp */
  LALForwardRealFFT( &status, stilde, signal, dataParams->fwdPlan );
  TEST_STATUS( &status );


  /* re-compute data normalization */
  memset( fcSegVec->data->segNorm->data, 0, 
      fcSegVec->data->segNorm->length * sizeof(REAL4) );
  segNormSum = 0;
  for ( k = 1; k < stilde->length; ++k )
  {
    REAL4 re = stilde->data[k].re;
    REAL4 im = stilde->data[k].im;
    REAL4 power = re * re + im * im;
    segNormSum += power * dataParams->wtildeVec->data[k].re;
    fcSegVec->data->segNorm->data[k] += segNormSum;
      
  }


  memset( filterInput->fcTmplt->data->data, 0,
      filterInput->fcTmplt->data->length
      * sizeof( *filterInput->fcTmplt->data->data ) );
  for ( k = 0; k < stilde->length; ++k )
  {
    filterInput->fcTmplt->data->data[k].re = stilde->data[k].re;
    filterInput->fcTmplt->data->data[k].im = stilde->data[k].im;
  }

  filterInput->fcTmplt->tmpltNorm  = 2 * tmplt.mu;
  filterInput->fcTmplt->tmpltNorm *= 2 * LAL_MRSUN_SI / ( 1e6 * LAL_PC_SI );
  filterInput->fcTmplt->tmpltNorm *= dataParams->dynRange;
  filterInput->fcTmplt->tmpltNorm *= filterInput->fcTmplt->tmpltNorm;


  for ( segment = 0; segment < fcSegVec->length; ++segment )
  {
    FILE *fp;
    SnglInspiralTable *event = NULL;

    filterInput->segment = fcSegVec->data + segment;

    LALFindChirpFilterSegment( &status, &event, filterInput, filterParams );
    TEST_STATUS( &status );

    fp = fopen( "td_rhosq.out", "w" );
    for ( i = 0; i < filterParams->rhosqVec->data->length; ++i )
      fprintf( fp, "%e\n", filterParams->rhosqVec->data->data[i] );
    fclose( fp );

    while ( event )
    {
      SnglInspiralTable *thisEvent = event;
      event = thisEvent->next;
      fprintf( stdout, "\nTD Observed Event:\n" );
      fprintf( stdout, "\tsnr = %e\n", thisEvent->snr );
      fprintf( stdout, "\tchisq = %e\n", thisEvent->chisq );
      fprintf( stdout, "\tsigmasq = %e\n", thisEvent->sigmasq );
      fprintf( stdout, "\td_Mpc = %e\n", thisEvent->eff_distance );
      LALFree( thisEvent );
    }
  }

  LALSDestroyVector( &status, &signal );
  TEST_STATUS( &status );

  LALCDestroyVector( &status, &stilde );
  TEST_STATUS( &status );

  return 0;
}




/*
 *
 * MakeData()
 *
 * Populate the dataSegVec with an injected inspiral and set the spectrum
 * and response vectors to unity.
 *
 */

int MakeData(
    DataSegmentVector *dataSegVec,
    REAL4 mass1,
    REAL4 mass2,
    REAL4 srate,
    REAL4 fmin,
    REAL4 fmax
    )
{
  InspiralTemplate tmplt;
  UINT4 i;
  UINT4 k;
  UINT4 n;
  FILE *fp;

  tmplt.mass1           = mass1;
  tmplt.mass2           = mass2;
  tmplt.totalMass       = tmplt.mass1 + tmplt.mass2;
  tmplt.mu              = tmplt.mass1 * tmplt.mass2 / tmplt.totalMass;
  tmplt.eta             = tmplt.mu / tmplt.totalMass;
  tmplt.ieta            = 1;
  tmplt.massChoice      = m1Andm2;
  tmplt.startTime       = 0;
  tmplt.startPhase      = 0;
  tmplt.fLower          = fmin;
  tmplt.fCutoff         = fmax;
  tmplt.tSampling       = srate;
  tmplt.signalAmplitude = 1;
  tmplt.nStartPad       = dataSegVec->data->chan->data->length / 2;
  tmplt.nEndPad         = 0;
  tmplt.order           = twoPN;
  tmplt.approximant     = TaylorT2;
  tmplt.massChoice      = m1Andm2;
  tmplt.OmegaS          = 0;
  tmplt.Theta           = 0;

  LALInspiralParameterCalc( &status, &tmplt );
  TEST_STATUS( &status );

  LALInspiralWaveLength( &status, &n, tmplt );
  TEST_STATUS( &status );
  if ( n > dataSegVec->data->chan->data->length )
  {
    fprintf( stderr, "Chirp is too long!\n" );
    exit( 1 );
  }

  memset( dataSegVec->data->chan->data->data, 0,
      dataSegVec->data->chan->data->length
      * sizeof( *dataSegVec->data->chan->data->data ) );
  LALInspiralWave( &status, dataSegVec->data->chan->data, &tmplt );
  TEST_STATUS( &status );

  dataSegVec->data->chan->deltaT = 1 / srate;
  dataSegVec->data->spec->deltaF = srate / dataSegVec->data->chan->data->length;

  dataSegVec->data->chan->epoch.gpsSeconds     = 0;
  dataSegVec->data->chan->epoch.gpsNanoSeconds = 0;

  for ( k = 0; k < dataSegVec->data->spec->data->length; ++k )
  {
    dataSegVec->data->spec->data->data[k]    = 1;
    dataSegVec->data->resp->data->data[k].re = 1;
    dataSegVec->data->resp->data->data[k].im = 0;
  }

  fp = fopen( "data.out", "w" );
  for ( i = 0; i < dataSegVec->data->chan->data->length; ++i )
    fprintf( fp, "%e\t%e\n", i * dataSegVec->data->chan->deltaT,
	dataSegVec->data->chan->data->data[i] );
  fclose( fp );

  return 0;
}
