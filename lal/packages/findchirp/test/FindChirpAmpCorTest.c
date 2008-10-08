/********************** <lalVerbatim file="COMPLETEFindChirpAmpCorTestCV">
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{FindChirpAmpCorlTest.c}}
\label{ss:FindChirpAmpCorTest.c}

Provides the necessary function to test the AmpCorPPN filter.


****************************************** </lalLaTeX><lalErrTable> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include "FindChirpTDTest.h"

#include <lal/LALRCSID.h>
NRCSID (FINDCHIRPAMPCORTESTC,"$Id$");

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
  const REAL4 mass2        = 10.;    /* solar masses */
  const REAL4 srate        = 16384;  /* Hz */
  const REAL4 fmin         = 100;    /* Hz */
  const REAL4 fmax         = 1000;   /* Hz */
  const REAL4 dynRange     = 1;
  int i;
  FILE *harm1, *harm2, *harm3;


  FindChirpInitParams initParams; /* need to populate this by hand */
  InspiralTemplate mytmplt;

  /* these are created by Start() */
  FindChirpFilterInput   *filterInput   = NULL;

  FindChirpFilterParams  *filterParams  = NULL;
  FindChirpSegmentVector *fcSegVec      = NULL;
  DataSegmentVector      *dataSegVec    = NULL;

  /* these are required for filtering; they are created by Init() */
  FindChirpTmpltParams *tmpltParams = NULL;
  FindChirpDataParams  *dataParams  = NULL;

  
  /* set initialization parameters */
  initParams.numSegments    = numSegments;
  initParams.numPoints      = numPoints;
  initParams.numChisqBins   = numChisqBins;
  initParams.approximant    = AmpCorPPN;
  initParams.createRhosqVec = 1;


  /* create objects needed by  filters */
  Start( &dataSegVec, &filterInput, &filterParams, &fcSegVec, &initParams );


  /* set filter parameters, e.g., thresholds for events */
  filterParams->deltaT         = 1 / srate;
  filterParams->rhosqThresh    = 1e-6;
  filterParams->chisqThresh    = 1e+6;


  /* create some fake data */
  MakeData( dataSegVec, mass1, mass2, srate, fmin, fmax );

  mytmplt.mass1           = mass1;
  mytmplt.mass2           = mass2;
  mytmplt.totalMass       = mytmplt.mass1 + mytmplt.mass2;
  mytmplt.mu              = mytmplt.mass1 * mytmplt.mass2 / mytmplt.totalMass;
  mytmplt.eta             = mytmplt.mu / mytmplt.totalMass;
  mytmplt.ieta            = 1;
  mytmplt.massChoice      = m1Andm2;
  mytmplt.startTime       = 0.;
  mytmplt.startPhase      = 0.;
  mytmplt.tSampling       = 1 / fcSegVec->data->deltaT;
  mytmplt.fLower          = fcSegVec->data->fLow;
  mytmplt.fCutoff         = fmax;
  mytmplt.signalAmplitude = 1;
  mytmplt.nStartPad       = 0;
  mytmplt.nEndPad         = 0;
  mytmplt.order           = twoPN;
  mytmplt.approximant     = AmpCorPPN;
  mytmplt.massChoice      = m1Andm2;
  mytmplt.OmegaS          = 0;
  mytmplt.Theta           = 0;
  mytmplt.inclination     = LAL_PI/2.0;




  /*
   * initialize specific parameters
   */

  initParams.approximant = AmpCorPPN;
  initParams.numPoints = dataSegVec->data->chan->data->length;

  Init( &tmpltParams, &dataParams, &initParams, srate, fmin, dynRange,								     invSpecTrunc );


  for(i = 0; i < dataSegVec->data->chan->data->length; i++)
    tmpltParams->PTFQ->data[i] = dataSegVec->data->chan->data->data[i];

  tmpltParams->taperTmplt = INSPIRAL_TAPER_NONE;

  fprintf( stderr, "Testing AmpCorTemplate...\n" );

  LALFindChirpAmpCorTemplate( &status, filterInput->fcTmplt, &mytmplt, 
                                                             tmpltParams  );

  harm1 = fopen("harm1.dat", "w");
  for(i=0; i < filterInput->fcTmplt->PTFQtilde->vectorLength; i++ )  
  fprintf( harm1, "%1.6e %1.6e\n", 
                         filterInput->fcTmplt->PTFQtilde->data[i].re,
                         filterInput->fcTmplt->PTFQtilde->data[i].im );
  fclose(harm1);

  harm2 = fopen("harm2.dat", "w");
  for(i = filterInput->fcTmplt->PTFQtilde->vectorLength; 
      i < 2 * filterInput->fcTmplt->PTFQtilde->vectorLength; i++ )  
  fprintf( harm2, "%1.6e %1.6e\n", 
                         filterInput->fcTmplt->PTFQtilde->data[i].re,
                         filterInput->fcTmplt->PTFQtilde->data[i].im );
  fclose(harm2);

  harm3 = fopen("harm3.dat", "w");
  for(i = 2*filterInput->fcTmplt->PTFQtilde->vectorLength; 
      i < 3 * filterInput->fcTmplt->PTFQtilde->vectorLength; i++ )  
  fprintf( harm3, "%1.6e %1.6e\n", 
                         filterInput->fcTmplt->PTFQtilde->data[i].re,
                         filterInput->fcTmplt->PTFQtilde->data[i].im );
  fclose(harm3);



  LALFindChirpTDData( &status, fcSegVec, dataSegVec, dataParams );

  filterInput->segment = fcSegVec->data;

  fprintf( stderr, "Testing AmpCorNormalize...\n" );
  LALFindChirpAmpCorNormalize( &status, filterInput->fcTmplt, 
                               filterInput->segment, dataParams );

  Stop( &dataSegVec, &filterInput, &filterParams, &fcSegVec, numChisqBins );
  Fini( &tmpltParams, &dataParams );

  initParams.approximant = GeneratePPN;




  fprintf( stderr, "Testing TDTemplate with same template...\n" );
  Start( &dataSegVec, &filterInput, &filterParams, &fcSegVec, &initParams );
  MakeData( dataSegVec, mass1, mass2, srate, fmin, fmax );
  Init( &tmpltParams, &dataParams, &initParams, srate, fmin, dynRange,								     invSpecTrunc );
  for(i = 0; i < dataSegVec->data->chan->data->length; i++)
    tmpltParams->xfacVec->data[i] = dataSegVec->data->chan->data->data[i];

  LALFindChirpTDTemplate( &status, filterInput->fcTmplt, &mytmplt,
                                                               tmpltParams );


  /* clean up memory and exit */


  Stop( &dataSegVec, &filterInput, &filterParams, &fcSegVec, numChisqBins );



  Fini( &tmpltParams, &dataParams );
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
 * Init(), Fini()
 *
 * Init() creates and initializes various structures needed for filtering.
 * Fini() destroys the allocated memory.
 *
 */

int Init(
    FindChirpTmpltParams **tmpltParams,
    FindChirpDataParams  **dataParams,
    FindChirpInitParams     *initParams,
    REAL4 srate,
    REAL4 fmin,
    REAL4 dynRange,
    UINT4 trunc
    )
{
  LALFindChirpTemplateInit( &status, tmpltParams, initParams );
  TEST_STATUS( &status );

  (*tmpltParams)->deltaT   = 1 / srate;
  (*tmpltParams)->fLow     = fmin;
  (*tmpltParams)->dynRange = dynRange;

  LALFindChirpDataInit( &status, dataParams, initParams );
  TEST_STATUS( &status );

  (*dataParams)->fLow         = fmin;
  (*dataParams)->dynRange     = dynRange;
  (*dataParams)->invSpecTrunc = trunc;


  return 0;
}

int Fini(
    FindChirpTmpltParams **tmpltParams,
    FindChirpDataParams  **dataParams
    )
{
  LALFindChirpTemplateFinalize( &status, tmpltParams );
  TEST_STATUS( &status );

  LALFindChirpDataFinalize( &status, dataParams );
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



