#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/RingSearch.h>

NRCSID( AVGSPECC, "$Id$" );

static int mycompare( const void *p, const void *q )
{
  const REAL4 *a = p;
  const REAL4 *b = q;
  int ans = 0;
  if ( *a < *b )
    ans = -1;
  else if ( *a > *b )
    ans = 1;
  return ans;
}

void
LALMedianSpectrum(
    LALStatus            *status,
    REAL4FrequencySeries *output,
    REAL4TimeSeries      *input,
    AvgSpecParams        *params
    )
{
  /* factor that converts median of exponential distribution to mean */
  const REAL4 corrfac = 1.0 / LAL_LN2;

  REAL4FrequencySeries *spectrum;
  REAL4TimeSeries       segment;
  REAL4Vector          *window = NULL;
  LALWindowParams       wpar;

  LALUnitPair unitPair;

  REAL4 *arr;
  REAL4  fac;

  UINT4 numseg;
  UINT4 segsz;
  UINT4 seg;
  UINT4 med;
  UINT4 k;


  INITSTATUS( status, "LALMedianSpectrum", AVGSPECC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( params, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( input, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );

  segsz = params->segsize;

  if ( ( 2 * input->data->length ) % segsz )
  {
    ABORT( status, RINGSEARCHH_ENSEG, RINGSEARCHH_MSGENSEG );
  }
  numseg = 2 * input->data->length / segsz - 1;
  if ( numseg < 2 )
  {
    ABORT( status, RINGSEARCHH_ELSEG, RINGSEARCHH_MSGELSEG );
  }

  /* set up metadata first */
  output->epoch    = input->epoch;
  output->deltaF   = 1.0 / ( segsz * input->deltaT );
  output->f0       = 0;
  unitPair.unitOne = &input->sampleUnits;
  unitPair.unitTwo = &input->sampleUnits;
  LALUnitMultiply( status->statusPtr, &output->sampleUnits, &unitPair );
  CHECKSTATUSPTR( status );
  unitPair.unitOne = &output->sampleUnits;
  unitPair.unitTwo = &lalSecondUnit;
  LALUnitMultiply( status->statusPtr, &output->sampleUnits, &unitPair );
  CHECKSTATUSPTR( status );


  /* allocate memory for spectra */
  spectrum = LALCalloc( numseg, sizeof( *spectrum ) );
  if ( ! spectrum )
  {
    ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
  }
  for ( seg = 0; seg < numseg; ++seg )
  {
    LALCreateVector( status->statusPtr, &spectrum[seg].data, segsz / 2 + 1 );
    BEGINFAIL( status )
    {
      while ( seg-- )
      {
        LALDestroyVector( status->statusPtr, &spectrum[seg].data );
      }
      LALFree( spectrum );
    }
    ENDFAIL( status );
  }

  /* allocate memory for temporary data segment */
  segment = *input;
  segment.data = NULL;
  LALCreateVector( status->statusPtr, &segment.data, segsz );
  BEGINFAIL( status )
  {
    for ( seg = 0; seg < numseg; ++seg )
    {
      LALDestroyVector( status->statusPtr, &spectrum[seg].data );
    }
    LALFree( spectrum );
  }
  ENDFAIL( status );

  /* allocate memory for window */
  LALCreateVector( status->statusPtr, &window, segsz );
  BEGINFAIL( status )
  {
    for ( seg = 0; seg < numseg; ++seg )
    {
      LALDestroyVector( status->statusPtr, &spectrum[seg].data );
    }
    LALFree( spectrum );
    LALDestroyVector( status->statusPtr, &segment.data );
  }
  ENDFAIL( status );

  /* compute window function */
  wpar.type = params->wintype;
  wpar.length = segsz;
  LALWindow( status->statusPtr, window, &wpar );
  BEGINFAIL( status )
  {
    for ( seg = 0; seg < numseg; ++seg )
    {
      LALDestroyVector( status->statusPtr, &spectrum[seg].data );
    }
    LALFree( spectrum );
    LALDestroyVector( status->statusPtr, &segment.data );
    LALDestroyVector( status->statusPtr, &window );
  }
  ENDFAIL( status );
  fac = sqrt( input->deltaT / wpar.sumofsquares );

  /* segment input data and compute spectra */
  for ( seg = 0; seg < numseg; ++seg )
  {
    UINT4 i;
    for ( i = 0; i < segment.data->length; ++i )
    {
      segment.data->data[i]  = input->data->data[i + seg * segsz / 2];
      segment.data->data[i] *= fac * window->data[i];
    }

    /* compute the power spectrum */
    LALRealPowerSpectrum( status->statusPtr, spectrum[seg].data, segment.data,
        params->fwdplan );
    BEGINFAIL( status )
    {
      for ( seg = 0; seg < numseg; ++seg )
      {
        LALDestroyVector( status->statusPtr, &spectrum[seg].data );
      }
      LALFree( spectrum );
      LALDestroyVector( status->statusPtr, &segment.data );
      LALDestroyVector( status->statusPtr, &window );
    }
    ENDFAIL( status );
    /* CHANGE BACK TO ORIGINAL NORMALIZATION -- JC */
    {
      REAL4Vector *myvector = spectrum[seg].data;
      UINT4 mybin;
      for ( mybin = 1; mybin < myvector->length - 1; ++mybin )
        myvector->data[mybin] *= 0.5;
    }
  }

  /* done with temporary segment vector and window */
  TRY( LALDestroyVector( status->statusPtr, &segment.data ), status );
  TRY( LALDestroyVector( status->statusPtr, &window ), status );

  /* allocate memory for sort array */
  arr = LALMalloc( numseg * sizeof( *arr ) );
  if ( ! arr )
  {
    for ( seg = 0; seg < numseg; ++seg )
    {
      LALDestroyVector( status->statusPtr, &spectrum[seg].data );
    }
    LALFree( spectrum );
    ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
  }

  /* loop over frequency bins finding median (numseg should be odd) */
  med = ( numseg + 1 ) / 2;
  for ( k = 0; k < segsz / 2 + 1; ++k )
  {
    for ( seg = 0; seg < numseg; ++seg )
    {
      arr[seg] = spectrum[seg].data->data[k];
    }
    qsort( arr, numseg, sizeof( *arr ), mycompare );
    /* correct by corrfac: appropriate factor for Gaussian noise */
    output->data->data[k] = 2 * corrfac * arr[med];
  }
  /* but DC shouldn't have a factor of two */
  output->data->data[0] /= 2;

  LALFree( arr );
  for ( seg = 0; seg < numseg; ++seg )
  {
    LALDestroyVector( status->statusPtr, &spectrum[seg].data );
  }
  LALFree( spectrum );

  output->f0 = 0;
  output->deltaF = 1.0 / ( segsz * input->deltaT );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
