#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>

#include "lalapps.h"
#include "spectrm.h"
#include "errutil.h"

RCSID( "$Id$" );


/* routine to compute an average spectrum from time series data */
/* FIXME: only uses median-mean method */
REAL4FrequencySeries *compute_average_spectrum(
    REAL4TimeSeries         *series,
    REAL8                    segmentDuration,
    REAL8                    strideDuration,
    REAL4FFTPlan            *fwdPlan,
    int                      whiteSpectrum
    )
{
  LALStatus status = blank_status;
  REAL4Window  *window  = NULL;
  REAL4FrequencySeries *spectrum;
  UINT4 segmentLength;
  UINT4 segmentStride;

  segmentLength  = floor( segmentDuration/series->deltaT + 0.5 );
  segmentStride  = floor( strideDuration/series->deltaT + 0.5 );

  spectrum       = LALCalloc( 1, sizeof( *spectrum ) );
  spectrum->data = XLALCreateREAL4Vector( segmentLength/2 + 1 );

  window = XLALCreateWelchREAL4Window( segmentLength );

  if ( whiteSpectrum ) /* just return a constant spectrum */
  {
    UINT4 k;
    REAL4 spec;
    spec = 2.0 * series->deltaT;
    verbose( "creating white spectrum with constant value %g\n", spec );
    for ( k = 1; k < spectrum->data->length - 1; ++k )
      spectrum->data->data[k] = spec;
    /* DC and Nyquist */
    spectrum->data->data[0] = 2.0 * spec;
    spectrum->data->data[spectrum->data->length - 1] = 2.0 * spec;
    spectrum->epoch  = series->epoch;
    spectrum->deltaF = 1.0/segmentDuration;
  }
  else /* compute average spectrum using median-mean method */
  {
    verbose( "estimating average spectrum using median-mean method\n" );
    XLALREAL4AverageSpectrumMedianMean( spectrum, series, segmentLength,
        segmentStride, window, fwdPlan );
  }

  LALSnprintf( spectrum->name, sizeof( spectrum->name ),
      "%s_SPEC", series->name );

  XLALDestroyREAL4Window( window );

  return spectrum;
}


/* routine to invert and truncate (to have compact time support) a spectrum */
int invert_spectrum(
    REAL4FrequencySeries *spectrum,
    REAL8                 dataSampleRate,
    REAL8                 strideDuration,
    REAL8                 truncateDuration,
    REAL8                 lowCutoffFrequency,
    REAL4FFTPlan         *fwdPlan,
    REAL4FFTPlan         *revPlan
    )
{
  REAL8 segmentDuration;
  UINT4 segmentLength;
  UINT4 segmentStride;
  UINT4 truncateLength;

  segmentDuration = 1.0/spectrum->deltaF;
  segmentLength = floor( segmentDuration * dataSampleRate + 0.5 );
  segmentStride = floor( strideDuration * dataSampleRate + 0.5 );
  if ( truncateDuration > 0.0 )
    truncateLength = floor( truncateDuration * dataSampleRate + 0.5 );
  else
    truncateLength = 0;

  verbose( "computing inverse spectrum with truncation length %d\n",
      truncateLength );

  XLALREAL4SpectrumInvertTruncate( spectrum, lowCutoffFrequency,
      segmentLength, truncateLength, fwdPlan, revPlan );

  LALSnprintf( spectrum->name, sizeof( spectrum->name ),
      "%s_INV", spectrum->name );

  return 0;
}


/* routine to scale a spectrum by the magnitude of the response function */
int calibrate_spectrum(
    REAL4FrequencySeries    *spectrum,
    COMPLEX8FrequencySeries *response,
    REAL8                    lowCutoffFrequency,
    int                      inverse
    )
{
  UINT4 cut;
  UINT4 k;

  if ( response )
  {
    /* compute low frequency cutoff */
    if ( lowCutoffFrequency > 0.0 )
      cut = lowCutoffFrequency / spectrum->deltaF;
    else
      cut = 0;

    /* apply the response function */
    if ( inverse ) /* divide by response */
    {
      for ( k = cut; k < spectrum->data->length; ++k )
      {
        REAL4 re = response->data->data[k].re;
        REAL4 im = response->data->data[k].im;
        spectrum->data->data[k] /= (re*re + im*im );
      }
      XLALUnitMultiply( &spectrum->sampleUnits, &spectrum->sampleUnits,
          &response->sampleUnits );
    }
    else /* multiply by response */
    {
      for ( k = cut; k < spectrum->data->length; ++k )
      {
        REAL4 re = response->data->data[k].re;
        REAL4 im = response->data->data[k].im;
        spectrum->data->data[k] *= (re*re + im*im );
      }
      XLALUnitDivide( &spectrum->sampleUnits, &spectrum->sampleUnits,
          &response->sampleUnits );
    }
    LALSnprintf( spectrum->name, sizeof( spectrum->name ),
        "%s_CAL", spectrum->name );
  }

  return 0;
}

