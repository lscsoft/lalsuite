/*
*  Copyright (C) 2007 Jolien Creighton, Lisa M. Goggin
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

#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/LIGOMetadataRingdownUtils.h>

#include "lalapps.h"
#include "spectrm.h"
#include "errutil.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* routine to compute an average spectrum from time series data */
REAL4FrequencySeries *compute_average_spectrum(
    REAL4TimeSeries         *series,
    int                      spectrumAlgthm,
    REAL8                    segmentDuration,
    REAL8                    strideDuration,
    REAL4FFTPlan            *fwdPlan,
    int                      whiteSpectrum
    )
{
  /*LALStatus status = blank_status;*/
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
  else /* compute average spectrum using either the median or the median-mean method */
  {
    switch ( spectrumAlgthm )
    {
      case LALRINGDOWN_SPECTRUM_MEDIAN:
        verbose( "estimating average spectrum using median method\n" );
        XLALREAL4AverageSpectrumMedian(spectrum, series, segmentLength,
          segmentStride, window, fwdPlan );
      break;
      case LALRINGDOWN_SPECTRUM_MEDIAN_MEAN:
        verbose( "estimating average spectrum using median-mean method\n" );
        XLALREAL4AverageSpectrumMedianMean( spectrum, series, segmentLength,
          segmentStride, window, fwdPlan );
      break;
      default:
        error( "unrecognized injection signal type\n" );
    }
  }

  snprintf( spectrum->name, sizeof( spectrum->name ),
      "%s_SPEC", series->name );

  XLALDestroyREAL4Window( window );

  return spectrum;
}

REAL4FrequencySeries *resample_psd(
  REAL4FrequencySeries *origspectrum,
  REAL8        sampleRate,
  REAL8        segmentDuration
  )
{
  REAL4FrequencySeries *rsmplspectrum;
  UINT4 segmentLength;
  REAL8 origDeltaF,rsmplDeltaF;
  REAL8 currF,belowF,aboveF,belowWeight,aboveWeight,currPower;
  UINT4 k;
  INT4 belowK,aboveK;

  segmentLength  = floor( segmentDuration*sampleRate + 0.5 );
  rsmplspectrum       = LALCalloc( 1, sizeof( *rsmplspectrum ) );
  rsmplspectrum->data = XLALCreateREAL4Vector( segmentLength/2 + 1 );

  origDeltaF = origspectrum->deltaF;
  rsmplDeltaF = 1.0/segmentDuration;
  rsmplspectrum->epoch = origspectrum->epoch;
  rsmplspectrum->f0 = origspectrum->f0;
  rsmplspectrum->deltaF = rsmplDeltaF;
  rsmplspectrum->sampleUnits = origspectrum->sampleUnits;
  snprintf( rsmplspectrum->name, sizeof( rsmplspectrum->name ),
      "%s_RSMPL", origspectrum->name);

  for (k=0; k < segmentLength/2 + 1; k++)
  {
    currF = k * rsmplDeltaF;
    belowK = floor( currF / origDeltaF);
    aboveK = ceil( currF / origDeltaF);
    belowF = belowK * origDeltaF;
    aboveF = aboveK * origDeltaF;
    if (belowK < 0)
      belowK = 0;
    if (aboveK < 0)
    {
      // Shouldn't this fail?
      aboveK = 0;
    }
    if (belowK > (INT4)origspectrum->data->length)
      belowK = origspectrum->data->length;
    if (aboveK > (INT4)origspectrum->data->length)
      aboveK = origspectrum->data->length;
    sanity_check( aboveK >= belowK );
    if (aboveK - belowK == 0)
    {
      belowWeight = 1;
      aboveWeight = 0;
    }
    else
    {
      belowWeight = (currF - belowF) / (aboveF - belowF);
      aboveWeight = (aboveF - currF) / (aboveF - belowF);
      sanity_check( aboveWeight >= 0);
      sanity_check( belowWeight >= 0);
    }
    currPower = origspectrum->data->data[belowK] * belowWeight;
    currPower += origspectrum->data->data[aboveK] * aboveWeight;
    rsmplspectrum->data->data[k] = currPower;
  }

  XLALDestroyREAL4Vector( origspectrum->data );
  LALFree( origspectrum );

  return rsmplspectrum;
}

REAL8FrequencySeries *generate_theoretical_psd(
    REAL4                    deltaT,
    REAL8                    segmentDuration,
    UINT4                    spectrumNumber,
    REAL8                    simScale
)
{
  REAL8FrequencySeries *spectrum;
  UINT4 segmentLength;

  segmentLength  = floor( segmentDuration/deltaT + 0.5 );

  spectrum       = LALCalloc( 1, sizeof( *spectrum ) );
  spectrum->data = XLALCreateREAL8Vector( segmentLength/2 + 1 );

  spectrum->deltaF = 1.0/segmentDuration;

  if (spectrumNumber == WHITE_PSD) /* just return a constant spectrum */
  {
    UINT4 k;
    REAL4 spec;
    spec = 2.0 * deltaT;
    verbose( "creating white PSD with constant value %g\n", spec * simScale * simScale );
    for ( k = 1; k < spectrum->data->length - 1; ++k )
      spectrum->data->data[k] = spec * simScale * simScale;
    /* DC and Nyquist are set to 0*/
    spectrum->data->data[0] = 0.0;
    spectrum->data->data[spectrum->data->length - 1] = 0.0;
    snprintf( spectrum->name, sizeof( spectrum->name ),
      "WHITE_NOISE_PSD" );
  }
  else if ( spectrumNumber == ILIGO_PSD )
  {
    verbose( "Creating initial LIGO PSD. PSD only generated from 30Hz \n" ); 
    /* FIXME  */
    REAL4 flow = 30;
    XLALSimNoisePSD(spectrum, flow, XLALSimNoisePSDiLIGOSRD);   
    snprintf( spectrum->name, sizeof( spectrum->name ),
      "iLIGO_PSD" );
  }
  else if ( spectrumNumber == ALIGO_PSD )
  {
    verbose( "Creating advanced LIGO high-power broad-band signal recycling "
             "PSD. PSD only generated from 10Hz \n ");
    REAL4 flow = 10;
    XLALSimNoisePSD(spectrum, flow, XLALSimNoisePSDaLIGOZeroDetHighPower);
    snprintf( spectrum->name, sizeof( spectrum->name ),
      "aLIGO_PSD" );
  }
  else
  {
    fprintf(stderr,"Spectrum number not valid. This message should not be seen."
            " Have you broken the code?? \n");
  }

  return spectrum;
}





/* routine to invert and truncate (to have compact time support) a spectrum */
int invert_spectrum(
    REAL4FrequencySeries *spectrum,
    REAL8                 dataSampleRate,
    REAL8                 UNUSED strideDuration,
    REAL8                 truncateDuration,
    REAL8                 lowCutoffFrequency,
    REAL4FFTPlan         *fwdPlan,
    REAL4FFTPlan         *revPlan
    )
{
  REAL8 segmentDuration;
  UINT4 segmentLength;
  UINT4 UNUSED segmentStride;
  UINT4 truncateLength;
  char name[LALNameLength];

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

  strncpy( name, spectrum->name, LALNameLength * sizeof(char) );
  snprintf( spectrum->name, sizeof( spectrum->name ),
      "%s_INV", name );

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
  char name[LALNameLength];

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
        REAL4 re = crealf(response->data->data[k]);
        REAL4 im = cimagf(response->data->data[k]);
        spectrum->data->data[k] /= (re*re + im*im );
      }
      XLALUnitMultiply( &spectrum->sampleUnits, &spectrum->sampleUnits,
          &response->sampleUnits );
    }
    else /* multiply by response */
    {
      for ( k = cut; k < spectrum->data->length; ++k )
      {
        REAL4 re = crealf(response->data->data[k]);
        REAL4 im = cimagf(response->data->data[k]);
        spectrum->data->data[k] *= (re*re + im*im );
      }
      XLALUnitDivide( &spectrum->sampleUnits, &spectrum->sampleUnits,
          &response->sampleUnits );
    }
  strncpy( name, spectrum->name, LALNameLength * sizeof(char) );
    snprintf( spectrum->name, sizeof( spectrum->name ),
        "%s_CAL", name );
  }

  return 0;
}

