/*
*  Copyright (C) 2007 Jolien Creighton
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
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>

#include "lalapps.h"
#include "segment.h"
#include "errutil.h"
#include "gpstime.h"

RCSID( "$Id$" );

/* routine to compute a single overwhitened data segment */
int compute_data_segment(
    COMPLEX8FrequencySeries  *segment,
    UINT4                     segmentNumber,
    REAL4TimeSeries          *series,
    REAL4FrequencySeries     *invspec,
    COMPLEX8FrequencySeries  *response,
    REAL8                     segmentDuration,
    REAL8                     strideDuration,
    REAL4FFTPlan             *fwdPlan
    )
{
  REAL4TimeSeries seg;
  REAL4Vector     vec;
  INT8            ns;
  UINT4           segmentLength;
  UINT4           segmentStride;

  segmentLength  = floor( segmentDuration/series->deltaT + 0.5 );
  segmentStride  = floor( strideDuration/series->deltaT + 0.5 );

  verbose( "computing overwhitened data segment %u\n", segmentNumber );

  /* name of this segment */
  snprintf( segment->name, sizeof( segment->name ), "%s_SEG_%u",
      series->name, segmentNumber );

  /* allocate memory for the data */
  segment->data = XLALCreateCOMPLEX8Vector( segmentLength/2 + 1 );

  /* create a time series that contains only the relevant data */
  seg        = *series;
  vec.length = segmentLength;
  vec.data   = seg.data->data + segmentNumber * segmentStride;
  seg.data   = &vec;
  ns  = epoch_to_ns( &seg.epoch );
  ns += sec_to_ns( segmentNumber * strideDuration );
  ns_to_epoch( &seg.epoch, ns );

  /* fft the data */
  XLALREAL4TimeFreqFFT( segment, &seg, fwdPlan );

  /* multiply by the response function */
  if ( response )
  {
    XLALCCVectorMultiply( segment->data, segment->data, response->data );
    XLALUnitMultiply( &segment->sampleUnits, &segment->sampleUnits,
        &response->sampleUnits );
  }

  /* multiply by the inverse spectrum */
  XLALSCVectorMultiply( segment->data, invspec->data, segment->data );
  XLALUnitMultiply( &segment->sampleUnits, &segment->sampleUnits,
      &invspec->sampleUnits );

  return 0;
}
