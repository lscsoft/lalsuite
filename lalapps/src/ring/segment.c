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

COMPLEX8FrequencySeries * compute_data_segments(
    UINT4                     numSegments,
    UINT4                    *segmentNumbers,
    REAL4TimeSeries          *series,
    REAL4FrequencySeries     *invspec,
    COMPLEX8FrequencySeries  *response,
    REAL8                     segmentDuration,
    REAL8                     strideDuration,
    REAL4FFTPlan             *fwdPlan
    )
{
  COMPLEX8FrequencySeries *segments;
  LALStatus       status = blank_status;
  REAL4TimeSeries seg;
  REAL4Vector     vec;
  INT8            ns;
  UINT4           segmentLength;
  UINT4           segmentStride;
  UINT4           sgmnt;

  segments = LALCalloc( numSegments, sizeof( *segments ) );

  segmentLength  = floor( segmentDuration/series->deltaT + 0.5 );
  segmentStride  = floor( strideDuration/series->deltaT + 0.5 );

  for ( sgmnt = 0; sgmnt < numSegments; ++sgmnt )
  {
    COMPLEX8FrequencySeries *thisSegment;
    UINT4 segmentNumber;

    thisSegment = segments + sgmnt;

    /* the number of this segment to compute */
    if ( segmentNumbers )
      segmentNumber = segmentNumbers[sgmnt];
    else
      segmentNumber = sgmnt;

    verbose( "computing overwhitened data segment %u\n", segmentNumber );

    /* name of this segment */
    LALSnprintf( thisSegment->name, sizeof( thisSegment->name ), "%s_SEG_%u",
        series->name, segmentNumber );

    /* allocate memory for the data */
    LAL_CALL( LALCCreateVector(&status, &thisSegment->data, segmentLength/2+1),
        &status );

    seg = *series;
    vec.length = segmentLength;
    vec.data   = seg.data->data + segmentNumber * segmentStride;
    seg.data   = &vec;

    ns  = epoch_to_ns( &seg.epoch );
    ns += sec_to_ns( segmentNumber * strideDuration );
    ns_to_epoch( &seg.epoch, ns );

    /* fft the data */
    LAL_CALL( LALTimeFreqRealFFT( &status, thisSegment, &seg, fwdPlan ),
        &status );

    /* multiply by the response function */
    if ( response )
    {
      LAL_CALL( LALCCVectorMultiply( &status, thisSegment->data,
            thisSegment->data, response->data ), &status );
      XLALUnitMultiply( &thisSegment->sampleUnits, &thisSegment->sampleUnits,
          &response->sampleUnits );
    }

    /* multiply by the inverse spectrum */
    LAL_CALL( LALSCVectorMultiply( &status, thisSegment->data, invspec->data,
          thisSegment->data ), &status );
    XLALUnitMultiply( &thisSegment->sampleUnits, &thisSegment->sampleUnits,
        &invspec->sampleUnits );
  }

  return segments;
}

