#ifndef SEGMENT_H
#define SEGMENT_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>

COMPLEX8FrequencySeries * compute_data_segments(
    UINT4                     numSegments,
    UINT4                    *segmentNumbers,
    REAL4TimeSeries          *series,
    REAL4FrequencySeries     *invspec,
    COMPLEX8FrequencySeries  *response,
    REAL8                     segmentDuration,
    REAL8                     strideDuration,
    REAL4FFTPlan             *fwdPlan
    );

#endif /* SEGMENT_H */
