#ifndef SEGMENT_H
#define SEGMENT_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>

int compute_data_segment(
    COMPLEX8FrequencySeries  *segment,
    UINT4                     segmentNumber,
    REAL4TimeSeries          *series,
    REAL4FrequencySeries     *invspec,
    COMPLEX8FrequencySeries  *response,
    REAL8                     segmentDuration,
    REAL8                     strideDuration,
    REAL4FFTPlan             *fwdPlan
    );

#endif /* SEGMENT_H */
