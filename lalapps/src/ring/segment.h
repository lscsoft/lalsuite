#ifndef SEGMENT_H
#define SEGMENT_H

/*
 *
 * Routine to create a single overwhitened data segment from a time series:
 * this performs the FFT to get the data into the frequency domain and then
 * multiplies it by the inverse power spectrum to overwhiten it.
 *
 */

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
