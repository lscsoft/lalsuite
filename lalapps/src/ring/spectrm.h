#ifndef INVSPEC_H
#define INVSPEC_H

#include <lal/LALDatatypes.h>

REAL4FrequencySeries *compute_average_spectrum(
    REAL4TimeSeries         *series,
    REAL8                    segmentDuration,
    REAL8                    strideDuration,
    REAL4FFTPlan            *fwdplan,
    int                      whiteSpectrum
    );

int invert_spectrum(
    REAL4FrequencySeries *spectrum,
    REAL8                 dataSampleRate,
    REAL8                 strideDuration,
    REAL8                 truncateDuration,
    REAL8                 lowCutoffFrequency,
    REAL4FFTPlan         *fwdplan,
    REAL4FFTPlan         *revplan
    );

int calibrate_spectrum(
    REAL4FrequencySeries    *spectrum,
    COMPLEX8FrequencySeries *response,
    REAL8                    lowCutoffFrequency,
    int                      inverse
    );

#endif /* INVSPEC_H */
