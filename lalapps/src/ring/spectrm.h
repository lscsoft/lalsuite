#ifndef INVSPEC_H
#define INVSPEC_H

/*
 *
 * Routine to compute the average spectrum from time series data.
 * Routine to invert and truncate (to have compact time support) a spectrum.
 * Routine to scale a spectrum by the magnitude of the response function.
 *
 */

#include <lal/LALDatatypes.h>


/* routine to compute an average spectrum from time series data */
REAL4FrequencySeries *compute_average_spectrum(
    REAL4TimeSeries         *series,
    REAL8                    segmentDuration,
    REAL8                    strideDuration,
    REAL4FFTPlan            *fwdplan,
    int                      whiteSpectrum
    );


/* routine to invert and truncate (to have compact time support) a spectrum */
int invert_spectrum(
    REAL4FrequencySeries *spectrum,
    REAL8                 dataSampleRate,
    REAL8                 strideDuration,
    REAL8                 truncateDuration,
    REAL8                 lowCutoffFrequency,
    REAL4FFTPlan         *fwdplan,
    REAL4FFTPlan         *revplan
    );


/* routine to scale a spectrum by the magnitude of the response function */
int calibrate_spectrum(
    REAL4FrequencySeries    *spectrum,
    COMPLEX8FrequencySeries *response,
    REAL8                    lowCutoffFrequency,
    int                      inverse
    );

#endif /* INVSPEC_H */
