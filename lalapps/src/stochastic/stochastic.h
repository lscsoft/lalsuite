/*
 * stochastic.h - SGWB Standalone Analysis Pipeline
 *              - header file
 *
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 * Tania Regimbau <Tania.Regimbau@astro.cf.ac.uk>
 *
 * $Id$
 */

#ifndef _STOCHASTIC_H
#define _STOCHASTIC_H

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID(STOCHASTICH, "$Id$");

/* helper functions */
void parse_options(INT4 argc, CHAR *argv[]);
REAL4TimeSeries *get_time_series(LALStatus *status, CHAR *ifo,
    CHAR *cacheFile, CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end,
    INT4 buffer);
REAL4TimeSeries *get_ligo_data(LALStatus *status, FrStream *stream,
    CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end);
REAL4TimeSeries *get_geo_data(LALStatus *status, FrStream *stream,
    CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end);
REAL8 delta_gps_to_float(LALStatus *status, LIGOTimeGPS end,
    LIGOTimeGPS start);
REAL4FrequencySeries *omega_gw(LALStatus *status, REAL4 alpha,
    REAL8 fRef, REAL4 omegaRef, UINT4 length, REAL8 f0, REAL8 deltaF,
    LIGOTimeGPS time);
REAL4FrequencySeries *overlap_reduction_function(LALStatus *status,
    UINT4 length, REAL8 f0, REAL8 deltaF, INT4 siteOne, INT4 siteTwo,
    LIGOTimeGPS time);
LIGOTimeGPS increment_gps(LALStatus *status, LIGOTimeGPS time,
    INT4 increment);
REAL4TimeSeries *cut_time_series(LALStatus *status,
    REAL4TimeSeries *input, LIGOTimeGPS start, LIGOTimeGPS end);
void write_ccspectra_frame(COMPLEX8FrequencySeries *series,
    CHAR *ifoOne, CHAR *ifoTwo, LIGOTimeGPS time, INT4 duration);
REAL4TimeSeries *rectangular_window(LALStatus *status, REAL8 deltaT,
    REAL8 f0, INT4 length);
REAL4TimeSeries *hann_window(LALStatus *status, REAL8 deltaT,
    REAL8 f0, INT4 length);
REAL4TimeSeries *data_window(LALStatus *status, REAL8 deltaT,
    REAL8 f0, INT4 length, INT4 hannDuration);
REAL4FrequencySeries *inverse_noise(LALStatus *status,
    REAL4FrequencySeries *psd, COMPLEX8FrequencySeries *response);
REAL4FrequencySeries *optimal_filter(LALStatus *status,
    REAL4FrequencySeries *overlap, REAL4FrequencySeries *omega,
    REAL4FrequencySeries *psdOne, REAL4FrequencySeries *psdTwo,
    REAL4WithUnits normalisation);

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */

/*
 * vim: et
 */
