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

typedef struct tagStreamPair {
	REAL4TimeSeries *streamOne;
	REAL4TimeSeries *streamTwo;
} StreamPair;

static void parse_options(INT4 argc, CHAR *argv[]);
static void readDataPair(LALStatus *status, StreamPair *streamPair,
		REAL4TimeSeries *seriesOne, REAL4TimeSeries *seriesTwo,
		LIGOTimeGPS start);

/* helper functions */
static REAL4TimeSeries *get_time_series(LALStatus *status, CHAR *ifo,
		CHAR *cacheFile, CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end,
		INT4 buffer);
static REAL4TimeSeries *get_ligo_data(LALStatus *status, FrStream *stream,
		CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end);
static REAL4TimeSeries *get_geo_data(LALStatus *status, FrStream *stream,
		CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end);
static REAL8 delta_gps_to_float(LALStatus *status, LIGOTimeGPS end,
		LIGOTimeGPS start);
static REAL4FrequencySeries *omega_gw(LALStatus *status, REAL4 alpha,
		REAL8 fRef, REAL4 omegaRef, UINT4 length, REAL8 f0, REAL8 deltaF,
    LIGOTimeGPS time);
static REAL4FrequencySeries *overlap_reduction_function(LALStatus *status,
    UINT4 length, REAL8 f0, REAL8 deltaF, INT4 siteOne, INT4 siteTwo,
    LIGOTimeGPS time);
static LIGOTimeGPS increment_gps(LALStatus *status, LIGOTimeGPS time,
		INT4 increment);
static REAL4TimeSeries *cut_time_series(LALStatus *status,
		REAL4TimeSeries *input, LIGOTimeGPS start, LIGOTimeGPS end);

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */
