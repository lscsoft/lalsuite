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

typedef struct tagReadDataPairParams {
	UINT8 start;
	INT4 duration;
	INT4 buffer;
} ReadDataPairParams;

static void parse_options(INT4 argc, CHAR *argv[]);
static void readDataPair(LALStatus *status, StreamPair *streamPair,
		ReadDataPairParams *params);
static void adam_readDataPair(LALStatus *status, StreamPair *streamPair,
		ReadDataPairParams *params);
static REAL4TimeSeries *get_time_series(LALStatus *status, CHAR *ifo,
		CHAR *cacheFile, CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end);
static REAL4TimeSeries *get_ligo_data(LALStatus *status, FrStream *stream,
		CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end);
static REAL4TimeSeries *get_geo_data(LALStatus *status, FrStream *stream,
		CHAR *channel, LIGOTimeGPS start, LIGOTimeGPS end);
static REAL8 DeltaGPStoFloat(LALStatus *status, LIGOTimeGPS *end,
		LIGOTimeGPS *start);

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */
