/*
 * stochastic_dev.h - SGWB Standalone Analysis Pipeline
 *                  - header file
 *                  - development branch
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

NRCSID (STOCHASTICDEVH, "$Id$" );

typedef struct tagStreamPair {
	REAL4TimeSeries *streamOne;
	REAL4TimeSeries *streamTwo;
} StreamPair;

typedef struct tagReadDataPairParams {
	UINT8 startTime;
	INT4 duration;
	CHAR *frameCacheOne;
	CHAR *frameCacheTwo;
	CHAR *ifoOne;
	CHAR *ifoTwo;
	CHAR *channelOne;
	CHAR *channelTwo;
	INT4 buffer;
	INT4 sampleRate;
	INT4 resampleRate;
} ReadDataPairParams;

typedef struct tagPSDEstimatorInput {
	REAL4TimeSeries *segmentA;
	REAL4TimeSeries *segmentC;
	COMPLEX8FrequencySeries *responseA;
	COMPLEX8FrequencySeries *responseC;
} PSDEstimatorInput;

typedef struct tagPSDEstimatorParams {
	INT4 psdTempLength;
	INT4 filterLength;
	INT4 numFMin;
	AverageSpectrumParams *psdParams;
} PSDEstimatorParams;

static void parseOptions(INT4 argc, CHAR *argv[]);
static void readDataPair(LALStatus *status, StreamPair *streamPair, \
		ReadDataPairParams *params);
static void psdEstimator(LALStatus *status, REAL4FrequencySeries *output, \
    PSDEstimatorInput input, PSDEstimatorParams params);

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */

/*
 * vim: et
 */
