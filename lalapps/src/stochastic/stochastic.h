/*
 * stochasticPipeline.h - SGWB Standalone Analysis Pipeline
 *                      - header file
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

NRCSID (STOCHASTICH, "$Id$" );

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

void parseOptions(INT4 argc, CHAR *argv[]);
void displayUsage(INT4 exitcode);
void readDataPair(LALStatus *status, StreamPair *streamPair,
		ReadDataPairParams *params);

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */
