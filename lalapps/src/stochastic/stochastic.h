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
	REAL4TimeSeries *stream1;
	REAL4TimeSeries *stream2;
} StreamPair;

typedef struct tagReadDataPairParams {
	UINT8 startTime;
	INT4 duration;
	CHAR *frameCache1;
	CHAR *frameCache2;
	CHAR *ifo1;
	CHAR *ifo2;
	CHAR *channel1;
	CHAR *channel2;
	INT4 buffer;
	INT4 sampleRate;
	INT4 resampleRate;
} ReadDataPairParams;


void readDataPair(LALStatus *status, StreamPair *streamPair,
		ReadDataPairParams *params);


#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */
