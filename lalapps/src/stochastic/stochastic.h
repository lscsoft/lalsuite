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
	UINT8 start;
	INT4 duration;
	INT4 buffer;
} ReadDataPairParams;

static void parseOptions(INT4 argc, CHAR *argv[]);
static void readDataPair(LALStatus *status, StreamPair *streamPair,
		ReadDataPairParams *params);

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTIC_H */
