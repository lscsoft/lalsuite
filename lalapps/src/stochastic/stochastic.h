/*
 * stochasticPipeline.h - SGWB Standalone Analysis Pipeline
 *                      - header file
 * Adam Mercer <ram@star.sr.bham.ac.uk>
 * Tania Regimbau <Tania.Regimbau@astro.cf.ac.uk>
 *
 * $Id$
 */

#ifndef _PIPELINE_H
#define _PIPELINE_H

#include <lal/Calibration.h>
#include <lal/DetectorSite.h>
#include <lal/LALStdlib.h>
#include <lal/SimulateSB.h>
#include <lal/Units.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (PIPELINEH, "$Id$" );

typedef struct tagStreamPair {
	REAL4TimeSeries *streamOne;
	REAL4TimeSeries *streamTwo;
} StreamPair;

typedef struct tagReadDataPairParams {
	UINT8 startTime;
	INT4 duration;
	CHAR *dataCacheOne;
	CHAR *dataCacheTwo;
	CHAR *ifoOne;
	CHAR *ifoTwo;
	CHAR *channelOne;
	CHAR *channelTwo;
	INT4 buffer;
	INT4 sampleRate;
	INT4 resampleRate;
} ReadDataPairParams;

typedef struct tagMonteCarloInput {
	CHAR *ifoOne;
	CHAR *ifoTwo;
	CHAR *calCacheOne;
	CHAR *calCacheTwo;
} MonteCarloInput;

typedef struct tagMonteCarloParams {
	UINT4 lengthSegment;
	UINT4 numSegment;
	UINT4 sampleRate;
	UINT8 startTime;
	INT4  seed;
	REAL8 fRef;
	REAL8 f0;
	REAL8 omegaRef;
	REAL8 alpha;
	INT4  siteOne;
	INT4  siteTwo;
} MonteCarloParams;

void parseOptions(INT4 argc, CHAR *argv[]);
void displayUsage(INT4 exitcode);
void readDataPair(LALStatus *status, StreamPair *streamPair,
		ReadDataPairParams *params);
void monteCarlo (LALStatus *status, SSSimStochBGOutput *MCoutput,
		MonteCarloInput  *MCinput, MonteCarloParams *MCparams);
void monteCarloSplice (LALStatus *status, SSSimStochBGOutput *MCoutput,
		MonteCarloInput  *MCinput, MonteCarloParams *MCparams);
void SinusoidalSplice(REAL4Vector **longData, REAL4Vector **shortData,
		REAL4Vector *output, UINT4 nSpliceSegs, UINT4 offset);

#ifdef  __cplusplus
}
#endif

#endif /* _PIPELINE_H */
