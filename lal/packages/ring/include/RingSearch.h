#ifndef _RINGSEARCH_H
#define _RINGSEARCH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Ring.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( RINGSEARCHH, "$Id$" );

#define RINGSEARCHH_ENULL 00001
#define RINGSEARCHH_ENNUL 00002
#define RINGSEARCHH_EALOC 00004
#define RINGSEARCHH_ESIZE 00010
#define RINGSEARCHH_ESZMM 00020
#define RINGSEARCHH_ENSEG 00040
#define RINGSEARCHH_EIOPT 00100
#define RINGSEARCHH_EFLOW 00200
#define RINGSEARCHH_EFREQ 00400
#define RINGSEARCHH_EQUAL 01000

#define RINGSEARCHH_MSGENULL "Null pointer"
#define RINGSEARCHH_MSGENNUL "Non-null pointer"
#define RINGSEARCHH_MSGEALOC "Memory allocation error"
#define RINGSEARCHH_MSGESIZE "Invalid segment size"
#define RINGSEARCHH_MSGESZMM "Size mismatch"
#define RINGSEARCHH_MSGENSEG "Non integer number of segments in data"
#define RINGSEARCHH_MSGEIOPT "Invalid option"
#define RINGSEARCHH_MSGEFLOW "Invalid low frequency cutoff"
#define RINGSEARCHH_MSGEFREQ "Invalid bank frequency range"
#define RINGSEARCHH_MSGEQUAL "Invalid bank quality range"

typedef struct
tagRingSearchParams
{
  UINT4                    segmentSize;
  UINT4                    numSegments;
  COMPLEX8FrequencySeries *dataSegment;
  REAL4FrequencySeries    *invSpectrum;
  RealFFTPlan             *forwardPlan;
  RealFFTPlan             *reversePlan;
  REAL4                    dynRangeFac;
  UINT4                    invSpecTrunc;
  REAL4                    lowFrequency;
  REAL4                    minFrequency;
  REAL4                    maxFrequency;
  REAL4                    minQuality;
  REAL4                    maxQuality;
  REAL4                    maxMismatch;
  REAL4                    sampleRate;
  RingTemplateBank        *templateBank;
  UINT4                    templatesSent;
  UINT4                    templatesDone;
  INT4                     searchMaster;
  INT4                     myProcNumber;
  INT4                     numSlaves;
  REAL4                    threshold;
  CHAR                     ifoName[3];
  INT4                     maximizeEvents;
  INT4                     keepResults;
  UINT4                    numResults;
  REAL4TimeSeries         *result;
}
RingSearchParams;

typedef struct
tagRingSearchData
{
  REAL4TimeSeries         *channel;
  REAL4FrequencySeries    *spectrum;
  COMPLEX8FrequencySeries *response;
}
RingSearchData;

typedef struct
tagRingEventList
{
  CHAR  ifoName[3];
  INT8  startTimeNS;
  REAL4 duration;
  REAL4 amplitude;
  REAL4 frequency;
  REAL4 quality;
  REAL4 mass;
  REAL4 snr;
  REAL4 confidence;
  struct tagRingEventList *next;
}
RingEventList;

typedef struct
tagRingSearchInput
{
  UINT4 startTemplate;
  UINT4 templatesToDo;
}
RingSearchInput;

void LALRingSearchInit(
    LALStatus         *status,
    RingSearchParams **searchParams,
    const CHAR       **argv,
    INT4               argc
    );

void
LALRingSearchFini(
    LALStatus         *status,
    RingSearchParams **searchParams
    );

void
LALRingSearchConditionData(
    LALStatus               *status,
    RingSearchParams        *params,
    RingSearchData          *data
    );

void
LALRingSearch(
    LALStatus         *status,
    RingEventList    **output,
    RingSearchInput   *input,
    RingSearchParams  *params
    );

#ifdef __cplusplus
#pragma {
}
#endif 

#endif /* _RINGSEARCH_H */
