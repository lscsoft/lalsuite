/********************************** <lalVerbatim file="ExcessPowerHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _EPSEARCH_H
#define _EPSEARCH_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/TFTransform.h>
#include <lal/ExcessPower.h>
#include <lal/BurstSearch.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/IIRFilter.h>
#include <lal/LALRCSID.h>
#include <lal/ResampleTimeSeries.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID (EPSEARCHH, "$Id$");

/******** <lalErrTable file="EPSearchHErrTab"> ********/
#define EPSEARCHH_ENULLP       1
#define EPSEARCHH_EPOSARG      2
#define EPSEARCHH_EALOC        3
#define EPSEARCHH_EPOW2        4
#define EPSEARCHH_EMALLOC      8
#define EPSEARCHH_EINCOMP      16
#define EPSEARCHH_EORDER       32
#define EPSEARCHH_ENONNULL     64
#define EPSEARCHH_ETILES       65
#define EPSEARCHH_EDELFT       128
#define EPSEARCHH_EARGS     129
#define EPSEARCHH_ENUMZ     130
#define EPSEARCHH_ESEGZ     131
#define EPSEARCHH_EOVLP     132
#define EPSEARCHH_EOVPF     133
#define EPSEARCHH_EMFBZ     134
#define EPSEARCHH_EMTBZ     135
#define EPSEARCHH_EFLOW     136
#define EPSEARCHH_EDELF     137
#define EPSEARCHH_ELTFZ     138
#define EPSEARCHH_ESIGM     139
#define EPSEARCHH_EALPH     140
#define EPSEARCHH_EFREE     141
#define EPSEARCHH_ENLAL     142
#define EPSEARCHH_ENDAT     143
#define EPSEARCHH_EDUTY     144
#define EPSEARCHH_EAMAX     145
#define EPSEARCHH_EE2MS     146
#define EPSEARCHH_ECHNL     147
#define EPSEARCHH_ESIM      148
#define EPSEARCHH_ESPEC     149
#define EPSEARCHH_EWIN      150
#define EPSEARCHH_EDATZ      151




#define EPSEARCHH_MSGENULLP    "Null pointer"
#define EPSEARCHH_MSGEPOSARG   "Arguments must be non-negative"
#define EPSEARCHH_MSGEALOC     "Memory allocation error"
#define EPSEARCHH_MSGEPOW2     "Length of supplied data must be a power of 2"
#define EPSEARCHH_MSGEMALLOC   "Malloc failure"
#define EPSEARCHH_MSGEINCOMP   "Incompatible arguments"
#define EPSEARCHH_MSGEORDER    "Routines called in illegal order"
#define EPSEARCHH_MSGENONNULL  "Null pointer expected"
#define EPSEARCHH_MSGETILES    "Malloc failed while assigning memory for a tile"
#define EPSEARCHH_MSGEDELFT    "Inconsistent deltaF in spectrum and data"
#define EPSEARCHH_MSGEARGS     "Wrong number of arguments"
#define EPSEARCHH_MSGENUMZ     "Data segment length is zero or negative"
#define EPSEARCHH_MSGESEGZ     "Number of data segments is zero or negative"
#define EPSEARCHH_MSGEOVLP     "Overlap of data segments is negative"
#define EPSEARCHH_MSGEOVPF     "Overlap of TF tiles is negative"
#define EPSEARCHH_MSGEMFBZ     "Smallest extent in freq of TF tile is zero" 
#define EPSEARCHH_MSGEMTBZ     "Smallest extent in time of TF tile is zero"
#define EPSEARCHH_MSGEFLOW     "Lowest frequency to be searched is <= 0"
#define EPSEARCHH_MSGEDELF     "Freq resolution of 1st time-freq plane is <=0"
#define EPSEARCHH_MSGELTFZ     "Length (Nf) of 1st TF plane (with Nt=1) is <=0"
#define EPSEARCHH_MSGESIGM     "Threshold number of sigma is <=1"
#define EPSEARCHH_MSGEALPH     "Default alpha value is out of range"
#define EPSEARCHH_MSGEFREE     "Memory free error"
#define EPSEARCHH_MSGENLAL     "Tried to allocate to non-null pointer"
#define EPSEARCHH_MSGENDAT     "No data read"
#define EPSEARCHH_MSGEDUTY     "Number of segments sent to slave is zero"
#define EPSEARCHH_MSGEAMAX     "The threshold value of alpha is negative"
#define EPSEARCHH_MSGEE2MS     "Number of events out of range[1..99]"
#define EPSEARCHH_MSGECHNL     "Channel name not valid"
#define EPSEARCHH_MSGESIM      "Invalid simulation type: 0, 1, or 2"
#define EPSEARCHH_MSGESPEC     "Invalid spectrum type: 0 or 1"
#define EPSEARCHH_MSGEWIN      "Invalid window type: 0 or 1"
#define EPSEARCHH_MSGEDATZ   "Got less data than expected"
/******** </lalErrTable> ********/


/* the names of the multidim data channels */
#define INPUTNAME_CHANNEL   "channel"
#define INPUTNAME_SPECTRUM  "spectrum"
#define INPUTNAME_RESPONSE  "response"



typedef struct
tagEPDataSegment
{
  REAL4TimeSeries              *data;
  REAL4FrequencySeries         *spec;
  COMPLEX8FrequencySeries      *resp;
  INT4                          endOfData;
  INT4                          newLock;
  INT4                          newCal;
  INT4                          number;
  REAL4IIRFilter               *preprocessing_filter;
}
EPDataSegment;

typedef struct
tagEPDataSegmentVector
{
  UINT4                         length;
  EPDataSegment                  *data;
}
EPDataSegmentVector;

typedef struct
tagEPInitParams
{
  UINT4                         numPoints;
  UINT4                         numSegments;
  UINT4                         segDutyCycle;
  AvgSpecMethod                 method;
}
EPInitParams;

typedef struct
tagEPSearchParams
{
  BOOLEAN                       searchMaster;
  BOOLEAN                       haveData;
  INT4                          cluster;
  UINT4                        *numSlaves;          
  UINT4                         simType;
  UINT4                         currentSegment;
  INT4                          ntotT;
  INT4                          numEvents;
  UINT4                         ovrlap;
  REAL8                         lambda;
  REAL8                         alphaThreshold;
  INT4                          events2Master;     /* Added Erik Kats */
  CHAR                         *channelName;       /* Added Erik Kats */
  EPInitParams                 *initParams;
  EPDataSegmentVector          *epSegVec;
  CreateTFTilingIn             *tfTilingInput;
  TFTiling                     *tfTiling;
  ComputeExcessPowerIn         *compEPInput;
  LALWindowParams               winParams;
  BOOLEAN                       printSpectrum;
}
EPSearchParams;

void
EPSearch (
            LALStatus               *status,
            REAL4TimeSeries         *tseries,
            EPSearchParams          *params,
            SnglBurstTable         **burstEvent
         );

void EPInitSearch(
        LALStatus             *status,
        void                 **searchParams,
        CHAR                  *argv[],
        INT4                   argc
        );

void
EPConditionData(
    LALStatus             *status,
    REAL4TimeSeries       *series,
    REAL4                  flow,
    REAL8                  resampledeltaT,
    ResampleTSFilter       resampleFiltType,
    void                  *searchParams
    );

void
EPFinalizeSearch(
    LALStatus             *status,
    void                 **searchParams
    );

void
LALTFTileToBurstEvent (
               LALStatus                            *status,
               SnglBurstTable                          *burstEvent,
               TFTile                               *event,
               INT8                                  tstart,
               EPSearchParams                       *params
               );

void LALWeighTFTileList (
        LALStatus         *status,
        TFTiling          *tfTiling,
        INT4               maxDOF
        );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif




