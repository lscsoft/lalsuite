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
#define EPSEARCHH_ENULLP    1
#define EPSEARCHH_EALOC     3
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
#define EPSEARCHH_EDUTY     144
#define EPSEARCHH_EAMAX     145
#define EPSEARCHH_EE2MS     146
#define EPSEARCHH_ECHNL     147
#define EPSEARCHH_ESIM      148
#define EPSEARCHH_ESPEC     149
#define EPSEARCHH_EWIN      150


#define EPSEARCHH_MSGENULLP    "Null pointer"
#define EPSEARCHH_MSGEALOC     "Memory allocation error"
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
#define EPSEARCHH_MSGEDUTY     "Number of segments sent to slave is zero"
#define EPSEARCHH_MSGEAMAX     "The threshold value of alpha is negative"
#define EPSEARCHH_MSGEE2MS     "Number of events out of range[1..99]"
#define EPSEARCHH_MSGECHNL     "Channel name not valid"
#define EPSEARCHH_MSGESIM      "Invalid simulation type: 0, 1, or 2"
#define EPSEARCHH_MSGESPEC     "Invalid spectrum type: 0 or 1"
#define EPSEARCHH_MSGEWIN      "Invalid window type: 0 or 1"
/******** </lalErrTable> ********/


/* the names of the multidim data channels */
#define INPUTNAME_CHANNEL   "channel"
#define INPUTNAME_SPECTRUM  "spectrum"
#define INPUTNAME_RESPONSE  "response"

typedef struct
tagEPSearchParams
{
	CHAR                  *channelName;
	BOOLEAN                cluster;
	BOOLEAN                printSpectrum;
	INT4                   eventLimit;
	UINT4                  windowLength;
	UINT4                  windowShift;
	REAL8                  alphaThreshold;
	AvgSpecMethod          method;
	CreateTFTilingIn      *tfTilingInput;
	ComputeExcessPowerIn  *compEPInput;
	WindowType             windowType;
}
EPSearchParams;

void
EPSearch (
            LALStatus               *status,
            REAL4TimeSeries         *tseries,
            EPSearchParams          *params,
            SnglBurstTable         **burstEvent
         );

void
EPConditionData(
    LALStatus             *status,
    REAL4TimeSeries       *series,
    REAL4                  flow,
    REAL8                  resampledeltaT,
    ResampleTSFilter       resampleFiltType,
    EPSearchParams        *params
    );

void
LALCountEPEvents (
               LALStatus                               *status,
               INT4                                 *numEvents,
               TFTiling                             *tfTiling,
               REAL8                                alphaThreshold
               );

void
LALTFTileToBurstEvent (
               LALStatus         *status,
               SnglBurstTable    *burstEvent,
               TFTile            *event,
               LIGOTimeGPS       *epoch,
               EPSearchParams    *params
               );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif




