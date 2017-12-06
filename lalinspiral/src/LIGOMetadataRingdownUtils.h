/*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef _LIGOMETADATARINGDOWNUTILS_H
#define _LIGOMETADATARINGDOWNUTILS_H

#ifdef  __cplusplus
extern "C" {
#endif

#include <lal/LIGOMetadataInspiralUtils.h>

/*
 *
 *  ringdown specific structures
 *
 */
typedef enum
{
  LALRINGDOWN_UNKNOWN_TEST,
  LALRINGDOWN_F_AND_Q,
  LALRINGDOWN_DS_SQ,
  LALRINGDOWN_DS_SQ_FQT
}
SnglRingdownParameterTest;

typedef enum
{
  LALRINGDOWN_RING_INJECT,
  LALRINGDOWN_IMR_INJECT,
  LALRINGDOWN_IMR_RING_INJECT,
  LALRINGDOWN_EOBNR_INJECT,
  LALRINGDOWN_PHENOM_INJECT
}
lalringdown_inject_type;

typedef enum
{
  LALRINGDOWN_SPECTRUM_MEDIAN,
  LALRINGDOWN_SPECTRUM_MEDIAN_MEAN
}
lalringdown_spectrum_type;

typedef enum
{
  LALRINGDOWN_DATATYPE_SIM,
  LALRINGDOWN_DATATYPE_ZERO,
  LALRINGDOWN_DATATYPE_UNCAL,
  LALRINGDOWN_DATATYPE_HT_REAL4,
  LALRINGDOWN_DATATYPE_HT_REAL8
}
lalringdown_data_type;

typedef struct
tagSnglRingdownAccuracy
{
  INT4        match;
  REAL4       kappa;
  INT8        dt;
  REAL4       df;
  REAL4       dQ;
  REAL4       ddeff;
  REAL4       ds_sq;
  SnglRingdownParameterTest test;
}
SnglRingdownAccuracy;


typedef struct
tagRingdownAccuracyList
{
  INT4                      match;
  SnglRingdownParameterTest test;
  SnglRingdownAccuracy      ifoAccuracy[LAL_NUM_IFO];
  INT8                      lightTravelTime[LAL_NUM_IFO][LAL_NUM_IFO];
  REAL8                     minimizerStep;
}
RingdownAccuracyList;

/*
 *  ringdown specific functions
 *
 *
 */

/* sngl ringdown */

void
LALFreeSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead
    );

int
XLALFreeSnglRingdown (
        SnglRingdownTable **eventHead
            );

void
LALSortSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

SnglRingdownTable *
XLALSortSnglRingdown (
    SnglRingdownTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

int
LALCompareSnglRingdownByTime (
    const void *a,
    const void *b
    );

void
LALCompareRingdowns (
    LALStatus                *status,
    SnglRingdownTable        *aPtr,
    SnglRingdownTable        *bPtr,
    RingdownAccuracyList     *params
    );

REAL8
XLALCompareRingdowns (
    SnglRingdownTable        *aPtr,
    SnglRingdownTable        *bPtr,
    RingdownAccuracyList     *params
    );

void
LALClusterSnglRingdownTable (
    LALStatus                  *status,
    SnglRingdownTable          *ringdownEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    );

SnglRingdownTable *
XLALVetoSingleRingdown (
    SnglRingdownTable          *eventHead,
    LALSegList                 *vetoSegs,
    CHAR                       *ifo
    );

int
XLALCoincRingdownIfosDiscard(
    CoincRingdownTable **coincHead,
    char                *ifos
    );

void
LALIfoCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    CHAR                       *ifo
    );

SnglRingdownTable *
XLALIfoCutSingleRingdown(
    SnglRingdownTable         **eventHead,
    char                       *ifo
    );

SnglRingdownTable *
XLALTimeCutSingleRingdown(
    SnglRingdownTable          *eventList,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    );

void
LALTimeCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    );

void
LALIfoCountSingleRingdown(
    LALStatus                  *status,
    UINT4                      *numTrigs,
    SnglRingdownTable          *input,
    InterferometerNumber        ifoNumber
    );


void
LALTimeSlideSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable          *triggerList,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime,
    LIGOTimeGPS                 slideTimes[LAL_NUM_IFO]
    );

SnglRingdownTable *
XLALPlayTestSingleRingdown(
    SnglRingdownTable          *eventHead,
    LALPlaygroundDataMask      *dataType
    );

void
LALPlayTestSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    );

int
XLALMaxSnglRingdownOverIntervals(
    SnglRingdownTable         **eventHead,
    INT8                       deltaT
    );

INT4
XLALCountSnglRingdown(
    SnglRingdownTable *head
    );

/* sim ringdown */
void
XLALSortSimRingdown(
    SimRingdownTable **head,
    int (*comparefunc)(const SimRingdownTable * const *,
      const SimRingdownTable * const *)
    );

int
XLALCompareSimRingdownByGeocentStartTime(
    const SimRingdownTable * const *a,
    const SimRingdownTable * const *b
    );

int
XLALFreeSimRingdown (
    SimRingdownTable **eventHead
    );

void
XLALPlayTestSimRingdown(
     SimRingdownTable          **eventHead,
     LALPlaygroundDataMask      *dataType
     );


int
XLALSimRingdownInSearchedData(
    SimRingdownTable         **eventHead,
    SearchSummaryTable       **summList
    );

INT8
XLALReturnSimRingdownStartTime (
    SimRingdownTable *event,
    CHAR             *ifo
    );

int
XLALSnglSimRingdownTest (
    SimRingdownTable  **simHead,
    SnglRingdownTable **eventHead,
    SimRingdownTable  **missedSimHead,
    SnglRingdownTable **missedSnglHead,
    INT8                injectWindowNS
    );

int
XLALCoincSimRingdownTest (
    SimRingdownTable   **simHead,
    CoincRingdownTable **coincHead,
    SimRingdownTable   **missedSimHead,
    CoincRingdownTable **missedCoincHead
    );

REAL8
XLAL2DRinca(
    SnglRingdownTable         *aPtr,
    SnglRingdownTable         *bPtr
    );

REAL8
XLAL3DRinca(
    SnglRingdownTable         *aPtr,
    SnglRingdownTable         *bPtr
    );

/* coinc ringdown */

void
LALCreateTwoIFORingdownCoincList(
    LALStatus                  *status,
    CoincRingdownTable        **coincOutput,
    SnglRingdownTable          *snglInput,
    RingdownAccuracyList       *accuracyParams
    );

void
LALCreateNIFORingdownCoincList(
    LALStatus                  *status,
    CoincRingdownTable        **coincHead,
    RingdownAccuracyList       *accuracyParams,
    INT4                        N
    );

void
LALRemoveRepeatedRingdownCoincs(
    LALStatus                  *status,
    CoincRingdownTable        **coincHead
    );

void
LALFreeCoincRingdown(
    LALStatus                  *status,
    CoincRingdownTable        **coincPtr
    );

void
XLALFreeCoincRingdown(
    CoincRingdownTable        **coincPtr
    );

void
LALAddSnglRingdownToCoinc(
    LALStatus                  *status,
    CoincRingdownTable        **coincPtr,
    SnglRingdownTable          *snglRingdown
    );

CoincRingdownTable *
XLALAddSnglRingdownToCoinc(
    CoincRingdownTable         *coincRingdown,
    SnglRingdownTable          *snglRingdown
    );

REAL8
XLALSnglRingdownCoincTest(
    CoincRingdownTable         *coincRingdown,
    SnglRingdownTable          *snglRingdown,
    RingdownAccuracyList       *accuracyParams
    );

CoincRingdownTable *
XLALRingdownDistanceCut(
    CoincRingdownTable        **coincRingdown,
    RingdownAccuracyList       *accuracyParams
    );

SnglRingdownTable *
XLALExtractSnglRingdownFromCoinc(
    CoincRingdownTable         *coincRingdown,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    );

void
LALExtractSnglRingdownFromCoinc(
    LALStatus                  *status,
    SnglRingdownTable         **snglPtr,
    CoincRingdownTable         *coincRingdown,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    );

int
XLALRecreateRingdownCoincFromSngls(
    CoincRingdownTable        **coincPtr,
    SnglRingdownTable          *snglRingdown
    );


INT8
XLALCoincRingdownTimeNS (
    const CoincRingdownTable         *coincRingdown
    );

REAL4
XLALCoincRingdownStat(
    CoincRingdownTable         *coincRingdown,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    );

int
XLALClusterCoincRingdownTable (
    CoincRingdownTable        **coincList,
    INT8                        dtimeNS,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams   *bittenLParams
    );

int
XLALCompareCoincRingdownByTime (
    const void *a,
    const void *b
    );

CoincRingdownTable *
XLALSortCoincRingdown (
    CoincRingdownTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

int
XLALCoincRingdownIfos (
    CoincRingdownTable  *coincRingdown,
    char                *ifos
    );

int
XLALCoincRingdownIfosCut(
    CoincRingdownTable **coincHead,
    char                *ifos
     );

UINT8
XLALCoincRingdownIdNumber (
    CoincRingdownTable  *coincRingdown
    );

CoincRingdownTable *
XLALCoincRingdownSlideCut(
    CoincRingdownTable **coincHead,
    int                  slideNum
    );

CoincRingdownTable *
XLALStatCutCoincRingdown (
    CoincRingdownTable         *eventHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       statCut
    );

SnglRingdownTable *
XLALCompleteCoincRingdown (
    CoincRingdownTable         *eventHead,
    int                         ifoList[LAL_NUM_IFO]
    );

CoincRingdownTable *
XLALPlayTestCoincRingdown(
    CoincRingdownTable         *eventHead,
    LALPlaygroundDataMask      *dataType
    );

INT4
XLALCountCoincRingdown(
    CoincRingdownTable *head
    );

void
LALRingdownH1H2Consistency(
    LALStatus                  *status,
    CoincRingdownTable        **coincRingdown,
    REAL4                       H2snrCutThreshold,
    LALSegList                 *vetoSegsH1,
    LALSegList                 *vetoSegsH2
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOMETADATARINGDOWNUTILS_H */
