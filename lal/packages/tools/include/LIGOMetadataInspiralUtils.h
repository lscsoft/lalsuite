/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Duncan Brown, Eirini Messaritaki, Gareth Jones, George Birthisel, Jolien Creighton, Kipp Cannon, Lisa M. Goggin, Patrick Brady, Peter Shawhan, Saikat Ray-Majumder, Stephen Fairhurst, Xavier Siemens, Craig Robinson , Thomas Cokelaer
*
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

#ifndef _LIGOMETADATAINSPIRALUTILS_H
#define _LIGOMETADATAINSPIRALUTILS_H

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

#include <lal/LIGOMetadataUtils.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>

/* SnglInspiralParameterTest contains an enum type for each of the
 * tests of mass parameters which are used */
typedef enum
{
  unspecified_test,
  no_test,
  m1_and_m2,
  psi0_and_psi3,
  mchirp_and_eta,
  tau0_and_tau3,
  ellipsoid
}
SnglInspiralParameterTest;

/*
 * The SnglInspiralAccuracy structure contains parameters used for
 * testing coincidence between two or more single inspiral tables.
 * These include a timing accuracy dt, five mass accuracies dm
 * (used for testing mass1 and mass2), deta, dmchirp, dpsi0, and
 * dpsi3. It also includes the parameters kappa and epsilon which
 * are used for testing consistency of effective distance. */
typedef struct
tagSnglInspiralAccuracy
{
  INT4 match;
  REAL4 epsilon;
  REAL4 kappa;
  INT8 dt;
  REAL4 dm;
  REAL4 deta;
  REAL4 dmchirp;
  REAL4 dmchirpHi;
  REAL4 highMass;
  REAL4 dpsi0;
  REAL4 dpsi3;
  REAL4 dtau0;
  REAL4 dtau3;
  SnglInspiralParameterTest test;
  INT4 exttrig;
}
SnglInspiralAccuracy;


/* The InspiralAccuracyList structure contains parameter accuracies
 * for each of the six global interferometers. These are stored in
 * the SnglInspiralAccuracy structure. The accuracies stored should
 * be the accuracy with which each instrument can determine the given
 * parameter. It also contains a \texttt{match} which is set to 1 to
 * signify that coincidence criteria are satisfied and 0 otherwise.
 * Finally, the SnglInspiralParameterTest must be specified. */
typedef struct
tagInspiralAccuracyList
{
  INT4 match;
  SnglInspiralParameterTest test;
  SnglInspiralAccuracy ifoAccuracy[LAL_NUM_IFO];
  INT8 lightTravelTime[LAL_NUM_IFO][LAL_NUM_IFO];
  REAL4 iotaCutH1H2;
  REAL4 iotaCutH1L1;
  REAL8 eMatch;
  INT4 exttrig;
}
InspiralAccuracyList;

/* The CoincInspiralStatParams structure contains the bitten L
 * parameter for each of the six global interferometers. These are
 * stored in the param_a and param_b structure */
typedef struct
tagCoincInspiralStatParams
{
  REAL4    param_a[LAL_NUM_IFO];
  REAL4    param_b[LAL_NUM_IFO];
  REAL4    eff_snr_denom_fac;
}
CoincInspiralStatParams;

/* The SnglInspiralClusterChoice provides three choices for clustering
 * a single inspiral table. The snr clustering returns the trigger
 * with the greatest signal to noise ratio; snr_and_chisq replaces the
 * existing trigger if the new trigger has both a greater snr and a
 * smaller chi squared value; snrsq_over_chisq selects the trigger
 * with the largest value of snr squared divided by the chi squared. */
typedef enum
{
  none,
  snr_and_chisq,
  snrsq_over_chisq,
  snr
}
SnglInspiralClusterChoice;

/* The CoincInspiralStatistic provides two choices for clustering a
 * single inspiral table. The snrsq clustering returns the trigger
 * with the greatest summed snr^2 from all instruments.
 * The snr_chi_stat replaces selects the trigger with the largest
 * value of the snr and chisq statistic and the bitten_l returns
 * the minimum among the summed snr^2 from all instruments and
 * the a * snr_i - b in each detector. The parameters a and b must
 * be provided by the user. */
typedef enum
{
  no_stat,
  snrsq,
  effective_snrsq,
  s3_snr_chi_stat,
  bitten_l,
  bitten_lsq,
  ifar
}
CoincInspiralStatistic;


/* The SnglInspiralBCVCalphafCut provides entries for cutting single
 * IFO triggers generated with the BCVC code. For each LSC IFO there
 * is a field lo and hi which corresponds to the area allowing
 * triggers. */
typedef struct
tagSnglInspiralBCVCalphafCut
{
  REAL4       h1_lo;
  REAL4       h1_hi;
  REAL4       h2_lo;
  REAL4       h2_hi;
  REAL4       l1_lo;
  REAL4       l1_hi;
  REAL4       psi0cut;
}
SnglInspiralBCVCalphafCut;

/*
 *
 *  ringdown specific structures
 *
 */
typedef enum
{
  unknown_test,
  f_and_Q,
  ds_sq,
  ds_sq_fQt
}
SnglRingdownParameterTest;

typedef enum
{
  ring_inject,
  imr_inject,
  imr_ring_inject,
  EOBNR_inject,
  Phenom_inject
}
inject_type;

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
  SnglInspiralParameterTest test;
  SnglRingdownAccuracy      ifoAccuracy[LAL_NUM_IFO];
  INT8                      lightTravelTime[LAL_NUM_IFO][LAL_NUM_IFO];
  REAL8                     minimizerStep;
}
RingdownAccuracyList;

/*
 *
 * inspiral specific functions
 *
 */

/* sngl inspiral */
void
LALFreeSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead
    );

int
XLALFreeSnglInspiral (
    SnglInspiralTable **eventHead
    );

void
LALSortSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

SnglInspiralTable *
XLALSortSnglInspiral (
    SnglInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

int
LALCompareSnglInspiralByMass (
    const void *a,
    const void *b
    );

int
LALCompareSnglInspiralByPsi (
    const void *a,
    const void *b
    );

int
LALCompareSnglInspiralByTime (
    const void *a,
    const void *b
    );

int
LALCompareSnglInspiralByID (
     const void *a,
     const void *b
     );

void
LALCompareSnglInspiral (
    LALStatus                *status,
    SnglInspiralTable        *aPtr,
    SnglInspiralTable        *bPtr,
    SnglInspiralAccuracy     *params
    );

int
XLALCompareInspirals (
    SnglInspiralTable        *aPtr,
    SnglInspiralTable        *bPtr,
    InspiralAccuracyList     *params
    );

void
LALCompareInspirals (
    LALStatus                *status,
    SnglInspiralTable        *aPtr,
    SnglInspiralTable        *bPtr,
    InspiralAccuracyList     *params
    );

void
LALClusterSnglInspiralTable (
    LALStatus                  *status,
    SnglInspiralTable         **inspiralEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    );

REAL4
XLALSnglInspiralStat(
    SnglInspiralTable         *snglInspiral,
    SnglInspiralClusterChoice  snglStat
    );

int
XLALClusterSnglInspiralTable (
    SnglInspiralTable         **inspiralList,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    );

SnglInspiralTable *
XLALTimeCutSingleInspiral(
    SnglInspiralTable          *eventList,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    );

void
LALTimeCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    );


void
LALSNRCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    REAL4                       snrCut
    );

SnglInspiralTable *
XLALSNRCutSingleInspiral(
    SnglInspiralTable          *eventHead,
    REAL4                       snrCut
    );

SnglInspiralTable *
XLALRsqCutSingleInspiral (
    SnglInspiralTable          *eventHead,
    REAL4                       rsqVetoTimeThresh,
    REAL4                       rsqMaxSnr,
    REAL4                       rsqAboveSnrCoeff,
    REAL4                       rsqAboveSnrPow
    );

SnglInspiralTable *
XLALVetoSingleInspiral (
    SnglInspiralTable          *eventHead,
    LALSegList                 *vetoSegs,
    const CHAR                 *ifo
    );

SnglInspiralTable *
XLALalphaFTemp (
    SnglInspiralTable          *eventHead
    );

void
LALBCVCVetoSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    SnglInspiralBCVCalphafCut   alphafParams
    );

void
LALalphaFCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    REAL4                       alphaFhi,
    REAL4                       alphaFlo
    );

void
LALIfoCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    CHAR                       *ifo
    );

SnglInspiralTable *
XLALIfoCutSingleInspiral(
    SnglInspiralTable         **eventHead,
    char                       *ifo
    );

void
LALIfoCountSingleInspiral(
    LALStatus                  *status,
    UINT4                      *numTrigs,
    SnglInspiralTable          *input,
    InterferometerNumber        ifoNumber
    );

int
XLALTimeSlideSegList(
    LALSegList                 *seglist,
    const LIGOTimeGPS          *ringStartTime,
    const LIGOTimeGPS          *ringEndTime,
    const LIGOTimeGPS          *slideTimes
    );

void
XLALTimeSlideSingleInspiral(
    SnglInspiralTable          *input,
    const LIGOTimeGPS          *ringStartTime,
    const LIGOTimeGPS          *ringEndTime,
    const LIGOTimeGPS           slideTimes[]
    );

SnglInspiralTable *
XLALPlayTestSingleInspiral(
    SnglInspiralTable          *eventHead,
    LALPlaygroundDataMask      *dataType
    );

void
LALPlayTestSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    );


void
LALCreateTrigBank(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    SnglInspiralParameterTest  *test
    );

void
LALIncaCoincidenceTest(
    LALStatus                  *status,
    SnglInspiralTable         **ifoAOutput,
    SnglInspiralTable         **ifoBOutput,
    SnglInspiralTable          *ifoAInput,
    SnglInspiralTable          *ifoBInput,
    SnglInspiralAccuracy       *errorParams
    );

void
LALTamaCoincidenceTest(
    LALStatus                  *status,
    SnglInspiralTable         **ifoAOutput,
    SnglInspiralTable         **ifoBOutput,
    SnglInspiralTable          *ifoAInput,
    SnglInspiralTable          *ifoBInput,
    SnglInspiralAccuracy       *errorParams,
    SnglInspiralClusterChoice   clusterchoice
    );

int
XLALMaxSnglInspiralOverIntervals(
    SnglInspiralTable         **eventHead,
    INT4                       deltaT
    );

INT4
XLALCountSnglInspiral(
    SnglInspiralTable *head
    );

INT4
XLALCountCoincInspiral(
    CoincInspiralTable *head
    );

/* coinc inspiral */
void
LALCreateTwoIFOCoincList(
    LALStatus                  *status,
    CoincInspiralTable        **coincOutput,
    SnglInspiralTable          *snglInput,
    InspiralAccuracyList       *accuracyParams
    );

void
LALCreateNIFOCoincList(
    LALStatus                  *status,
    CoincInspiralTable        **coincHead,
    InspiralAccuracyList       *accuracyParams,
    INT4                        N
    );

void
LALRemoveRepeatedCoincs(
    LALStatus                  *status,
    CoincInspiralTable        **coincHead
    );

void
LALFreeCoincInspiral(
    LALStatus                  *status,
    CoincInspiralTable        **coincPtr
    );

int
XLALFreeCoincInspiral(
    CoincInspiralTable        **coincPtr
    );


void
LALAddSnglInspiralToCoinc(
    LALStatus                  *status,
    CoincInspiralTable        **coincPtr,
    SnglInspiralTable          *snglInspiral
    );

CoincInspiralTable *
XLALAddSnglInspiralToCoinc(
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral
    );

void
XLALSnglInspiralCoincTest(
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    );

void
LALSnglInspiralCoincTest(
    LALStatus                  *status,
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    );

void
XLALInspiralPsi0Psi3CutBCVC(
    CoincInspiralTable        **coincInspiral
    );

void
XLALInspiralIotaCutBCVC(
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams
    );

void
LALInspiralDistanceCutCleaning(
    LALStatus                  *status,
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams,
    REAL4 			snrThreshold,
    SummValueTable             *summValueList,
    LALSegList                 *vetoSegsH1,
    LALSegList                 *vetoSegsH2
    );

void
XLALInspiralDistanceCutBCVC(
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams
    );

void
XLALInspiralDistanceCut(
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams
    );

void
XLALInspiralSNRCutBCV2(
    CoincInspiralTable        **coincInspiral
    );

SnglInspiralTable *
XLALExtractSnglInspiralFromCoinc(
    CoincInspiralTable         *coincInspiral,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    );

void
LALExtractSnglInspiralFromCoinc(
    LALStatus                  *status,
    SnglInspiralTable         **snglPtr,
    CoincInspiralTable         *coincInspiral,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    );

int
XLALCreateCoincSlideTable(
    CoincInspiralSlideTable   **slideTableHead,
    INT4                        numSlides
    );

REAL4
XLALSetupCoincSlideTable(
    CoincInspiralSlideTable    *slideTableHead,
    CoincInspiralTable         *coincSlideHead,
    char                       *timeAnalyzedFileName,
    REAL4                       timeModifier,
    INT4                        numSlides
    );

int
XLALRecreateCoincFromSngls(
    CoincInspiralTable        **coincPtr,
    SnglInspiralTable         **snglInspiral
    );

void
LALCoincCutSnglInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **snglPtr
    );

int
XLALGenerateCoherentBank(
    SnglInspiralTable         **coherentBank,
    CoincInspiralTable         *coincInput,
    CHAR                       *ifos
    );

INT8
XLALCoincInspiralTimeNS (
    const CoincInspiralTable         *coincInspiral
    );

REAL4
XLALCoincInspiralStat(
    const CoincInspiralTable   *coincInspiral,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    );

int XLALComputeAndStoreEffectiveSNR(   	CoincInspiralTable *head,
					CoincInspiralStatistic *stat,
					CoincInspiralStatParams *par
    );

int
XLALClusterCoincInspiralTable (
    CoincInspiralTable        **coincList,
    INT8                        dtimeNS,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams *bittenLParams
    );

int
XLALCompareCoincInspiralByTime (
    const void *a,
    const void *b
    );

int
XLALCompareCoincInspiralByEffectiveSnr (
    const void *a,
    const void *b
    );

CoincInspiralTable *
XLALSortCoincInspiral (
    CoincInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

CoincInspiralTable *
XLALSortCoincInspiralByStat (
    CoincInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *),
    CoincInspiralStatParams *statParams,
    CoincInspiralStatistic *stat
    );


int
XLALCoincInspiralIfos (
    CoincInspiralTable  *coincInspiral,
    char                *ifos
    );

int
XLALCoincInspiralIfosCut(
    CoincInspiralTable **coincHead,
    char                *ifos
    );

int
XLALCoincInspiralIfosDiscard(
    CoincInspiralTable **coincHead,
    char                *ifos
    );

UINT8
XLALCoincInspiralIdNumber (
    CoincInspiralTable  *coincInspiral
    );

CoincInspiralTable *
XLALCoincInspiralSlideCut(
    CoincInspiralTable **coincHead,
    int                  slideNum
    );

CoincInspiralTable *
XLALStatCutCoincInspiral (
    CoincInspiralTable         *eventHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       statCut
    );

int
XLALCalcExpFitNLoudestBackground (
    CoincInspiralTable         *coincSlideHead,
    int                         fitNum,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                      *fitStat,
    REAL4                      *fitA,
    REAL4                      *fitB
    );

REAL4
XLALRateCalcCoincInspiral (
    CoincInspiralTable         *eventZeroHead,
    CoincInspiralTable         *eventSlideHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       timeAnalyzed,
    REAL4                       fitStat,
    REAL4                       fitA,
    REAL4                       fitB
    );

REAL4
XLALRateErrorCalcCoincInspiral (
    CoincInspiralTable         *eventZeroHead,
    CoincInspiralSlideTable    *eventSlideHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    int                         numSlides,
    REAL4                       timeAnalyzed,
    REAL4                       fitStat,
    REAL4                       fitA,
    REAL4                       fitB
    );

CoincInspiralTable *
XLALRateStatCutCoincInspiral (
    CoincInspiralTable         *eventZeroHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       statCut,
    REAL4                       rateCut
    );

SnglInspiralTable *
XLALCompleteCoincInspiral (
    CoincInspiralTable         *eventHead,
    int                         ifoList[LAL_NUM_IFO]
    );

CoincInspiralTable *
XLALPlayTestCoincInspiral(
    CoincInspiralTable         *eventHead,
    LALPlaygroundDataMask      *dataType
    );

CoincInspiralTable *
XLALMeanMassCut(
    CoincInspiralTable         *eventHead,
    char                       *massCut,
    REAL4                       massRangeLow,
    REAL4                       massRangeHigh,
    REAL4                       mass2RangeLow,
    REAL4                       mass2RangeHigh
    );

SnglInspiralTable *
XLALMassCut(
    SnglInspiralTable          *eventHead,
    char                       *massCut,
    REAL4                       massRangeLow,
    REAL4                       massRangeHigh,
    REAL4                       mass2RangeLow,
    REAL4                       mass2RangeHigh
    );

/* sim inspiral */

void
LALGalacticInspiralParamsToSimInspiralTable(
    LALStatus                  *status,
    SimInspiralTable           *output,
    GalacticInspiralParamStruc *input,
    RandomParams               *params
    );

void
LALInspiralSiteTimeAndDist(
    LALStatus         *status,
    SimInspiralTable  *output,
    LALDetector       *detector,
    LIGOTimeGPS       *endTime,
    REAL4             *effDist,
    SkyPosition       *skyPos
    );

void
LALPopulateSimInspiralSiteInfo(
    LALStatus                  *status,
    SimInspiralTable           *output
    );

void
XLALSortSimInspiral(
    SimInspiralTable **head,
    int (*comparefunc)(const SimInspiralTable * const *,
      const SimInspiralTable * const *)
    );

int
XLALCompareSimInspiralByGeocentEndTime(
	const SimInspiralTable * const *a,
	const SimInspiralTable * const *b
    );

int
XLALFreeSimInspiral (
    SimInspiralTable **eventHead
    );

void
XLALPlayTestSimInspiral(
    SimInspiralTable          **eventHead,
    LALPlaygroundDataMask      *dataType
    );

int
XLALSimInspiralInSearchedData(
    SimInspiralTable         **eventHead,
    SearchSummaryTable       **summList
    );

int
XLALSimInspiralChirpMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minChirpMass,
    REAL4                maxChirpMass
    );

int
XLALSimInspiralCompMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minCompMass,
    REAL4                maxCompMass,
    REAL4                minCompMass2,
    REAL4                maxCompMass2
    );

int
XLALSimInspiralTotalMassCut(
    SimInspiralTable   **eventHead,
    REAL4                minTotalMass,
    REAL4                maxTotalMass
    );

INT8
XLALReturnSimInspiralEndTime (
    SimInspiralTable *event,
    CHAR             *ifo
    );

int
XLALSnglSimInspiralTest (
    SimInspiralTable  **simHead,
    SnglInspiralTable **eventHead,
    SimInspiralTable  **missedSimHead,
    SnglInspiralTable **missedSnglHead,
    INT8                injectWindowNS
    );

int
XLALCoincSimInspiralTest (
    SimInspiralTable   **simHead,
    CoincInspiralTable **coincHead,
    SimInspiralTable   **missedSimHead,
    CoincInspiralTable **missedCoincHead
   );

/* multi inspiral */

void
LALFreeMultiInspiral (
    LALStatus          *status,
    MultiInspiralTable **eventHead
    );

int
XLALFreeMultiInspiral (
    MultiInspiralTable **eventHead
    );

MultiInspiralTable *
XLALSortMultiInspiral (
    MultiInspiralTable *eventHead,
    int(*comparfunc)   (const void *, const void *)
    );

REAL4
XLALMultiInspiralStat(
    MultiInspiralTable         *multiInspiral,
    SnglInspiralClusterChoice  multiStat
    );

int
XLALClusterMultiInspiralTable (
    MultiInspiralTable         **inspiralList,
    INT8                         dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    );

MultiInspiralTable *
XLALTimeCutMultiInspiral(
    MultiInspiralTable          *eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    );

MultiInspiralTable *
XLALSNRCutMultiInspiral (
    MultiInspiralTable          *eventHead,
    REAL4                       snrCut
    );

MultiInspiralTable *
XLALPlayTestMultiInspiral(
    MultiInspiralTable          *eventHead,
    LALPlaygroundDataMask      *dataType
    );

INT4 XLALCountMultiInspiral( MultiInspiralTable *head );

int
LALCompareMultiInspiralByTime (
    const void *a,
    const void *b
    );

int
XLALMultiSimInspiralTest (
    SimInspiralTable  **simHead,
    MultiInspiralTable **eventHead,
    SimInspiralTable  **missedSimHead,
    MultiInspiralTable **missedMultiHead,
    INT8                injectWindowNS
    );

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
XLALReturnSimRingdownEndTime (
    SimRingdownTable *event,
    CHAR             *ifo
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

void
LALSnglRingdownCoincTest(
    LALStatus                  *status,
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
    LALSegList                 *vetoSegsH2
    );


/* inspiral param accuracy */

void
XLALPopulateAccuracyParams(
       InspiralAccuracyList  *accuracyParams
);

void
XLALPopulateAccuracyParamsExt(
       InspiralAccuracyList  *accuracyParams,
       const LIGOTimeGPS     *gpstime,
       const REAL8            ra_deg,
       const REAL8            dec_deg
);

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LIGOMETADATAINSPIRALUTILS_H */
