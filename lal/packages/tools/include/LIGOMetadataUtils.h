/*----------------------------------------------------------------------- 
 * 
 * File Name: LIGOMetadataUtils.h
 *
 * Author: Brown, D. A. and Fairhurst, S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LIGOMetadataUtilsHV">
Author: Brown, D. A. and Fairhurst, S.
$Id$
</lalVerbatim> 
<lalLaTeX>
\section{Header \texttt{LIGOMetadataUtils.h}}
\label{s:LIGOMetadataItils.h}

Provides functions for manipulating the LAL structures that correspond
to the LIGO metadata database tables defined in \texttt{LIGOMetadataTables.h}.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LIGOMetadataUtils.h>
\end{verbatim}

\noindent This header provides prototypes for routines that perform processing
on the LAL structures that correspond to the LIGO metadata database tables
defined in \texttt{LIGOMetadataTables.h}, such as sorting and eliminating 
duplictaes. The functions specific to a particular metadata table (e.g. 
\texttt{sngl\_inspiral}, \texttt{sngl\_burst}, etc.) are all prototyped in 
this header.

\subsection*{Types}

\noindent None.

</lalLaTeX>
#endif

#ifndef _LIGOMETADATAUTILS_H
#define _LIGOMETADATAUTILS_H

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

#include <lal/LIGOMetadataTables.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/Segments.h>

NRCSID( LIGOMETADATAUTILSH, "$Id$" );

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define LIGOMETADATAUTILSH_ENULL 1
#define LIGOMETADATAUTILSH_ENNUL 2
#define LIGOMETADATAUTILSH_ETIME 3
#define LIGOMETADATAUTILSH_ECOOR 4
#define LIGOMETADATAUTILSH_ESGAP 5
#define LIGOMETADATAUTILSH_ESDUB 6
#define LIGOMETADATAUTILSH_ETEST 7
#define LIGOMETADATAUTILSH_EDET 8
#define LIGOMETADATAUTILSH_MSGENULL "Null pointer"
#define LIGOMETADATAUTILSH_MSGENNUL "Non-null pointer"
#define LIGOMETADATAUTILSH_MSGETIME "Invalid GPS Time"
#define LIGOMETADATAUTILSH_MSGECOOR "Invalid Coordinate System"
#define LIGOMETADATAUTILSH_MSGESGAP "Gap in Search Summary Input"
#define LIGOMETADATAUTILSH_MSGESDUB "Repeated data in Search Summary Input"
#define LIGOMETADATAUTILSH_MSGETEST "Unknown parameter test for sorting events"
#define LIGOMETADATAUTILSH_MSGEDET "Unknown detector"


/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
\idx[Type]{LALPlaygroundDataMask}
\idx[Type]{SnglInspiralParameterTest}
\idx[Type]{SnglInspiralAccuracy}
\idx[Type]{SnglInspiralClusterChoice}
\idx[Type]{SnglBurstAccuracy}





\subsubsection*{Type \texttt{LALPlaygroundDataMask}}
#endif
/* <lalVerbatim> */
typedef enum 
{
  unspecified_data_type, 
  playground_only,
  exclude_play, 
  all_data
}
LALPlaygroundDataMask;
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{LALPlaygroundDataMask} contains an enum type for describing the
subset of data to be used, \texttt{playground\_only}, \texttt{exclude\_play}
and \texttt{all\_data}.
\subsubsection*{Type \texttt{LALPlaygroundDataMask}}
</lalLaTeX>
#endif


/*
 *
 * inspiral specific structures 
 *
 */

#if 0
<lalLaTeX>
\subsubsection*{Type \texttt{SnglInspiralParameterTest}}
</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef enum 
{ 
  no_test,
  m1_and_m2, 
  psi0_and_psi3, 
  mchirp_and_eta, 
  mchirp_and_eta_ext 
} 
SnglInspiralParameterTest;
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{SnglInspiralParameterTest} contains an enum type for each of the
tests of mass parameters which are used.
\subsubsection*{Type \texttt{SnglInspiralAccuracy}}
</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagSnglInspiralAccuracy
{
  INT4        match;
  REAL4       epsilon;
  REAL4       kappa;
  INT8        dt;
  REAL4       dm;
  REAL4       deta;
  REAL4       dmchirp;
  REAL4	      dmchirpHi;
  REAL4       highMass;
  REAL4       dpsi0;
  REAL4       dpsi3;
  SnglInspiralParameterTest test;
}
SnglInspiralAccuracy;
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{SnglInspiralAccuracy} structure contains parameters used for 
testing coincidence between two or more single inspiral tables.  These include
a timing accuracy \texttt{dt}, five mass accuracies \texttt{dm} (used for 
testing \texttt{mass1} and \texttt{mass2}), \texttt{deta}, \texttt{dmchirp},
\texttt{dpsi0} and \texttt{dpsi3}.  It also includes the parameters 
\texttt{kappa} and \texttt{epsilon} which are used for testing consistency of
effective distance.
\subsubsection*{Type \texttt{SnglInspiralClusterChoice}}
</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagInspiralAccuracyList
{
  INT4                      match;
  SnglInspiralParameterTest test;
  SnglInspiralAccuracy      ifoAccuracy[LAL_NUM_IFO];
  INT8                      lightTravelTime[LAL_NUM_IFO][LAL_NUM_IFO];
  REAL4                     iotaCutH1H2;
  REAL4                     iotaCutH1L1;
}
InspiralAccuracyList;
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{InspiralAccuracyList} structure contains parameter accuracies for
each of the six global interferometers.  These are stored in the
\texttt{SnglInspiralAccuracy} structure.  The accuracies stored should be the
accuracy with which each instrument can determine the given parameter.  It also
contains a \texttt{match} which is set to 1 to signify that coincidence
criteria are satisfied and 0 otherwise.  Finally, the
\texttt{SnglInspiralParameterTest} must be specified.
\subsubsection*{Type \texttt{SnglInspiralClusterChoice}}
</lalLaTeX>
#endif

/* <lalVerbatim> */
typedef struct
tagCoincInspiralBittenLParams
{
  REAL4    param_a[LAL_NUM_IFO];
  REAL4    param_b[LAL_NUM_IFO];
}
CoincInspiralBittenLParams;
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{CoincInspiralBittenLParams} structure contains the bitten L parameter for
each of the six global interferometers.  These are stored in the
\texttt{param\_a} and \texttt{param\_b} structure. 
</lalLaTeX>
#endif

/* <lalVerbatim> */
typedef enum 
{ 
  none,
  snr_and_chisq, 
  snrsq_over_chisq, 
  snr 
} 
SnglInspiralClusterChoice;
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{SnglInspiralClusterChoice} provides three choices for clustering
a single inspiral table.  The\texttt{snr} clustering returns the trigger
with the greatest signal to noise ratio; \texttt{snr\_and\_chisq} replaces
the existing trigger if the new trigger has \textit{both} a greater snr and
a smaller chi squared value; \texttt{snrsq\_over\_chisq} selects the trigger
with the largest value of snr squared divided by the chi squared.
</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef enum 
{ 
  no_stat,
  snrsq,
  effective_snrsq,
  s3_snr_chi_stat,
  bitten_l 
} 
CoincInspiralStatistic;
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{CoincInspiralStatistic} provides two choices for clustering
a single inspiral table.  The\texttt{snrsq} clustering returns the trigger
with the greatest summed snr$^{2}$ from all instruments.  The 
\texttt{snr\_chi\_stat} replaces selects the trigger
with the largest value of the snr and chisq statistic and the \texttt{bitten\_l}
returns the minimum among the summed snr$^{2}$ from all instruments and the
$a\times snr_i - b$ in each detector. The parameters  $a$ and $b$ must be
provided by the user.
</lalLaTeX>
#endif
/* <lalVerbatim> */
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
/*</lalVerbatim> */
#if 0
<lalLaTeX>
The \texttt{SnglInspiralBCVCalphafCut} provides entries for cutting single IFO triggers generated with the BCVC code. For each LSC IFO there is a field \texttt{lo} and \texttt{hi} which corresponds to the area allowing triggers.
</lalLaTeX>
#endif
/*
 *
 * burst specific structures 
 *
 */

#if 0
<lalLaTeX>
\subsubsection*{Type \texttt{SnglSnglBurstAccuracy}}
</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagSnglBurstAccuracy
{
  INT4  difference;
  REAL4 dRhoPlus;
  REAL4 dRhoMinus;
  INT8  dtime;
  REAL4 dm;
}
SnglBurstAccuracy;
/*</lalVerbatim> */



/*
 * 
 *  ringdown specific structures
 *    
 */
typedef enum
{
  unknown_test,
  f_and_Q
}
SnglRingdownParameterTest;

typedef struct
tagSnglRingdownAccuracy
{  
  INT4        match;
  INT8        dt;
  REAL4       df;
  REAL4       dQ;
  REAL4       ddeff;
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
}
RingdownAccuracyList;
  


/*
 *
 * general manipulation functions
 *
 */

int 
XLALCountProcessTable(
    ProcessTable *head
    );

int 
XLALCountProcessParamsTable(
    ProcessParamsTable *head
    );

int 
XLALCountMultiInspiralTable(
    MultiInspiralTable *head
    );

int 
XLALIFONumber( 
    const char *ifo 
    );

void 
XLALReturnIFO( 
    char                *ifo,
    InterferometerNumber IFONumber 
    );

void
XLALReturnDetector(
    LALDetector           *det,
    InterferometerNumber   IFONumber 
    );


int
XLALPlaygroundInSearchSummary (
    SearchSummaryTable *ssTable,
    LIGOTimeGPS        *inPlayTime,
    LIGOTimeGPS        *outPlayTime
    );

void
LALPlaygroundInSearchSummary (
    LALStatus          *status,
    SearchSummaryTable *ssTable,
    LIGOTimeGPS        *inPlayTime,
    LIGOTimeGPS        *outPlayTime
    );


void
LALTimeCheckSearchSummary (
    LALStatus          *status,
    SearchSummaryTable *ssTable,
    LIGOTimeGPS        *startTime,
    LIGOTimeGPS        *endTime
    );
 
int
LALCompareSearchSummaryByInTime (
    const void *a,
    const void *b
    );

int
LALCompareSearchSummaryByOutTime (
    const void *a,
    const void *b
    );

int
XLALTimeSortSearchSummary(
    SearchSummaryTable  **summHead,
    int(*comparfunc)    (const void *, const void *)
    );

void
LALTimeSortSearchSummary (
    LALStatus            *status,
    SearchSummaryTable  **summHead,
    int(*comparfunc)    (const void *, const void *)
    );

void
LALIfoScanSearchSummary(
    LALStatus                  *status,
    SearchSummaryTable        **output,
    SearchSummaryTable         *input,
    CHAR                       *ifo
    );

void
LALDistanceScanSummValue (
    LALStatus            *status,
    SummValueTable       *summList,
    LIGOTimeGPS          gps,
    CHAR                 *ifo,
    REAL4                 distance
    );

void
LALCheckOutTimeFromSearchSummary (
    LALStatus            *status,
    SearchSummaryTable   *summList,
    CHAR                 *ifo,
    LIGOTimeGPS          *startTime,
    LIGOTimeGPS          *endTime
    );

SearchSummaryTable *
XLALIfoScanSearchSummary(
    SearchSummaryTable         *input,
    CHAR                       *ifos
    );

void
LALIfoScanSummValue(
    LALStatus                  *status,
    SummValueTable            **output,
    SummValueTable             *input,
    CHAR                       *ifo
    );

int
LALCompareSummValueByTime (
    const void *a,
    const void *b
    );

int
XLALTimeSortSummValue (
    SummValueTable      **summHead,
    int(*comparfunc)    (const void *, const void *)
    );

void
LALTimeSortSummValue (
    LALStatus            *status,
    SummValueTable      **summHead,
    int(*comparfunc)    (const void *, const void *)
    );

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

void
LALCompareSnglInspiral (
    LALStatus                *status,
    SnglInspiralTable        *aPtr,
    SnglInspiralTable        *bPtr,
    SnglInspiralAccuracy     *params
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
    REAL4                       rsqMaxSnr
    );

SnglInspiralTable *
XLALVetoSingleInspiral (
    SnglInspiralTable          *eventHead,
    LALSegList                 *vetoSegs, 
    CHAR 		       *ifo
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

void
LALTimeSlideSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable          *input,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime,
    LIGOTimeGPS                 slideTimes[]
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
    LALStatus *status, 
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams,
    REAL4 			snrThreshold,
    SummValueTable            **summValueList,
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
XLALRecreateCoincFromSngls(
    CoincInspiralTable        **coincPtr,
    SnglInspiralTable          *snglInspiral
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
    CoincInspiralTable         *coincInspiral
    );

REAL4
XLALCoincInspiralStat(
    CoincInspiralTable         *coincInspiral,
    CoincInspiralStatistic      coincStat,
    CoincInspiralBittenLParams *bittenLParams
    );

int 
XLALClusterCoincInspiralTable (
    CoincInspiralTable        **coincList,
    INT8                        dtimeNS,
    CoincInspiralStatistic      coincStat,
    CoincInspiralBittenLParams *bittenLParams
    );

int
XLALCompareCoincInspiralByTime (
    const void *a,
    const void *b
    );

CoincInspiralTable *
XLALSortCoincInspiral (
    CoincInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
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
    CoincInspiralBittenLParams *bittenLParams,
    REAL4                       statCut
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
    SearchSummaryTable        *summList
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


/*
 *
 * burst specific functions
 *
 */

SnglBurstTable *
XLALVetoSnglBurst(
    SnglBurstTable             *eventHead,
    LALSegList                 *vetoSegs 
    );

void
XLALFreeSnglBurst (
    SnglBurstTable *event
);

INT4 
XLALCountSnglBurst( 
    SnglBurstTable *head 
);

void
XLALSnglBurstAssignIDs(
	SnglBurstTable *head
);

void
XLALSortSnglBurst(
	SnglBurstTable **head,
	int (*comparefunc)(const SnglBurstTable * const *, const SnglBurstTable * const *)
);

void
LALSortSnglBurst(
	LALStatus *status,
	SnglBurstTable **head,
	int (*comparefunc)(const SnglBurstTable * const *, const SnglBurstTable * const *)
);

int
XLALCompareSnglBurstByStartTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByPeakTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByExactPeakTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareStringBurstByTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareStringBurstByTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareStringBurstByAmplitude(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);


int
XLALCompareSnglBurstBySNR(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByPeakTimeAndSNR(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);


int
XLALCompareSnglBurstByTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstSnglInspiralByTime(
	const SnglBurstTable * const *a,
	const SnglInspiralTable * const *b
);

int
XLALCompareSnglBurstByLowFreq(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByFreq(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByStartTimeAndLowFreq(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByPeakTimeAndFreq(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByIFOPeakTimeAndFreq(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByIFOTimeAndFreq(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstByIFO(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurst(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSnglBurstSnglInspiral(
	const SnglBurstTable * const *a,
	const SnglInspiralTable * const *b
);

void
LALCompareSnglBurstSnglInspiral(
	LALStatus *status,
	const SnglBurstTable *a,
	const SnglInspiralTable *b,
	int *difference,
	INT8 deltaT
);

int
XLALCompareSimBurstSnglInspiral(
	const SimBurstTable * const *a,
	const SnglInspiralTable * const *b
);

int
XLALCompareSimBurstAndSnglBurstByTimeandFreq(
	const SimBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSimBurstAndSnglBurstByTime(
	const SimBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareCoincBurstByStartTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b,
	const SnglBurstTable * const *c,
	const SnglBurstTable * const *d
);

int
XLALCompareCoincBurstByPeakTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b,
	const SnglBurstTable * const *c,
	const SnglBurstTable * const *d
);

int
XLALCompareSimInspiralAndSnglBurst(
	const SimInspiralTable * const *a,
	const SnglBurstTable * const *b
);

void
LALCompareSimInspiralAndSnglBurst(
	LALStatus *status,
	const SimInspiralTable *a,
	const SnglBurstTable *b,
	int *match
);

void
XLALSnglBurstCluster(
	SnglBurstTable *a,
	const SnglBurstTable *b
);

void
XLALStringBurstCluster(
	SnglBurstTable *a,
	const SnglBurstTable *b
);

void
XLALCoincBurstCluster(
	SnglBurstTable *a,
	SnglBurstTable *b,
	SnglBurstTable *c,
	SnglBurstTable *d
);

void
XLALClusterSnglBurstTable(
	SnglBurstTable  **list,
	int (*bailoutfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	int (*testfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	void (*clusterfunc)(SnglBurstTable *, const SnglBurstTable *)
);

void
XLALClusterTWOCoincSnglBurstTable(
	SnglBurstTable  **list,
	int (*bailoutfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	int (*testfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	void (*clusterfunc)(SnglBurstTable *, const SnglBurstTable *)
);

void
XLALClusterCoincSnglBurstTable(
	SnglBurstTable  **list,
	int (*bailoutfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	int (*testfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *, const SnglBurstTable * const *, const SnglBurstTable * const * ),
	void (*clusterfunc)(SnglBurstTable *, SnglBurstTable *, SnglBurstTable *, SnglBurstTable *)
);

void
XLALTimeSlideSnglBurst (
	SnglBurstTable *triggerlist,
	INT8 startTime,
	INT8 stopTime,
	INT8 slideTime
);

void
XLALIfoCutSnglBurst (
	SnglBurstTable  **eventHead,
	CHAR            *ifo
);


/*
 *  ringdown specific functions
 * 
 *      */

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

int
LALCompareSnglRingdownByTime (
    const void *a,
    const void *b
    );

void
LALClusterSnglRingdownTable (
    LALStatus                  *status,
    SnglRingdownTable          *ringdownEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    );

void
LALIfoCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    CHAR                       *ifo
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

void
LALSnglRingdownCoincTest(
    LALStatus                  *status,
    CoincRingdownTable         *coincRingdown,
    SnglRingdownTable          *snglRingdown,
    RingdownAccuracyList       *accuracyParams
    );

void
LALExtractSnglRingdownFromCoinc(
    LALStatus                  *status,
    SnglRingdownTable         **snglPtr,
    CoincRingdownTable         *coincRingdown,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    );


#if 0
<lalLaTeX>
\vfill{\footnotesize\input{LIGOMetadataUtilsHV}}

\newpage\input{LIGOMetadataUtilsC}
\newpage\input{SnglInspiralUtilsC}
\newpage\input{CoincInspiralUtilsC}
\newpage\input{SimInspiralUtilsC}
</lalLaTeX>
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LIGOMETADATAUTILS_H */

