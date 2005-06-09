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
  mchirp_and_eta 
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
 * general manipulation functions
 *
 */

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
LALCheckOutTimeFromSearchSummary (
    LALStatus            *status,
    SearchSummaryTable   *summList,
    CHAR                 *ifo,
    LIGOTimeGPS          *startTime,
    LIGOTimeGPS          *endTime
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
    SnglInspiralTable          *inspiralEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    );

void
LALTimeCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    );

int
XLALTimeCutSingleInspiral(
    SnglInspiralTable         **eventHead,
    INT8                        startTimeNS,
    INT8                        endTimeNS
    );

void
LALalphaFCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    REAL4                       alphaFcut
    );

void
LALIfoCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    CHAR                       *ifo
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
LALExtractSnglInspiralFromCoinc(
    LALStatus                  *status,
    SnglInspiralTable         **snglPtr,
    CoincInspiralTable         *coincInspiral,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    );

void
LALCoincCutSnglInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **snglPtr
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




/*
 *
 * burst specific functions
 *
 */

INT4 
XLALCountSnglBurst( 
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
XLALCompareStringBurstByTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

void XLALClusterStringBurstTable(
	SnglBurstTable  **list,
	int (*bailoutfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	int (*testfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *)
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
XLALCompareSnglBurst(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
);

void
LALCompareSnglBurst(
	LALStatus *status,
	const SnglBurstTable *a,
	const SnglBurstTable *b,
	int *difference
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
XLALCompareSimBurstAndSnglBurstByTimeandFreq(
	const SimBurstTable * const *a,
	const SnglBurstTable * const *b
);

int
XLALCompareSimBurstAndSnglBurstByTime(
	const SimBurstTable * const *a,
	const SnglBurstTable * const *b
);

void
LALCompareSimBurstAndSnglBurst(
	LALStatus *status,
	const SimBurstTable *a,
	const SnglBurstTable *b,
	int (*testfunc)(const SimBurstTable * const *, const SnglBurstTable * const *),
	int *match
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
XLALClusterSnglBurstTable(
	SnglBurstTable  **list,
	int (*bailoutfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	int (*testfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *)
);

void
LALClusterSnglBurstTable(
	LALStatus       *status,
	SnglBurstTable  **list,
	int (*bailoutfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *),
	int (*testfunc)(const SnglBurstTable * const *, const SnglBurstTable * const *)
);

void
XLALTimeSlideSnglBurst (
	SnglBurstTable *triggerlist,
	INT8 startTime,
	INT8 stopTime,
	INT8 slideTime
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

