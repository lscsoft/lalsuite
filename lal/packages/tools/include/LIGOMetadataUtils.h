/*----------------------------------------------------------------------- 
 * 
 * File Name: LIGOMetadataUtils.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LIGOMetadataUtilsHV">
Author: Brown, D. A.
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
#define LIGOMETADATAUTILSH_ESSGAP 5
#define LIGOMETADATAUTILSH_ESSDUB 6
#define LIGOMETADATAUTILSH_ETEST 7
#define LIGOMETADATAUTILSH_EDET 8
#define LIGOMETADATAUTILSH_MSGENULL "Null pointer"
#define LIGOMETADATAUTILSH_MSGENNUL "Non-null pointer"
#define LIGOMETADATAUTILSH_MSGETIME "Invalid GPS Time"
#define LIGOMETADATAUTILSH_MSGECOOR "Invalid Coordinate System"
#define LIGOMETADATAUTILSH_MSGESSGAP "Gap in Search Summary Input"
#define LIGOMETADATAUTILSH_MSGESSDUB "Repeated data in Search Summary Input"
#define LIGOMETADATAUTILSH_MSGETEST "Unknown parameter test for sorting events"
#define LIGOMETADATAUTILSH_MSGEDET "Unknown detector"


/* </lalErrTable> */

#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif


/*
 *
 * inspiral specific structures 
 *
 */

typedef enum 
{ 
  no_test,
  m1_and_m2, 
  psi0_and_psi3, 
  mchirp_and_eta 
} 
SnglInspiralParameterTest;


typedef enum 
{
  unspecified_data_type, 
  playground_only,
  exclude_play, 
  all_data
}
LALPlaygroundDataMask;


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


typedef enum 
{ 
  none,
  snr_and_chisq, 
  snrsq_over_chisq, 
  snr 
} 
SnglInspiralClusterChoice;


/*
 *
 * burst specific structures 
 *
 */


typedef struct
tagSnglBurstAccuracy
{
  INT4  match;
  REAL4 dRhoPlus;
  REAL4 dRhoMinus;
  INT8  dtime;
  REAL4 dm;
}
SnglBurstAccuracy;


/*
 *
 * general manipulation functions
 *
 */


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

void
LALIfoScanSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **output,
    SnglInspiralTable          *input,
    CHAR                       *ifo
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


/* coinc inspiral */
void
LALCreateTwoIFOCoincList(
    LALStatus                  *status,
    CoincInspiralTable        **coincOutput,
    SnglInspiralTable          *snglInput,
    SnglInspiralAccuracy       *errorParams
    );


void
LALSnglInspiralLookup(
    LALStatus                  *status,
    SnglInspiralTable         **snglInspiralPtr,
    CoincInspiralTable         *coincInspiral,
    char                       *ifo 
    );


void
LALAddSnglInspiralToCoinc(
    LALStatus                  *status,
    CoincInspiralTable        **coincPtr,
    SnglInspiralTable          *snglInspiral
    );


void
LALCreateNewCoinc(
    LALStatus                  *status,
    CoincInspiralTable        **coincPtr,
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral
    );


void
LALSnglInspiralCoincTest(
    LALStatus                  *status,
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    SnglInspiralAccuracy       *errorParams
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


/*
 *
 * burst specific functions
 *
 */


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
XLALCompareSnglBurstByTime(
	const SnglBurstTable * const *a,
	const SnglBurstTable * const *b
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
	int *match
);

int
XLALCompareSimBurstAndSnglBurst(
	const SimBurstTable * const *a,
	const SnglBurstTable * const *b
);

void
LALCompareSimBurstAndSnglBurst(
	LALStatus *status,
	const SimBurstTable *a,
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

