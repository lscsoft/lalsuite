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
#define LIGOMETADATAUTILSH_MSGENULL "Null pointer"
#define LIGOMETADATAUTILSH_MSGENNUL "Non-null pointer"
#define LIGOMETADATAUTILSH_MSGETIME "Invalid GPS Time"
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


/*
 *
 * inspiral specific functions
 *
 */


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
LALCompareSnglInspiralByTime (
    const void *a,
    const void *b
    );

int
LALCompareSnglInspiralByPsi (
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
LALCompareSnglInspiralBCV (
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

/*
 *
 * burst specific functions
 *
 */


void
LALSortSnglBurst(
    LALStatus          *status,
    SnglBurstTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    );

int
LALCompareSnglBurstByTime(
    const void *a,
    const void *b
    );

int
LALCompareSnglBurstByTimeAndFreq(
    const void *a,
    const void *b
    );

void
LALCompareSnglBurst(
    LALStatus             *status,
    SnglBurstTable        *aPtr,
    SnglBurstTable        *bPtr,
    SnglBurstAccuracy     *params
    );

void
LALCompareSimBurstAndSnglBurst(
    LALStatus             *status,
    SimBurstTable         *aPtr,
    SnglBurstTable        *bPtr,
    SnglBurstAccuracy     *params
    );

void
LALClusterSnglBurstTable (
	      LALStatus        *status,
              SnglBurstTable   *burstEvent,
	      INT4             *nevents
	      );

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{LIGOMetadataUtilsHV}}

\newpage\input{LIGOMetadataUtilsC}
\newpage\input{SnglInspiralUtilsC}
</lalLaTeX>
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LIGOMETADATAUTILS_H */

