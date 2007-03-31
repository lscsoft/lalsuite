/*-----------------------------------------------------------------------
 *
 * File Name: NullStatistic.h
 *
 * Author: Messaritaki, E.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="NullStatisticHV">
Author: Messaritaki, E.
$Id$
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{NullStatistic.h}}
\label{s:NullStatistic.h}

\noindent Provides core prototypes, structures and functions to calculate
the null statistic for the binary inspiral searches.

\subsection*{Calculation of the null statistic}

The null statistic will be defined here.

\subsubsection*{Synopsis}

\begin{verbatim}
#include <lal/NullStatistic.h>
\end{verbatim}

\input{NullStatisticHDoc}

</lalLaTeX>
#endif

#ifndef _NULLSTATISTICH_H
#define _NULLSTATISTICH_H

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALDatatypes.h>
#include <lal/DataBuffer.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/LALInspiralBank.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID (NULLSTATISTICH, "$Id$");

#if 0
<lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
#endif
/* <lalErrTable> */
#define NULLSTATH_ENULL 1
#define NULLSTATH_ENNUL 2
#define NULLSTATH_EALOC 3
#define NULLSTATH_ENUMZ 4
#define NULLSTATH_ESEGZ 5
#define NULLSTATH_EDTZO 6
#define NULLSTATH_EFREE 7
#define NULLSTATH_ETHRN 8
#define NULLSTATH_ESMSM 9
#define NULLSTATH_ENDET 10
#define NULLSTATH_MSGENULL "Null pointer"
#define NULLSTATH_MSGENNUL "Non-null pointer"
#define NULLSTATH_MSGEALOC "Memory allocation error"
#define NULLSTATH_MSGENUMZ "Invalid number of points in segment"
#define NULLSTATH_MSGESEGZ "Invalid number of segments"
#define NULLSTATH_MSGEDTZO "deltaT is zero or negative"
#define NULLSTATH_MSGEFREE "Error freeing memory"
#define NULLSTATH_MSGETHRN "null statistic threshold is negative"
#define NULLSTATH_MSGESMSM "Size mismatch between vectors"
#define NULLSTATH_MSGENDET "Number of detectors should be 2; H1 and H2"
/* </lalErrTable> */


typedef struct
tagNullStatInitParams
{
  UINT4   numDetectors;
  UINT4   numSegments;
  UINT4   numPoints;
  UINT4   nullStatOut;
}
NullStatInitParams;


typedef struct
tagNullStatInputParams
{
  DetectorVector             detVector;
  InspiralTemplate          *tmplt;
  CoherentInspiralCVector   *multiCData;
}
NullStatInputParams;
  
typedef struct
tagNullStatParams
{
  INT4              numTmplts;
  UINT4             maxOverChirp;
  UINT4             numDetectors;
  UINT4             numSegments;
  UINT4             numPoints;
  REAL4             fLow;
  REAL4             deltaT;
  REAL4             nullStatThresh;
  REAL8             sigmasq[4];
  REAL4             templateNorm;
  INT4              segmentLength;
  UINT4             nullStatOut;
  REAL4TimeSeries   nullStatVec;
}
NullStatParams;


typedef struct
tagNullStatCVector
{
  UINT4                numDetectors;
  COMPLEX8TimeSeries   *cData;
}
NullStatCVector;

