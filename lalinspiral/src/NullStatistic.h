/*
*  Copyright (C) 2007 Eirini Messaritaki
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
tagChanNames
{
  CHAR  *chanNameH1[LALNameLength];
  CHAR  *chanNameH2[LALNameLength];
}
ChanNames;


typedef struct
tagCVector
{
  UINT4                 numDetectors;
  COMPLEX8TimeSeries   *cData[LAL_NUM_IFO];
}
CVector;


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
  InspiralTemplate      *tmplt; /*probably not necessary */
  DetectorVector         detVector;
  CVector               *CData;
}
NullStatInputParams;


typedef struct
tagNullStatParams
{
  INT4               numTmplts;
  UINT4              maxOverChirp; /* probably not necessary */
  UINT4              numDetectors;
  UINT4              numSegments;
  UINT4              numPoints;
  REAL4              nullStatThresh;
  REAL4              eventNullStat;
  REAL4              sigmasq[6];
  INT4               segmentLength;
  DetectorVector    *detVector;
  UINT4              nullStatOut;
  REAL4TimeSeries   *nullStatVec;
}
NullStatParams;


int
XLALNullStatisticInputInit (
   NullStatInputParams    **input,
   NullStatInitParams      *init
   );

int
XLALNullStatisticParamsInit(
   NullStatParams      **output,
   NullStatInitParams   *init
   );

int
XLALComputeNullStatistic (
  MultiInspiralTable    **eventList,
  NullStatInputParams    *input,
  NullStatParams         *params
  );

int
XLALNullStatisticParamsFinal(
   NullStatParams      **output
   );

int
XLALNullStatisticInputFinal (
   NullStatInputParams    **input
   );


#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _NULLSTATISTICH_H */

