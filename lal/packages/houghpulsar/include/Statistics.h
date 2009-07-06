/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes
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

/**
 *  \file Statistics.h
 *  \ingroup moduleHoughPulsar
 *  \brief Functions prototypes and structures for calculating statistical properties of Hough maps
 *  \author Krishnan, B.
 *
 *  $Id$
 *
 *  Created by Badri Krishnan on July 09, 2003


\par Description

Given a total Hough map, the function LALHoughStatistics() calculates
the maximum number count, minimum number count, average and
standard deviation of the number counts in a given Hough map.  The
function LALHoughHistogram() produces a histogram of the number counts.

*/


/* ***********************************<lalVerbatim file="StatisticsHV">
Author: Krishnan, B., Sintes, A.M.
$Id$
*************************************</lalVerbatim> */

/* <lalLaTeX> *************************************************
\section{Header \texttt{Statistics.h}}
\label{s:Statistics.h}
Computes statistics of the Hough maps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Synopsis}

\begin{verbatim}
#include <lal/Statistics.h>
\end{verbatim}

Given a total Hough map, this calculates the maximum number count, minimum
number count, average and standard deviation and produces a histogram of the
number counts.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Error conditions}
\vspace{0.1in}
\input{StatisticsHErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection*{Structures}
\begin{verbatim}
struct HoughStats
\end{verbatim}
\index{\texttt{HoughStats}}
\noindent This structure stores the statistics of a Hough map.  The fields are:

\begin{description}
\item[\texttt{UINT4    maxCount}]
\item[\texttt{UINT4    maxIndex[2]}]
\item[\texttt{UINT4    minCount}]
\item[\texttt{UINT4    minIndex[2]}]
\item[\texttt{REAL8    avgCount }]
\item[\texttt{REAL8    stdDev }]
\end{description}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{StatisticsHV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newpage\input{StatisticsC}
%%%%%%%%%%Test program. %%
\newpage\input{TestStatisticsC}

*************************************************** </lalLaTeX> */

/* double inclusion protection */
#ifndef _STATISTICS_H
#define _STATISTICS_H

/* *************
 *    Includes. This header may include others; if so, they go immediately
 *    after include-loop protection. Includes should appear in the following
 *    order:
 *    a. Standard library includes
 *    b. LDAS includes
 *    c. LAL includes
 */
#include<lal/Date.h>
#include<lal/LALDatatypes.h>
#include<lal/HoughMap.h>
#include<lal/LALStdlib.h>
#include<lal/LALConstants.h>

/*  ****************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif

/* ***************************************
 *   Assignment of Id string using NRCSID()
 */

NRCSID( STATISTICSH, "$Id$");

/* ***************************************
 *   Error codes and messages. This must be auto-extracted for
 *    inclusion in the documentation.
 */

/* <lalErrTable file="StatisticsHErrorTable"> */
#define STATISTICSH_ENULL 1
#define STATISTICSH_EVAL 2
#define STATISTICSH_MSGENULL "Null Pointer"
#define STATISTICSH_MSGEVAL "Invalid Value"
/* </lalErrTable>  */

/* ******************************************************
 *   Macros. But, note that macros are deprecated.
 *    They could be moved to the modules where are needed
 */


/* *******************************************************
 *  Constant Declarations. (discouraged)
 */

/* *****************************************************
 *   Structure, enum, union, etc., typdefs.
 */

/** Structure for storing statistical information about a Hough map */
typedef struct tagHoughStats {
  HoughTT    maxCount;    /**< maximum number count */
  UINT2      maxIndex[2]; /**< loctaion of maximum number count */
  HoughTT    minCount;    /**< minimum number count */
  UINT2      minIndex[2]; /**< location of minimum number count */
  REAL8      avgCount;    /**< average number count */
  REAL8      stdDev;      /**< standard deviation of number counts */
} HoughStats;

/*
 *  Extern Global variables. (discouraged)
 */


/* ***************************************************
 *  Functions Declarations (i.e., prototypes).
 */
/** Calculates max, min, average and standard deviation of Hough number counts */
void LALHoughStatistics(LALStatus      *status,
		        HoughStats     *out, /**< output containing statistics */
		        HOUGHMapTotal  *in); /**< hough map */

/** Calculates number count histogram */
void LALHoughHistogram(LALStatus       *status,
		       UINT8Vector     *out,  /**< histogram */
		       HOUGHMapTotal   *in);  /**< hough map*/

void LALHoughmapMeanVariance( LALStatus     *status,
			      REAL8         *mean,
			      REAL8         *variance,
			      HOUGHMapTotal *in);

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif  /* end of double inclusion protection */













