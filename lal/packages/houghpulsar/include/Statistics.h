/* ***************************************************
 *  File Name: Statistics.h
 *
 *  Author: Krishnan, B.
 *
 *  $Id$ 
 *
 *  Created by Badri Krishnan on July 09, 2003
***************************************************** */


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

typedef struct tagHoughStats {
  UINT4    maxCount;
  UINT2    maxIndex[2];
  UINT4    minCount; 
  UINT2    minIndex[2];
  REAL8    avgCount;
  REAL8    stdDev;
} HoughStats;

/*
 *  Extern Global variables. (discouraged) 
 */
  

/* ***************************************************
 *  Functions Declarations (i.e., prototypes).
 */
 
void LALHoughStatistics(LALStatus      *status,
		        HoughStats     *out, /* output containing statistics */ 
		        HOUGHMapTotal  *in); /* hough map */

void LALHoughHistogram(LALStatus       *status, 
		       UINT4Vector     *out,  /* histogram */ 
		       HOUGHMapTotal   *in);  /* hough map*/

/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif  /* end of double inclusion protection */













