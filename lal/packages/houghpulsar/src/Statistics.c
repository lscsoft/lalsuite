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
 *
 *  \file Statistics.c
 *  \author Krishnan, B
 *  \ingroup pulsarHough
 *  \brief Module for calculating statistics and histogram of a total Hough map.
 *
 *  $Id$
 *
 *  Bug in variance corrected -- Feb 01, 2004
 *  Generalized to handle REAL4 number counts -- Nov 01, 2005

\par Description

The fuction LALHoughStatistics() calculates the maximum number count, minimum
number count, average and standard deviation of a given total Hough map.
The input HOUGHMapTotal *in is  a total Hough map and the output is a
structure HoughStats *out.

LALHoughHistogram@ produces a histogram of the number counts in a total Hough map.
The input is of type HOUGHMapTotal *in and the output UINT4Vector *out.

*/

/************************************<lalVerbatim file="StatisticsCV">
Authors: Krishnan, B., Sintes, A.M.
$Id$
*************************************</lalVerbatim> */

/* <lalLaTeX>  **********************************************

\subsection{Module \texttt{Statistics.c}}
\label{ss:Statistics.c}
Calculation of statistical quantities of the Hough map and a histogram of the
number counts.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StatisticsD}
\index{\verb&LALHoughStatistics()&}
\index{\verb&LALHoughHistogram()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}
The fuction \verb@LALHoughStatistics@ calculates the maximum number count, minimum
   number count, average and standard deviation of a given total Hough map.
The input \verb@HOUGHMapTotal *in@ is  a total Hough map
and the output is a structure \verb@HoughStats *out@.

The \verb@LALHoughHistogram@ produces a histogram of
   the number counts in a total Hough map. The input is of type \verb@HOUGHMapTotal *in@
and the output \verb@UINT4Vector *out@.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
%\begin{verbatim}
%LALZDestroyVector()
%\end{verbatim}
%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%
\vfill{\footnotesize\input{StatisticsCV}}
********************************************** </lalLaTeX> */


#include <lal/Statistics.h>
/* #include "./Statistics.h" */

NRCSID (STATISTICSC, "$Id$");

/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/** Given a total hough map, this calculates the maximum number count, minimum
   number count, average and standard deviation */

/* *******************************  <lalVerbatim file="StatisticsD"> */
void LALHoughStatistics( LALStatus     *status,
			 HoughStats    *out,
		         HOUGHMapTotal *in)
{ /*-------------------------------------------------</lalVerbatim> */

  INT4   i,j, xSide, ySide, mObsCoh;
  INT4   maxIndexX, maxIndexY, minIndexX, minIndexY;
  REAL8  average, nPixel, variance, ep, temp, sum;
  HoughTT max, min;
  /*--------------------------------------------------------------- */

  INITSTATUS (status, "LALHoughStatistics", STATISTICSC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are not null */
  ASSERT (in, status, STATISTICSH_ENULL, STATISTICSH_MSGENULL);
  ASSERT (out, status, STATISTICSH_ENULL, STATISTICSH_MSGENULL);

  /* make sure input hough map is ok*/
  ASSERT (in->xSide > 0, status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);
  ASSERT (in->ySide > 0, status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);
  ASSERT (in->mObsCoh > 0, status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);

  /* read input parameters */
  xSide = in->xSide;
  ySide = in->ySide;
  mObsCoh = in->mObsCoh;

  /* first find maximum, minimum and average */
  maxIndexX = maxIndexY = minIndexX = minIndexY = 0;
  sum = 0;
  max = min = in->map[0];
  /* loop over hough map and calculate statistics */
  for (i = 0; i < ySide; i++){
    for (j = 0; j < xSide; j++){
      /* read the current number count */
      temp = (REAL4) in->map[i*xSide + j];

      /* add number count to sum */
      sum += temp;

      /* look for maximum */
      if (temp > (REAL4)max){
	max = in->map[i*xSide + j];
        maxIndexX = j;
	maxIndexY = i;
      }

      /* look for minimum */
      if (temp < (REAL4)min){
	min = in->map[i*xSide + j];
	minIndexX = j;
	minIndexY = i;
      }
    }
  }

  nPixel = 1.0 * xSide * ySide;
  average = sum / nPixel;

  /* now do another loop to find the variance and standard deviation */
  /* look at numerical recipes in C for calculation of variance */
  variance = 0.0;
  ep = 0.0;
  for (i = 0; i < ySide; i++){
    for (j = 0; j < xSide; j++){
      temp = (REAL4) (in->map[i*xSide + j]) - average;
      ep += temp;
      variance += temp*temp;
    }
  }
  /* the ep is supposed to reduce the rounding off errors */
  variance = (variance - ep*ep/nPixel)/(nPixel-1);


  /* copy results to output structure */
  out->maxCount = max;
  out->maxIndex[0] = maxIndexX;
  out->maxIndex[1] = maxIndexY;
  out->minCount = min;
  out->minIndex[0] = minIndexX;
  out->minIndex[1] = minIndexY;
  out->avgCount = average;
  out->stdDev = sqrt(variance);


  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}



void LALHoughmapMeanVariance( LALStatus     *status,
			      REAL8         *mean,
			      REAL8         *variance,
			      HOUGHMapTotal *in)
{ /*-------------------------------------------------</lalVerbatim> */

  INT4   i,j, xSide, ySide, mObsCoh, nPixel;
  REAL8  sum, tempVariance, tempMean, ep, temp;
  /*--------------------------------------------------------------- */

  INITSTATUS (status, "LALHoughmapMeanVariance", STATISTICSC);
  ATTATCHSTATUSPTR (status);


  xSide = in->xSide;
  ySide = in->ySide;
  nPixel = 1.0 * xSide * ySide;


  sum = 0.0;
  /* loop over hough map and calculate statistics */
  for (i = 0; i < ySide; i++){
    for (j = 0; j < xSide; j++){
      /* read the current number count */
      sum += in->map[i*xSide + j];
    }
  }

  tempMean = sum/nPixel;

  tempVariance = 0.0;
  ep = 0.0;
  for (i = 0; i < ySide; i++){
    for (j = 0; j < xSide; j++){
      temp = (REAL4) (in->map[i*xSide + j]) - tempMean;
      ep += temp;
      tempVariance += temp*temp;
    }
  }
  /* the ep is supposed to reduce the rounding off errors */
  tempVariance = (tempVariance - ep*ep/nPixel)/(nPixel-1);

  *mean = tempMean;
  *variance = tempVariance;

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/** given a total hough map, this function produces a histogram of
   the number counts */
/* *******************************  <lalVerbatim file="StatisticsD"> */
void LALHoughHistogram(LALStatus      *status,
		       UINT8Vector    *out,
		       HOUGHMapTotal  *in)
{ /* *********************************************  </lalVerbatim> */


  INT4   i, j, length, xSide, ySide, temp;

  INITSTATUS (status, "LALHoughHistogram", STATISTICSC);
  ATTATCHSTATUSPTR (status);

  /* make sure arguments are not null */
  ASSERT (in, status, STATISTICSH_ENULL, STATISTICSH_MSGENULL);
  ASSERT (in->map, status, STATISTICSH_ENULL, STATISTICSH_MSGENULL);
  ASSERT (out, status, STATISTICSH_ENULL, STATISTICSH_MSGENULL);

  /* make sure input hough map is ok*/
  ASSERT (in->xSide > 0, status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);
  ASSERT (in->ySide > 0, status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);
  ASSERT (in->mObsCoh > 0, status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);

  /* make sure output length is same as mObsCoh+1 */
  ASSERT (out->length == in->mObsCoh+1, status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);

  length = out->length;
  xSide = in->xSide;
  ySide = in->ySide;

  /* initialize histogram vector*/
  for (i=0; i < length; i++) out->data[i] = 0;

  /* loop over hough map and find histogram */
  for (i = 0; i < ySide; i++){
    for (j = 0; j < xSide; j++){
      /* read the current number count and round it off to the
         nearest integer -- useful when the number counts are
         floats as when we use weights */
      temp = (INT4)(in->map[i*xSide + j]);

      if ( (temp > length) || (temp < 0) ) {
	ABORT ( status, STATISTICSH_EVAL, STATISTICSH_MSGEVAL);
      }
      /* add to relevant entry in histogram */
      out->data[temp] += 1;
    }
  }

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}





