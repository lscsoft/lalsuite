/*
*  Copyright (C) 2007 Badri Krishnan, Jolien Creighton
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
 * File Name: RngMedBias.c
 *
 * Authors: Krishnan, B.  Itoh, Y.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="RngMedBiasCV">
Author: Krishnan, B., Itoh, Y.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{MOdule \texttt{RngMedBias.c}}
\label{ss:RngMedBias.c}
Routine for finding bias in median for exponential distribution
To be used with any code which uses the running median to estimate PSD.

For the exponential distribution with unit mean and variance, the value of the
median is $\log(2.0)$ in the limit of infinite sample size. Thus, if we are
using the running median to estimate the PSD, there is a correction factor
of $\log(2.0)$.  However, for finite sample sizes (i.e. for finite block size
values), there is a bias in the estimator of the median and the correction
factor is different.  This program returns the correct normalization factor
for block sizes from 1 to 1000.  For larger values it returns $\log(2.0)$ and
returns and error for smaller values.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{RngMedBiasD}
\index{\verb&LALRngMedBias()&}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHO()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{RngMedBiasCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*********************************************** </lalLaTeX> */


#include <lal/RngMedBias.h>

NRCSID (RNGMEDBIASC, "$Id$");

/*
 * The functions that make up the guts of this module
 */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="RngMedBiasD"> */
void LALRngMedBias (LALStatus   *status,
		 REAL8       *biasFactor,
		 INT4        blkSize
                 )
{/*   *********************************************  </lalVerbatim> */

  REAL8 temp;
  INT4 plusminus, count;


  /* --------------------------------------------- */
  INITSTATUS (status, "LALRngMedBias", RNGMEDBIASC);
  ATTATCHSTATUSPTR (status);

  /* check arguments are not null and block size is positive*/
  ASSERT (biasFactor, status, RNGMEDBIASH_ENULL, RNGMEDBIASH_MSGENULL);
  ASSERT (blkSize > 0, status,  RNGMEDBIASH_EVAL, RNGMEDBIASH_MSGEVAL);

  /* if blkSize is even, reduce it by one */
  if ( (blkSize % 2) == 0)
    blkSize -= 1;

  /* now sum alternating harmonic series upto blkSize */
  temp = 0.0;
  plusminus = 1;
  for (count = 0; count < blkSize; count++)
    {
      temp += plusminus / (count + 1.0);
      plusminus *= -1;
    }

  *biasFactor =  temp;


  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}








