/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Steve Berukoff
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

/*********************************** <lalVerbatim file="ComputeSkyCV">
Author:  Berukoff, S.J., Papa, M.A.
$Id$
************************************ </lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{ComputeSky.c}}\label{ss:ComputeSky.c}
Computes the phase model coefficients necessary for a successful demodulation.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ComputeSkyCP}
\index{\texttt{ComputeSky()}}

\subsubsection*{Description}
Given an input index which refers to the sky patch under consideration, this
routine returns the phase model coefficients $A_{s\alpha}$ and $B_{s\alpha}$
which are needed to correctly account for the phase variance of a signal over
time.  The \verb@CSParams@ parameter structure contains relevant information
for this routine to properly run.  In particular, it contains an array of
timestamps in \verb@LIGOTimeGPS@ format, which are the GPS times of the first
data from each SFT.  The \verb@input@ is an \verb@INT4@ variable
\verb@iSkyCoh@, which is the index of the sky patch under consideration.  For
each sky patch, this code needs to be run once; the necessary phase model
coefficients are calculated, and can then be applied to the relevant spindown
parameter sets one is using in their search.

\subsubsection*{Algorithm}
The routine uses a simplistic nested for-loop structure.  The outer loop is
over the number of SFTs in the observation timescale; this accounts for the
temporal variability of the phase model coefficients.  The inner loop is over
the number of spindown parameters in one set.  Inside the inner loop, the
values are calculated using the analytical formulae given in the
\verb@ComputeSky.h@ documentation.

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

The reference-time, at which the pulsar spin-parameters are defined, is
taken to be the start-time *INTERPRETED* as an SSB time (i.e. no translation
is done, the times are numerically equal!).

\vfill{\footnotesize\input{ComputeSkyCV}}

</lalLaTeX> */

#include <math.h>
#include <lal/LALConstants.h>
#include "ComputeSky.h"
#include "LALBarycenter.h"

NRCSID (COMPUTESKYC, "$Id ComputeSky.c $");

static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps);
static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f);

/* <lalVerbatim file="ComputeSkyCP"> */
void
LALComputeSky (LALStatus *status,
	       REAL8 *skyConst,
	       INT8 iSkyCoh,
	       CSParams	*params)
{  /* </lalVerbatim> */

  INT4	m, n;
  REAL8	t;
  REAL8	basedTbary;

  REAL8	dTbary;
  REAL8	tBary;
  REAL8	tB0;
  INITSTATUS (status, "LALComputeSky", COMPUTESKYC);
  ATTATCHSTATUSPTR(status);

  /* Check for non-negativity of sky positions in SkyCoh[] */
  ASSERT(iSkyCoh>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);

  /* Check to make sure sky positions are loaded */
  ASSERT(params->skyPos!=NULL, status, COMPUTESKYH_ENULL, COMPUTESKYH_MSGENULL);
  ASSERT(params->skyPos!=NULL, status, COMPUTESKYH_ENULL, COMPUTESKYH_MSGENULL);

  /* Check to make sure parameters are loaded and reasonable */
  ASSERT(params->spinDwnOrder>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);
  ASSERT(params->mObsSFT>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);
  ASSERT(params->tSFT>=0, status, COMPUTESKYH_ENEGA, COMPUTESKYH_MSGENEGA);

  for(n=0;n<params->mObsSFT;n++)
    {
      ASSERT(params->tGPS[n].gpsSeconds>=0, status, COMPUTESKYH_ENEGA, 	COMPUTESKYH_MSGENEGA);
    }

  /* Check to make sure pointer to output is not NULL */
  ASSERT(skyConst!=NULL, status, COMPUTESKYH_ENNUL, COMPUTESKYH_MSGENNUL);

  params->baryinput->alpha=params->skyPos[iSkyCoh];
  params->baryinput->delta=params->skyPos[iSkyCoh+1];

  /* NOTE: we DO NOT translate the start-time into the SSB frame,
   * as this would result in a source-position dependent reference-time.
   * Instead we simply take the GPS start-time and interpret it as the
   * SSB reference-time
   */
  /*
    params->baryinput->tgps.gpsSeconds=params->tGPS[0].gpsSeconds;
    params->baryinput->tgps.gpsNanoSeconds=params->tGPS[0].gpsNanoSeconds;

    LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat);
    LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth);
    TimeToFloat(&tB0, &(params->emit->te));
  */
  TimeToFloat(&tB0, &(params->tGPS[0]));


  for (n=0; n<params->mObsSFT; n++)
   {
     t=(REAL8)(params->tGPS[n].gpsSeconds)+(REAL8)(params->tGPS[n].gpsNanoSeconds)*1.0E-9+0.5*params->tSFT;

     FloatToTime(&(params->baryinput->tgps), &t);

     TRY( LALBarycenterEarth(status->statusPtr, params->earth, &(params->baryinput->tgps), params->edat),
	  status);
     TRY( LALBarycenter(status->statusPtr, params->emit, params->baryinput, params->earth), status);

     TimeToFloat(&tBary, &(params->emit->te));

     dTbary = tBary-tB0;

     for (m=0; m<params->spinDwnOrder+1; m++)
       {
	 basedTbary = pow(dTbary, (REAL8)m);
	 skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m]=1.0/((REAL8)m+1.0)*basedTbary*dTbary-0.5*params->tSFT*params->emit->tDot*basedTbary;
	 skyConst[2*n*(params->spinDwnOrder+1)+2*(INT4)m+1]= params->tSFT*params->emit->tDot*basedTbary;
       }
   }
  /* Normal Exit */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/* Internal routines */
static void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps)
{
  INT4 x, y;

  x=tgps->gpsSeconds;
  y=tgps->gpsNanoSeconds;
  *f=(REAL8)x+(REAL8)y*1.e-9;
}

static void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f)
{
  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;

  temp0 = floor(*f);     /* this is tgps.S */
  temp1 = (*f) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2);
  temp4 = (temp2-temp3) * 0.1;

  tgps->gpsSeconds = (INT4)temp0;
  tgps->gpsNanoSeconds = (INT4)temp4;
}
