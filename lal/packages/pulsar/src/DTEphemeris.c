/*
*  Copyright (C) 2007 David M. Whitbeck, Ian Jones, Reinhard Prix
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

/******************************* <lalVerbatim file="DTEphemerisCV">
Author: Jones, D. I.,   Owen, B. J.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{DTEphemeris.c}}
\label{ss:DTEphemeris.c}

Computes the barycentric arrival time of an incoming wavefront using
accurate ephemeris-based data files of the Sun and Earth's motions.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DTEphemerisCP}
\idx{LALDTEphemeris()}
\idx{LALTEphemeris()}

\subsubsection*{Description}


These routines compute the barycentric time transformation and its
derivatives.  That is, if a signal originating from a right ascension
$\alpha$ and declination $\delta$ on the sky and arrives at the
detector at a time $t$, then it will pass the centre of the solar
system at a time $t_b(t,\alpha,\delta)$.

The input/output features of this function are nearly identical to
those of \texttt{DTBaryPtolemaic()}, whose documentation should be
consulted  for the details. One important difference in calling this
function is that the user has to supply the initialised ephemeris-data
in the \verb+PulsarTimesParamStruc->ephemeris+ and the detector-data
in \verb+PulsarTimesParamStruc->site+.

\texttt{DTBaryPtolemaic()} uses the Ptolemaic approximation to model
the Earth/Sun system, while \texttt{DTEphemeris()} uses accurate
ephemeris data read in from files in the calling function, and passed
into \texttt{DTEphemeris()} using the \texttt{EphemerisData}
structure, which is a member of the \texttt{PulsarTimesParamStruc}.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel                LALBarycenterEarth()
LALBarycenter()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DTEphemerisCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/PtoleMetric.h>
#include <lal/StackMetric.h>
#include <lal/Date.h>


NRCSID(DTEPHEMERISC,"$Id$");

/* <lalVerbatim file="DTEphemerisCP"> */
void
LALDTEphemeris( LALStatus             *status,
	        REAL8Vector           *drv,
	        REAL8Vector           *var,
	        PulsarTimesParamStruc *tev )
{ /* </lalVerbatim> */
  LIGOTimeGPS tGPS;           /* Input structure to BartcenterEarth()  */
  const EphemerisData *eph;         /* Input structure to BarycenterEarth()  */
  EarthState earth;           /* Output structure of BarycenterEarth() */
                              /* and input structure to Barycenter()   */
  BarycenterInput baryin;     /* Input structure for Barycenter()      */
  EmissionTime emit;          /* Output structure of Barycenter()      */
  UINT4 numDeriv; 	 /* number of derivatives to compute */

  REAL8 upper, lower;         /* Quantities for finite differnecing */
  REAL8 d_alpha, d_delta;

  INITSTATUS(status,"DTEphemeris",DTEPHEMERISC);
  ATTATCHSTATUSPTR( status );

  /* Make sure parameter structures and their fields exist. */
  ASSERT(drv,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(drv->data,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(var,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(var->data,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(tev,status, PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  /* Make sure array sizes are consistent. */
  /* need at least [ t, alpha, delta ] */
  ASSERT( var->length >= 3, status, PULSARTIMESH_EBAD, PULSARTIMESH_MSGEBAD);
  /* result has to be at least [T(t), and optionally dT/dt, dT/dalpha, dT/ddelta] */
  ASSERT( drv->length >= 1, status, PULSARTIMESH_EBAD, PULSARTIMESH_MSGEBAD);

  /* Make sure ephermis and detector data have been passed */
  ASSERT (tev->ephemeris != NULL, status, PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT (tev->site != NULL, status, PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);

  numDeriv = drv->length - 1; /* number of derivatives of T to compute..*/

  /* First compute the location, velocity, etc... of the Earth:  */

  /* Set the GPS time: */
  tGPS = tev->epoch;
  XLALGPSAdd ( &tGPS, var->data[0] );	/* time relative to epoch */

  /* Set the ephemeris data: */
  eph = tev->ephemeris;

  TRY( LALBarycenterEarth( status->statusPtr, &earth, &tGPS, eph ), status );
  /* Now "earth" contains position of center of Earth. */

  /* Now do the barycentering.  Set the input parameters: */

  /* Get time delay for detector vertex. */
  baryin.tgps = tGPS;

  /* Set the detector site...*/
  baryin.site = *(tev->site);

  /* ...remembering to divide the coordinates by the speed of light: */
  baryin.site.location[0] /= LAL_C_SI;
  baryin.site.location[1] /= LAL_C_SI;
  baryin.site.location[2] /= LAL_C_SI;

  baryin.alpha = var->data[1];
  baryin.delta = var->data[2];

  /* Set 1/distance to zero: */
  baryin.dInv = 0.e0;

  TRY( LALBarycenter( status->statusPtr, &emit, &baryin, &earth ), status );
  /* Now "emit" contains detector position, velocity, time, tdot. */

  /* Now assemble output: */

  /* Subtract off epoch: => barycentered time since epoch */
  drv->data[0] = XLALGPSDiff( &(emit.te), &(tev->epoch) );

  /* ---------- calculate only the requested derivatives ---------- */

  /* ----- derivative dT/dt ----- */
  if ( numDeriv >= 1 )
    {
      drv->data[1] = emit.tDot; /* dtb/dt */
    }
  d_alpha = d_delta = 0.001;	      /* Set finite difference step sizes: */
  /* ----- dT/dalpha ----- */
  if ( numDeriv >= 2 )
    {
      /* Need to finite difference to get d(tb)/d(alpha), d(tb)/d(delta) */

      /* Get dtb/da by finite differencing.   */
      /* Default upper and lower alpha values: */
      upper = var->data[1] + d_alpha;
      lower = var->data[1] - d_alpha;
      /* Overwrite if alpha is too close to zero or 2 PI: */
      if(var->data[1] < d_alpha)
	lower = var->data[1];
      if(var->data[1] > (LAL_TWOPI-d_alpha))
	upper = var->data[1];
      /* Evaluate emit at upper value: */
      baryin.alpha = upper;
      TRY( LALBarycenter( status->statusPtr, &emit, &baryin, &earth ), status );
      drv->data[2] = emit.te.gpsSeconds + 1e-9*emit.te.gpsNanoSeconds;
      /* Evaluate emit at lower value: */
      baryin.alpha = lower;
      TRY( LALBarycenter( status->statusPtr, &emit, &baryin, &earth ), status );
      drv->data[2] -= emit.te.gpsSeconds + 1e-9*emit.te.gpsNanoSeconds;
      /* Divide by alpha interval: */
      drv->data[2] /= (upper-lower);
      baryin.alpha = var->data[1];
    } /* if numDeriv >= 2 */
  /* ----- dT/ddelta ----- */
  if ( numDeriv >= 3 )
    {
      /* Get dtb/dd by finite differencing.   */
      /* Default upper and lower alpha values: */
      upper = var->data[2] + d_delta;
      lower = var->data[2] - d_delta;
      /* Overwrite if delta is too close to PI/2 or -PI/2: */
      if(var->data[2] < (-LAL_PI_2+d_alpha))
	lower = var->data[2];
      if(var->data[2] > (LAL_PI_2-d_alpha))
	upper = var->data[2];
      /* Evaluate emit at upper value: */
      baryin.delta = upper;
      TRY( LALBarycenter( status->statusPtr, &emit, &baryin, &earth ), status );
      drv->data[3] = emit.te.gpsSeconds + 1e-9*emit.te.gpsNanoSeconds;
      /* Evaluate emit at lower value: */
      baryin.delta = lower;
      TRY( LALBarycenter( status->statusPtr, &emit, &baryin, &earth ), status );
      drv->data[3] -= emit.te.gpsSeconds + 1e-9*emit.te.gpsNanoSeconds;
      /* Divide by delta interval: */
      drv->data[3] /= (upper-lower);
      baryin.alpha = var->data[2];
    } /* if numDeriv >= 3 */

  /* Go home */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*Computes the barycentric time using Ephemeris data.*/

/* <lalVerbatim file="DTEphemerisCP"> */
void
LALTEphemeris( LALStatus   *status,
	       REAL8 *tBary,
	       REAL8Vector *var,
	       PulsarTimesParamStruc *tev )
{ /* </lalVerbatim>*/
  LIGOTimeGPS tGPS; /* Input structure to BarycenterEarth() */
  const EphemerisData *eph; /* Input structure to BarycenterEarth() */
  EarthState earth; /* Output structure of BarycenterEarth() */
  BarycenterInput baryin; /* Input structure for Barycenter() */
  EmissionTime emit; /*Output structure of Barycenter() */

  INITSTATUS(status,"TEphemeris",DTEPHEMERISC);
  ATTATCHSTATUSPTR( status );

  /*Make sure that param structs and their fields exist. */
  ASSERT(tBary,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(var,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(var->data,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(tev,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  /*Make sure that array sizes are consistent. */
  ASSERT(var->length>2,status,PULSARTIMESH_EBAD,PULSARTIMESH_MSGEBAD);
  /*Make sure ephemeris and detector data have been passed. */
  ASSERT(tev->ephemeris!=NULL,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  ASSERT(tev->site!= NULL,status,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);

  /* set the GPS time */
  tGPS.gpsSeconds=floor(var->data[0])+tev->epoch.gpsSeconds;
  tGPS.gpsNanoSeconds=1e9*fmod(var->data[0],1.0)+tev->epoch.gpsNanoSeconds;

  /* set the ephemeris data */
  eph=tev->ephemeris;

  /* journey to the center of the Earth */
  TRY(LALBarycenterEarth(status->statusPtr,&earth,&tGPS,eph),status);

  /* time delay for detector vertex */
  baryin.tgps.gpsSeconds=tGPS.gpsSeconds;
  baryin.tgps.gpsNanoSeconds=tGPS.gpsNanoSeconds;

  /* set detector site */
  baryin.site=*(tev->site);

  /*divide by speed of light */
  baryin.site.location[0]/=LAL_C_SI;
  baryin.site.location[1]/=LAL_C_SI;
  baryin.site.location[2]/=LAL_C_SI;

  /* sky positions */
  baryin.alpha=var->data[1];
  baryin.delta=var->data[2];

  /* set 1/distance to zero */
  baryin.dInv=0.e0;

  /* this gives emit the position, velocity, time etc */
  TRY(LALBarycenter(status->statusPtr,&emit,&baryin,&earth),status);

  *tBary = emit.te.gpsSeconds+1.0e-9*emit.te.gpsNanoSeconds;
  *tBary -= tev->epoch.gpsSeconds+1.0e-9*tev->epoch.gpsNanoSeconds;

  RETURN(status);
}
