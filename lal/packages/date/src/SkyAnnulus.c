/*
*  Copyright (C) 2007 Alexander Dietz, Jolien Creighton, Robert Adam Mercer
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
 * File Name: SkyAnnulus.c
 *
 * Author: Dietz, A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SkyAnnulusCV">
Author: Dietz, A.
$Id$
</lalVerbatim>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetectorSite.h>
#include <lal/AVFactories.h>
#include <lal/SkyAnnulus.h>

NRCSID( SKYANNULUSSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SkyAnnulus.c}}

\noindent
Provides methods to calculate the annulus given two or three detectors with a certain time difference.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SkyAnnulusCP}
\idx{LALComputeSingleAnnulus(
		   LALStatus          *status,
		   INT2Vector         *sky,
		   LALPlaceAndGPS     *placetime1,
		   LALPlaceAndGPS     *placetime2,
		   REAL4              dtErr,
		   INT2               nx,
		   INT2               ny)}
\idx{LALComputeTripleAnnuli( LALStatus          *status,
		       INT2Vector         **sky,
		       LALPlaceAndGPS     *placetime1,
		       LALPlaceAndGPS     *placetime2,
		       LALPlaceAndGPS     *placetime3,
		       REAL4              dtErr,
		       INT2               nx,
		       INT2               ny)}
\idx{LALExistTripleAnnuli( LALStatus
			   *status,
		     INT2               *exist,
		     LALPlaceAndGPS     *placetime1,
		     LALPlaceAndGPS     *placetime2,
		     LALPlaceAndGPS     *placetime3,
		     REAL4              dtErr,
		     INT2               nx,
		     INT2               ny)}

\subsubsection*{Description}

\noident
\texttt{LALComputeSingleAnnulus} calculates the single annulus when two detectors and two times are given in \texttt{placetime1} and \texttt{placetime2} assuming an error in time resolution of \texttt{dtErr}. The result is a vector \texttt{sky} with dimensions \texttt{nx}x\texttt{ny} representing the whole sky. The index of the vector ($\mathrm{index} = i + j\cdot nx$) represents the coordinates ($\alpha=360^\circ \cdot i / nx$ and $\delta= 180^\circ \cdot i / ny -90$.
If the desired position lies within the annulus the element of the vector is 1, 0 else.
\texttt{LALComputeTripleAnnuli} calculates the annuli when {\em three} detectors and times are given. The result is a vector of four vectors, with the first vector (index 0) containing the complete overlap of all annuli and the next three vectors (indixes 1,2,3) containing the annuli of two of the time intervals each.
\texttt{LALExistTripleAnnuli} just checks if there is at least one position in the sky overlapped by all three annuli created with \texttt{LALComputeTripleAnnuli}.
If that is the case, \texttt{exist} will contain the value 1, 0 else.


\subsubsection*{Algorithm}

\noindent

\subsubsection*{Uses}

\noindent For each module it is not needed to allocate some memory for the structure \option{sky}, this is done in the code.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{SkyAnnulusCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="SkyAnnulusCP"> */
void
LALComputeSingleAnnulus (
		   LALStatus          *status,
		   INT2Vector         *sky,
		   LALPlaceAndGPS     *placetime1,
		   LALPlaceAndGPS     *placetime2,
		   REAL4              dtErr,
		   INT2               nx,
		   INT2               ny
		  )
     /* </lalVerbatim> */
{
  INT4 i,j, myindex;
  REAL4 alpha, delta;
  REAL8 dt, dtLow, dtUpp;
  REAL8 deltat;

  SkyPosition coords;
  TwoDetsTimeAndASource tdtaas;


  INITSTATUS( status, "LALComputeSingleAnnulus", SKYANNULUSSC );

  /* create vector */
  LALI2CreateVector( status, &sky, nx*ny);

  /* calculate allowed area of time difference */
  dt = XLALGPSDiff(placetime2->p_gps, placetime1->p_gps);
  dtLow=dt-dtErr;
  dtUpp=dt+dtErr;

  /* loop over sky positions */
  for ( i=0;i<nx;i++ ) {
    for ( j=0;j<ny;j++ ) {

      /* create coordinates in units of degrees */
      alpha= ((REAL4) i)*360./((REAL4) nx);
      delta= ((REAL4) j)*180.0/((REAL4) ny)-90.0;

      coords.longitude=alpha*LAL_PI_180;
      coords.latitude=delta*LAL_PI_180;
      coords.system=COORDINATESYSTEM_EQUATORIAL;
      tdtaas.p_det_and_time1=placetime1;
      tdtaas.p_det_and_time2=placetime2;
      tdtaas.p_source=&coords;

      /* check if this coordinate is in position */
      LALTimeDelay(status, &deltat, &tdtaas);

      /* check if deltat is in range of dt */
      myindex=i+j*nx;
      if ( dtLow<deltat && deltat<dtUpp ) {
	sky->data[myindex]=1;
      } else {
	sky->data[myindex]=0;
      }

    }
  }

  RETURN( status );

}


																																																		/* <lalVerbatim file="SkyAnnulusCP">  */
void
LALComputeTripleAnnuli (
		       LALStatus          *status,
		       INT2Vector         **sky,
		       LALPlaceAndGPS     *placetime1,
		       LALPlaceAndGPS     *placetime2,
		       LALPlaceAndGPS     *placetime3,
		       REAL4              dtErr,
		       INT2               nx,
		       INT2               ny
		  )
  /* </lalVerbatim>  */
{
  INT4 i,j, myindex;

  INITSTATUS( status, "LALComputeTripleAnnuli", SKYANNULUSSC );

  LALI2CreateVector( status, &sky[0], nx*ny);
  LALI2CreateVector( status, &sky[1], nx*ny);
  LALI2CreateVector( status, &sky[2], nx*ny);
  LALI2CreateVector( status, &sky[3], nx*ny);

  /* calculate the three sky anulli */
  LALComputeSingleAnnulus( status, sky[1], placetime1, placetime2, dtErr, nx, ny);
  LALComputeSingleAnnulus( status, sky[2], placetime2, placetime3, dtErr, nx, ny);
  LALComputeSingleAnnulus( status, sky[3], placetime1, placetime3, dtErr, nx, ny);

  /* compute triple coincidence */
  for ( i=0;i<nx;i++ ) {
    for ( j=0;j<ny;j++ ) {
      myindex=i+j*nx;

      if ( sky[1]->data[myindex]==1 && sky[2]->data[myindex]==1 &&  sky[3]->data[myindex]==1 ) {
	sky[0]->data[myindex]=1;
      } else {
	sky[0]->data[myindex]=0;
      }

    }
  }

  RETURN( status );

}



/* <lalVerbatim file="SkyAnnulusCP"> */
void
LALExistTripleAnnuli(
		     LALStatus          *status,
		     INT2               *exist,
		     LALPlaceAndGPS     *placetime1,
		     LALPlaceAndGPS     *placetime2,
		     LALPlaceAndGPS     *placetime3,
		     REAL4              dtErr,
		     INT2               nx,
		     INT2               ny
		     )
     /* </lalVerbatim> */
{
  INT4 i,j,myindex;
  static INT2Vector* sky[4];

  INITSTATUS( status, "LALExistTripleAnnuli", SKYANNULUSSC );

  /* compute triple anulli */
  LALComputeTripleAnnuli( status, sky, placetime1, placetime2, placetime3,
			 dtErr, nx, ny);

  *exist=0;

  /* check result */
  for ( i=0;i<nx;i++ ) {
    for ( j=0;j<ny;j++ ) {
      myindex=i+j*nx;

      if (sky[0]->data[myindex]==1 ) {
	*exist=1;
	continue;
      }

    }
  }

  RETURN( status );
}
