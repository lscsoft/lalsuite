/*
*  Copyright (C) 2007 Jolien Creighton, Sean Seader
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
 * File Name: CoherentInspiralInput.c
 *
 * Author: Seader, S. E.  Brown, D.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <time.h>
#include <math.h>

#include <lal/LALRCSID.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>
#include <lal/DetectorSite.h>
#include <lal/Random.h>
#include <lal/LALInspiral.h>
#include <lal/CoherentInspiral.h>
#include <lal/LALStatusMacros.h>

NRCSID( COHERENTINSPIRALINPUTC, "$Id$");

#define rint(x) (floor((x)+0.5))
double modf( double value, double *integerPart );

void
LALFindChirpCreateCoherentInput(
     LALStatus                  *status,
     COMPLEX8TimeSeries         **coherentInputData,
     COMPLEX8TimeSeries         *input,
     SnglInspiralTable          *templt,
     REAL4                      coherentSegmentLength,/*in seconds*/
     INT4                       corruptedDataLength /*in timepoints */
     )
{
  COMPLEX8TimeSeries      *cohInputData = NULL;
  LIGOTimeGPS              end_time;
  UINT8                    eventID = 0;
  INT4                     numPoints = 0;
  REAL4                    cohSegLength = 0.0;
  INT4                     inputEpochSeconds = 0;
  INT4                     inputEpochNanoSeconds = 0;
  REAL8                    deltaT = 0.0;
  REAL8                    intpart = 0.0;
  REAL8                    fracpart = 0.0;
  REAL8                    tempTime = 0.0;
  INT4                     crupDataLength = 0;
  INT4                     cohSegEnd = 0;
  INT4                     cohSegStart = 0;
  INT4                     nonCorruptEnd = 0;
  INT4                     nonCorruptStart = 0;
  /*INT4                     overlap         = 0;*/
  INT4                     fullCohSegLength = 0;
  INT4                     eventTimePoint = 0;

  INITSTATUS( status, "LALFindChirpCreateCoherentInput",
	      COHERENTINSPIRALINPUTC );
  ATTATCHSTATUSPTR( status );

  /* Ensure that arguments are reasonable */

  /* check that the output handle exists, but points to a null pointer */
  ASSERT( coherentInputData, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*coherentInputData, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* check that a valid COMPLEX8TimeSeries is input */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->data->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->deltaT > 0, status, FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( input->epoch.gpsSeconds > 0, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->epoch.gpsNanoSeconds >= 0, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->data->length, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that a valid snglInspiralTable is input */
  ASSERT( templt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( templt->end_time.gpsSeconds > 0, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( templt->end_time.gpsNanoSeconds >= 0, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check for valid lengths */
  ASSERT( coherentSegmentLength > 0, status, FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ );
  ASSERT( corruptedDataLength >= 0, status, FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ );


  /* Get necessary info from input structures */

  end_time = templt->end_time;
  eventID  = templt->event_id->id;
  numPoints = input->data->length;
  cohSegLength = coherentSegmentLength;
  inputEpochSeconds = input->epoch.gpsSeconds;
  inputEpochNanoSeconds = input->epoch.gpsNanoSeconds;
  deltaT = input->deltaT;
  crupDataLength = corruptedDataLength;

  /* Now determine if the event lies in corrupted data */

  nonCorruptStart = crupDataLength;
  nonCorruptEnd = numPoints - crupDataLength;
  cohSegStart = (INT4) rint( ((end_time.gpsSeconds + end_time.gpsNanoSeconds*1.0e-9) - (inputEpochSeconds + inputEpochNanoSeconds*1.0e-9)-cohSegLength)/deltaT - 1.0);
  cohSegEnd = cohSegStart + 2* (INT4)cohSegLength/deltaT + 1;
  eventTimePoint = (INT4) rint( ((end_time.gpsSeconds + end_time.gpsNanoSeconds*1.0e-9) - (inputEpochSeconds + inputEpochNanoSeconds*1.0e-9))/deltaT - 1.0);

  if( eventTimePoint < nonCorruptEnd && eventTimePoint >= nonCorruptStart )
    {
      /* The coherent chunk is outside of the corrupted data */

      cohInputData = *coherentInputData = (COMPLEX8TimeSeries *)
	LALCalloc(1, sizeof(COMPLEX8TimeSeries) );

      fullCohSegLength = 2*cohSegLength/deltaT;

      LALCCreateVector(status->statusPtr, &(cohInputData->data), fullCohSegLength);
      CHECKSTATUSPTR( status );
      /* store C-data snippet and its associated info */
      memcpy(cohInputData->name, input->name, LALNameLength * sizeof(CHAR) );
      memcpy(cohInputData->data->data, &(input->data->data[cohSegStart]), fullCohSegLength * sizeof(COMPLEX8) );
      cohInputData->deltaT = deltaT;
      /*CHECK: Below is a temporary fix for communicating sigmasq
	from inspiral to coherent_inspiral via the f0 member of
	the COMPLEX8TimeSeries structure*/
      cohInputData->f0 = (REAL8) templt->sigmasq;
      tempTime = inputEpochSeconds + inputEpochNanoSeconds*1.0e-9 + cohSegStart * deltaT;
      fracpart = modf(tempTime, &intpart);
      cohInputData->epoch.gpsSeconds = (INT4) intpart;
      cohInputData->epoch.gpsNanoSeconds = (INT4) (fracpart*1.0e9);
    }
  else
    {
      /* return a null pointer cohInputData */
    }
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
