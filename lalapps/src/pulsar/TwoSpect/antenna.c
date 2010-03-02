/*
*  Copyright (C) 2010 Evan Goetz
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

#include <math.h>
#include <lal/Date.h>
#include "antenna.h"


EphemerisData * new_Ephemeris(CHAR *earth_ephemeris, CHAR *sun_ephemeris)
{
   
   EphemerisData *ephemdata = (EphemerisData*)XLALMalloc(sizeof(EphemerisData));
   ephemdata->ephiles.earthEphemeris = earth_ephemeris;
   ephemdata->ephiles.sunEphemeris = sun_ephemeris;
   initEphemeris(ephemdata);
   
   return ephemdata;
   
}


void free_Ephemeris(EphemerisData *ephemdata)
{
   
   XLALFree((EphemerisData*)ephemdata);
   
}


void initEphemeris(EphemerisData *ephemdata)
{
   
   LALStatus status;
   status.statusPtr = NULL;
   
   LALInitBarycenter(&status, ephemdata);
   
}


INT4Vector * CompBinShifts(REAL4 freq, REAL4Vector *velocities, REAL4 Tcoh, REAL4 dopplerMultiplier)
{
   
   INT4Vector *binshifts = XLALCreateINT4Vector(velocities->length);
   INT4 ii;
   
   for (ii=0; ii<(INT4)binshifts->length; ii++) binshifts->data[ii] = (INT4)roundf(dopplerMultiplier*freq*velocities->data[ii]*Tcoh);
   
   return binshifts;
   
}


REAL4Vector * CompAntennaPatternWeights(REAL4 ra, REAL4 dec, REAL8 t0, REAL4 Tcoh, REAL8 Tobs, LALDetector det)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(2*(Tobs/Tcoh)-1);    //Number of FFTs
   REAL8 fplus, fcross;
   REAL4Vector *antptrnweights = XLALCreateREAL4Vector((UINT4)numffts);
   
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0+(ii+1)*Tcoh*0.5);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+(ii+1)*Tcoh*0.5 - floor(t0+(ii+1)*Tcoh*0.5))*1e9);
      REAL8 gmst = XLALGreenwichMeanSiderealTime(&gpstime);
      XLALComputeDetAMResponse(&fplus, &fcross, det.response, ra, dec, 0.0, gmst);
      antptrnweights->data[ii] = (REAL4)(fplus*fplus + fcross*fcross);
   }
   
   return antptrnweights;

}



REAL4Vector * CompAntennaVelocity(REAL4 ra, REAL4 dec, REAL8 t0, REAL4 Tcoh, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(2*(Tobs/Tcoh)-1);    //Number of FFTs
   LALStatus status;
   status.statusPtr = NULL;
   
   REAL4Vector *velocities = XLALCreateREAL4Vector((UINT4)numffts);
   
   REAL8 detvel[3];
   for (ii=0; ii<(INT4)velocities->length; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0+(ii+1)*Tcoh*0.5);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+(ii+1)*Tcoh*0.5 - floor(t0+(ii+1)*Tcoh*0.5))*1e9);
   
      LALDetectorVel(&status, detvel, &gpstime, det, edat);
      velocities->data[ii] = (REAL4)(detvel[0]*cosf(ra)*cosf(dec) + detvel[1]*sinf(ra)*cosf(dec) + detvel[2]*sinf(dec));
   }
   
   return velocities;
   
}


REAL4 CompDetectorDeltaVmax(REAL8 t0, REAL4 Tcoh, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(2*(Tobs/Tcoh)-1);    //Number of FFTs
   LALStatus status;
   status.statusPtr = NULL;
   
   REAL8 detvel[3];
   REAL8 detvel0[3];
   REAL8 dv[3];
   REAL4 deltaVmax = 0.0;
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0+(ii+1)*Tcoh*0.5);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+(ii+1)*Tcoh*0.5 - floor(t0+(ii+1)*Tcoh*0.5))*1e9);
   
      if (ii==0) LALDetectorVel(&status, detvel0, &gpstime, det, edat);
      else {
         LALDetectorVel(&status, detvel, &gpstime, det, edat);
         dv[0] = detvel[0] - detvel0[0];
         dv[1] = detvel[1] - detvel0[1];
         dv[2] = detvel[2] - detvel0[2];
         REAL4 deltaV = (REAL4)sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
         if (deltaV > deltaVmax) deltaVmax = deltaV;
      }
   }
   
   return deltaVmax;
   
}


