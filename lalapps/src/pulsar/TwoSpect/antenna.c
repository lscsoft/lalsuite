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
#include "antenna.h"


EphemerisData * new_Ephemeris(CHAR *earth_ephemeris, CHAR *sun_ephemeris)
{
   
   LALStatus status;
   status.statusPtr = NULL;
   EphemerisData *ephemdata = (EphemerisData*)XLALMalloc(sizeof(EphemerisData));
   ephemdata->ephiles.earthEphemeris = earth_ephemeris;
   ephemdata->ephiles.sunEphemeris = sun_ephemeris;
   initEphemeris(ephemdata);
   
   return ephemdata;
   
}



void free_Ephemeris(EphemerisData *ephemdata)
{
   
   XLALFree((CHAR*)ephemdata->ephiles.earthEphemeris);
   XLALFree((CHAR*)ephemdata->ephiles.sunEphemeris);
   XLALFree((PosVelAcc*)ephemdata->ephemE);
   XLALFree((PosVelAcc*)ephemdata->ephemS);
   XLALFree((EphemerisData*)ephemdata);
   
}



void initEphemeris(EphemerisData *ephemdata)
{
   
   LALStatus status;
   status.statusPtr = NULL;
   
   LALInitBarycenter(&status, ephemdata);
   
}



void CompBinShifts(INT4Vector *output, REAL8 freq, REAL4Vector *velocities, REAL8 Tcoh, REAL4 dopplerMultiplier)
{
   
   INT4 ii;
   
   for (ii=0; ii<(INT4)velocities->length; ii++) output->data[ii] = (INT4)round(dopplerMultiplier*freq*velocities->data[ii]*Tcoh);
   
}



void CompAntennaPatternWeights(REAL4Vector *output, REAL4 ra, REAL4 dec, REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   REAL8 fplus, fcross;
   
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
      REAL8 gmst = XLALGreenwichMeanSiderealTime(&gpstime);
      XLALComputeDetAMResponse(&fplus, &fcross, det.response, ra, dec, 0.0, gmst);
      output->data[ii] = (REAL4)(fplus*fplus + fcross*fcross);
   }

}



void CompAntennaVelocity(REAL4Vector *output, REAL4 ra, REAL4 dec, REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   LALStatus status;
   status.statusPtr = NULL;
   
   REAL8 detvel[3];
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
   
      LALDetectorVel(&status, detvel, &gpstime, det, edat);
      output->data[ii] = (REAL4)(detvel[0]*cos(ra)*cos(dec) + detvel[1]*sin(ra)*cos(dec) + detvel[2]*sin(dec));
   }
   
}



REAL4 CompDetectorDeltaVmax(REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   LALStatus status;
   status.statusPtr = NULL;
   
   REAL8 detvel[3];
   REAL8 detvel0[3];
   REAL8 dv[3];
   REAL4 deltaVmax = 0.0;
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
   
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


