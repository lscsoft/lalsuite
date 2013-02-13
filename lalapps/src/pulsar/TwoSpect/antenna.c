/*
*  Copyright (C) 2010, 2012, 2013 Evan Goetz
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


static const LALStatus empty_status;

//Compute the number of integer bin shifts per SFT
// bin shift = f0*v*Tcoh
// where f0 is frequency, v is velocity in units of c, and Tcoh is the SFT coherence length
// an optional dopplerMultiplier value could be multiplied if desired (default value is 1.0)
void CompBinShifts(INT4Vector *output, REAL8 freq, REAL4Vector *velocities, REAL8 Tcoh, REAL4 dopplerMultiplier)
{
   
   INT4 ii;
   
   for (ii=0; ii<(INT4)velocities->length; ii++) output->data[ii] = (INT4)round(dopplerMultiplier*freq*velocities->data[ii]*Tcoh);
   
} /* CompBinShifts() */


//Compute the antenna pattern weights
//If linPolOn = 1, then the output weights are Fplus*Fplus for the given polarization angle
//If linPolOn = 0, then the output weights are Fplus*Fplus + Fcross*Fcross
void CompAntennaPatternWeights(REAL4Vector *output, REAL4 ra, REAL4 dec, REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, INT4 linPolOn, REAL8 polAngle, LALDetector det)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   REAL8 fplus, fcross;
   
   for (ii=0; ii<numffts; ii++) {
      
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
      REAL8 gmst = XLALGreenwichMeanSiderealTime(&gpstime);
      if (XLAL_IS_REAL8_FAIL_NAN(gmst)) {
         fprintf(stderr,"%s: XLALGreenwichMeanSiderealTime(%.9d.%.9d) failed.\n", __func__, gpstime.gpsSeconds, gpstime.gpsNanoSeconds);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      XLALComputeDetAMResponse(&fplus, &fcross, det.response, ra, dec, polAngle, gmst);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: XLALComputeDetAMResponse() failed.\n", __func__);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      if (!linPolOn) output->data[ii] = (REAL4)(fplus*fplus + fcross*fcross);
      else output->data[ii] = (REAL4)(fplus*fplus);
      
   } /* for ii < numffts */

} /* CompAntennaPatternWeights() */


//Compute the antenna velocity
void CompAntennaVelocity(REAL4Vector *output, REAL4 ra, REAL4 dec, REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   LALStatus status = empty_status;
   
   REAL8 detvel[3];
   for (ii=0; ii<numffts; ii++) {
      
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
   
      LALDetectorVel(&status, detvel, &gpstime, det, edat);
      if (status.statusCode!=0) {
         fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", __func__, status.statusCode);
         XLAL_ERROR_VOID(XLAL_EFUNC);
      }
      
      output->data[ii] = (REAL4)(detvel[0]*cos(ra)*cos(dec) + detvel[1]*sin(ra)*cos(dec) + detvel[2]*sin(dec));
      
   } /* for ii < numffts */
   
} /* CompAntennaVelocity() */


//Determine the maximum change in velocity
REAL4 CompDetectorDeltaVmax(REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   LALStatus status = empty_status;
   
   REAL8 detvel[3];
   REAL8 detvel0[3];
   REAL8 dv[3];
   REAL4 deltaVmax = 0.0;
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
   
      if (ii==0) {
         LALDetectorVel(&status, detvel0, &gpstime, det, edat);
         if (status.statusCode!=0) {
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", __func__, status.statusCode);
            XLAL_ERROR_REAL4(XLAL_EFUNC);
         }
      } else {
         LALDetectorVel(&status, detvel, &gpstime, det, edat);
         if (status.statusCode!=0) {
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", __func__, status.statusCode);
            XLAL_ERROR_REAL4(XLAL_EFUNC);
         }
         dv[0] = detvel[0] - detvel0[0];
         dv[1] = detvel[1] - detvel0[1];
         dv[2] = detvel[2] - detvel0[2];
         REAL4 deltaV = (REAL4)sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
         if (deltaV > deltaVmax) deltaVmax = deltaV;
      } /* if ii==0 else ... */
   } /* for ii < numffts */
   
   return deltaVmax;
   
} /* CompDetectorDeltaVmax() */


//From a given t0 start time, determine the maximum velocity over the observation time
REAL4 CompDetectorVmax(REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   LALStatus status = empty_status;
   
   REAL8 detvel[3];
   REAL8 detvel0[3];
   REAL4 Vmax = 0.0;
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
      
      if (ii==0) {
         LALDetectorVel(&status, detvel0, &gpstime, det, edat);
         if (status.statusCode!=0) {
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", __func__, status.statusCode);
            XLAL_ERROR_REAL4(XLAL_EFUNC);
         }
         Vmax = (REAL4)sqrt(detvel[0]*detvel[0] + detvel[1]*detvel[1] + detvel[2]*detvel[2]);
      } else {
         LALDetectorVel(&status, detvel, &gpstime, det, edat);
         if (status.statusCode!=0) {
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", __func__, status.statusCode);
            XLAL_ERROR_REAL4(XLAL_EFUNC);
         }
         REAL4 V = (REAL4)sqrt(detvel[0]*detvel[0] + detvel[1]*detvel[1] + detvel[2]*detvel[2]);
         if (V > Vmax) Vmax = V;
      } /* if ii==0 else ... */
   } /* for ii < numffts */
   
   return Vmax;
   
} /* CompDetectorVmax() */


