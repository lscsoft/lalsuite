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



void CompBinShifts(INT4Vector *output, REAL8 freq, REAL4Vector *velocities, REAL8 Tcoh, REAL4 dopplerMultiplier)
{
   
   INT4 ii;
   
   for (ii=0; ii<(INT4)velocities->length; ii++) output->data[ii] = (INT4)round(dopplerMultiplier*freq*velocities->data[ii]*Tcoh);
   
} /* CompBinShifts() */



void CompAntennaPatternWeights(REAL4Vector *output, REAL4 ra, REAL4 dec, REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   REAL8 fplus, fcross;
   
   for (ii=0; ii<numffts; ii++) {
      
      LIGOTimeGPS gpstime = {0,0};
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tcoh-SFToverlap) + 0.5*Tcoh);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh - floor(t0+ii*(Tcoh-SFToverlap)+0.5*Tcoh))*1e9);
      REAL8 gmst = XLALGreenwichMeanSiderealTime(&gpstime);
      if (XLAL_IS_REAL8_FAIL_NAN(gmst)) {
         fprintf(stderr,"%s: XLALGreenwichMeanSiderealTime(%.9d.%.9d) failed.\n", fn, gpstime.gpsSeconds, gpstime.gpsNanoSeconds);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      XLALComputeDetAMResponse(&fplus, &fcross, det.response, ra, dec, 0.0, gmst);
      if (xlalErrno!=0) {
         fprintf(stderr,"%s: XLALComputeDetAMResponse() failed.\n", fn);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      output->data[ii] = (REAL4)(fplus*fplus + fcross*fcross);
      
   } /* for ii < numffts */

} /* CompAntennaPatternWeights() */



void CompAntennaVelocity(REAL4Vector *output, REAL4 ra, REAL4 dec, REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   const CHAR *fn = __func__;
   
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
      if (status.statusCode!=0) {
         fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", fn, status.statusCode);
         XLAL_ERROR_VOID(fn, XLAL_EFUNC);
      }
      
      output->data[ii] = (REAL4)(detvel[0]*cos(ra)*cos(dec) + detvel[1]*sin(ra)*cos(dec) + detvel[2]*sin(dec));
      
   } /* for ii < numffts */
   
} /* CompAntennaVelocity() */



REAL4 CompDetectorDeltaVmax(REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   const CHAR *fn = __func__;
   
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
   
      if (ii==0) {
         LALDetectorVel(&status, detvel0, &gpstime, det, edat);
         if (status.statusCode!=0) {
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", fn, status.statusCode);
            XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
         }
      } else {
         LALDetectorVel(&status, detvel, &gpstime, det, edat);
         if (status.statusCode!=0) {
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", fn, status.statusCode);
            XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
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


REAL4 CompDetectorVmax(REAL8 t0, REAL8 Tcoh, REAL8 SFToverlap, REAL8 Tobs, LALDetector det, EphemerisData *edat)
{
   
   const CHAR *fn = __func__;
   
   INT4 ii;
   INT4 numffts = (INT4)floor(Tobs/(Tcoh-SFToverlap)-1);    //Number of FFTs
   LALStatus status;
   status.statusPtr = NULL;
   
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
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", fn, status.statusCode);
            XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
         }
      } else {
         LALDetectorVel(&status, detvel, &gpstime, det, edat);
         if (status.statusCode!=0) {
            fprintf(stderr,"%s: LALDetectorVel() failed with error code %d.\n", fn, status.statusCode);
            XLAL_ERROR_REAL4(fn, XLAL_EFUNC);
         }
         REAL4 V = (REAL4)sqrt(detvel[0]*detvel[0] + detvel[1]*detvel[1] + detvel[2]*detvel[2]);
         if (V > Vmax) Vmax = V;
      } /* if ii==0 else ... */
   } /* for ii < numffts */
   
   return Vmax;
   
} /* CompDetectorVmax() */


