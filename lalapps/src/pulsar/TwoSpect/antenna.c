
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
   
   for (ii=0; ii<binshifts->length; ii++) binshifts->data[ii] = (INT4)roundf(dopplerMultiplier*freq*velocities->data[ii]*Tcoh);
   
   return binshifts;
   
}


REAL4Vector * CompAntennaPatternWeights(REAL4 ra, REAL4 dec, REAL8 t0, REAL4 Tcoh, REAL8 Tobs, LALDetector det)
{
   
   INT4 ii;
   INT4 numffts = (INT4)floor(2*(Tobs/Tcoh)-1);    //Number of FFTs
   REAL8 fplus, fcross;
   REAL4Vector *antptrnweights = XLALCreateREAL4Vector((UINT4)numffts);
   
   for (ii=0; ii<numffts; ii++) {
      LIGOTimeGPS time = {0,0};
      time.gpsSeconds = (INT4)floor(t0+(ii+1)*Tcoh*0.5);
      time.gpsNanoSeconds = (INT4)floor((t0+(ii+1)*Tcoh*0.5 - floor(t0+(ii+1)*Tcoh*0.5))*1e9);
      REAL8 gmst = XLALGreenwichMeanSiderealTime(&time);
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
   for (ii=0; ii<velocities->length; ii++) {
      LIGOTimeGPS time = {0,0};
      time.gpsSeconds = (INT4)floor(t0+(ii+1)*Tcoh*0.5);
      time.gpsNanoSeconds = (INT4)floor((t0+(ii+1)*Tcoh*0.5 - floor(t0+(ii+1)*Tcoh*0.5))*1e9);
   
      LALDetectorVel(&status, detvel, &time, det, edat);
      velocities->data[ii] = (REAL4)(detvel[0]*cosf(ra)*cosf(dec) + detvel[1]*sinf(ra)*cosf(dec) + detvel[2]*sinf(dec));
   }
   
   return velocities;
   
}


