/*
*  Copyright (C) 2010, 2012, 2013, 2014 Evan Goetz
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
#include <lal/SSBtimes.h>
#include <lal/DetResponse.h>
#include <lal/Velocity.h>
#include "antenna.h"

/**
 * \brief Compute the number of integer bin shifts per SFT
 *
 * bin shift = f0*v*Tsft, where f0 is frequency, v is velocity in units of c, and Tsft is the SFT coherence length.
 * An optional dopplerMultiplier value could be multiplied if desired (default value is 1.0)
 * \param [out] output            Pointer to INT4Vector of bin shift values
 * \param [in]  ssbTimes          Pointer to the interferometer-specific SSBtimes
 * \param [in]  freq              Frequency from which to compute the bin shifts
 * \param [in]  Tsft              Coherence length of the SFTs
 * \param [in]  dopplerMultiplier Multiplicative factor to increase or decrease the bin shifts (standard physics = 1.0)
 * \return Status value
 */
INT4 CompBinShifts(INT4Vector *output, const SSBtimes *ssbTimes, const REAL8 freq, const REAL8 Tsft, const REAL4 dopplerMultiplier)
{
   for (UINT4 ii=0; ii<output->length; ii++) output->data[ii] = (INT4)round(dopplerMultiplier*(ssbTimes->Tdot->data[ii]-1.0)*freq*Tsft);
   return XLAL_SUCCESS;
}

/**
 * \brief Compute the antenna pattern weights
 *
 * If linPolOn = 0, then the output weights are Fplus*Fplus + Fcross*Fcross,
 * or if linPolOn = 1, then the output weights are Fplus*Fplus for the given polarization angle
 * \param [out] output     Pointer to REAL4Vector of antenna pattern weights
 * \param [in]  skypos     Sky position measured in RA and Dec in radians
 * \param [in]  t0         Start time (GPS seconds)
 * \param [in]  Tsft       Coherence length of the SFTs (in seconds)
 * \param [in]  SFToverlap Overlap of the SFTs (in seconds)
 * \param [in]  Tobs       Observation time (in seconds)
 * \param [in]  linPolOn   Flag for linear polarizations (linPolOn = 0 is circular polarization weights, linPolOn = 1 is linear polarization weights)
 * \param [in]  polAngle   Only used for when linPolOn = 1. Polarization angle from which to compute the antenna pattern weights
 * \param [in]  det        A LALDetector struct
 * \return Status value
 */
INT4 CompAntennaPatternWeights(REAL4VectorAligned *output, const SkyPosition skypos, const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const BOOLEAN linPolOn, const REAL8 polAngle, const LALDetector det)
{

   XLAL_CHECK( output != NULL, XLAL_EINVAL );

   INT4 numffts = (INT4)floor(Tobs/(Tsft-SFToverlap)-1);    //Number of FFTs
   REAL8 fplus, fcross;

   for (INT4 ii=0; ii<numffts; ii++) {

      LIGOTimeGPS gpstime = LIGOTIMEGPSZERO;
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tsft-SFToverlap) + 0.5*Tsft);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tsft-SFToverlap)+0.5*Tsft - floor(t0+ii*(Tsft-SFToverlap)+0.5*Tsft))*1e9);
      REAL8 gmst = XLALGreenwichMeanSiderealTime(&gpstime);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

      XLALComputeDetAMResponse(&fplus, &fcross, (const REAL4(*)[3])det.response, skypos.longitude, skypos.latitude, polAngle, gmst);
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

      if (!linPolOn) output->data[ii] = (REAL4)(fplus*fplus + fcross*fcross);
      else output->data[ii] = (REAL4)(fplus*fplus);

   } /* for ii < numffts */

   return XLAL_SUCCESS;

} /* CompAntennaPatternWeights() */

INT4 CompAntennaPatternWeights2(REAL4VectorAligned *output, const SkyPosition skypos, const LIGOTimeGPSVector *timestamps, const LALDetector det, const REAL8 *cosi, const REAL8 *psi)
{

   XLAL_CHECK( output != NULL, XLAL_EINVAL );

   REAL8 onePlusCosiSqOver2Sq = 1.0, cosiSq = 1.0;
   if (cosi!=NULL) {
      onePlusCosiSqOver2Sq = 0.25*(1.0 + (*cosi)*(*cosi))*(1.0 + (*cosi)*(*cosi));
      cosiSq = (*cosi)*(*cosi);
   }

   REAL8 fplus, fcross;
   for (UINT4 ii=0; ii<timestamps->length; ii++) {
      REAL8 gmst = XLALGreenwichMeanSiderealTime(&(timestamps->data[ii]));
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

      if (psi!=NULL) {
         XLALComputeDetAMResponse(&fplus, &fcross, det.response, skypos.longitude, skypos.latitude, *psi, gmst);
         XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
         output->data[ii] = fplus*fplus*onePlusCosiSqOver2Sq + fcross*fcross*cosiSq;
      } else {
         output->data[ii] = 0.0;
         for (UINT4 jj=0; jj<16; jj++) {
            XLALComputeDetAMResponse(&fplus, &fcross, det.response, skypos.longitude, skypos.latitude, 0.0625*LAL_PI*jj, gmst);
            XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
            output->data[ii] += fplus*fplus*onePlusCosiSqOver2Sq + fcross*fcross*cosiSq;
         }
         output->data[ii] *= 0.0625;
      }
   }

   return XLAL_SUCCESS;

} // CompAntennaPatternWeights()

/**
 * Compute the antenna velocity
 * \param [out] output     Pointer to REAL4Vector of antenna velocities
 * \param [in]  skypos     Sky position measured in RA and Dec in radians
 * \param [in]  t0         Start time (GPS seconds)
 * \param [in]  Tsft       Coherence length of the SFTs (in seconds)
 * \param [in]  SFToverlap Overlap of the SFTs (in seconds)
 * \param [in]  Tobs       Observation time (in seconds)
 * \param [in]  det        A LALDetector struct
 * \param [in]  edat       Pointer to EphemerisData
 * \return Status value
 */
INT4 CompAntennaVelocity(REAL4VectorAligned *output, const SkyPosition skypos, const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const LALDetector det, EphemerisData *edat)
{

   XLAL_CHECK( output != NULL && edat != NULL, XLAL_EINVAL );

   INT4 numffts = (INT4)floor(Tobs/(Tsft-SFToverlap)-1);    //Number of FFTs
   LALStatus XLAL_INIT_DECL(status);

   REAL8 detvel[3];
   for (INT4 ii=0; ii<numffts; ii++) {

      LIGOTimeGPS gpstime = LIGOTIMEGPSZERO;
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tsft-SFToverlap) + 0.5*Tsft);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tsft-SFToverlap)+0.5*Tsft - floor(t0+ii*(Tsft-SFToverlap)+0.5*Tsft))*1e9);

      LALDetectorVel(&status, detvel, &gpstime, det, edat);
      XLAL_CHECK( status.statusCode == 0, XLAL_EFUNC );

      output->data[ii] = (REAL4)(detvel[0]*cos(skypos.longitude)*cos(skypos.latitude) + detvel[1]*sin(skypos.longitude)*cos(skypos.latitude) + detvel[2]*sin(skypos.latitude));

   } /* for ii < numffts */

   return XLAL_SUCCESS;

} /* CompAntennaVelocity() */


/**
 * Compute the maximum change in antenna velocity
 * \param [in]  t0         Start time (GPS seconds)
 * \param [in]  Tsft       Coherence length of the SFTs (in seconds)
 * \param [in]  SFToverlap Overlap of the SFTs (in seconds)
 * \param [in]  Tobs       Observation time (in seconds)
 * \param [in]  det        A LALDetector struct
 * \param [in]  edat       Pointer to EphemerisData
 * \return Maximum change in antenna velocity
 */
REAL4 CompDetectorDeltaVmax(const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const LALDetector det, EphemerisData *edat)
{

   XLAL_CHECK( edat != NULL, XLAL_EINVAL );

   INT4 numffts = (INT4)floor(Tobs/(Tsft-SFToverlap)-1);    //Number of FFTs
   LALStatus XLAL_INIT_DECL(status);

   REAL8 detvel[3], detvel0[3], dv[3];
   REAL4 deltaVmax = 0.0;
   for (INT4 ii=0; ii<numffts; ii++) {
      LIGOTimeGPS gpstime = LIGOTIMEGPSZERO;
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tsft-SFToverlap) + 0.5*Tsft);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tsft-SFToverlap)+0.5*Tsft - floor(t0+ii*(Tsft-SFToverlap)+0.5*Tsft))*1e9);

      if (ii==0) {
         LALDetectorVel(&status, detvel0, &gpstime, det, edat);
         XLAL_CHECK_REAL4( status.statusCode == 0, XLAL_EFUNC );
      } else {
         LALDetectorVel(&status, detvel, &gpstime, det, edat);
         XLAL_CHECK_REAL4( status.statusCode == 0, XLAL_EFUNC );

         dv[0] = detvel[0] - detvel0[0];
         dv[1] = detvel[1] - detvel0[1];
         dv[2] = detvel[2] - detvel0[2];
         REAL4 deltaV = (REAL4)sqrt(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
         if (deltaV > deltaVmax) deltaVmax = deltaV;
      } /* if ii==0 else ... */
   } /* for ii < numffts */

   return deltaVmax;

} /* CompDetectorDeltaVmax() */


/**
 * Compute the largest magnitude of antenna velocity
 * \param [in]  t0         Start time (GPS seconds)
 * \param [in]  Tsft       Coherence length of the SFTs (in seconds)
 * \param [in]  SFToverlap Overlap of the SFTs (in seconds)
 * \param [in]  Tobs       Observation time (in seconds)
 * \param [in]  det        A LALDetector struct
 * \param [in]  edat       Pointer to EphemerisData
 * \return Maximum magnitude of antenna velocity
 */
REAL4 CompDetectorVmax(const REAL8 t0, const REAL8 Tsft, const REAL8 SFToverlap, const REAL8 Tobs, const LALDetector det, EphemerisData *edat)
{

   XLAL_CHECK( edat != NULL, XLAL_EINVAL );

   INT4 numffts = (INT4)floor(Tobs/(Tsft-SFToverlap)-1);    //Number of FFTs
   LALStatus XLAL_INIT_DECL(status);

   LIGOTimeGPS gpstime = LIGOTIMEGPSZERO;
   gpstime.gpsSeconds = (INT4)floor(t0 + 0.5*Tsft);
   gpstime.gpsNanoSeconds = (INT4)floor((t0+0.5*Tsft - floor(t0+0.5*Tsft))*1e9);

   REAL8 detvel[3];
   LALDetectorVel(&status, detvel, &gpstime, det, edat);
   XLAL_CHECK_REAL4( status.statusCode == 0, XLAL_EFUNC );
   REAL4 Vmax = (REAL4)sqrt(detvel[0]*detvel[0] + detvel[1]*detvel[1] + detvel[2]*detvel[2]);

   for (INT4 ii=1; ii<numffts; ii++) {
      gpstime.gpsSeconds = (INT4)floor(t0 + ii*(Tsft-SFToverlap) + 0.5*Tsft);
      gpstime.gpsNanoSeconds = (INT4)floor((t0+ii*(Tsft-SFToverlap)+0.5*Tsft - floor(t0+ii*(Tsft-SFToverlap)+0.5*Tsft))*1e9);

      LALDetectorVel(&status, detvel, &gpstime, det, edat);
      XLAL_CHECK_REAL4( status.statusCode == 0, XLAL_EFUNC );
      REAL4 V = (REAL4)sqrt(detvel[0]*detvel[0] + detvel[1]*detvel[1] + detvel[2]*detvel[2]);
      if (V > Vmax) Vmax = V;
   } /* for ii < numffts */

   return Vmax;

} /* CompDetectorVmax() */
REAL4 CompDetectorVmax2(const LIGOTimeGPSVector *timestamps, const LALDetector det, EphemerisData *edat)
{
   XLAL_CHECK( timestamps!=NULL && edat != NULL, XLAL_EINVAL );

   LALStatus XLAL_INIT_DECL(status);

   LIGOTimeGPS gpstime = timestamps->data[0];

   REAL8 detvel[3];
   LALDetectorVel(&status, detvel, &gpstime, det, edat);
   XLAL_CHECK_REAL4( status.statusCode == 0, XLAL_EFUNC );
   REAL4 Vmax = (REAL4)sqrt(detvel[0]*detvel[0] + detvel[1]*detvel[1] + detvel[2]*detvel[2]);

   for (UINT4 ii=1; ii<timestamps->length; ii++) {
      gpstime = timestamps->data[ii];

      LALDetectorVel(&status, detvel, &gpstime, det, edat);
      XLAL_CHECK_REAL4( status.statusCode == 0, XLAL_EFUNC );
      REAL4 V = (REAL4)sqrt(detvel[0]*detvel[0] + detvel[1]*detvel[1] + detvel[2]*detvel[2]);
      if (V > Vmax) Vmax = V;
   } /* for ii < numffts */

   return Vmax;
}
