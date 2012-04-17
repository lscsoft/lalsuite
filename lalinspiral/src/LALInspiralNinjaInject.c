/*
*  Copyright (C) 2007 Patrick Brady
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
 * File Name: LALSimNinjaInject.c
 *
 * Author: Pekowsky, L.,  Harry, I.
 *
 *
 *-----------------------------------------------------------------------
*/


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>

#include <lal/LALConstants.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrameStream.h>
#include <lal/LALSimulation.h>
#include <lal/NRWaveInject.h>
#include <lal/LALInspiral.h>


void AddNumRelStrainModesREAL8(
        LALStatus *status,
        REAL8TimeSeries **seriesPlus,
        REAL8TimeSeries **seriesCross,
        SimInspiralTable *thisinj
);

REAL8TimeSeries *XLALNRInjectionStrain(
        const char *ifo,
        SimInspiralTable *inj
);

void XLALSimInjectNinjaSignals(
        REAL4TimeSeries* chan,
        const char *ifo,
        REAL8 dynRange,
        SimInspiralTable* events
);


void AddNumRelStrainModesREAL8(LALStatus      *status,  /**< pointer to LALStatus structure */
                            REAL8TimeSeries   **seriesPlus, /**< [out]  h+, hx data    */
                            REAL8TimeSeries   **seriesCross, /**< [out]  h+, hx data    */
                            SimInspiralTable  *thisinj  /**< [in]   injection data */)
{
  INT4 modeL, modeM, modeLlo, modeLhi;
  INT4 len, lenPlus, lenCross, k;
  CHAR *channel_name_plus;
  CHAR *channel_name_cross;
  FrStream  *frStream = NULL;
  FrCache frCache;
  LIGOTimeGPS epoch;
  REAL8TimeSeries  *modePlus=NULL;
  REAL8TimeSeries  *modeCross=NULL;
  REAL8 massMpc, timeStep;
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  modeLlo = thisinj->numrel_mode_min;
  modeLhi = thisinj->numrel_mode_max;

  /* create a frame cache and open the frame stream */
  frCache.numFrameFiles     = 1;
  frCache.frameFiles        = LALCalloc(1, sizeof(frCache.frameFiles[0]));
  frCache.frameFiles[0].url = thisinj->numrel_data;
  frStream                  = XLALFrCacheOpen( &frCache );

  /* the total mass of the binary in Mpc */
  massMpc = (thisinj->mass1 + thisinj->mass2) * LAL_MRSUN_SI / ( LAL_PC_SI * 1.0e6);

  /* Time step in dimensionful units */
  timeStep = (thisinj->mass1 + thisinj->mass2) * LAL_MTSUN_SI;

  /* start time of waveform -- set it to something */
  epoch.gpsSeconds     = thisinj->geocent_end_time.gpsSeconds;
  epoch.gpsNanoSeconds = thisinj->geocent_end_time.gpsNanoSeconds;

  /* loop over l values */
  for ( modeL = modeLlo; modeL <= modeLhi; modeL++ ) {

    /* loop over m values */
    for ( modeM = -modeL; modeM <= modeL; modeM++ ) {
      /* read numrel waveform */
      /* first the plus polarization */
      channel_name_plus = XLALGetNinjaChannelName("plus", modeL, modeM);
      /*get number of data points */
      lenPlus = XLALFrGetVectorLength ( channel_name_plus, frStream );

      /* now the cross polarization */
      channel_name_cross = XLALGetNinjaChannelName("cross", modeL, modeM);
      /*get number of data points */
      lenCross = XLALFrGetVectorLength ( channel_name_cross, frStream );

      /* skip on to next mode if mode doesn't exist */
      if ( (lenPlus <= 0) || (lenCross <= 0) || (lenPlus != lenCross) ) {
        XLALClearErrno();
        LALFree(channel_name_plus);
        LALFree(channel_name_cross);
        continue;
      }

      /* note: lenPlus and lenCross must be equal if we got this far*/
      len = lenPlus;

      /* allocate and read the plus/cross time series */
      modePlus = XLALCreateREAL8TimeSeries ( channel_name_plus, &epoch, 0, 0, &lalDimensionlessUnit, len);
      memset(modePlus->data->data, 0, modePlus->data->length*sizeof(REAL8));
      XLALFrGetREAL8TimeSeries ( modePlus, frStream );
      XLALFrRewind( frStream );
      LALFree(channel_name_plus);

      modeCross = XLALCreateREAL8TimeSeries ( channel_name_cross, &epoch, 0, 0, &lalDimensionlessUnit, len);
      memset(modeCross->data->data, 0, modeCross->data->length*sizeof(REAL8));
      XLALFrGetREAL8TimeSeries ( modeCross, frStream );
      XLALFrRewind( frStream );
      LALFree(channel_name_cross);

      /* scale and add */
      if (*seriesPlus == NULL) {
          *seriesPlus = XLALCreateREAL8TimeSeries ( "hplus", &epoch, 0, 0, &lalDimensionlessUnit, len);
          memset((*seriesPlus)->data->data, 0, (*seriesPlus)->data->length*sizeof(REAL8));
          (*seriesPlus)->deltaT = modePlus->deltaT;
      }

      if (*seriesCross == NULL) {
          *seriesCross = XLALCreateREAL8TimeSeries ( "hcross", &epoch, 0, 0, &lalDimensionlessUnit, len);
          memset((*seriesCross)->data->data, 0, (*seriesCross)->data->length*sizeof(REAL8));
          (*seriesCross)->deltaT = modeCross->deltaT;
      }

      XLALOrientNRWaveTimeSeriesREAL8( modePlus, modeCross, modeL, modeM, thisinj->inclination, thisinj->coa_phase );

      for (k = 0; k < len; k++) {
        (*seriesPlus)->data->data[k]  += massMpc * modePlus->data->data[k];
        (*seriesCross)->data->data[k] += massMpc * modeCross->data->data[k];
      }

      /* we are done with seriesPlus and Cross for this iteration */
      XLALDestroyREAL8TimeSeries (modePlus);
      XLALDestroyREAL8TimeSeries (modeCross);
    } /* end loop over modeM values */
  } /* end loop over modeL values */
  (*seriesPlus)->deltaT  *= timeStep;
  (*seriesCross)->deltaT *= timeStep;
  XLALFrClose( frStream );
  LALFree(frCache.frameFiles);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

REAL8TimeSeries *
XLALNRInjectionStrain(const char *ifo, SimInspiralTable *inj)
{
  static LALStatus status;

  REAL8TimeSeries *plus      = NULL;
  REAL8TimeSeries *plus_int  = NULL;
  REAL8TimeSeries *cross     = NULL;
  REAL8TimeSeries *cross_int = NULL;
  REAL8TimeSeries *strain    = NULL;

  INT4 sampleRate                = 16384;
  InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;
  LALDetector det;

  /* Add the modes together */
  AddNumRelStrainModesREAL8(&status, &plus, &cross, inj);
  /* Place at distance */
  for (uint j = 0; j < plus->data->length; j++)
  {
    plus->data->data[j]  /= inj->distance;
    cross->data->data[j] /= inj->distance;
  }

  plus->sampleUnits  = lalADCCountUnit;
  cross->sampleUnits = lalADCCountUnit;

  /* Interpolate to desired sample rate */
  plus_int  = XLALInterpolateNRWaveREAL8(plus, sampleRate);
  cross_int = XLALInterpolateNRWaveREAL8(cross, sampleRate);

  /* look up detector */
  memset( &det, 0, sizeof(LALDetector) );
  ifoNumber = XLALIFONumber( ifo );
  XLALReturnDetector( &det, ifoNumber );

  /* Use Jolien's method to place on the sky */
  strain = XLALSimDetectorStrainREAL8TimeSeries(plus_int, cross_int,
           inj->longitude, inj->latitude, inj->polarization, &det);

  /* We want the end time to be the time of largest amplitude */
  REAL8 offset = 0;
  XLALFindNRCoalescencePlusCrossREAL8(&offset, plus_int, cross_int);
  XLALGPSAdd( &(strain->epoch), -offset);

  XLALDestroyREAL8TimeSeries (plus);
  XLALDestroyREAL8TimeSeries (plus_int);
  XLALDestroyREAL8TimeSeries (cross);
  XLALDestroyREAL8TimeSeries (cross_int);

  return strain;
}

void XLALSimInjectNinjaSignals(
  REAL4TimeSeries* chan,
  const char *ifo,
  REAL8 dynRange,
  SimInspiralTable* events
)
{
  /* New REAL8, NINJA-2 code */
  UINT4 j;
  SimInspiralTable *thisInj = NULL;
  
  REAL8TimeSeries *tempStrain = NULL;
  REAL8TimeSeries *tempChan   = NULL;

  /* Make a REAL8 version of the channel data    */
  /* so we can call Jolien's new inject function */
  tempChan = XLALCreateREAL8TimeSeries(
             chan->name,
             &(chan->epoch),
             chan->f0,
             chan->deltaT,
             &(chan->sampleUnits),
             chan->data->length);
  for ( j = 0 ; j < tempChan->data->length ; ++j )
  {
    tempChan->data->data[j] = (REAL8) ( chan->data->data[j] );
  }

  /* loop over injections */
  for ( thisInj = events; thisInj; thisInj = thisInj->next )
  {
    tempStrain = XLALNRInjectionStrain(ifo, thisInj);
    for ( j = 0 ; j < tempStrain->data->length ; ++j )
    {
      tempStrain->data->data[j] *= dynRange;
    }

    XLALSimAddInjectionREAL8TimeSeries( tempChan, tempStrain, NULL);
    XLALDestroyREAL8TimeSeries(tempStrain);
  } /* loop over injections */

  /* Back to REAL4 */
  for ( j = 0 ; j < tempChan->data->length ; ++j )
  {
    chan->data->data[j] = (REAL4) ( tempChan->data->data[j] );
  }

  XLALDestroyREAL8TimeSeries(tempChan);
}

