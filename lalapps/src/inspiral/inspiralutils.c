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
 * File Name: inspiralutils.c
 *
 * Author: Brown, D. A., Krishnan, B, Vitale S.
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

#include <lalapps.h>
#include <series.h>
#include <processtable.h>
#include <lalappsfrutils.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/PrintFTSeries.h>
#include <lal/LALFrStream.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Calibration.h>
#include <lal/FrameCalibration.h>
#include <lal/Window.h>
#include <lal/TimeFreqFFT.h>
#include <lal/IIRFilter.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>
#include <lal/LALFrameL.h>

#include <lal/LALSimulation.h>
#include <lal/LALSimNoise.h>

#include "inspiral.h"

/*
#include <lal/AVFactories.h>
#include <lal/FrequencySeries.h>
#include <lal/LALNoiseModels.h>
#include <lal/ConfigFile.h>
#include <lal/Units.h>
*/

REAL4 compute_candle_distance(REAL4 candleM1, REAL4 candleM2,
    REAL4 thissnr, REAL8 chanDeltaT, INT4 nPoints,
    REAL8FrequencySeries *spec, UINT4 cut)
{
  UINT4 k;
  REAL8 sigmaSqSum = 0;
  REAL8 distance = 0;
  REAL8 negativeSevenOverThree = -7.0/3.0;
  REAL8 totalMass = candleM1 + candleM2;
  REAL8 mu = candleM1 * candleM2 / totalMass;
  REAL8 distNorm = 2.0 * LAL_MRSUN_SI / (1.0e6 * LAL_PC_SI );
  REAL8 a = sqrt( (5.0 * mu) / 96.0 ) *
    pow( totalMass / ( LAL_PI * LAL_PI ), 1.0/3.0 ) *
    pow( LAL_MTSUN_SI / chanDeltaT, -1.0/6.0 );
  REAL8 sigmaSq = 4.0 * ( chanDeltaT / (REAL8) nPoints ) *
    distNorm * distNorm * a * a;
  REAL8 f_max = 1.0 / (6.0 * sqrt(6.0) * LAL_PI * totalMass * LAL_MTSUN_SI);
  REAL8 f = 0;

  for ( k = cut, f = spec->deltaF * cut;
      k < spec->data->length && f < f_max;
      ++k, f = spec->deltaF * k )
  {
    sigmaSqSum +=
      pow( (REAL8) k / (REAL8) nPoints, negativeSevenOverThree )
      / spec->data->data[k];
  }

  sigmaSq *= sigmaSqSum;

  distance = sqrt( sigmaSq ) / thissnr;

  return distance;
}

REAL4 XLALCandleDistanceTD(
    Approximant approximant,
    REAL4 candleM1,
    REAL4 candleM2,
    REAL4 candlesnr,
    REAL8 chanDeltaT,
    INT4 nPoints,
    REAL8FrequencySeries *spec,
    UINT4 cut)
{

  LALStatus      status   = blank_status;

  InspiralTemplate  tmplt;
  REAL4Vector    *waveform = NULL;
  COMPLEX8Vector *waveFFT  = NULL;
  REAL4FFTPlan   *fwdPlan  = NULL;

  REAL8          sigmaSq;
  REAL8          distance;
  UINT4          i;

  memset( &tmplt, 0, sizeof(tmplt) );

  /* Create storage for TD and FD template */
  waveform = XLALCreateREAL4Vector( nPoints );
  waveFFT  = XLALCreateCOMPLEX8Vector( spec->data->length );
  fwdPlan  = XLALCreateForwardREAL4FFTPlan( nPoints, 0 );

  /* Populate the template parameters */
  tmplt.mass1 = candleM1;
  tmplt.mass2 = candleM2;
  tmplt.ieta  = 1;
  tmplt.approximant = approximant;
  tmplt.tSampling   = 1.0/chanDeltaT;
  tmplt.order       = LAL_PNORDER_PSEUDO_FOUR; /* Hardcode for EOBNR for now */
  tmplt.fLower      = spec->deltaF *cut;
  tmplt.distance    = 1.0e6 * LAL_PC_SI; /* Mpc */
  tmplt.massChoice  = m1Andm2;
  tmplt.fCutoff     = tmplt.tSampling / 2.0 - spec->deltaF;

  /* From this, calculate the other parameters */
  LAL_CALL( LALInspiralParameterCalc( &status, &tmplt ), &status );

  /* Generate the waveform */
  LAL_CALL( LALInspiralWave( &status, waveform, &tmplt ), &status );

  XLALREAL4ForwardFFT( waveFFT, waveform, fwdPlan );

  sigmaSq = 0.0;
  for ( i = cut; i < waveFFT->length; i++ )
  {
    sigmaSq += ( crealf(waveFFT->data[i]) * crealf(waveFFT->data[i])
            + cimagf(waveFFT->data[i]) * cimagf(waveFFT->data[i]) )/ spec->data->data[i];
  }

  sigmaSq *= 4.0 * chanDeltaT / (REAL8)nPoints;

  /* Now calculate the distance */
  distance = sqrt( sigmaSq ) / (REAL8)candlesnr;

  /* Clean up! */
  XLALDestroyREAL4Vector( waveform );
  XLALDestroyCOMPLEX8Vector( waveFFT );
  XLALDestroyREAL4FFTPlan( fwdPlan );

  return (REAL4)distance;
}

SummValueTable **add_summvalue_table(SummValueTable **newTable,
    LIGOTimeGPS gpsStartTime, LIGOTimeGPS gpsEndTime,
    const CHAR *programName, const CHAR *ifoName,
    const CHAR *summValueName, const CHAR *comment, REAL8 value
    )
{
  /* add value to summ_value table */
  *newTable = (SummValueTable *) LALCalloc( 1, sizeof(SummValueTable) );
  snprintf( (*newTable)->program, LIGOMETA_PROGRAM_MAX,
      "%s", programName );
  (*newTable)->version = 0;
  (*newTable)->start_time = gpsStartTime;
  (*newTable)->end_time = gpsEndTime;
  snprintf( (*newTable)->ifo, LIGOMETA_IFO_MAX, "%s", ifoName );
  snprintf( (*newTable)->name, LIGOMETA_SUMMVALUE_NAME_MAX,
      "%s", summValueName );
  snprintf( (*newTable)->comment, LIGOMETA_SUMMVALUE_COMM_MAX,
      "%s", comment );
  (*newTable)->value = value;

  return (newTable);
}



void AddNumRelStrainModes(  LALStatus              *status,     /**< pointer to LALStatus structure */
                            REAL4TimeVectorSeries  **outStrain, /**< [out]  h+, hx data    */
                            SimInspiralTable *thisinj           /**< [in]   injection data */)
{
  INT4 modeL, modeM, modeLlo, modeLhi;
  INT4 len, lenPlus, lenCross, k;
  CHAR *channel_name_plus;
  CHAR *channel_name_cross;
  LALFrStream  *frStream = NULL;
  LALCache frCache;
  LIGOTimeGPS epoch;
  REAL4TimeSeries  *seriesPlus=NULL;
  REAL4TimeSeries  *seriesCross=NULL;
  REAL8 massMpc;
  REAL4TimeVectorSeries *sumStrain=NULL;
  REAL4TimeVectorSeries *tempStrain=NULL;
  /*   NRWaveMetaData thisMetaData; */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  modeLlo = thisinj->numrel_mode_min;
  modeLhi = thisinj->numrel_mode_max;

  /* create a frame cache and open the frame stream */
  frCache.length = 1;
  frCache.list = LALCalloc(1, sizeof(frCache.list[0]));
  frCache.list[0].url = thisinj->numrel_data;
  frStream = XLALFrStreamCacheOpen( &frCache );

  /* the total mass of the binary in Mpc */
  massMpc = (thisinj->mass1 + thisinj->mass2) * LAL_MRSUN_SI / ( LAL_PC_SI * 1.0e6);

  /* start time of waveform -- set it to something */
  epoch.gpsSeconds     = 0;
  epoch.gpsNanoSeconds = 0;

  /* loop over l values */
  for ( modeL = modeLlo; modeL <= modeLhi; modeL++ ) {

    /* loop over m values */
    for ( modeM = -modeL; modeM <= modeL; modeM++ ) {
      /* read numrel waveform */
      /* first the plus polarization */
      channel_name_plus = XLALGetNinjaChannelName("plus", modeL, modeM);
      /*get number of data points */
      lenPlus = XLALFrStreamGetVectorLength ( channel_name_plus, frStream );

      /* now the cross polarization */
      channel_name_cross = XLALGetNinjaChannelName("cross", modeL, modeM);
      /*get number of data points */
      lenCross = XLALFrStreamGetVectorLength ( channel_name_cross, frStream );

      /* skip on to next mode if mode doesn't exist */
      if ( (lenPlus <= 0) || (lenCross <= 0) || (lenPlus != lenCross) ) {
        XLALClearErrno();
        LALFree(channel_name_plus);
        LALFree(channel_name_cross);
        continue;
      }

      /* note: lenPlus and lenCross must be equal if we got this far*/
      /* len = lenPlus; */
      len = lenPlus;

      /* allocate and read the plus/cross time series */
      seriesPlus = XLALCreateREAL4TimeSeries ( channel_name_plus, &epoch, 0, 0, &lalDimensionlessUnit, len);
      memset(seriesPlus->data->data, 0, seriesPlus->data->length*sizeof(REAL4));
      XLALFrStreamGetREAL4TimeSeries ( seriesPlus, frStream );
      XLALFrStreamRewind( frStream );
      LALFree(channel_name_plus);

      seriesCross = XLALCreateREAL4TimeSeries ( channel_name_cross, &epoch, 0, 0, &lalDimensionlessUnit, len);
      memset(seriesCross->data->data, 0, seriesCross->data->length*sizeof(REAL4));
      XLALFrStreamGetREAL4TimeSeries ( seriesCross, frStream );
      XLALFrStreamRewind( frStream );
      LALFree(channel_name_cross);

      /* allocate memory for tempStrain */
      tempStrain = LALCalloc(1, sizeof(*tempStrain));
      tempStrain->data = XLALCreateREAL4VectorSequence(2, len);
      tempStrain->deltaT = LAL_MTSUN_SI * (thisinj->mass1 + thisinj->mass2) * seriesPlus->deltaT ;
      tempStrain->f0 = seriesPlus->f0;
      tempStrain->sampleUnits = seriesPlus->sampleUnits;
      memset(tempStrain->data->data, 0, tempStrain->data->length * tempStrain->data->vectorLength*sizeof(REAL4));

      /* now copy the data and scale amplitude corresponding to distance of 1Mpc*/
      for (k = 0; k < len; k++) {
        tempStrain->data->data[k] = massMpc * seriesPlus->data->data[k];
        tempStrain->data->data[len + k] = massMpc * seriesCross->data->data[k];            }

        /* we are done with seriesPlus and Cross for this iteration */
        XLALDestroyREAL4TimeSeries (seriesPlus);
        XLALDestroyREAL4TimeSeries (seriesCross);
        seriesPlus = NULL;
        seriesCross = NULL;

        /* compute the h+ and hx for given inclination and coalescence phase*/
        XLALOrientNRWave( tempStrain, modeL, modeM, thisinj->inclination, thisinj->coa_phase );
        if (sumStrain == NULL) {

          sumStrain = LALCalloc(1, sizeof(*sumStrain));
          sumStrain->data =  XLALCreateREAL4VectorSequence(2, tempStrain->data->vectorLength);
          sumStrain->deltaT = tempStrain->deltaT;
          sumStrain->f0 = tempStrain->f0;
          sumStrain->sampleUnits = tempStrain->sampleUnits;

          memset(sumStrain->data->data,0.0,2*tempStrain->data->vectorLength*sizeof(REAL4));

          sumStrain = XLALSumStrain( sumStrain, tempStrain );

        } else {
          sumStrain = XLALSumStrain( sumStrain, tempStrain );
        }

        /* clear memory for strain */
        if (tempStrain->data != NULL) {
          XLALDestroyREAL4VectorSequence ( tempStrain->data );
          LALFree( tempStrain );
          tempStrain = NULL;
      }
    } /* end loop over modeM values */
  } /* end loop over modeL values */

  XLALFrStreamClose( frStream );
  LALFree(frCache.list);
  *outStrain = sumStrain;
  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/** Main function for injecting numetrical relativity waveforms.
    Takes as input a list of injections, and adds h(t) to a given
    timeseries for a specified ifo and a dynamic range factor.
    Updated/generalized version of InjectNumRelWaveforms that allows
    arbitrary LIGO and Virgo PSDs and integration starting frequencies
*/
void InjectNumRelWaveformsUsingPSDREAL8(LALStatus *status,         /**< pointer to LALStatus structure */
                            REAL8TimeSeries      *chan,         /**< [out] the output time series */
                            SimInspiralTable     *injections,   /**< [in] list of injections */
                            CHAR                 ifo[3],        /**< [in] 2 char code for interferometer */
                            REAL8                freqLowCutoff, /**< [in] Lower cutoff frequency */
                            REAL8                snrLow,        /**< [in] lower cutoff value of snr */
                            REAL8                snrHigh,       /**< TO BE DOCUMENTED */
                            REAL8FrequencySeries *ligoPSD,        /**< [in] PSD to use for LIGO SNRs.  If NULL, use initial PSD */
                            REAL8                ligoSnrLowFreq,  /**< [in] Frequency at which to start integration for LIGO SNRs */
                            REAL8FrequencySeries *virgoPSD,       /**< [in] PSD to use for Virgo SNRs.  If NULL, use initial PSD */
                            REAL8                virgoSnrLowFreq, /**< [in] Frequency at which to start integration for Virgo SNRs */
                            CHAR                 *fname)       /**< [in] higher cutoff value of snr */
{
  SimInspiralTable *thisInj = NULL;
  REAL8 startFreq, startFreqHz, massTotal;
  REAL8 thisSNR;
  SimInspiralTable *simTableOut=NULL;
  SimInspiralTable *thisInjOut=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( chan, status, INSPIRALH_ENULL, INSPIRALH_MSGENULL );
  ASSERT( ifo, status, INSPIRALH_ENULL, INSPIRALH_MSGENULL );


  /* loop over injections */
  for ( thisInj = injections; thisInj; thisInj = thisInj->next )
    {

      startFreq = start_freq_from_frame_url(thisInj->numrel_data);
      massTotal = (thisInj->mass1 + thisInj->mass2) * LAL_MTSUN_SI;
      startFreqHz = startFreq / ( LAL_TWOPI * massTotal);

      if (startFreqHz < freqLowCutoff)
        {
          REAL8TimeSeries *strain = NULL;
          strain  = XLALNRInjectionStrain(ifo, thisInj);

          if (ifo[0] == 'V')
            thisSNR = calculate_snr_from_strain_and_psd_real8( strain, virgoPSD, virgoSnrLowFreq, ifo );
          else
            thisSNR = calculate_snr_from_strain_and_psd_real8( strain, ligoPSD, ligoSnrLowFreq, ifo );

           /* set channel name */
           snprintf( chan->name, LALNameLength * sizeof( CHAR ),
                    "%s:STRAIN", ifo );

          printf("Injection at %d.%d in ifo %s has SNR %f\n",
                   thisInj->geocent_end_time.gpsSeconds,
                   thisInj->geocent_end_time.gpsNanoSeconds,
                   ifo,
                   thisSNR);

          if ((thisSNR < snrHigh) && (thisSNR > snrLow))
            {
              /* simTableOut will be null only the first time */
              if ( simTableOut == NULL) {
                simTableOut = (SimInspiralTable *)LALCalloc( 1, sizeof(SimInspiralTable) );
                memcpy(simTableOut, thisInj, sizeof(*thisInj));
                simTableOut->next = NULL;
                thisInjOut = simTableOut;
              }
              else {
                thisInjOut->next = (SimInspiralTable *)LALCalloc( 1, sizeof(SimInspiralTable) );
                memcpy(thisInjOut->next, thisInj, sizeof(*thisInj));
                thisInjOut->next->next = NULL;
                thisInjOut = thisInjOut->next;
              }

              XLALSimAddInjectionREAL8TimeSeries( chan, strain, NULL);
            }

          XLALDestroyREAL8TimeSeries (strain);
        }
      else
        {
           fprintf( stderr, "Skipping injection at %d because it turns on at %f Hz, "
                            "but the low frequency cutoff is %f\n",
                            thisInj->geocent_end_time.gpsSeconds, startFreqHz, freqLowCutoff);
        }
    } /* loop over injectionsj */


  /* write and free the output simInspiral table */
  if ( simTableOut ) {

    LIGOLwXMLStream xmlfp;
    MetadataTable dummyTable;
    dummyTable.simInspiralTable = simTableOut;

    /* write the xml table of actual injections */
    if (fname) {
      memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
      TRY( LALOpenLIGOLwXMLFile( status->statusPtr, &xmlfp, fname ), status );
      TRY( LALBeginLIGOLwXMLTable( status->statusPtr, &xmlfp, sim_inspiral_table ),
                status );

      TRY( LALWriteLIGOLwXMLTable( status->statusPtr, &xmlfp, dummyTable,
                                        sim_inspiral_table ), status );

      TRY( LALEndLIGOLwXMLTable ( status->statusPtr, &xmlfp ), status );
      TRY( LALCloseLIGOLwXMLFile ( status->statusPtr, &xmlfp ), status );
    }
  }

  while (simTableOut) {
    thisInjOut = simTableOut;
    simTableOut = simTableOut->next;
    LALFree(thisInjOut);
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


/** Main function for injecting numetrical relativity waveforms.
    Takes as input a list of injections, and adds h(t) to a given
    timeseries for a specified ifo and a dynamic range factor.
*/
void InjectNumRelWaveforms (LALStatus           *status,       /**< pointer to LALStatus structure */
                            REAL4TimeSeries     *chan,         /**< [out] the output time series */
                            SimInspiralTable    *injections,   /**< [in] list of injections */
                            CHAR                ifo[3],        /**< [in] 2 char code for interferometer */
                            REAL8               dynRange,      /**< [in] dynamic range factor for scaling time series */
                            REAL8               freqLowCutoff, /**< [in] Lower cutoff frequency */
                            REAL8               snrLow,        /**< [in] lower cutoff value of snr */
                            REAL8               snrHigh,       /**< TO BE DOCUMENTED */
                            CHAR                *fname)        /**< [in] higher cutoff value of snr */
{

  REAL4TimeVectorSeries *tempStrain=NULL;
  SimInspiralTable    *thisInj = NULL;
  REAL8 startFreq, startFreqHz, massTotal;
  REAL8 thisSNR;
  SimInspiralTable *simTableOut=NULL;
  SimInspiralTable *thisInjOut=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( chan, status, INSPIRALH_ENULL, INSPIRALH_MSGENULL );
  ASSERT( ifo, status, INSPIRALH_ENULL, INSPIRALH_MSGENULL );


  /* loop over injections */
  for ( thisInj = injections; thisInj; thisInj = thisInj->next )
    {

      startFreq = start_freq_from_frame_url(thisInj->numrel_data);
      massTotal = (thisInj->mass1 + thisInj->mass2) * LAL_MTSUN_SI;
      startFreqHz = startFreq / ( LAL_TWOPI * massTotal);

       if (startFreqHz < freqLowCutoff)
       {
          TRY( AddNumRelStrainModes( status->statusPtr, &tempStrain, thisInj),
               status);
          thisSNR = calculate_ligo_snr_from_strain( tempStrain, thisInj, ifo);

          /* set channel name */
          snprintf( chan->name, LALNameLength * sizeof( CHAR ),
                    "%s:STRAIN", ifo );

          if ((thisSNR < snrHigh) && (thisSNR > snrLow))
            {
              /* simTableOut will be null only the first time */
              if ( simTableOut == NULL) {
                simTableOut = (SimInspiralTable *)LALCalloc( 1, sizeof(SimInspiralTable) );
                memcpy(simTableOut, thisInj, sizeof(*thisInj));
                simTableOut->next = NULL;
                thisInjOut = simTableOut;
              }
              else {
                thisInjOut->next = (SimInspiralTable *)LALCalloc( 1, sizeof(SimInspiralTable) );
                memcpy(thisInjOut->next, thisInj, sizeof(*thisInj));
                thisInjOut->next->next = NULL;
                thisInjOut = thisInjOut->next;
              }

              TRY( LALInjectStrainGW( status->statusPtr, chan, tempStrain, thisInj,
                                      ifo, dynRange), status);
            }

          XLALDestroyREAL4VectorSequence ( tempStrain->data);
          tempStrain->data = NULL;
          LALFree(tempStrain);
          tempStrain = NULL;
        }

    } /* loop over injectionsj */


  /* write and free the output simInspiral table */
  if ( simTableOut ) {

    LIGOLwXMLStream xmlfp;
    MetadataTable dummyTable;
      dummyTable.simInspiralTable = simTableOut;

    /* write the xml table of actual injections */
    memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
    TRY( LALOpenLIGOLwXMLFile( status->statusPtr, &xmlfp, fname ), status );
    TRY( LALBeginLIGOLwXMLTable( status->statusPtr, &xmlfp, sim_inspiral_table ),
              status );

    TRY( LALWriteLIGOLwXMLTable( status->statusPtr, &xmlfp, dummyTable,
                                      sim_inspiral_table ), status );

    TRY( LALEndLIGOLwXMLTable ( status->statusPtr, &xmlfp ), status );
    TRY( LALCloseLIGOLwXMLFile ( status->statusPtr, &xmlfp ), status );

  }

  while (simTableOut) {
    thisInjOut = simTableOut;
    simTableOut = simTableOut->next;
    LALFree(thisInjOut);
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/** Main function for injecting numetrical relativity waveforms.
    Takes as input a list of injections, and adds h(t) to a given
    timeseries for a specified ifo and a dynamic range factor.
*/
void InjectNumRelWaveformsREAL8 (LALStatus      *status,       /**< pointer to LALStatus structure */
                            REAL8TimeSeries     *chan,         /**< [out] the output time series */
                            SimInspiralTable    *injections,   /**< [in] list of injections */
                            CHAR                ifo[3],        /**< [in] 2 char code for interferometer */
                            REAL8               freqLowCutoff, /**< [in] Lower cutoff frequency */
                            REAL8               snrLow,        /**< [in] lower cutoff value of snr */
                            REAL8               snrHigh,       /**< TO BE DOCUMENTED */
                            CHAR                *fname)       /**< [in] higher cutoff value of snr */
{
  SimInspiralTable *thisInj = NULL;
  REAL8 startFreq, startFreqHz, massTotal;
  REAL8 thisSNR;
  SimInspiralTable *simTableOut=NULL;
  SimInspiralTable *thisInjOut=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT( chan, status, INSPIRALH_ENULL, INSPIRALH_MSGENULL );
  ASSERT( ifo, status, INSPIRALH_ENULL, INSPIRALH_MSGENULL );


  /* loop over injections */
  for ( thisInj = injections; thisInj; thisInj = thisInj->next )
    {

      startFreq = start_freq_from_frame_url(thisInj->numrel_data);
      massTotal = (thisInj->mass1 + thisInj->mass2) * LAL_MTSUN_SI;
      startFreqHz = startFreq / ( LAL_TWOPI * massTotal);

      if (startFreqHz < freqLowCutoff)
        {
          REAL8TimeSeries *strain = NULL;
          strain  = XLALNRInjectionStrain(ifo, thisInj);
          thisSNR = calculate_ligo_snr_from_strain_real8(strain, ifo);

           /* set channel name */
           snprintf( chan->name, LALNameLength * sizeof( CHAR ),
                    "%s:STRAIN", ifo );

          if ((thisSNR < snrHigh) && (thisSNR > snrLow))
            {
              /* simTableOut will be null only the first time */
              if ( simTableOut == NULL) {
                simTableOut = (SimInspiralTable *)LALCalloc( 1, sizeof(SimInspiralTable) );
                memcpy(simTableOut, thisInj, sizeof(*thisInj));
                simTableOut->next = NULL;
                thisInjOut = simTableOut;
              }
              else {
                thisInjOut->next = (SimInspiralTable *)LALCalloc( 1, sizeof(SimInspiralTable) );
                memcpy(thisInjOut->next, thisInj, sizeof(*thisInj));
                thisInjOut->next->next = NULL;
                thisInjOut = thisInjOut->next;
              }

              XLALSimAddInjectionREAL8TimeSeries( chan, strain, NULL);
            }

          XLALDestroyREAL8TimeSeries (strain);
        }
    } /* loop over injectionsj */


  /* write and free the output simInspiral table */
  if ( simTableOut ) {

    LIGOLwXMLStream xmlfp;
    MetadataTable dummyTable;
    dummyTable.simInspiralTable = simTableOut;

    /* write the xml table of actual injections */
    memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
    TRY( LALOpenLIGOLwXMLFile( status->statusPtr, &xmlfp, fname ), status );
    TRY( LALBeginLIGOLwXMLTable( status->statusPtr, &xmlfp, sim_inspiral_table ),
              status );

    TRY( LALWriteLIGOLwXMLTable( status->statusPtr, &xmlfp, dummyTable,
                                      sim_inspiral_table ), status );

    TRY( LALEndLIGOLwXMLTable ( status->statusPtr, &xmlfp ), status );
    TRY( LALCloseLIGOLwXMLFile ( status->statusPtr, &xmlfp ), status );

  }

  while (simTableOut) {
    thisInjOut = simTableOut;
    simTableOut = simTableOut->next;
    LALFree(thisInjOut);
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

REAL8 start_freq_from_frame_url(CHAR  *url)
{

  FrameH *frame=NULL;
  FrFile *frFile=NULL;
  FrHistory *frHist=NULL;
  FrHistory *thisHist;
  CHAR *comment=NULL;
  CHAR *token=NULL;
  REAL8 ret=0;
  CHAR *path;

  /* convert url to path by skipping protocol part of protocol:path */
  path = strchr(url, ':');
  if (path == NULL)
    path = url;
  else
    path++; /* skip the ':' -- now on the path */

  frFile =  FrFileINew( path );
  frame = FrameRead (frFile);
  frHist = frame->history;
  thisHist = frHist;
  while (thisHist) {

    /* get history comment string and parse it */
    comment = LALCalloc(1, (strlen(thisHist->comment)+1)*sizeof(CHAR));
    strcpy(comment, thisHist->comment);

    token = strtok(comment,":");

    if (strstr(token,"freqStart22") || strstr(token,"freq_start_22")) {
      token = strtok(NULL,":");
      ret = atof(token);
    }

    LALFree(comment);
    comment = NULL;

    thisHist = thisHist->next;
  }

  FrFileIEnd( frFile );
  return ret;
}


REAL8 calculate_ligo_snr_from_strain_real8(  REAL8TimeSeries *strain,
                                       const CHAR            ifo[3])
{

  REAL8 ret = -1, snrSq, freq, psdValue;
  REAL8 deltaF;
  REAL8FFTPlan *pfwd;
  COMPLEX16FrequencySeries *fftData;
  UINT4 k;

  /* create the time series */
  deltaF  = strain->deltaT * strain->data->length;
  fftData = XLALCreateCOMPLEX16FrequencySeries( strain->name,  &(strain->epoch),
                                               0, deltaF, &lalDimensionlessUnit,
                                               strain->data->length/2 + 1 );

  /* perform the fft */
  pfwd = XLALCreateForwardREAL8FFTPlan( strain->data->length, 0 );
  XLALREAL8TimeFreqFFT( fftData, strain, pfwd );

  /* compute the SNR for initial LIGO at design */
  for ( snrSq = 0, k = 0; k < fftData->data->length; k++ )
    {
      freq = fftData->deltaF * k;

      if ( ifo[0] == 'V' )
        {
          if (freq < 35)
            continue;

          LALVIRGOPsd( NULL, &psdValue, freq );
          psdValue /= 9e-46;
        }
      else
        {
          if (freq < 40)
            continue;

          LALLIGOIPsd( NULL, &psdValue, freq );
        }

      fftData->data->data[k] /= 3e-23;
      snrSq += creal(fftData->data->data[k]) * creal(fftData->data->data[k]) / psdValue;
      snrSq += cimag(fftData->data->data[k]) * cimag(fftData->data->data[k]) / psdValue;
    }
  snrSq *= 4*fftData->deltaF;

  XLALDestroyREAL8FFTPlan( pfwd );
  XLALDestroyCOMPLEX16FrequencySeries( fftData );

  ret = sqrt(snrSq);
  return ret;
}


REAL8 calculate_snr_from_strain_and_psd_real8(  REAL8TimeSeries *strain,
                                       REAL8FrequencySeries  *psd,
                                       REAL8                 startFreq,
                                       const CHAR            ifo[3])
{

  REAL8 ret = -1, snrSq, freq, psdValue;
  REAL8 deltaF;
  REAL8FFTPlan *pfwd;
  COMPLEX16FrequencySeries *fftData;
  UINT4 k;

  /* create the time series */
  deltaF  = strain->deltaT * strain->data->length;
  fftData = XLALCreateCOMPLEX16FrequencySeries( strain->name,  &(strain->epoch),
                                               0, deltaF, &lalDimensionlessUnit,
                                               strain->data->length/2 + 1 );

  /* perform the fft */
  pfwd = XLALCreateForwardREAL8FFTPlan( strain->data->length, 0 );
  XLALREAL8TimeFreqFFT( fftData, strain, pfwd );

  /* The PSD, if provided, comes in as it was in the original file  */
  /* since we don't know deltaF until we get here.  Interpolate now */
  if ( psd )
  {
    psd = XLALInterpolatePSD(psd, 1.0 / deltaF);
  }

  /* compute the SNR for initial LIGO at design */
  for ( snrSq = 0, k = 0; k < fftData->data->length; k++ )
  {
    freq = fftData->deltaF * k;

    if ( psd )
    {
      if ( freq < startFreq || k > psd->data->length )
        continue;
      psdValue  = psd->data->data[k];
      psdValue /= 9e-46;
    }
    else if ( ifo[0] == 'V' )
    {
      if (freq < 35)
        continue;
      LALVIRGOPsd( NULL, &psdValue, freq );
      psdValue /= 9e-46;
    }
    else
    {
      if (freq < 40)
        continue;
      LALLIGOIPsd( NULL, &psdValue, freq );
    }

    fftData->data->data[k] /= 3e-23;

    snrSq += creal(fftData->data->data[k]) * creal(fftData->data->data[k]) / psdValue;
    snrSq += cimag(fftData->data->data[k]) * cimag(fftData->data->data[k]) / psdValue;
  }

  snrSq *= 4*fftData->deltaF;

  XLALDestroyREAL8FFTPlan( pfwd );
  XLALDestroyCOMPLEX16FrequencySeries( fftData );

  if ( psd )
    XLALDestroyREAL8FrequencySeries( psd );

  ret = sqrt(snrSq);

  printf("Obtained snr=%f\n", ret);
  return ret;
}


REAL8 calculate_ligo_snr_from_strain(  REAL4TimeVectorSeries *strain,
                                       SimInspiralTable      *thisInj,
                                       const CHAR            ifo[3])
{

  REAL8 ret = -1, snrSq, freq, psdValue;
  REAL8 sampleRate = 4096, deltaF;
  REAL4TimeSeries *chan = NULL;
  REAL4FFTPlan *pfwd;
  COMPLEX8FrequencySeries *fftData;
  UINT4 k;

  /* create the time series */
  chan    = XLALCalculateNRStrain( strain, thisInj, ifo, sampleRate );
  deltaF  = chan->deltaT * strain->data->vectorLength;
  fftData = XLALCreateCOMPLEX8FrequencySeries( chan->name,  &(chan->epoch),
                                               0, deltaF, &lalDimensionlessUnit,
                                               chan->data->length/2 + 1 );

  /* perform the fft */
  pfwd = XLALCreateForwardREAL4FFTPlan( chan->data->length, 0 );
  XLALREAL4TimeFreqFFT( fftData, chan, pfwd );

  /* compute the SNR for initial LIGO at design */
  for ( snrSq = 0, k = 0; k < fftData->data->length; k++ )
    {
      freq = fftData->deltaF * k;

      if ( ifo[0] == 'V' )
        {
          if (freq < 35)
            continue;

          LALVIRGOPsd( NULL, &psdValue, freq );
          psdValue /= 9e-46;
        }
      else
        {
          if (freq < 40)
            continue;

          LALLIGOIPsd( NULL, &psdValue, freq );
        }

      fftData->data->data[k] /= 3e-23;
      snrSq += crealf(fftData->data->data[k]) * crealf(fftData->data->data[k]) / psdValue;
      snrSq += cimagf(fftData->data->data[k]) * cimagf(fftData->data->data[k]) / psdValue;
    }

  snrSq *= 4*fftData->deltaF;

  XLALDestroyREAL4FFTPlan( pfwd );
  XLALDestroyCOMPLEX8FrequencySeries( fftData );

  XLALDestroyREAL4Vector ( chan->data);
  LALFree(chan);

  ret = sqrt(snrSq);
  return ret;
}

int
XLALPsdFromFile(REAL8FrequencySeries **psd,  /**< [out] The PSD */
                const CHAR *filename)        /**< [in] name of the file to be read */
{
  REAL8FrequencySeries *ret;
  LALParsedDataFile *cfgdata=NULL;
  LIGOTimeGPS stubEpoch;
  UINT4 length, k, r;
  REAL8 freq, value;
  REAL8 step1=0, deltaF=0;
  int retval;

  /* XLALParseDataFile checks that filename is not null for us */
  retval = XLALParseDataFile(&cfgdata, filename);
  if ( retval != XLAL_SUCCESS ) {
    XLAL_ERROR ( retval );
  }

  /*number of data points */
  length = cfgdata->lines->nTokens;

  /* allocate memory */
  ret = XLALCreateREAL8FrequencySeries("PSD", &stubEpoch, 0, 0, &lalHertzUnit, length);
  if (ret == NULL) {
    XLALDestroyParsedDataFile( cfgdata );
    XLALPrintError ("%s: XLALPsdFromFile() failed.\n", __func__ );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* now get the data */
  for (k = 0; k < length; k++) {
    r = sscanf(cfgdata->lines->tokens[k], "%lf%lf", &freq, &value);

    /* Check the data file format */
    if ( r != 2 ) {
      XLALDestroyParsedDataFile( cfgdata );
      XLALPrintError ("%s: XLALPsdFromFile() failed on bad line in psd file.\n", __func__ );
      XLAL_ERROR ( XLAL_EFUNC );
    }

    if (deltaF == 0) {
      if (step1 == 0) {
        step1 = freq;
      } else {
        deltaF = (freq - step1);
        ret->deltaF = deltaF;
      }
    }

    ret->data->data[k] = value;
  }

  (*psd) = ret;

  XLALDestroyParsedDataFile( cfgdata );

  return XLAL_SUCCESS;
}


/** Function for interpolating PSD to a given sample rate
  */
REAL8FrequencySeries *
XLALInterpolatePSD( REAL8FrequencySeries *in,      /**< input strain time series */
                    REAL8                deltaFout /**< sample rate of time series */)
{
  REAL8FrequencySeries *ret=NULL;
  REAL8 deltaFin, r, y_1, y_2;
  UINT4 k, lo, numPoints;

  deltaFin = in->deltaF;

  /* length of output vector */
  numPoints = (UINT4) (in->data->length * deltaFin / deltaFout);

  /* allocate memory */
  ret = LALCalloc(1, sizeof(*ret));
  if (!ret)
  {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  ret->data = XLALCreateREAL8Vector( numPoints );
  if (! ret->data)
  {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  ret->deltaF = deltaFout;

  /* copy values from in which should be the same */
  ret->epoch       = in->epoch;
  ret->f0          = in->f0;
  ret->sampleUnits = in->sampleUnits;
  strcpy(ret->name, in->name);

  /* go over points of output vector and interpolate linearly
     using closest points of input */
  for (k = 0; k < numPoints; k++) {
    lo = (UINT4)( k*deltaFout / deltaFin);

    /* y_1 and y_2 are the input values at x1 and x2 */
    /* here we need to make sure that we don't exceed
       bounds of input vector */
    if ( lo < in->data->length - 1) {
      y_1 = in->data->data[lo];
      y_2 = in->data->data[lo+1];

      /* we want to calculate y_2*r + y_1*(1-r) where
         r = (x-x1)/(x2-x1) */
      r = k*deltaFout / deltaFin - lo;

      ret->data->data[k] = y_2 * r + y_1 * (1 - r);
    }
    else {
      ret->data->data[k] = 0.0;
    }
  }

  return ret;
}


void get_FakePsdFromString(REAL8FrequencySeries* PsdFreqSeries,char* FakePsdName, REAL8 StartFreq)
{
  /* Call XLALSimNoisePSD to fill the REAL8FrequencySeries PsdFreqSeries (must been already allocated by callers). FakePsdName contains the label of the fake PSD */
  if (!strcmp("LALAdVirgo",FakePsdName))
  {
    XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDAdvVirgo);
  }
  else if (!strcmp("LALVirgo",FakePsdName))
  {
    XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDVirgo);
  }
  else if(!strcmp("LALAdLIGO",FakePsdName))
  {
    XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDaLIGOZeroDetHighPower);
  }
  else if(!strcmp("LALLIGO",FakePsdName))
  {
    XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDiLIGOSRD);
  }
  else
  {
    fprintf(stderr,"Unknown fake PSD %s. Known types are LALLIGO, LALAdLIGO, LALAdVirgo, LALVirgo. Exiting...\n",FakePsdName);
    exit(1);
  }
}


REAL8 calculate_lalsim_snr(SimInspiralTable *inj, char *IFOname, REAL8FrequencySeries *psd, REAL8 start_freq)
{
  /* Calculate and return the single IFO SNR
   *
   * Required options:
   *
   * inj:     SimInspiralTable entry for which the SNR has to be calculated
   * IFOname: The canonical name (e.g. H1, L1, V1) name of the IFO for which the SNR must be calculated
   * PSD:     PSD curve to be used for the overlap integrap
   * start_freq: lower cutoff of the overlap integral
   *
   * */

  int ret=0;
  INT4 errnum=0;
  UINT4 j=0;
  /* Fill detector site info */
  LALDetector*  detector=NULL;
  detector=calloc(1,sizeof(LALDetector));
  if(!strcmp(IFOname,"H1"))
    memcpy(detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
  if(!strcmp(IFOname,"H2"))
    memcpy(detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
  if(!strcmp(IFOname,"LLO")||!strcmp(IFOname,"L1"))
    memcpy(detector,&lalCachedDetectors[LALDetectorIndexLLODIFF],sizeof(LALDetector));
  if(!strcmp(IFOname,"V1")||!strcmp(IFOname,"VIRGO"))
    memcpy(detector,&lalCachedDetectors[LALDetectorIndexVIRGODIFF],sizeof(LALDetector));

  Approximant approx=TaylorF2;
  approx=XLALGetApproximantFromString(inj->waveform);
  LALSimulationDomain modelDomain;

  switch(approx)
  {
    case GeneratePPN:
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case EOBNRv2HM:
    case SpinTaylor:
    case SpinTaylorT4:
    case SpinQuadTaylor:
    case SpinTaylorFrameless:
    case PhenSpinTaylorRD:
    case NumRel:
      modelDomain=LAL_SIM_DOMAIN_TIME;
      break;
    case TaylorF1:
    case TaylorF2:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case IMRPhenomA:
    case IMRPhenomB:
      modelDomain=LAL_SIM_DOMAIN_FREQUENCY;
      break;
    default:
      fprintf(stderr,"ERROR. Unknown approximant number %i. Unable to choose time or frequency domain model.",approx);
      exit(1);
      break;
  }

  REAL8 m1,m2, s1x,s1y,s1z,s2x,s2y,s2z,phi0,f_min,f_max,iota,polarization,

  /* No tidal PN terms until injtable is able to get them */
  lambda1=0.0,lambda2=0.0;

  LALSimInspiralWaveformFlags *waveFlags= XLALSimInspiralCreateWaveformFlags();

  /* Spin and tidal interactions at the highest level (default) until injtable stores them.
   * When spinO and tideO are added to injtable we can un-comment those lines and should be ok
   *
  int spinO = inj->spinO;
  int tideO = inj->tideO;
  XLALSimInspiralSetSpinOrder(waveFlags, *(LALSimInspiralSpinOrder*) spinO);
  XLALSimInspiralSetTidalOrder(waveFlags, *(LALSimInspiralTidalOrder*) tideO);
  */

  /* When nonGR terms stored in the table, we can add them here. (If they are only phase deformations, they won't change the SNR though) */
  LALSimInspiralTestGRParam *nonGRparams=NULL;
  /* Linked list of non-GR parameters. Pass in NULL (or None in python) for standard GR waveforms */

  LALPNOrder  order;   /* Phase order of the model   */
  INT4        amporder=0;
  order = XLALGetOrderFromString(inj->waveform);
  amporder = inj->amp_order;
  /* Read parameters */
  m1=inj->mass1*LAL_MSUN_SI;
  m2=inj->mass2*LAL_MSUN_SI;
  s1x=inj->spin1x;
  s1y=inj->spin1y;
  s1z=inj->spin1z;
  s2x=inj->spin2x;
  s2y=inj->spin2y;
  s2z=inj->spin2z;
  iota=inj->inclination;
  f_min=start_freq;
  phi0=inj->coa_phase;
  polarization=inj->polarization;
  REAL8 latitude=inj->latitude;
  REAL8 longitude=inj->longitude;

  LIGOTimeGPS epoch;
  memcpy(&epoch,&(inj->geocent_end_time),sizeof(LIGOTimeGPS));

  /* Hardcoded values of srate and segment length. If changed here they must also be changed in inspinj.c */
  REAL8 srate=4096.0;
  const CHAR *WF=inj->waveform;
  /* Increase srate for EOB WFs */
  if (strstr(WF,"EOB"))
    srate=8192.0;
  REAL8 segment=64.0;

  f_max=(srate/2.0-(1.0/segment));
  size_t seglen=(size_t) segment*srate;
  REAL8 deltaF=1.0/segment;
  REAL8 deltaT=1.0/srate;

  /* Frequency domain h+ and hx. They are going to be filled either by a FD WF or by the FFT of a TD WF*/
  COMPLEX16FrequencySeries *freqHplus;
  COMPLEX16FrequencySeries *freqHcross;
  freqHplus=  XLALCreateCOMPLEX16FrequencySeries("fhplus",
    &epoch,
    0.0,
    deltaF,
    &lalDimensionlessUnit,
    seglen/2+1
  );

  freqHcross=XLALCreateCOMPLEX16FrequencySeries("fhcross",
    &epoch,
    0.0,
    deltaF,
    &lalDimensionlessUnit,
    seglen/2+1
  );

  /* If the approximant is on the FD call XLALSimInspiralChooseFDWaveform */
  if (modelDomain == LAL_SIM_DOMAIN_FREQUENCY)
  {

    COMPLEX16FrequencySeries *hptilde=NULL;
    COMPLEX16FrequencySeries *hctilde=NULL;
    XLAL_TRY(ret=XLALSimInspiralChooseFDWaveform(&hptilde,&hctilde, phi0, deltaF, m1, m2,
      s1x, s1y, s1z, s2x, s2y, s2z, f_min, 0.0, LAL_PC_SI * 1.0e6,
      iota, lambda1, lambda2, waveFlags, nonGRparams,
      amporder, order, approx),errnum
    );

    if(!hptilde|| hptilde->data==NULL || hptilde->data->data==NULL ||!hctilde|| hctilde->data==NULL || hctilde->data->data==NULL)
    {
      XLALPrintError(" ERROR in XLALSimInspiralChooseFDWaveform(): error generating waveform. errnum=%d. Exiting...\n",errnum );
      exit(1);
    }

    COMPLEX16 *dataPtr = hptilde->data->data;
    for (j=0; j<(UINT4) freqHplus->data->length; ++j)
    {
      if(j < hptilde->data->length)
      {
        freqHplus->data->data[j] = dataPtr[j];
      }
      else
      {
        freqHplus->data->data[j]=0.0 + I*0.0;
      }
    }
    dataPtr = hctilde->data->data;
    for (j=0; j<(UINT4) freqHplus->data->length; ++j)
    {
      if(j < hctilde->data->length)
      {
        freqHcross->data->data[j] = dataPtr[j];
      }
      else
      {
        freqHcross->data->data[j]=0.0+0.0*I;
      }
    }
    /* Clean */
    if(hptilde) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    if(hctilde) XLALDestroyCOMPLEX16FrequencySeries(hctilde);

  }
  else
  {

    /* Otherwise use XLALSimInspiralChooseTDWaveform */
    REAL8FFTPlan *timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 0 );
    REAL8TimeSeries *hplus=NULL;
    REAL8TimeSeries *hcross=NULL;
    REAL8TimeSeries *timeHplus=NULL;
    REAL8TimeSeries *timeHcross=NULL;

    timeHcross=XLALCreateREAL8TimeSeries("timeModelhCross",
      &epoch,
      0.0,
      deltaT,
      &lalDimensionlessUnit,
      seglen
    );
    timeHplus=XLALCreateREAL8TimeSeries("timeModelhplus",
      &epoch,
      0.0,
      deltaT,
      &lalDimensionlessUnit,
      seglen
    );

    XLAL_TRY(ret=XLALSimInspiralChooseTDWaveform(&hplus, &hcross, phi0, deltaT,
        m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, f_min, 0., LAL_PC_SI*1.0e6,
        iota, lambda1, lambda2, waveFlags, nonGRparams, amporder, order, approx),
        errnum);

    if (ret == XLAL_FAILURE || hplus == NULL || hcross == NULL)
    {
      XLALPrintError(" ERROR in XLALSimInspiralChooseTDWaveform(): error generating waveform. errnum=%d. Exiting...\n",errnum );
      exit(1);
    }

    memset(timeHplus->data->data, 0, sizeof (REAL8)*timeHplus->data->length);
    memset(timeHcross->data->data, 0, sizeof (REAL8)*timeHcross->data->length);
    memcpy(timeHplus->data->data, hplus->data->data,hplus->data->length*sizeof(REAL8));
    memcpy(timeHcross->data->data, hcross->data->data ,hplus->data->length*sizeof(REAL8));

    for (j=0; j<(UINT4) freqHplus->data->length; ++j)
    {
      freqHplus->data->data[j]=0.0+I*0.0;
      freqHcross->data->data[j]=0.0+I*0.0;
    }

    /* FFT into freqHplus and freqHcross */
    XLALREAL8TimeFreqFFT(freqHplus,timeHplus,timeToFreqFFTPlan);
    XLALREAL8TimeFreqFFT(freqHcross,timeHcross,timeToFreqFFTPlan);

    /* Clean... */
    if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
    if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
    if ( timeHplus ) XLALDestroyREAL8TimeSeries(timeHplus);
    if ( timeHcross ) XLALDestroyREAL8TimeSeries(timeHcross);
    if (timeToFreqFFTPlan) LALFree(timeToFreqFFTPlan);

  }

  /* The WF has been generated and is in freqHplus/cross. Now project into the IFO frame */
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double HSquared;
  double GPSdouble=(REAL8) inj->geocent_end_time.gpsSeconds+ (REAL8) inj->geocent_end_time.gpsNanoSeconds*1.0e-9;
  double gmst;
  LIGOTimeGPS GPSlal;
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  /* Fill Fplus and Fcross*/
  XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4 (*)[3])detector->response,longitude, latitude, polarization, gmst);
  /* And take the distance into account */
  FplusScaled  = Fplus  / (inj->distance);
  FcrossScaled = Fcross / (inj->distance);

  REAL8 timedelay = XLALTimeDelayFromEarthCenter(detector->location,longitude, latitude, &GPSlal);
  REAL8 timeshift =  timedelay;
  REAL8 twopit    = LAL_TWOPI * timeshift;

  UINT4 lower = (UINT4)ceil(f_min / deltaF);
  UINT4 upper = (UINT4)floor(f_max / deltaF);
  REAL8 re = cos(twopit*deltaF*lower);
  REAL8 im = -sin(twopit*deltaF*lower);

  /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
  REAL8 dim = -sin(twopit*deltaF);
  REAL8 dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);
  REAL8 TwoDeltaToverN = 2.0 *deltaT / ((double) seglen);

  REAL8 plainTemplateReal,  plainTemplateImag,templateReal,templateImag;
  REAL8 newRe, newIm,temp;
  REAL8 this_snr=0.0;
  for (j=lower; j<=(UINT4) upper; ++j)
  {
    /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
    plainTemplateReal = FplusScaled * creal(freqHplus->data->data[j])
                        +  FcrossScaled *creal(freqHcross->data->data[j]);
    plainTemplateImag = FplusScaled * cimag(freqHplus->data->data[j])
                        +  FcrossScaled * cimag(freqHcross->data->data[j]);

    /* do time-shifting...             */
    /* (also un-do 1/deltaT scaling): */
    templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
    templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
    HSquared  = templateReal*templateReal + templateImag*templateImag ;
    temp = ((TwoDeltaToverN * HSquared) / psd->data->data[j]);
    this_snr  += temp;
    /* Now update re and im for the next iteration. */
    newRe = re + re*dre - im*dim;
    newIm = im + re*dim + im*dre;

    re = newRe;
    im = newIm;
  }

  /* Clean */
  if (freqHcross) XLALDestroyCOMPLEX16FrequencySeries(freqHcross);
  if (freqHplus) XLALDestroyCOMPLEX16FrequencySeries(freqHplus);
  if (waveFlags) XLALSimInspiralDestroyWaveformFlags(waveFlags);
  if (nonGRparams) XLALSimInspiralDestroyTestGRParam(nonGRparams);
  if (detector) free(detector);

  return sqrt(this_snr*2.0);

}
