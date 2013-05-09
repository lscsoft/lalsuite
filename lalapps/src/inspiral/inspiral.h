/*
*  Copyright (C) 2007 Duncan Brown, Patrick Brady
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

#ifndef _LALAPPS_INSPIRAL_H
#define _LALAPPS_INSPIRAL_H

#include <lalapps.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LALSimInspiral.h>
#include <lal/AVFactories.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Inject.h>
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>
#include <lal/FrameStream.h>


#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


#define INSPIRALH_ENULL   1
#define INSPIRALH_EFILE   2
#define INSPIRALH_ENONULL 3
#define INSPIRALH_ENOMEM  4
#define INSPIRALH_EVAL    5

#define INSPIRALH_MSGENULL    "Null pointer"
#define INSPIRALH_MSGEFILE    "Error in file-IO"
#define INSPIRALH_MSGENONULL  "Not a Null pointer"
#define INSPIRALH_MSGENOMEM   "Memory ellocation error"
#define INSPIRALH_MSGEVAL     "Invalid value"


REAL4 compute_candle_distance(
    REAL4 candleM1,
    REAL4 candleM2,
    REAL4 snr,
    REAL8 chanDeltaT,
    INT4 nPoints,
    REAL8FrequencySeries *spec,
    UINT4 cut);

REAL4 XLALCandleDistanceTD(
    Approximant approximant,
    REAL4 candleM1,
    REAL4 candleM2,
    REAL4 candlesnr,
    REAL8 chanDeltaT,
    INT4 nPoints,
    REAL8FrequencySeries *spec,
    UINT4 cut);

SummValueTable **add_summvalue_table(
    SummValueTable **newTable,
    LIGOTimeGPS gpsStartTime,
    LIGOTimeGPS gpsEndTime,
    const CHAR *programName,
    const CHAR *ifoName,
    const CHAR *summValueName,
    const CHAR *comment,
    REAL8 value
    );

void AddNumRelStrainModes( LALStatus              *status,
                           REAL4TimeVectorSeries  **outStrain,
                           SimInspiralTable *thisinj);

void AddNumRelStrainModesREAL8( LALStatus   *status,
                           REAL8TimeSeries  **outPlus,
                           REAL8TimeSeries  **outCross,
                           SimInspiralTable *thisinj);

void InjectNumRelWaveforms (LALStatus           *status,
                            REAL4TimeSeries     *chan,
                            SimInspiralTable    *injections,
                            CHAR                ifo[3],
                            REAL8               dynRange,
                            REAL8               freqLowCutoff,
                            REAL8               snrLow,
                            REAL8               snrHigh,
                            CHAR                *fnameOutXML);

void InjectNumRelWaveformsREAL8 (LALStatus      *status,
                            REAL8TimeSeries     *chan,
                            SimInspiralTable    *injections,
                            CHAR                ifo[3],
                            REAL8               freqLowCutoff,
                            REAL8               snrLow,
                            REAL8               snrHigh,
                            CHAR                *fname);

void InjectNumRelWaveformsUsingPSDREAL8(LALStatus *status,
                            REAL8TimeSeries      *chan,
                            SimInspiralTable     *injections,
                            CHAR                 ifo[3],
                            REAL8                freqLowCutoff,
                            REAL8                snrLow,
                            REAL8                snrHigh,
                            REAL8FrequencySeries *ligoPSD,
                            REAL8                ligoSnrLowFreq,
                            REAL8FrequencySeries *virgoPSD,
                            REAL8                virgoSnrLowFreq,
                            CHAR                 *fname);


REAL8 start_freq_from_frame_url(CHAR  *url);

REAL8 calculate_ligo_snr_from_strain(REAL4TimeVectorSeries *strain,
                                     SimInspiralTable      *thisInj,
                                     const CHAR            ifo[3]);

REAL8 calculate_ligo_snr_from_strain_real8(REAL8TimeSeries *strain,
                                       const CHAR            ifo[3]);

REAL8 calculate_snr_from_strain_and_psd_real8(REAL8TimeSeries *strain,
                                       REAL8FrequencySeries  *psd,
                                       REAL8                 startFreq,
                                       const CHAR            ifo[3]);

REAL8TimeSeries *XLALNRInjectionStrain(const char *ifo,
                                       SimInspiralTable *inj);

int XLALPsdFromFile (REAL8FrequencySeries **psd, const CHAR *filename);

REAL8FrequencySeries *XLALInterpolatePSD(
              REAL8FrequencySeries *in,
              REAL8 deltaFout
            );

REAL8 calculate_lalsim_snr(SimInspiralTable *inj, char *IFOname, REAL8FrequencySeries *psd, REAL8 start_freq);
void get_FakePsdFromString(REAL8FrequencySeries* PsdFreqSeries,char* FakePsdName, REAL8 StartFreq);


#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif           /* Close double-include protection */
