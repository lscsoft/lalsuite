/*
*  Copyright (C) 2015 Reinhard Prix
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

#include <lal/XLALError.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LogPrintf.h>
#include <lal/CWMakeFakeData.h>
#include <lal/LALConstants.h>
#include <lal/ExtrapolatePulsarSpins.h>
#include <lal/ComputeFstat.h>
#include <lal/DetectorStates.h>
#include <lal/LFTandTSutils.h>
#include <lal/LALString.h>
#include <lal/UserInput.h>

// benchmark ComputeFstat() functions for performance and memory usage

typedef struct
{
  CHAR *FstatMethod;		//!< select which method/algorithm to use to compute the F-statistic
  REAL8 Freq;
  REAL8 f1dot;
  REAL8Vector *FreqResolution;
  INT4Vector *numFreqBins;
  REAL8 Tseg;
  INT4 numSegments;
  LALStringVector *IFOs;
  CHAR *outputInfo;
  INT4 numTrials;

  // ----- developer options
  REAL8 Tsft;
  BOOLEAN reuseInput;   // only useful for checking workspace management
} UserInput_t;

// ---------- main ----------
int
main ( int argc, char *argv[] )
{
  // ---------- handle user input ----------
  UserInput_t XLAL_INIT_DECL(uvar_s);
  UserInput_t *uvar = &uvar_s;

  uvar->FstatMethod = XLALStringDuplicate("ResampBest");
  uvar->Freq = 100;
  uvar->f1dot = -3e-9;
  uvar->FreqResolution = XLALCreateREAL8Vector ( 2 );
  uvar->FreqResolution->data[0] = 1;
  uvar->FreqResolution->data[1] = 10;
  uvar->numFreqBins = XLALCreateINT4Vector ( 2 );
  uvar->numFreqBins->data[0] = 1000;
  uvar->numFreqBins->data[1] = 100000;
  uvar->Tseg = 60 * 3600;
  uvar->numSegments = 90;
  uvar->numTrials = 1;

  uvar->Tsft = 1800;
  uvar->reuseInput = 1;

  XLAL_CHECK ( (uvar->IFOs = XLALCreateStringVector ( "H1", NULL )) != NULL, XLAL_EFUNC );
  uvar->outputInfo = NULL;

  XLAL_CHECK ( XLALRegisterUvarMember ( FstatMethod,    STRING,         0, OPTIONAL,  "F-statistic method to use. Available methods: %s", XLALFstatMethodHelpString() ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( Freq,           REAL8,          0, OPTIONAL,  "Search frequency in Hz" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( f1dot,          REAL8,          0, OPTIONAL,  "Search spindown f1dot in Hz/s" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( FreqResolution, REAL8Vector,    0, OPTIONAL,  "Range of frequency resolution factor 'r' (st dFreq = 1/(r*T)) [2-number range input]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( Tseg,           REAL8,          0, OPTIONAL,  "Coherent segment length" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( numSegments,    INT4,           0, OPTIONAL,  "Number of semi-coherent segment" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( numFreqBins,    INT4Vector,     0, OPTIONAL,  "Range of number of frequency bins to search [2-number range input]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( IFOs,    	STRINGVector,   0, OPTIONAL,  "IFOs to use" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( numTrials,    	INT4,           0, OPTIONAL,  "Number of repeated trials to run (with potentially randomized parameters)" ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALRegisterUvarMember ( outputInfo,     STRING,         0, OPTIONAL, "Append Resampling internal info into this file") == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALRegisterUvarMember ( Tsft,           REAL8,          0, DEVELOPER, "SFT length" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember ( reuseInput,     BOOLEAN,        0, DEVELOPER, "Re-use FstatInput from previous setups (only useful for checking workspace management)" ) == XLAL_SUCCESS, XLAL_EFUNC );

  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    return EXIT_FAILURE;
  }
  // check user input
  XLAL_CHECK ( uvar->numSegments >= 1, XLAL_EINVAL );
  XLAL_CHECK ( (uvar->FreqResolution->length == 1) || (uvar->FreqResolution->length == 2), XLAL_EINVAL );
  XLAL_CHECK ( uvar->FreqResolution->data[0] > 0, XLAL_EINVAL );
  REAL8 FreqResolutionMin, FreqResolutionMax;
  FreqResolutionMin = FreqResolutionMax = uvar->FreqResolution->data[0];
  if ( uvar->FreqResolution->length == 2 )
    {
      XLAL_CHECK ( uvar->FreqResolution->data[1] > 0, XLAL_EINVAL );
      XLAL_CHECK ( uvar->FreqResolution->data[1] > uvar->FreqResolution->data[0], XLAL_EINVAL );
      FreqResolutionMax = uvar->FreqResolution->data[1];
    }

  XLAL_CHECK ( uvar->Freq > 0, XLAL_EINVAL );
  XLAL_CHECK ( uvar->Tseg > uvar->Tsft, XLAL_EINVAL );
  XLAL_CHECK ( uvar->Tsft > 1, XLAL_EINVAL );
  XLAL_CHECK ( (uvar->numFreqBins->length == 1) || (uvar->numFreqBins->length == 2), XLAL_EINVAL );
  XLAL_CHECK ( uvar->numFreqBins->data[0] > 0, XLAL_EINVAL );
  UINT4 numFreqBinsMax, numFreqBinsMin;
  numFreqBinsMin = numFreqBinsMax = uvar->numFreqBins->data[0];
  if ( uvar->numFreqBins->length == 2 )
    {
      XLAL_CHECK ( uvar->numFreqBins->data[1] > 0, XLAL_EINVAL );
      XLAL_CHECK ( uvar->numFreqBins->data[1] > uvar->numFreqBins->data[0], XLAL_EINVAL );
      numFreqBinsMax = uvar->numFreqBins->data[1];
    }

  XLAL_CHECK ( uvar->numTrials >= 1, XLAL_EINVAL );
  // ---------- end: handle user input ----------

  // common setup over repeated trials
  FstatMethodType FstatMethod;
  XLAL_CHECK ( XLALParseFstatMethodString ( &FstatMethod, uvar->FstatMethod ) == XLAL_SUCCESS, XLAL_EFUNC );

  EphemerisData *ephem;
  XLAL_CHECK ( (ephem = XLALInitBarycenter ( TEST_DATA_DIR "earth00-19-DE405.dat.gz", TEST_DATA_DIR "sun00-19-DE405.dat.gz" )) != NULL, XLAL_EFUNC );
  REAL8 memBase = XLALGetPeakHeapUsageMB();

  UINT4 numDetectors = uvar->IFOs->length;
  // ----- setup injection and data parameters
  LIGOTimeGPS startTime = {711595934, 0};
  LIGOTimeGPS startTime_l = startTime;
  LIGOTimeGPS endTime_l;
  SFTCatalog **catalogs;
  XLAL_CHECK ( (catalogs = XLALCalloc ( uvar->numSegments, sizeof( catalogs[0] ))) != NULL, XLAL_ENOMEM );

  for ( INT4 l = 0; l < uvar->numSegments; l ++ )
    {
      endTime_l = startTime_l;
      XLALGPSAdd( &endTime_l, uvar->Tseg );
      MultiLIGOTimeGPSVector *multiTimestamps;
      XLAL_CHECK ( (multiTimestamps = XLALMakeMultiTimestamps ( startTime_l, uvar->Tseg, uvar->Tsft, 0, numDetectors )) != NULL, XLAL_EFUNC );
      XLAL_CHECK ( (catalogs[l] = XLALMultiAddToFakeSFTCatalog ( NULL, uvar->IFOs, multiTimestamps )) != NULL, XLAL_EFUNC );
      XLALDestroyMultiTimestamps ( multiTimestamps );
      startTime_l = endTime_l;
    } // for l < numSegments
  LIGOTimeGPS endTime = endTime_l;
  UINT4 numSFTsPerSeg = catalogs[0]->length;

  PulsarSpinRange XLAL_INIT_DECL(spinRange);
  LIGOTimeGPS refTime = { startTime.gpsSeconds - 2.3 * uvar->Tseg, 0 };
  spinRange.refTime = refTime;
  spinRange.fkdot[0] = uvar->Freq;
  spinRange.fkdot[1] = uvar->f1dot;
  spinRange.fkdotBand[1] = 0;
  REAL8 asini = 0, Period = 0, ecc = 0;
  REAL8 minCoverFreq, maxCoverFreq;

  PulsarDopplerParams XLAL_INIT_DECL(Doppler);
  Doppler.refTime = refTime;
  Doppler.Alpha = 0.5;
  Doppler.Delta = 0.5;
  memcpy ( &Doppler.fkdot, &spinRange.fkdot, sizeof(Doppler.fkdot) );;
  Doppler.period = Period;
  Doppler.ecc = ecc;
  Doppler.asini = asini;

  // ----- setup optional Fstat arguments
  FstatOptionalArgs optionalArgs = FstatOptionalArgsDefaults;
  MultiNoiseFloor XLAL_INIT_DECL(injectSqrtSX);
  injectSqrtSX.length = numDetectors;
  for ( UINT4 X=0; X < numDetectors; X ++ ) {
    injectSqrtSX.sqrtSn[X] = 1;
  }
  optionalArgs.injectSqrtSX = &injectSqrtSX;
  optionalArgs.FstatMethod = FstatMethod;
  optionalArgs.collectTiming = 1;

  FILE *timingLogFILE = NULL;
  if ( uvar->outputInfo != NULL )
    {
      XLAL_CHECK ( (timingLogFILE = fopen (uvar->outputInfo, "ab")) != NULL, XLAL_ESYS, "Failed to open '%s' for appending\n", uvar->outputInfo );
    }

  FstatInputVector *inputs;
  FstatQuantities whatToCompute = (FSTATQ_2F | FSTATQ_2F_PER_DET);
  FstatResults *results = NULL;
  REAL8 tauF1NoBuf = 0;
  REAL8 tauF1Buf = 0;
  // ---------- main loop over repeated trials ----------
  for ( INT4 i = 0; i < uvar->numTrials; i ++ )
    {
      // randomize numFreqBins
      UINT4 numFreqBins_i = numFreqBinsMin + (UINT4)round ( 1.0 * (numFreqBinsMax - numFreqBinsMin) * rand() / RAND_MAX );
      // randomize FreqResolution
      REAL8 FreqResolution_i = FreqResolutionMin + 1.0 * ( FreqResolutionMax - FreqResolutionMin ) * rand() / RAND_MAX;

      XLAL_CHECK ( (inputs = XLALCreateFstatInputVector ( uvar->numSegments )) != NULL, XLAL_EFUNC );

      REAL8 dFreq = 1.0 / ( FreqResolution_i * uvar->Tseg );

      REAL8 FreqBand = numFreqBins_i * dFreq;
      fprintf ( stderr, "trial %d/%d: Tseg = %.1f d, numSegments = %d, Freq = %.1f Hz, f1dot = %.1e Hz/s, FreqResolution r = %f, numFreqBins = %d [dFreq = %.2e Hz, FreqBand = %.2e Hz]\n",
                i+1, uvar->numTrials, uvar->Tseg / 86400.0, uvar->numSegments, uvar->Freq, uvar->f1dot, FreqResolution_i, numFreqBins_i, dFreq, FreqBand );

      spinRange.fkdotBand[0] = FreqBand;
      XLAL_CHECK ( XLALCWSignalCoveringBand ( &minCoverFreq, &maxCoverFreq, &startTime, &endTime, &spinRange, asini, Period, ecc ) == XLAL_SUCCESS, XLAL_EFUNC );

      UINT4 numBinsSFT = ceil ( (maxCoverFreq - minCoverFreq) * uvar->Tsft + 2 * 8 );
      REAL8 memSFTs = uvar->numSegments * numSFTsPerSeg * ( sizeof(SFTtype) + numBinsSFT * sizeof(COMPLEX8)) / 1e6;

      // create per-segment input structs
      for ( INT4 l = 0; l < uvar->numSegments; l ++ )
        {
          if ( uvar->reuseInput && l > 0 ) {
            optionalArgs.prevInput = inputs->data[0];
          } else {
            optionalArgs.prevInput = NULL;
          }
          XLAL_CHECK ( (inputs->data[l] = XLALCreateFstatInput ( catalogs[l], minCoverFreq, maxCoverFreq, dFreq, ephem, &optionalArgs )) != NULL, XLAL_EFUNC );
        }

      // ----- compute Fstatistics over segments
      REAL8 tauF1NoBuf_i = 0;
      REAL8 tauF1Buf_i = 0;
      for ( INT4 l = 0; l < uvar->numSegments; l ++ )
        {
          XLAL_CHECK ( XLALComputeFstat ( &results, inputs->data[l], &Doppler, numFreqBins_i, whatToCompute ) == XLAL_SUCCESS, XLAL_EFUNC );

          REAL8 Fstat_tauF1Buf, Fstat_tauF1NoBuf;
          XLAL_CHECK ( XLALGetFstatTiming ( inputs->data[l], &Fstat_tauF1Buf, &Fstat_tauF1NoBuf ) == XLAL_SUCCESS, XLAL_EFUNC );
          tauF1NoBuf_i += Fstat_tauF1NoBuf;
          tauF1Buf_i   += Fstat_tauF1Buf;
          // ----- output timing details to file if requested
          if ( timingLogFILE != NULL ) {
            XLAL_CHECK ( AppendFstatTimingInfo2File ( inputs->data[l], timingLogFILE ) == XLAL_SUCCESS, XLAL_EFUNC );
          }

        } // for l < numSegments

      tauF1NoBuf_i /= uvar->numSegments;
      tauF1Buf_i   /= uvar->numSegments;
      REAL8 memMaxCompute = XLALGetPeakHeapUsageMB() - memBase;

      tauF1Buf   += tauF1Buf_i;
      tauF1NoBuf += tauF1NoBuf_i;

      fprintf (stderr, "%-15s: tauF1Buf = %.2g s, tauF1NoBuf = %.2g s, memSFTs = %.1f MB, memMaxCompute = %.1f MB\n",
               XLALGetFstatInputMethodName ( inputs->data[0] ), tauF1Buf_i, tauF1NoBuf_i, memSFTs, memMaxCompute  );

      XLALDestroyFstatInputVector ( inputs );
    } // for i < numTrials

  tauF1Buf   /= uvar->numTrials;
  tauF1NoBuf /= uvar->numTrials;

  fprintf (stderr, "\nAveraged timings: <tauF1Buf> = %.2g s, <tauF1NoBuf> = %.2g s\n", tauF1Buf, tauF1NoBuf );

  // ----- free memory ----------
  if ( timingLogFILE != NULL ) {
    fclose ( timingLogFILE );
  }

  for ( INT4 l = 0; l < uvar->numSegments; l ++ )
    {
      XLALDestroySFTCatalog ( catalogs[l] );
    }
  XLALFree ( catalogs );
  XLALDestroyFstatResults ( results );
  XLALDestroyUserVars();
  XLALDestroyEphemerisData ( ephem );

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()
