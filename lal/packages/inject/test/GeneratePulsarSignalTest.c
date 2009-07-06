/*
*  Copyright (C) 2007 Gregory Mendell, Jolien Creighton, Reinhard Prix
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

/************************* <lalVerbatim file="GeneratePulsarSignalTestCV">
Author: Mendell, G.
$Id$
**************************************************** </lalVerbatim> */

/* NOTES: */
/* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */
/* 10/12/04 gam; include INCLUDE_SEQUENTIAL_MISMATCH to make make pulsar parameters more varied. */
/* 10/12/04 gam; update error conditions. */
/* 02/02/05 gam; if NOT pSFTandSignalParams->resTrig > 0 should not create trigArg etc... */
/* 02/02/05 gam; fix warnings about setting ephemeris files pointers equal to constant strings and unused variables. */
/* 09/07/05 gam; Add Dterms parameter to LALFastGeneratePulsarSFTs; use this to fill in SFT bins with fake data as per LALDemod else fill in bin with zero */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{GeneratePulsarSignalTest.c}}
\label{ss:GeneratePulsarSignalTest.c}

\subsubsection*{Usage}
\begin{verbatim}
GeneratePulsarSignalTest
\end{verbatim}

No command line options are currently supported. However, preprocessor flags
can be set to print output for debugging purposes.

\subsubsection*{Description}

This test program calls and compares the output from LALGeneratePulsarSignal
and LALSignalToSFTs with the output from LALComputeSkyAndZeroPsiAMResponse
and LALFastGeneratePulsarSFTs for a variety of signal parameters.

The current code only compares the modulus of the output SFTs,
not the phases.

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define GENERATEPULSARSIGNALTESTC_ENORM  0
#define GENERATEPULSARSIGNALTESTC_EIFO   1
#define GENERATEPULSARSIGNALTESTC_EMOD   2
#define GENERATEPULSARSIGNALTESTC_EBIN   3
#define GENERATEPULSARSIGNALTESTC_EBINS  4

#define GENERATEPULSARSIGNALTESTC_MSGENORM  "Normal exit"
#define GENERATEPULSARSIGNALTESTC_MSGEIFO   "IFO not supported"
#define GENERATEPULSARSIGNALTESTC_MSGEMOD   "SFT max power from LALSignalToSFTs and LALFastGeneratePulsarSFTs differs"
#define GENERATEPULSARSIGNALTESTC_MSGEBIN  "SFT freq with max power from LALSignalToSFTs and LALFastGeneratePulsarSFTs differs by more than 1 bin"
#define GENERATEPULSARSIGNALTESTC_MSGEBINS "SFTs freq with max power from LALSignalToSFTs and LALFastGeneratePulsarSFTs differs too often"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Notes}

See the pulsar search code in lalapps for more
example uses of the functions tested by this code.

\vfill{\footnotesize\input{GeneratePulsarSignalTestCV}}

******************************************************* </lalLaTeX> */

/* preprocessor flags that control output */
/* First two are used with PRINT_OUTPUTSFT
   and PRINT_FASTOUTPUTSFT to specify which
   SFT from which test to print to stdout */
#define TESTNUMBER_TO_PRINT 1
#define SFTINDEX_TO_PRINT 0
/* #define PRINT_OUTPUTSFT */
/* #define PRINT_FASTOUTPUTSFT */
/* #define PRINT_MAXSFTPOWER */
/* #define PRINT_MAXDIFFSFTPOWER */
/* #define PRINT_DIFFATMAXPOWER */
/* #define PRINT_ERRORATMAXPOWER */
/* #define PRINT_OVERALLMAXDIFFSFTPOWER */
#define INCLUDE_SEQUENTIAL_MISMATCH
/* #define INCLUDE_RANDVAL_MISMATCH */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StreamInput.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/VectorOps.h>
#include <lal/DetectorSite.h>
#include <lal/Units.h>
#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Date.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/Random.h>

/* prototype for function that runs the tests */
/* void RunGeneratePulsarSignalTest(LALStatus *status, int argc, char **argv); */ /* 02/02/05 gam */
void RunGeneratePulsarSignalTest(LALStatus *status);

NRCSID( GENERATEPULSARSIGNALTESTC, "$Id$" );

/*********************************************************/
/*                                                       */
/* START FUNCTION: main                                  */
/*                                                       */
/*********************************************************/
int main(int argc, char **argv)
{

  int mainArgc = 0; /* 02/02/05 gam; */
  static LALStatus status; /* status structure */

  lalDebugLevel = 0;

  /* initialize status */
  status.statusCode = 0;
  status.statusPtr = NULL;
  status.statusDescription = NULL;

  mainArgc = argc; /* 02/02/05 gam; get rid of warning about unused argc */

  /* Actual test are performed by this function */
  /* RunGeneratePulsarSignalTest(&status,argc,argv); */ /* 02/02/05 gam */
  RunGeneratePulsarSignalTest(&status);

  if (status.statusCode) {
     fprintf(stderr,"Error: statusCode = %i statusDescription = %s \n", status.statusCode, status.statusDescription);
     return 1;
  }

  LALCheckMemoryLeaks();

  if ( lalDebugLevel & LALINFO ) {
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"
       "        %s\n", *argv, __FILE__, __LINE__,
       GENERATEPULSARSIGNALTESTC, (GENERATEPULSARSIGNALTESTC_MSGENORM) );
  }

  return GENERATEPULSARSIGNALTESTC_ENORM;
}
/*********************************************************/
/*                                                       */
/* END FUNCTION: main                                    */
/*                                                       */
/*********************************************************/

/*********************************************************/
/*                                                       */
/* START FUNCTION: RunGeneratePulsarSignalTest           */
/*                                                       */
/*********************************************************/
/* void RunGeneratePulsarSignalTest(LALStatus *status, int argc, char **argv) */ /* 02/02/05 gam */
void RunGeneratePulsarSignalTest(LALStatus *status)
{
  INT4    i, j, k;         /* all purpose indices */
  INT4    iSky = 0;        /* for loop over sky positions */
  INT4    jDeriv = 0;      /* for loop over spindowns */
  INT4    testNumber = 0;  /* which test is being run */

  /* The input and output structs used by the functions being tested */
  PulsarSignalParams *pPulsarSignalParams = NULL;
  REAL4TimeSeries *signalvec = NULL;
  SFTParams *pSFTParams = NULL;
  SkyConstAndZeroPsiAMResponse *pSkyConstAndZeroPsiAMResponse;
  SFTandSignalParams *pSFTandSignalParams;
  SFTVector *outputSFTs = NULL;
  SFTVector *fastOutputSFTs = NULL;
  REAL4 renorm; /* to renormalize SFTs to account for different effective sample rates */

  /* containers for detector and ephemeris data */
  LALDetector cachedDetector;
  CHAR IFO[6] = "LHO";
  EphemerisData *edat = NULL;
  CHAR sunFile[] = "sun00-04.dat";     /* 02/02/05 gam */
  CHAR earthFile[] = "earth00-04.dat"; /* 02/02/05 gam */
  INT4 leap; /* 2nd arg to LALLeapSecFormatAndAcc is INT4 while edat->leap is INT2. */

  /* containers for sky position and spindown data */
  REAL8 **skyPosData;
  REAL8 **freqDerivData;
  INT4 numSkyPosTotal = 8;
  INT4 numSpinDown = 1;
  INT4 numFreqDerivTotal = 2;
  INT4 numFreqDerivIncludingNoSpinDown; /* set below */
  INT4 iFreq = 0;
  REAL8 cosTmpDEC;
  REAL8 tmpDeltaRA;
  REAL8 DeltaRA = 0.01;
  REAL8 DeltaDec = 0.01;
  REAL8 DeltaFDeriv1 = -1.0e-10;
  REAL8 DeltaFDeriv2 = 0.0;
  REAL8 DeltaFDeriv3 = 0.0;
  REAL8 DeltaFDeriv4 = 0.0;
  REAL8 DeltaFDeriv5 = 0.0;

  /* containers for number of SFTs, times to generate SFTs, reference time, and gap */
  INT4 numSFTs = 4;
  REAL8 f0SFT = 943.12;
  REAL8 bandSFT = 0.5;
  UINT4 gpsStartTimeSec = 765432109; /* GPS start time of data requested seconds; example = Apr 08 2004 04:01:36 UTC */
  UINT4 gpsStartTimeNan = 0;          /* GPS start time of data requested nanoseconds */
  LIGOTimeGPSVector *timeStamps;
  LIGOTimeGPS GPSin;   /* holder for reference-time for pulsar parameters at the detector; will convert to SSB! */
  REAL8 sftGap = 73.0; /* extra artificial gap between SFTs */
  REAL8 tSFT   = 1800.0;
  REAL8 duration = tSFT + (numSFTs - 1)*(tSFT + sftGap);
  INT4 nBinsSFT = (INT4)(bandSFT*tSFT + 0.5);

  /* additional parameters that determine what signals to test; note f0SGNL and bandSGNL must be compatible with f0SFT and bandSFT */
  REAL8 f0SGNL = f0SFT + bandSFT/2.0;
  REAL8 dfSGNL = 1.0/tSFT;
  REAL8 bandSGNL = 0.005;
  INT4  nBinsSGNL = (INT4)(bandSGNL*tSFT + 0.5);
  INT4  nsamples = 18000; /* nsamples from SFT header; 2 x this would be the effective number of time samples used to create an SFT from raw data */
  INT4  Dterms = 3;       /* 09/07/05 gam; use Dterms to fill in SFT bins with fake data as per LALDemod else fill in bin with zero */
  REAL8 h_0 = 7.0e-22;    /* Source amplitude; use arbitrary small number for default */
  REAL8 cosIota;        /* cosine of inclination angle iota of the source */

  /* variables for comparing differences */
  REAL4 maxDiffSFTMod, diffAtMaxPower;
  REAL4 overallMaxDiffSFTMod;
  REAL4 maxMod, fastMaxMod;
  INT4  jMaxMod, jFastMaxMod,jMaxDiff,jOverallMaxDiff;
  INT4  iOverallMaxDiffSFTMod;
  REAL4 tmpDiffSFTMod, sftMod, fastSFTMod;
  REAL4 smallMod = 1.e-30;
  REAL4 epsDiffMod;
  REAL4 epsBinErrorRate;   /* 10/12/04 gam; Allowed bin error rate */
  INT4  binErrorCount = 0; /* 10/12/04 gam; Count number of bin errors  */

  /* randval is always set to a default value or given a random value to generate certain signal parameters or mismatch */
  REAL4 randval;
  #ifdef INCLUDE_RANDVAL_MISMATCH
   INT4 seed=0;
   INT4 rndCount;
   RandomParams *randPar=NULL;
   FILE *fpRandom;
  #endif

  INITSTATUS( status, "RunGeneratePulsarSignalTest", GENERATEPULSARSIGNALTESTC );
  ATTATCHSTATUSPTR(status);

  /* generate timeStamps */
  timeStamps = (LIGOTimeGPSVector *)LALMalloc(sizeof(LIGOTimeGPSVector));
  timeStamps->data =(LIGOTimeGPS *)LALMalloc (numSFTs*sizeof(LIGOTimeGPS));
  timeStamps->length = numSFTs;
  for (i = 0; i < numSFTs; i++) {
      timeStamps->data[i].gpsSeconds = gpsStartTimeSec + (UINT4)(i*(tSFT + sftGap));
      timeStamps->data[i].gpsNanoSeconds = 0;
  } /* for i < numSFTs */

  /* generate skyPosData */
  skyPosData=(REAL8 **)LALMalloc(numSkyPosTotal*sizeof(REAL8 *));
  for(iSky=0;iSky<numSkyPosTotal;iSky++)
  {
        skyPosData[iSky] = (REAL8 *)LALMalloc(2*sizeof(REAL8));
        /* Try fairly random sky positions skyPosData[iSky][0] = RA, skyPosData[iSky][1] = DEC */
        if (iSky == 0) {
          skyPosData[iSky][0] = 0.02*LAL_TWOPI;
          skyPosData[iSky][1] = 0.03*LAL_PI/2.0;
        } else if (iSky == 1) {
          skyPosData[iSky][0] = 0.23*LAL_TWOPI;
          skyPosData[iSky][1] = -0.57*LAL_PI/2.0;
        } else if (iSky == 2) {
          skyPosData[iSky][0] = 0.47*LAL_TWOPI;
          skyPosData[iSky][1] = 0.86*LAL_PI/2.0;
        } else if (iSky == 3) {
          skyPosData[iSky][0] = 0.38*LAL_TWOPI;
          skyPosData[iSky][1] = -0.07*LAL_PI/2.0;
        } else if (iSky == 4) {
          skyPosData[iSky][0] = 0.65*LAL_TWOPI;
          skyPosData[iSky][1] = 0.99*LAL_PI/2.0;
        } else if (iSky == 5) {
          skyPosData[iSky][0] = 0.72*LAL_TWOPI;
          skyPosData[iSky][1] = -0.99*LAL_PI/2.0;
        } else if (iSky == 6) {
          skyPosData[iSky][0] = 0.81*LAL_TWOPI;
          skyPosData[iSky][1] = 0.19*LAL_PI/2.0;
        } else if (iSky == 7) {
          skyPosData[iSky][0] = 0.99*LAL_TWOPI;
          skyPosData[iSky][1] = 0.01*LAL_PI/2.0;
        } else {
          skyPosData[iSky][0] = 0.0;
          skyPosData[iSky][1] = 0.0;
        } /* END if (k == 0) ELSE ... */
  } /* END for(iSky=0;iSky<numSkyPosTotal;iSky++) */

  freqDerivData = NULL; /* 02/02/05 gam */
  if (numSpinDown > 0) {
    freqDerivData=(REAL8 **)LALMalloc(numFreqDerivTotal*sizeof(REAL8 *));
    for(jDeriv=0;jDeriv<numFreqDerivTotal;jDeriv++) {
        freqDerivData[jDeriv] = (REAL8 *)LALMalloc(numSpinDown*sizeof(REAL8));
        if (jDeriv == 0) {
          for(k=0;k<numSpinDown;k++) {
            freqDerivData[jDeriv][k] = 0.0;
          }
        } else if (jDeriv == 1) {
          for(k=0;k<numSpinDown;k++) {
            freqDerivData[jDeriv][k] = 1.e-9; /* REALLY THIS IS ONLY GOOD FOR numSpinDown = 1; */
          }
        } else {
          for(k=0;k<numSpinDown;k++) {
            freqDerivData[jDeriv][k] = 0.0;
          }
        } /* END if (k == 0) ELSE ... */
    } /* END for(jDeriv=0;jDeriv<numFreqDerivTotal;jDeriv++) */
    numFreqDerivIncludingNoSpinDown = numFreqDerivTotal;
  } else {
    numFreqDerivIncludingNoSpinDown = 1;  /* Even if numSpinDown = 0 still need to count case of zero spindown. */
  }  /* END if (numSpinDown > 0) ELSE ... */

  /* Initialize ephemeris data */
  edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
  /* edat->ephiles.sunEphemeris = "../../pulsar/test/sun00-04.dat";
  edat->ephiles.earthEphemeris = "../../pulsar/test/earth00-04.dat"; */
  /* 07/30/04 gam; added this line to Makefile.am: TESTS_ENVIRONMENT = LAL_DATA_PATH=$(top_srcdir)/packages/pulsar/test */
  /* edat->ephiles.sunEphemeris = "sun00-04.dat";
  edat->ephiles.earthEphemeris = "earth00-04.dat"; */
  edat->ephiles.sunEphemeris = sunFile;     /* 02/02/05 gam */
  edat->ephiles.earthEphemeris = earthFile; /* 02/02/05 gam */
  leap = XLALLeapSeconds ( timeStamps->data[0].gpsSeconds ); /**< [In] Seconds relative to GPS epoch.*/
  edat->leap = (INT2)leap;
  LALInitBarycenter(status->statusPtr, edat);
  CHECKSTATUSPTR (status);

  /* Allocate memory for PulsarSignalParams and initialize */
  pPulsarSignalParams = (PulsarSignalParams *)LALMalloc(sizeof(PulsarSignalParams));
  pPulsarSignalParams->pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;
  pPulsarSignalParams->pulsar.spindown = NULL;
  if (numSpinDown > 0) {
    LALDCreateVector(status->statusPtr, &(pPulsarSignalParams->pulsar.spindown),((UINT4)numSpinDown));
    CHECKSTATUSPTR (status);
  }
  pPulsarSignalParams->orbit = NULL;
  pPulsarSignalParams->transfer = NULL;
  /* Set up pulsar site */
  if (strstr(IFO, "LHO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  } else if (strstr(IFO, "LLO")) {
       cachedDetector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  } else if (strstr(IFO, "GEO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  } else if (strstr(IFO, "VIRGO")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  } else if (strstr(IFO, "TAMA")) {
      cachedDetector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  } else {
      /* "Invalid or null IFO" */
      ABORT( status, GENERATEPULSARSIGNALTESTC_EIFO, GENERATEPULSARSIGNALTESTC_MSGEIFO);
  }
  pPulsarSignalParams->site = &cachedDetector;
  pPulsarSignalParams->ephemerides = edat;
  pPulsarSignalParams->startTimeGPS.gpsSeconds = (INT4)gpsStartTimeSec;
  pPulsarSignalParams->startTimeGPS.gpsNanoSeconds = (INT4)gpsStartTimeNan;
  pPulsarSignalParams->duration = (UINT4)duration;
  pPulsarSignalParams->samplingRate = (REAL8)ceil(2.0*bandSFT); /* Make sampleRate an integer so that T*samplingRate = integer for integer T */
  pPulsarSignalParams->fHeterodyne = f0SFT;

  GPSin.gpsSeconds = timeStamps->data[0].gpsSeconds;
  GPSin.gpsNanoSeconds = timeStamps->data[0].gpsNanoSeconds;

  /* Allocate memory for SFTParams and initialize */
  pSFTParams = (SFTParams *)LALCalloc(1, sizeof(SFTParams));
  pSFTParams->Tsft = tSFT;
  pSFTParams->timestamps = timeStamps;
  pSFTParams->noiseSFTs = NULL;
  pSFTParams->make_v2SFTs = 1;

  #ifdef INCLUDE_RANDVAL_MISMATCH
    /* Initial seed and randPar to use LALUniformDeviate to generate random mismatch during Monte Carlo. */
    fpRandom = fopen("/dev/urandom","r");
    rndCount = fread(&seed, sizeof(INT4),1, fpRandom);
    fclose(fpRandom);
    /* seed = 1234; */ /* Test value */
    LALCreateRandomParams(status->statusPtr, &randPar, seed);
    CHECKSTATUSPTR (status);
  #endif

  /* allocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  pSkyConstAndZeroPsiAMResponse = (SkyConstAndZeroPsiAMResponse *)LALMalloc(sizeof(SkyConstAndZeroPsiAMResponse));
  pSkyConstAndZeroPsiAMResponse->skyConst = (REAL8 *)LALMalloc((2*numSpinDown*(numSFTs+1)+2*numSFTs+3)*sizeof(REAL8));
  pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi = (REAL4 *)LALMalloc(numSFTs*sizeof(REAL4));
  pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi = (REAL4 *)LALMalloc(numSFTs*sizeof(REAL4));
  pSFTandSignalParams = (SFTandSignalParams *)LALMalloc(sizeof(SFTandSignalParams));
  /* create lookup table (LUT) values for doing trig */
  /* pSFTandSignalParams->resTrig = 64; */ /* length sinVal and cosVal; resolution of trig functions = 2pi/resTrig */
  /* pSFTandSignalParams->resTrig = 128; */ /* 10/08/04 gam; length sinVal and cosVal; domain = -2pi to 2pi inclusive; resolution = 4pi/resTrig */
  pSFTandSignalParams->resTrig = 0; /* 10/12/04 gam; turn off using LUTs since this is more typical. */
  /* 02/02/05 gam; if NOT pSFTandSignalParams->resTrig > 0 should not create trigArg etc... */
  if (pSFTandSignalParams->resTrig > 0) {
    pSFTandSignalParams->trigArg = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
    pSFTandSignalParams->sinVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
    pSFTandSignalParams->cosVal  = (REAL8 *)LALMalloc((pSFTandSignalParams->resTrig+1)*sizeof(REAL8));
    for (k=0; k<=pSFTandSignalParams->resTrig; k++) {
       /* pSFTandSignalParams->trigArg[k]= ((REAL8)LAL_TWOPI) * ((REAL8)k) / ((REAL8)pSFTandSignalParams->resTrig); */ /* 10/08/04 gam */
       pSFTandSignalParams->trigArg[k]= -1.0*((REAL8)LAL_TWOPI) + 2.0 * ((REAL8)LAL_TWOPI) * ((REAL8)k) / ((REAL8)pSFTandSignalParams->resTrig);
       pSFTandSignalParams->sinVal[k]=sin( pSFTandSignalParams->trigArg[k] );
       pSFTandSignalParams->cosVal[k]=cos( pSFTandSignalParams->trigArg[k] );
    }
  }
  pSFTandSignalParams->pSigParams = pPulsarSignalParams;
  pSFTandSignalParams->pSFTParams = pSFTParams;
  pSFTandSignalParams->nSamples = nsamples;
  pSFTandSignalParams->Dterms = Dterms;  /* 09/07/05 gam */

  /*********************************************************/
  /*                                                       */
  /* START SECTION: LOOP OVER SKY POSITIONS                */
  /*                                                       */
  /*********************************************************/
  for(iSky=0;iSky<numSkyPosTotal;iSky++) {

    /* set source sky position declination (DEC) */
    randval = 0.5; /* Gives default value */
    #ifdef INCLUDE_SEQUENTIAL_MISMATCH
      randval = ( (REAL4)(iSky) )/( (REAL4)(numSkyPosTotal) );
    #endif
    #ifdef INCLUDE_RANDVAL_MISMATCH
       LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
    #endif
    pPulsarSignalParams->pulsar.position.latitude = skyPosData[iSky][1] + (((REAL8)randval) - 0.5)*DeltaDec;
    cosTmpDEC = cos(skyPosData[iSky][1]);
    if (cosTmpDEC != 0.0) {
          tmpDeltaRA = DeltaRA/cosTmpDEC;
    } else {
          tmpDeltaRA = 0.0; /* We are at a celestial pole */
    }

    /* set source sky position right ascension (RA) */
    randval = 0.5; /* Gives default value */
    #ifdef INCLUDE_SEQUENTIAL_MISMATCH
      randval = ( (REAL4)(iSky) )/( (REAL4)(numSkyPosTotal) );
    #endif
    #ifdef INCLUDE_RANDVAL_MISMATCH
      LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
    #endif
    pPulsarSignalParams->pulsar.position.longitude = skyPosData[iSky][0] + (((REAL8)randval) - 0.5)*tmpDeltaRA;

    /* Find reference time in SSB for this sky positions */
    LALConvertGPS2SSB(status->statusPtr,&(pPulsarSignalParams->pulsar.refTime), GPSin, pPulsarSignalParams);
    CHECKSTATUSPTR (status);

    /* one per sky position fill in SkyConstAndZeroPsiAMResponse for use with LALFastGeneratePulsarSFTs */
    LALComputeSkyAndZeroPsiAMResponse (status->statusPtr, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
    CHECKSTATUSPTR (status);

    /*********************************************************/
    /*                                                       */
    /* START SECTION: LOOP OVER SPINDOWN                     */
    /*                                                       */
    /*********************************************************/
    for(jDeriv=0;jDeriv<numFreqDerivIncludingNoSpinDown;jDeriv++) {
     /* source spindown parameters */
     if (numSpinDown > 0) {
       for(k=0;k<numSpinDown;k++) {
         randval = 0.5; /* Gives default value */
         #ifdef INCLUDE_SEQUENTIAL_MISMATCH
            if (freqDerivData[jDeriv][k] < 0.0) {
               randval = ( (REAL4)(iSky) )/( (REAL4)(numSkyPosTotal) );
            } else {
               randval = 0.5; /* If derivative is not negative (i.e., it is zero) then do not add in mismatch; keep it zero. */
            }
         #endif
         #ifdef INCLUDE_RANDVAL_MISMATCH
           LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
         #endif
         if (k == 0) {
          pPulsarSignalParams->pulsar.spindown->data[k] = freqDerivData[jDeriv][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv1;
         } else if (k == 1) {
          pPulsarSignalParams->pulsar.spindown->data[k] = freqDerivData[jDeriv][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv2;
         } else if (k == 2) {
          pPulsarSignalParams->pulsar.spindown->data[k] = freqDerivData[jDeriv][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv3;
         } else if (k == 3) {
          pPulsarSignalParams->pulsar.spindown->data[k] = freqDerivData[jDeriv][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv4;
         } else if (k == 4) {
          pPulsarSignalParams->pulsar.spindown->data[k] = freqDerivData[jDeriv][k] + (((REAL8)randval) - 0.5)*DeltaFDeriv5;
         } /* END if (k == 0) ELSE ... */
       }
     }

     /****************************************************/
     /*                                                  */
     /* START SECTION: LOOP OVER FREQUENCIES             */
     /*                                                  */
     /****************************************************/
     for(iFreq=0;iFreq<nBinsSGNL;iFreq++) {

       /* set source orientation psi */
       randval = 0.5; /* Gives default value */
       #ifdef INCLUDE_SEQUENTIAL_MISMATCH
            randval = ( (REAL4)(iFreq) )/( (REAL4)(nBinsSGNL) );
       #endif
       #ifdef INCLUDE_RANDVAL_MISMATCH
          LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
       #endif
       pPulsarSignalParams->pulsar.psi = (randval - 0.5) * ((REAL4)LAL_PI_2);

       /* set angle between source spin axis and direction from source to SSB, cosIota */
       randval = 1.0; /* Gives default value */
       #ifdef INCLUDE_SEQUENTIAL_MISMATCH
            randval = ( (REAL4)(iFreq) )/( (REAL4)(nBinsSGNL) );
       #endif
       #ifdef INCLUDE_RANDVAL_MISMATCH
          LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
       #endif
       cosIota = 2.0*((REAL8)randval) - 1.0;

       /* h_0 is fixed above; get A_+ and A_x from h_0 and cosIota */
       pPulsarSignalParams->pulsar.aPlus = (REAL4)(0.5*h_0*(1.0 + cosIota*cosIota));
       pPulsarSignalParams->pulsar.aCross = (REAL4)(h_0*cosIota);

       /* get random value for phi0 */
       randval = 0.125; /* Gives default value pi/4*/
       #ifdef INCLUDE_SEQUENTIAL_MISMATCH
            randval = ( (REAL4)(iFreq) )/( (REAL4)(nBinsSGNL) );
       #endif
       #ifdef INCLUDE_RANDVAL_MISMATCH
          LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
       #endif
       pPulsarSignalParams->pulsar.phi0 = ((REAL8)randval) * ((REAL8)LAL_TWOPI);

       /* randval steps through various mismatches to frequencies centered on a bin */
       /* Note that iFreq = nBinsSGNL/2 gives randval = 0.5 gives no mismatch from frequency centered on a bin */
       randval = ( (REAL4)(iFreq) )/( (REAL4)(nBinsSGNL) );
       #ifdef INCLUDE_RANDVAL_MISMATCH
          LALUniformDeviate(status->statusPtr, &randval, randPar); CHECKSTATUSPTR (status);
       #endif
       pPulsarSignalParams->pulsar.f0 = f0SGNL + iFreq*dfSGNL + (((REAL8)randval) - 0.5)*dfSGNL;

       testNumber++; /* Update count of which test we about to do. */

       /* FIRST: Use LALGeneratePulsarSignal and LALSignalToSFTs to generate outputSFTs */
       signalvec = NULL;
       LALGeneratePulsarSignal(status->statusPtr, &signalvec, pPulsarSignalParams);
       CHECKSTATUSPTR (status);
       outputSFTs = NULL;
       LALSignalToSFTs(status->statusPtr, &outputSFTs, signalvec, pSFTParams);
       CHECKSTATUSPTR (status);

       #ifdef PRINT_OUTPUTSFT
           if (testNumber == TESTNUMBER_TO_PRINT) {
              i=SFTINDEX_TO_PRINT; /* index of which outputSFT to output */
              fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
              fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
              fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
              fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
              fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",iFreq,f0SGNL + iFreq*dfSGNL,pPulsarSignalParams->pulsar.f0);
              fprintf(stdout,"outputSFTs->data[%i].data->length = %i \n",i,outputSFTs->data[i].data->length);
              renorm = ((REAL4)nsamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
              for(j=0;j<nBinsSFT;j++) {
                /* fprintf(stdout,"%i %g %g \n", j, renorm*outputSFTs->data[i].data->data[j].re, renorm*outputSFTs->data[i].data->data[j].im); */
                fprintf(stdout,"%i %g \n",j,renorm*renorm*outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + renorm*renorm*outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im);
                fflush(stdout);
              }
           }
       #endif

       /* SECOND: Use LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs to generate outputSFTs */
       LALFastGeneratePulsarSFTs (status->statusPtr, &fastOutputSFTs, pSkyConstAndZeroPsiAMResponse, pSFTandSignalParams);
       CHECKSTATUSPTR (status);

       #ifdef PRINT_FASTOUTPUTSFT
          if (testNumber == TESTNUMBER_TO_PRINT) {
            REAL4  fPlus;
            REAL4  fCross;
            i=SFTINDEX_TO_PRINT; /* index of which outputSFT to output */
            fPlus = pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi[i]*cos(2.0*pPulsarSignalParams->pulsar.psi) + pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi[i]*sin(2.0*pPulsarSignalParams->pulsar.psi);
            fCross = pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi[i]*cos(2.0*pPulsarSignalParams->pulsar.psi) - pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi[i]*sin(2.0*pPulsarSignalParams->pulsar.psi);
            fprintf(stdout,"iFreq = %i, inject h_0 = %23.10e \n",iFreq,h_0);
            fprintf(stdout,"iFreq = %i, inject cosIota = %23.10e, A_+ = %23.10e, A_x = %23.10e \n",iFreq,cosIota,pPulsarSignalParams->pulsar.aPlus,pPulsarSignalParams->pulsar.aCross);
            fprintf(stdout,"iFreq = %i, inject psi = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.psi);
            fprintf(stdout,"iFreq = %i, fPlus, fCross = %23.10e,  %23.10e \n",iFreq,fPlus,fCross);
            fprintf(stdout,"iFreq = %i, inject phi0 = %23.10e \n",iFreq,pPulsarSignalParams->pulsar.phi0);
            fprintf(stdout,"iFreq = %i, search f0 = %23.10e, inject f0 = %23.10e \n",iFreq,f0SGNL + iFreq*dfSGNL,pPulsarSignalParams->pulsar.f0);
            fprintf(stdout,"fastOutputSFTs->data[%i].data->length = %i \n",i,fastOutputSFTs->data[i].data->length);
            fflush(stdout);
            for(j=0;j<nBinsSFT;j++) {
               /* fprintf(stdout,"%i %g %g \n",j,fastOutputSFTs->data[i].data->data[j].re,fastOutputSFTs->data[i].data->data[j].im); */
               fprintf(stdout,"%i %g \n",j,fastOutputSFTs->data[i].data->data[j].re*fastOutputSFTs->data[i].data->data[j].re + fastOutputSFTs->data[i].data->data[j].im*fastOutputSFTs->data[i].data->data[j].im);
               fflush(stdout);
            }
          }
       #endif

       /* find maximum difference in power */
       epsDiffMod = 0.20; /* maximum allowed percent difference */ /* 10/12/04 gam */
       overallMaxDiffSFTMod = 0.0;
       iOverallMaxDiffSFTMod = -1;
       jOverallMaxDiff = -1;
       for (i = 0; i < numSFTs; i++) {
          renorm = ((REAL4)nsamples)/((REAL4)(outputSFTs->data[i].data->length - 1));
          maxDiffSFTMod = 0.0;
          diffAtMaxPower = 0.0;
          maxMod = 0.0;
          fastMaxMod = 0.0;
          jMaxMod = -1;
          jFastMaxMod = -1;
          jMaxDiff = -1;
          /* Since doppler shifts can move the signal by an unknown number of bins search the whole band for max modulus: */
          for(j=0;j<nBinsSFT;j++) {
               sftMod = renorm*renorm*outputSFTs->data[i].data->data[j].re*outputSFTs->data[i].data->data[j].re + renorm*renorm*outputSFTs->data[i].data->data[j].im*outputSFTs->data[i].data->data[j].im;
               sftMod = sqrt(sftMod);
               fastSFTMod = fastOutputSFTs->data[i].data->data[j].re*fastOutputSFTs->data[i].data->data[j].re + fastOutputSFTs->data[i].data->data[j].im*fastOutputSFTs->data[i].data->data[j].im;
               fastSFTMod = sqrt(fastSFTMod);
               if (fabs(sftMod) > smallMod) {
                   tmpDiffSFTMod = fabs((sftMod - fastSFTMod)/sftMod);
                   if (tmpDiffSFTMod > maxDiffSFTMod) {
                       maxDiffSFTMod = tmpDiffSFTMod;
                       jMaxDiff = j;
                   }
                   if (tmpDiffSFTMod > overallMaxDiffSFTMod) {
                       overallMaxDiffSFTMod = tmpDiffSFTMod;
                       iOverallMaxDiffSFTMod = i;
                       jOverallMaxDiff = j;
                   }
                   if (sftMod > maxMod) {
                       maxMod = sftMod;
                       jMaxMod = j;
                   }
                   if (fastSFTMod > fastMaxMod) {
                       fastMaxMod = fastSFTMod;
                       jFastMaxMod = j;
                   }
               }
          }
          if (fabs(maxMod) > smallMod) {
             diffAtMaxPower = fabs((maxMod - fastMaxMod)/maxMod);
          }
          #ifdef PRINT_MAXSFTPOWER
            fprintf(stdout,"maxSFTMod, testNumber %i, SFT %i, bin %i = %g \n",testNumber,i,jMaxMod,maxMod);
            fprintf(stdout,"maxFastSFTMod, testNumber %i, SFT %i, bin %i = %g \n",testNumber,i,jFastMaxMod,fastMaxMod);
            fflush(stdout);
          #endif
          #ifdef PRINT_MAXDIFFSFTPOWER
            fprintf(stdout,"maxDiffSFTMod, testNumber %i, SFT %i, bin %i = %g \n",testNumber,i,jMaxDiff,maxDiffSFTMod);
            fflush(stdout);
          #endif
          #ifdef PRINT_DIFFATMAXPOWER
            fprintf(stdout,"diffAtMaxPower, testNumber %i, SFT %i, bins %i and %i = %g \n",testNumber,i,jMaxMod,jFastMaxMod,diffAtMaxPower);
            fflush(stdout);
          #endif
          #ifdef PRINT_ERRORATMAXPOWER
            if (diffAtMaxPower > epsDiffMod) {
              fprintf(stdout,"diffAtMaxPower, testNumber %i, SFT %i, bins %i and %i = %g exceeded epsDiffMod = %g \n",testNumber,i,jMaxMod,jFastMaxMod,diffAtMaxPower,epsDiffMod);
              fflush(stdout);
              /* break; */ /* only report 1 error per test */
            }
            if (jMaxMod != jFastMaxMod) {
              fprintf(stdout,"MaxPower occurred in different bins: testNumber %i, SFT %i, bins %i and %i\n",testNumber,i,jMaxMod,jFastMaxMod);
              fflush(stdout);
              /* break; */ /* only report 1 error per test */
            }
          #endif
          if (jMaxMod != jFastMaxMod) {
              binErrorCount++; /* 10/12/04 gam; count up bin errors; if too ABORT at bottom of code */
          }
          if ( diffAtMaxPower > epsDiffMod ) {
            ABORT( status, GENERATEPULSARSIGNALTESTC_EMOD, GENERATEPULSARSIGNALTESTC_MSGEMOD);
          }
          /* 10/12/04 gam; turn on test above and add test below */
          if ( fabs(((REAL8)(jMaxMod - jFastMaxMod))) >  1.1 ) {
            ABORT( status, GENERATEPULSARSIGNALTESTC_EBIN,   GENERATEPULSARSIGNALTESTC_MSGEBIN);
          }
       } /* END for(i = 0; i < numSFTs; i++) */
       #ifdef PRINT_OVERALLMAXDIFFSFTPOWER
         fprintf(stdout,"overallMaxDiffSFTMod, testNumber = %i, SFT %i, bin %i = %g \n",testNumber,iOverallMaxDiffSFTMod,jOverallMaxDiff,overallMaxDiffSFTMod);
         fflush(stdout);
       #endif

       /* 09/07/05 gam; Initialize fastOutputSFTs since only 2*Dterms bins are changed by LALFastGeneratePulsarSFTs */
       for (i = 0; i < numSFTs; i++) {
          for(j=0;j<nBinsSFT;j++) {
             fastOutputSFTs->data[i].data->data[j].re = 0.0;
             fastOutputSFTs->data[i].data->data[j].im = 0.0;
          }
       }

       LALDestroySFTVector(status->statusPtr, &outputSFTs);
       CHECKSTATUSPTR (status);

       LALFree(signalvec->data->data);
       LALFree(signalvec->data);
       LALFree(signalvec);

     } /* END for(iFreq=0;iFreq<nBinsSGNL;iFreq++) */
     /****************************************************/
     /*                                                  */
     /* END SECTION: LOOP OVER FREQUENCIES               */
     /*                                                  */
     /****************************************************/
   } /* END for(jDeriv=0;jDeriv<numFreqDerivIncludingNoSpinDown;jDeriv++) */
  /*********************************************************/
  /*                                                       */
  /* END SECTION: LOOP OVER SPINDOWN                       */
  /*                                                       */
  /*********************************************************/

  } /* END for(iSky=0;iSky<numSkyPosTotal;iSky++) */
  /*********************************************************/
  /*                                                       */
  /* END SECTION: LOOP OVER SKY POSITIONS                  */
  /*                                                       */
  /*********************************************************/

  /* 10/12/04 gam; check if too many bin errors */
  epsBinErrorRate = 0.20;  /* 10/12/04 gam; maximum allowed bin errors */
  if ( (((REAL4)binErrorCount)/((REAL4)testNumber)) > epsBinErrorRate ) {
            ABORT( status, GENERATEPULSARSIGNALTESTC_EBINS, GENERATEPULSARSIGNALTESTC_MSGEBINS);
  }

  #ifdef INCLUDE_RANDVAL_MISMATCH
    LALDestroyRandomParams(status->statusPtr, &randPar);
    CHECKSTATUSPTR (status);
  #endif

  /* fprintf(stdout,"Total number of tests completed = %i. \n", testNumber);
  fflush(stdout); */

  LALFree(pSFTParams);
  if (numSpinDown > 0) {
    LALDDestroyVector(status->statusPtr, &(pPulsarSignalParams->pulsar.spindown));
    CHECKSTATUSPTR (status);
  }
  LALFree(pPulsarSignalParams);

  /* deallocate memory for structs needed by LALComputeSkyAndZeroPsiAMResponse and LALFastGeneratePulsarSFTs */
  LALDestroySFTVector(status->statusPtr, &fastOutputSFTs);
  CHECKSTATUSPTR (status);
  LALFree(pSkyConstAndZeroPsiAMResponse->fCrossZeroPsi);
  LALFree(pSkyConstAndZeroPsiAMResponse->fPlusZeroPsi);
  LALFree(pSkyConstAndZeroPsiAMResponse->skyConst);
  LALFree(pSkyConstAndZeroPsiAMResponse);
  /* 02/02/05 gam; if NOT pSFTandSignalParams->resTrig > 0 should not create trigArg etc... */
  if (pSFTandSignalParams->resTrig > 0) {
    LALFree(pSFTandSignalParams->trigArg);
    LALFree(pSFTandSignalParams->sinVal);
    LALFree(pSFTandSignalParams->cosVal);
  }
  LALFree(pSFTandSignalParams);

  /* deallocate skyPosData */
  for(i=0;i<numSkyPosTotal;i++)
  {
      LALFree(skyPosData[i]);
  }
  LALFree(skyPosData);

  if (numSpinDown > 0) {
    /* deallocate freqDerivData */
    for(i=0;i<numFreqDerivTotal;i++)
    {
        LALFree(freqDerivData[i]);
    }
    LALFree(freqDerivData);
  }

  LALFree(timeStamps->data);
  LALFree(timeStamps);

  LALFree ( edat->ephemE );
  LALFree ( edat->ephemS );
  LALFree ( edat);

  CHECKSTATUSPTR (status);
  DETATCHSTATUSPTR (status);
}
/*********************************************************/
/*                                                       */
/* END FUNCTION: RunGeneratePulsarSignalTest             */
/*                                                       */
/*********************************************************/
