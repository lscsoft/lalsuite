/*
*  Copyright (C) 2007 Tania Regimbau
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

/*
 * stochastic_preprocess.c 
 *
 * Tania Regimbau <regimbau@obs-nice.fr>  
 *
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#include <unistd.h>
#include <getopt.h>

#include <lal/LALFrameL.h>

#include <lal/AVFactories.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/ComplexFFT.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/IIRFilter.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/PrintFTSeries.h>
#include <lal/PrintVector.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/SimulateSB.h>
#include <lal/StreamInput.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/VectorOps.h>
#include <lalapps.h>

#include "stochastic_virgo.h"


/* cvs info */
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_stochastic_preprocess"

/* variables for getopt options parsing */
extern char *optarg;
extern int optind;

/* flags for getopt_long */
static int double_flag = 0;
static int double_high_pass_flag = 0;
static int recenter_flag = 0;
static int apply_mask_flag = 0;
static int high_pass_flag = 0;
static int overlap_hann_flag = 0;
static int verbose_flag = 0;
static int test_flag = 0;


/* parameters for the stochastic search */
UINT8 startTime = 700000100;
UINT8 stopTime =  700000300;

CHAR frameCache1 [200]= "H1.cache";
CHAR frameCache2[200] = "H2.cache";

CHAR channel1[LALNameLength]= "H1:STRAIN";
CHAR channel2[LALNameLength]= "H2/STRAIN";
CHAR ifo1[LALNameLength] = "H1";
CHAR ifo2[LALNameLength] = "H2";
INT4 site1 = 1;
INT4 site2 = 1;

/* data */
INT4 segmentDuration = 60;
INT4 numSegments = 3;
INT4 sampleRate1 = 16384.;
INT4 sampleRate2 = 16384.;
INT4 resampleRate1 = 1024.;
INT4 resampleRate2 = 1024.;
INT4 padData=1;
REAL4 calibFactor=1.e18;

/* frequency band */
REAL8 deltaF = 0.25;
INT4 fMin = 50;
INT4 fMax = 500;

/* window parameters */
/* 60 s for pure Hann, 1 s for Tukey, 0s for rectangular window */
INT4 hannDuration = 1;

/* high pass filtering parameters */
REAL4 highPassFreq = 40.;
REAL4 highPassAt ;
INT4  highPassOrder = 6;

/* double to single precision high pass filter */
REAL4 DhighPassFreq = 50.;
INT4  DhighPassOrder = 8;
REAL4 DhighPassAt = 0.9;

/* number of bins for frequency masking */ 
INT4 maskBin = 1;

/* arguments associated with test flag */
INT4 testInter = 0;
INT4 testSeg = 0;


INT4 main(INT4 argc, CHAR *argv[])
 {

  /* variable declarations */

  /* status pointer */
  LALStatus status;

  /* counters */
  INT4 i;
  INT4 lInter, lSeg;
  
  /* frame variables */
  LALCache *frCache1,*frCache2;
  LALFrStream *frStream1,*frStream2;
  FrChanIn frChanIn1, frChanIn2;
  FrOutPar frPSD1 = { ifo1, "PSD", ProcDataChannel, 1, 0, 0 };
  FrOutPar frPSD2 = { ifo2, "PSD", ProcDataChannel, 1, 0, 0 };
  FrOutPar frCC12 = { ifo1, "CC", ProcDataChannel, 1, 0, 0 };

  /* input data segment */
  LIGOTimeGPS gpsStartTime,gpsStartPadTime,gpsMidSegTime;
  INT4 duration, durationEff, extrasec; /*intervalDuration;*/
  INT4 numIntervals, midSegment;
  INT4 segmentLength1,segmentLength2;
  INT4 segmentPadDuration,segmentPadLength1,segmentPadLength2;
#if 0
  INT4 intervalLength,intervalLength1,intervalLength2;
  INT4 segmentShift;
#endif
  ResampleTSParams resampleParams1, resampleParams2; 
  REAL8TimeSeries segmentPadD1,segmentPadD2;
  REAL4TimeSeries segmentPad1,segmentPad2,segment1,segment2;
  REAL4Vector *seg1[numSegments],*seg2[numSegments];
  /*REAL4Vector *segBuffer1, *segBuffer2;*/
  LALUnit dataUnit = {-18,{0,0,0,0,0,1,0},{0,0,0,0,0,0,0}};

  /* high pass filtering */
  PassBandParamStruc highpassParam,DhighpassParam;

  /* window for segment data streams */
  REAL4Window *dataWindow1, *dataWindow2;

  /* hann window */
  REAL4Window *hannWindow1 = NULL, *hannWindow2 = NULL;
  
  /* data structures for PSDs */
  INT4 overlapPSDLength1,overlapPSDLength2;
  INT4 PSDTempLength1,PSDTempLength2;
  INT4 windowPSDLength1,windowPSDLength2;
  INT4 filterLength;
  INT4 numPointInf, numPointSup;
  AverageSpectrumParams specparPSD1,specparPSD2;
  REAL4FrequencySeries PSDTemp1,PSDTemp2,avPSD1,avPSD2;
  REAL4Vector *PSD1[numSegments], *PSD2[numSegments];
  LALUnit PSDUnits = {-36,{0,0,1,0,0,2,0},{0,0,0,0,0,0,0}};

  /* zeropad and fft structures */
  SZeroPadAndFFTParameters zeroPadParams1, zeroPadParams2;
  RealFFTPlan *fftDataPlan1 = NULL;
  RealFFTPlan *fftDataPlan2 = NULL;
  COMPLEX8FrequencySeries hBarTildeTemp1,hBarTildeTemp2,hBarTilde1,hBarTilde2;
  COMPLEX8FrequencySeries h1StarH2Temp,h1StarH2;
  FrequencySamplingParams freqParams;
  UINT4 zeroPadLength1, zeroPadLength2;
  UINT4 fftDataLength1, fftDataLength2;
  UINT4 hBarLength;
  LALUnit hBarUnit = {-18,{0,0,1,0,0,1,0},{0,0,0,0,0,0,0}};  
  LALUnit h2BarUnit = {-36,{0,0,2,0,0,2,0},{0,0,0,0,0,0,0}};
  
  /* error handler */
  status.statusPtr = NULL;

  lal_errhandler = LAL_ERR_EXIT;

  /* parse command line options */
  parseOptions(argc, argv);

  /* get durations etc...*/
  duration = stopTime - startTime;
#if 0
  intervalDuration = numSegments * segmentDuration;
#endif
  midSegment = (INT4)((numSegments-1)/2); 
  numIntervals = (INT4)((duration-2*padData)/segmentDuration)-numSegments+1;

#if 0
  if (overlap_hann_flag)
   segmentShift = segmentDuration / 2;
  else
   segmentShift = segmentDuration;
#endif

  /* recenter */
  if (recenter_flag)
   {
    durationEff = (INT4)((duration-2*padData)/segmentDuration)*segmentDuration;
    extrasec = duration - durationEff;
    startTime = startTime + (INT4)(extrasec/2);
    stopTime = startTime + durationEff; 
   }

  if (verbose_flag)
    {fprintf(stdout, "%d\n",numIntervals);}


  /* set channels */
  frChanIn1.name = channel1;
  frChanIn2.name = channel2;
  frChanIn1.type = SimDataChannel;
  frChanIn2.type = SimDataChannel;

 
  if (verbose_flag)
   fprintf(stdout, "Opening first frame cache...\n");
  /* open first frame cache */
  frCache1=NULL;frStream1=NULL;
  frCache1 = XLALCacheImport( frameCache1 );
  LAL_CALL( LALFrCacheOpen( &status, &frStream1, frCache1), &status );

  if (verbose_flag)
   fprintf(stdout, "Opening second frame cache...\n");

  /* open second frame cache */
  frCache2=NULL;frStream2=NULL;
  frCache2 = XLALCacheImport( frameCache2 );
  LAL_CALL( LALFrCacheOpen( &status, &frStream2, frCache2), &status);

   /* set the mode of the frame stream to fail on gaps or time errors */
  frStream1->mode = LAL_FR_STREAM_VERBOSE_MODE;
  frStream2->mode = LAL_FR_STREAM_VERBOSE_MODE;

  
  /* set resample parameters */
  if (resampleRate1 != sampleRate1)
   {
    resampleParams1.deltaT = 1.0 / (REAL8)resampleRate1;
    resampleParams1.filterType = defaultButterworth;
   }
  if (resampleRate2 != sampleRate2)
   {
    resampleParams2.deltaT = 1.0 / (REAL8)resampleRate2;
    resampleParams2.filterType = defaultButterworth;
   }

  /* initialize gps time structure */  
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0.;
  gpsStartPadTime.gpsSeconds = startTime - padData;
  gpsStartPadTime.gpsNanoSeconds = 0.;
  gpsMidSegTime.gpsSeconds = startTime + segmentDuration;
  gpsMidSegTime.gpsNanoSeconds = 0.;

  /* set length for data segments */
#if 0
  intervalLength = (intervalDuration + 2 * padData);
  intervalLength1 = (UINT4)(intervalLength * resampleRate1);
  intervalLength2 = (UINT4)(intervalLength * resampleRate2);
#endif
  segmentPadDuration = (segmentDuration + 2 * padData);
  segmentPadLength1 = (UINT4)(segmentPadDuration * sampleRate1);
  segmentPadLength2 = (UINT4)(segmentPadDuration * sampleRate2);
  segmentLength1 = (UINT4)(segmentDuration * resampleRate1);
  segmentLength2 = (UINT4)(segmentDuration * resampleRate2);

  /* set metadata fields for data segments */ 
  strncpy(segmentPad1.name, "segmentPad1", LALNameLength);
  strncpy(segmentPad2.name, "segmentPad2", LALNameLength);
  segmentPad1.sampleUnits = segmentPad2.sampleUnits = lalStrainUnit; 
  segmentPad1.epoch = segmentPad2.epoch = gpsStartPadTime;
  segmentPad1.deltaT = 1./(REAL8)sampleRate1;
  segmentPad2.deltaT = 1./(REAL8)sampleRate2;
  segmentPad1.f0 = segmentPad2.f0 = 0;

  if(double_flag)
   {
    strncpy(segmentPadD1.name, "segmentPadD1", LALNameLength);
    strncpy(segmentPadD2.name, "segmentPadD2", LALNameLength);
    segmentPadD1.sampleUnits = segmentPadD2.sampleUnits = lalStrainUnit; 
    segmentPadD1.epoch = segmentPadD2.epoch = gpsStartPadTime;
    segmentPadD1.deltaT = 1./(REAL8)sampleRate1;
    segmentPadD2.deltaT = 1./(REAL8)sampleRate2;
    segmentPadD1.f0 = segmentPadD2.f0 = 0;
   }

  strncpy(segment1.name, "segment1", LALNameLength);
  strncpy(segment2.name, "segment2", LALNameLength);
  segment1.sampleUnits = segment2.sampleUnits = dataUnit;
  segment1.epoch = segment2.epoch = gpsStartTime;
  segment1.deltaT = 1./(REAL8)resampleRate1;
  segment2.deltaT = 1./(REAL8)resampleRate2;
  segment1.f0 = segment2.f0 = 0;

  if (verbose_flag)
   fprintf(stdout, "Allocating memory for data segments...\n");

  /* allocate memory for data segments */

  segment1.data = segment2.data = NULL;
  LAL_CALL(LALSCreateVector(&status,&(segment1.data),segmentLength1), 
           &status);
  LAL_CALL( LALSCreateVector(&status,&(segment2.data),segmentLength2), 
            &status ); 
  memset(segment1.data->data, 0,
         segment1.data->length * sizeof(*segment1.data->data));
  memset(segment2.data->data, 0,
         segment2.data->length * sizeof(*segment2.data->data));
   

  for (i=0;i<numSegments;i++)
   {
    seg1[i]= seg2[i] = NULL;
    LAL_CALL(LALCreateVector(&status,&(seg1[i]),segmentLength1), 
             &status);
    LAL_CALL(LALCreateVector(&status,&(seg2[i]),segmentLength2), 
             &status);
    memset(seg1[i]->data, 0,seg1[i]->length * sizeof(*seg1[i]->data));
    memset(seg2[i]->data, 0,seg2[i]->length * sizeof(*seg2[i]->data));
   }
 

  /* structure for high pass filtering */
  if (high_pass_flag)
   {
     highpassParam.nMax = highPassOrder;
     highpassParam.f1 = -1.0;
     highpassParam.f2 = highPassFreq;
     highpassParam.a1 = -1.0;
     highpassParam.a2 = highPassAt;
   }

  if (double_high_pass_flag)
   {
     DhighpassParam.nMax = DhighPassOrder;
     DhighpassParam.f1 = -1.0;
     DhighpassParam.f2 = DhighPassFreq;
     DhighpassParam.a1 = -1.0;
     DhighpassParam.a2 = DhighPassAt;
   }

  /* generate windows */
  if (verbose_flag)
   { fprintf(stdout, "Generating data segment windows...\n");}

  if (overlap_hann_flag)
   { hannDuration = segmentDuration;}

  dataWindow1 = XLALCreateTukeyREAL4Window(segmentLength1, hannDuration * resampleRate1 / segmentLength1);
  dataWindow2 = XLALCreateTukeyREAL4Window(segmentLength2, hannDuration * resampleRate2 / segmentLength2);

  /* PSDs */
  /* set parameters for PSD estimation */
  windowPSDLength1 = (UINT4)(resampleRate1 / deltaF);
  windowPSDLength2 = (UINT4)(resampleRate2 / deltaF);
  overlapPSDLength1 = windowPSDLength1 / 2;
  overlapPSDLength2 = windowPSDLength2 / 2;
  PSDTempLength1 = (windowPSDLength1 / 2) + 1;
  PSDTempLength2 = (windowPSDLength2 / 2) + 1;
  numPointInf = (UINT4)(fMin / deltaF);
  numPointSup = (UINT4)(fMax / deltaF);
  filterLength = numPointSup - numPointInf + 1;


  /* create fft plan */
  if (verbose_flag)
   fprintf(stdout, "Creating FFT plan for PSD estimation...\n");

  specparPSD1.method = specparPSD2.method = useMean;
  specparPSD1.plan = specparPSD2.plan = NULL;
  specparPSD1.window =  specparPSD2.window = NULL;
  specparPSD1.overlap = overlapPSDLength1;
  specparPSD2.overlap = overlapPSDLength2;
  LAL_CALL(LALCreateForwardRealFFTPlan(&status,&specparPSD1.plan,
             windowPSDLength1,0), &status );
  LAL_CALL(LALCreateForwardRealFFTPlan(&status,&specparPSD2.plan,
             windowPSDLength2,0), &status );


  /* create window for PSD estimation */

  if (verbose_flag)
   fprintf(stdout, "Creating window for PSD estimation...\n");
  specparPSD1.window = XLALCreateHannREAL4Window(windowPSDLength1);
  specparPSD2.window = XLALCreateHannREAL4Window(windowPSDLength2);

  /* set metadata fields for PSDs */
  strncpy(PSDTemp1.name, "PSDTemp1", LALNameLength);
  strncpy(PSDTemp2.name, "PSDTemp2", LALNameLength);
  PSDTemp1.sampleUnits = PSDTemp2.sampleUnits = PSDUnits;
  PSDTemp1.epoch = PSDTemp2.epoch = gpsStartTime;
  PSDTemp1.deltaF = PSDTemp2.deltaF = deltaF;
  PSDTemp1.f0 = PSDTemp2.f0 = 0;

  if (verbose_flag)
   { fprintf(stdout, "Allocating memory for PSDs...\n");}

  /* allocate memory for PSDs */
  PSDTemp1.data = PSDTemp2.data = NULL;
  LAL_CALL(LALCreateVector(&status,&(PSDTemp1.data),PSDTempLength1),
	    &status );
  LAL_CALL(LALCreateVector(&status,&(PSDTemp2.data),PSDTempLength2), 
	    &status );
  memset(PSDTemp1.data->data, 0, 
         PSDTemp1.data->length * sizeof(*PSDTemp1.data->data));
  memset(PSDTemp2.data->data, 0, 
         PSDTemp2.data->length * sizeof(*PSDTemp2.data->data));

  /* reduced frequency band average PSDs */
  /* set metadata fields for reduced frequency band average PSDs */
  strncpy(avPSD1.name, "avPSD1", LALNameLength);
  strncpy(avPSD2.name, "avPSD2", LALNameLength);
  avPSD1.sampleUnits = avPSD2.sampleUnits = PSDUnits;
  avPSD1.epoch = avPSD2.epoch = gpsMidSegTime; 
  avPSD1.deltaF = avPSD2.deltaF = deltaF;
  avPSD1.f0 = avPSD2.f0 = fMin;

  if (verbose_flag)
   fprintf(stdout, "Allocating memory for reduced frequency band PSDs...\n");

  /* allocate memory for reduced frequency band average PSDs */
  avPSD1.data = avPSD2.data = NULL;
  LAL_CALL(LALCreateVector(&status,&avPSD1.data, filterLength), &status);
  LAL_CALL(LALCreateVector(&status,&avPSD2.data, filterLength), &status);
  memset(avPSD1.data->data,0,avPSD1.data->length * sizeof(*avPSD1.data->data));
  memset(avPSD2.data->data,0,avPSD2.data->length * sizeof(*avPSD2.data->data));

  for (i=0;i<numSegments;i++)
   {
    PSD1[i]= PSD2[i] = NULL;
    LAL_CALL(LALCreateVector(&status,&(PSD1[i]),filterLength), 
             &status);
    LAL_CALL(LALCreateVector(&status,&(PSD2[i]),filterLength), 
             &status);
    memset(PSD1[i]->data, 0, PSD1[i]->length * sizeof(*PSD1[i]->data));
    memset(PSD2[i]->data, 0, PSD2[i]->length * sizeof(*PSD2[i]->data));
   }


  /* zeropad and fft */
  /* zeropad lengths */
  zeroPadLength1 = 2 * segmentLength1;
  zeroPadLength2 = 2 * segmentLength2;
  fftDataLength1 = segmentLength1 + 1;
  fftDataLength2 = segmentLength2 + 1;
  
  if(fftDataLength1<=fftDataLength2)
   hBarLength  = fftDataLength1;
  else hBarLength  = fftDataLength2;

  freqParams.length = filterLength;
  freqParams.f0 = fMin;
  freqParams.deltaF = deltaF;
  
  
  /* create fft plan */
  LAL_CALL(LALCreateForwardRealFFTPlan(&status,&fftDataPlan1,zeroPadLength1,0),&status);
  LAL_CALL(LALCreateForwardRealFFTPlan(&status,&fftDataPlan2,zeroPadLength2,0),&status);

  /* set metadata fields for zeropad ffts */
  strncpy(hBarTildeTemp1.name, "hBarTildeTemp1", LALNameLength);
  strncpy(hBarTildeTemp2.name, "hBarTildeTemp2", LALNameLength);
  hBarTildeTemp1.epoch = hBarTildeTemp2.epoch = gpsMidSegTime;

  hBarTildeTemp1.f0 = hBarTildeTemp2.f0 = 0.;
  hBarTildeTemp1.deltaF = hBarTildeTemp2.deltaF = 1./(2.*(REAL8)segmentDuration); 
  hBarTildeTemp1.sampleUnits = hBarTildeTemp2.sampleUnits= hBarUnit; 
 

 /* allocate memory for zeropad */
 if (verbose_flag)
  fprintf(stdout, "Allocating memory for zeropad...\n");

  hBarTildeTemp1.data = hBarTildeTemp2.data = NULL;
  LAL_CALL( LALCCreateVector(&status, &(hBarTildeTemp1.data), fftDataLength1), 
            &status );
  LAL_CALL( LALCCreateVector(&status, &(hBarTildeTemp2.data), fftDataLength2), 
            &status );
  memset( hBarTildeTemp1.data->data, 0, 
          hBarTildeTemp1.data->length * sizeof(*hBarTildeTemp1.data->data));
  memset( hBarTildeTemp2.data->data, 0, 
          hBarTildeTemp2.data->length * sizeof(*hBarTildeTemp2.data->data));


  strncpy(hBarTilde1.name, "hBarTilde1", LALNameLength);
  strncpy(hBarTilde2.name, "hBarTilde2", LALNameLength);
  hBarTilde1.epoch = hBarTilde2.epoch = hBarTildeTemp1.epoch;
  hBarTilde1.f0 = hBarTilde2.f0  = 0.;
  hBarTilde1.deltaF = hBarTilde2.deltaF = hBarTildeTemp1.deltaF; 
  hBarTilde1.sampleUnits = hBarTilde2.sampleUnits = hBarTildeTemp1.sampleUnits; 

  hBarTilde1.data = hBarTilde2.data =  NULL;
  LAL_CALL( LALCCreateVector(&status, &(hBarTilde1.data), hBarLength), 
            &status );
  LAL_CALL( LALCCreateVector(&status, &(hBarTilde2.data), hBarLength), 
            &status );

  memset( hBarTilde1.data->data, 0, 
          hBarTilde1.data->length * sizeof(*hBarTilde1.data->data));
  memset( hBarTilde2.data->data, 0, 
          hBarTilde2.data->length * sizeof(*hBarTilde2.data->data)); 
  
  strncpy(h1StarH2.name, "h1StarH2Temp", LALNameLength);
  strncpy(h1StarH2.name, "h1StarH2", LALNameLength);
  h1StarH2Temp.epoch = h1StarH2.epoch = hBarTilde1.epoch;
  h1StarH2Temp.f0 = 0.; h1StarH2.f0 =fMin;
  h1StarH2Temp.deltaF = hBarTilde1.deltaF; h1StarH2.deltaF = deltaF; 
  h1StarH2Temp.sampleUnits = h1StarH2.sampleUnits = h2BarUnit; 
  h1StarH2Temp.data = h1StarH2.data = NULL;
  LAL_CALL( LALCCreateVector(&status, &(h1StarH2Temp.data), hBarLength), 
            &status );
  LAL_CALL( LALCCreateVector(&status, &(h1StarH2.data), filterLength), 
            &status );
  memset( h1StarH2Temp.data->data, 0, 
          h1StarH2Temp.data->length * sizeof(*h1StarH2Temp.data->data));
  memset( h1StarH2.data->data, 0, 
          h1StarH2.data->length * sizeof(*h1StarH2.data->data));

  /* set zeropad parameters */
  zeroPadParams1.fftPlan = fftDataPlan1;
  zeroPadParams2.fftPlan = fftDataPlan2;
  zeroPadParams1.window = dataWindow1;
  zeroPadParams2.window = dataWindow2;
  zeroPadParams1.length = zeroPadLength1;
  zeroPadParams2.length = zeroPadLength2;
 	   
  
  /* START HERE */


  /* loop over interval */
 
  for (lInter=0;lInter<numIntervals;lInter++)
   { 
    
    /* initialise average PSD vector */
    for (i=0;i<filterLength;i++)
     {
      avPSD1.data->data[i]=0.;
      avPSD2.data->data[i]=0.;
     }

    /* firstpass: get first interval */
    if(lInter==0)
     { 
      for (lSeg=0;lSeg<numSegments;lSeg++)
       {
              
        /* read data */ 

        /* set metadata fields for data segments */ 
        segmentPad1.epoch = segmentPad2.epoch = gpsStartPadTime;
        segmentPad1.deltaT = 1./(REAL8)sampleRate1;
        segmentPad2.deltaT = 1./(REAL8)sampleRate2;   

        /* allocate memory for data segments */
        segmentPad1.data = segmentPad2.data = NULL;
        LAL_CALL(LALSCreateVector(&status,&(segmentPad1.data),segmentPadLength1),&status );
        LAL_CALL(LALSCreateVector(&status,&(segmentPad2.data),segmentPadLength2), &status );
        memset(segmentPad1.data->data, 0,
               segmentPad1.data->length * sizeof(*segmentPad1.data->data));
        memset(segmentPad2.data->data, 0,
               segmentPad2.data->length * sizeof(*segmentPad2.data->data));

	if (double_flag)
	 {
          /* set metadata fields for data segments */ 
          segmentPadD1.epoch = segmentPadD2.epoch = gpsStartPadTime;
          segmentPadD1.deltaT = 1./(REAL8)sampleRate1;
          segmentPadD2.deltaT = 1./(REAL8)sampleRate2;   

          /* allocate memory for data segments */
          segmentPadD1.data = segmentPadD2.data = NULL;
          LAL_CALL(LALDCreateVector(&status,&(segmentPadD1.data),segmentPadLength1),&status );
          LAL_CALL(LALDCreateVector(&status,&(segmentPadD2.data),segmentPadLength2), &status );
          memset(segmentPadD1.data->data, 0,
                 segmentPadD1.data->length * sizeof(*segmentPadD1.data->data));
          memset(segmentPadD2.data->data, 0,
                 segmentPadD2.data->length * sizeof(*segmentPadD2.data->data));
          
          if (verbose_flag)
           fprintf(stdout, "request data at GPS time %d\n",
                   gpsStartPadTime.gpsSeconds);

          /* read first channel */	
          if (verbose_flag)
           fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);
          LAL_CALL(LALFrSeek( &status, &(gpsStartPadTime), frStream1), &status );
          LAL_CALL(LALFrGetREAL8TimeSeries(&status,&segmentPadD1,&frChanIn1,frStream1),&status );
        
          /* read second channel */	
          if (verbose_flag)
           fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);
          LAL_CALL(LALFrSeek(&status,&gpsStartPadTime,frStream2),&status);
          LAL_CALL(LALFrGetREAL8TimeSeries(&status,&segmentPadD2,&frChanIn2,frStream2),&status);
       
         /* print */
	  if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
          {
           LALDPrintTimeSeries(&segmentPadD1, "segmentPadD1.dat");
  	   LALDPrintTimeSeries(&segmentPadD2, "segmentPadD2.dat");
          }

          if(double_high_pass_flag)
	   {

            LAL_CALL(LALButterworthREAL8TimeSeries(&status,&segmentPadD1, 
	 					 &DhighpassParam),&status);
            LAL_CALL(LALButterworthREAL8TimeSeries(&status,&segmentPadD2, 
                                                 &DhighpassParam),&status);
	   }
      
          /* cast to REAL4  */
          for (i=0;i<(INT4)segmentPadD1.data->length;i++)           
           segmentPad1.data->data[i]= (REAL4)segmentPadD1.data->data[i];
          for (i=0;i<(INT4)segmentPadD2.data->length;i++)           
           segmentPad2.data->data[i]= (REAL4)segmentPadD2.data->data[i];

          /* clean up */
          LAL_CALL(LALDDestroyVector(&status, &(segmentPadD1.data)),&status );
          LAL_CALL(LALDDestroyVector(&status, &(segmentPadD2.data)),&status );

         }
        else
	 {
          if (verbose_flag)
           fprintf(stdout, "request data at GPS time %d\n",
                   gpsStartPadTime.gpsSeconds);
          /* read first channel */	
          if (verbose_flag)
           fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);
          LAL_CALL(LALFrSeek( &status, &(gpsStartPadTime), frStream1), &status );
          LAL_CALL(LALFrGetREAL4TimeSeries(&status,&segmentPad1,&frChanIn1,frStream1),&status );
        
          /* read second channel */	
          if (verbose_flag)
           fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);
          LAL_CALL(LALFrSeek(&status,&gpsStartPadTime,frStream2),&status);
          LAL_CALL(LALFrGetREAL4TimeSeries(&status,&segmentPad2,&frChanIn2,frStream2),&status);
       
         /* print */
	 if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
          {
           LALSPrintTimeSeries(&segmentPad1, "segmentPad1.dat");
  	   LALSPrintTimeSeries(&segmentPad2, "segmentPad2.dat");
          }
	 }
        /* resample */
 
        /* resample first data stream */
        if (resampleRate1 != sampleRate1)
         {
          if (verbose_flag)
           fprintf(stdout, "Resampling first data stream to %d Hz...\n",resampleRate1);
          LAL_CALL(LALResampleREAL4TimeSeries(&status,&segmentPad1,&resampleParams1),&status);	
         }

        /* resample second data stream */
        if (resampleRate2 != sampleRate2)
         {
          if (verbose_flag)
           fprintf(stdout, "Resampling second data stream to %d Hz...\n",resampleRate2);
         LAL_CALL(LALResampleREAL4TimeSeries(&status,&segmentPad2,&resampleParams2),&status);
         }
     
        /* high pass filtering */
        if (high_pass_flag)
         {               
          if (verbose_flag)
           fprintf(stdout, "High pass filter data streams...\n");
          LAL_CALL(LALButterworthREAL4TimeSeries(&status,&segmentPad1,&highpassParam),&status);
          LAL_CALL(LALButterworthREAL4TimeSeries(&status,&segmentPad2,&highpassParam),&status);
         } 
            
        /* throw away pad data */      
        segment1.epoch = segment2.epoch = gpsStartTime;
        for (i=0;i<segmentLength1;i++)
         segment1.data->data[i]=calibFactor*segmentPad1.data->data[i+padData*resampleRate1];
        for (i=0;i<segmentLength2;i++)
         segment2.data->data[i]=calibFactor*segmentPad2.data->data[i+padData*resampleRate2];
        /* clean up */
        LAL_CALL( LALDestroyVector(&status, &(segmentPad1.data)), &status );
        LAL_CALL( LALDestroyVector(&status, &(segmentPad2.data)), &status );
      
        /* print */
        if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
         {
          LALSPrintTimeSeries(&segment1, "segment1.dat");
  	  LALSPrintTimeSeries(&segment2, "segment2.dat");
         }
       
        /* store in memory */
        for (i=0;i<segmentLength1;i++)
         seg1[lSeg]->data[i] = segment1.data->data[i];
        for (i=0;i<segmentLength2;i++)
         seg2[lSeg]->data[i] = segment2.data->data[i];
       

        /* PSD estimation */
        if (verbose_flag)
         fprintf(stdout, "Estimating PSDs...\n");

        /* compute PSDs */
        PSDTemp1.epoch = PSDTemp2.epoch = gpsStartTime;
        LAL_CALL(LALREAL4AverageSpectrum(&status,&PSDTemp1,&segment1, 
                 &specparPSD1),&status);
        LAL_CALL(LALREAL4AverageSpectrum(&status,&PSDTemp2,&segment2, 
                 &specparPSD2),&status);
   
        /* print */
        if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
         {
          LALSPrintFrequencySeries(&PSDTemp1, "PSDTemp1.dat");
          LALSPrintFrequencySeries(&PSDTemp2, "PSDTemp2.dat");
         }

        /* reduce to the optimal filter frequency range and store*/
        if (verbose_flag)
         fprintf(stdout, "Reduce to optimal filter range...\n");
        
        for (i=0;i<filterLength;i++)
         {
	  PSD1[lSeg]->data[i] =  PSDTemp1.data->data[i + numPointInf];
	  PSD2[lSeg]->data[i] =  PSDTemp2.data->data[i + numPointInf];
         }
     

        /* compute average PSD, excluding middle segment*/
        if (verbose_flag)
         fprintf(stdout, "sum PSDs...\n");
       
        if(lSeg==midSegment)
	 { 
          if (verbose_flag)
	    fprintf(stdout, "skip middle segment...\n");}
         else
	   {
           for (i=0;i<filterLength;i++)
            {
             avPSD1.data->data[i]=avPSD1.data->data[i]+PSD1[lSeg]->data[i];
             avPSD2.data->data[i]=avPSD2.data->data[i]+PSD2[lSeg]->data[i];
	    }
           
	  }
        /* increment epoch */
        gpsStartPadTime.gpsSeconds=gpsStartPadTime.gpsSeconds+segmentDuration;
        gpsStartTime.gpsSeconds=gpsStartTime.gpsSeconds+segmentDuration;  
       }
      if (verbose_flag)
        fprintf(stdout, "Compute average PSD...\n");

      avPSD1.epoch = avPSD2.epoch = gpsMidSegTime;
      for (i=0;i<filterLength;i++)
       {
        avPSD1.data->data[i]=avPSD1.data->data[i]/(numSegments-1);
        avPSD2.data->data[i]=avPSD2.data->data[i]/(numSegments-1);
       }

      /* print */
      if ((test_flag)&&(lInter==testInter))
       {
        LALSPrintFrequencySeries(&avPSD1, "avPSD1.dat");
        LALSPrintFrequencySeries(&avPSD2, "avPSD2.dat");
       }
       }
    else
     {
      /* shift segments */
      if (verbose_flag)
       fprintf(stdout, "Shift Segment\n interval %d out of %d\n",lInter, numIntervals);
    
      for (lSeg=0;lSeg<numSegments-1;lSeg++)
       {
        for (i=0;i<segmentLength1;i++)
	 seg1[lSeg]->data[i] = seg1[lSeg+1]->data[i];
        for (i=0;i<segmentLength2;i++)
	 seg2[lSeg]->data[i] = seg2[lSeg+1]->data[i];
       }
 
      /* read extra segment */

      /* read data */ 

      /* set metadata fields for data segments */ 
      segmentPad1.epoch = segmentPad2.epoch = gpsStartPadTime;
      segmentPad1.deltaT = 1./(REAL8)sampleRate1;
      segmentPad2.deltaT = 1./(REAL8)sampleRate2;   

      /* allocate memory for data segments */
      segmentPad1.data = segmentPad2.data = NULL;
      LAL_CALL(LALSCreateVector(&status,&(segmentPad1.data),segmentPadLength1),&status );
      LAL_CALL(LALSCreateVector(&status,&(segmentPad2.data),segmentPadLength2), &status );
      memset(segmentPad1.data->data, 0,
             segmentPad1.data->length * sizeof(*segmentPad1.data->data));
      memset(segmentPad2.data->data, 0,
             segmentPad2.data->length * sizeof(*segmentPad2.data->data));

      if (double_flag)
       {
        /* set metadata fields for data segments */ 
        segmentPadD1.epoch = segmentPadD2.epoch = gpsStartPadTime;
        segmentPadD1.deltaT = 1./(REAL8)sampleRate1;
        segmentPadD2.deltaT = 1./(REAL8)sampleRate2;   

        /* allocate memory for data segments */
        segmentPadD1.data = segmentPadD2.data = NULL;
        LAL_CALL(LALDCreateVector(&status,&(segmentPadD1.data),segmentPadLength1),&status );
        LAL_CALL(LALDCreateVector(&status,&(segmentPadD2.data),segmentPadLength2), &status );
        memset(segmentPadD1.data->data, 0,
               segmentPadD1.data->length * sizeof(*segmentPadD1.data->data));
        memset(segmentPadD2.data->data, 0,
               segmentPadD2.data->length * sizeof(*segmentPadD2.data->data));
          
        if (verbose_flag)
         fprintf(stdout, "request data at GPS time %d\n",
                 gpsStartPadTime.gpsSeconds);
        /* read first channel */	
        if (verbose_flag)
         fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);
        LAL_CALL(LALFrSeek( &status, &(gpsStartPadTime), frStream1), &status );
        LAL_CALL(LALFrGetREAL8TimeSeries(&status,&segmentPadD1,&frChanIn1,frStream1),&status );
        
        /* read second channel */	
        if (verbose_flag)
         fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);
        LAL_CALL(LALFrSeek(&status,&gpsStartPadTime,frStream2),&status);
        LAL_CALL(LALFrGetREAL8TimeSeries(&status,&segmentPadD2,&frChanIn2,frStream2),&status);
       
        /* print */
        if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
         {
          LALDPrintTimeSeries(&segmentPadD1, "segmentPadD1.dat");
          LALDPrintTimeSeries(&segmentPadD2, "segmentPadD2.dat");
         }
         if(double_high_pass_flag)
	  {
           /* high pass the data  */                
           LAL_CALL(LALButterworthREAL8TimeSeries(&status,&segmentPadD1, 
	 					  &DhighpassParam),&status);
           LAL_CALL(LALButterworthREAL8TimeSeries(&status,&segmentPadD2, 
                                                  &DhighpassParam),&status);
	  }
        
        /* cast to REAL4  */
        for (i=0;i<(INT4)segmentPadD1.data->length;i++)           
         segmentPad1.data->data[i]= (REAL4)segmentPadD1.data->data[i];
        for (i=0;i<(INT4)segmentPadD2.data->length;i++)           
         segmentPad2.data->data[i]= (REAL4)segmentPadD2.data->data[i];

        /* clean up */
        LAL_CALL(LALDDestroyVector(&status,&(segmentPadD1.data)),&status );
        LAL_CALL(LALDDestroyVector(&status,&(segmentPadD2.data)),&status );	    
       }
      else
       {
        if (verbose_flag)
         fprintf(stdout, "request data at GPS time %d\n",
                 gpsStartPadTime.gpsSeconds);
        /* read first channel */	
        if (verbose_flag)
         fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);
        LAL_CALL(LALFrSeek( &status, &(gpsStartPadTime), frStream1), &status );
        LAL_CALL(LALFrGetREAL4TimeSeries(&status,&segmentPad1,&frChanIn1,frStream1),&status );
        
        /* read second channel */	
        if (verbose_flag)
         fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);
        LAL_CALL(LALFrSeek(&status,&gpsStartPadTime,frStream2),&status);
        LAL_CALL(LALFrGetREAL4TimeSeries(&status,&segmentPad2,&frChanIn2,frStream2),&status);
       
       /* print */
       if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
        {
         LALSPrintTimeSeries(&segmentPad1, "segmentPad1.dat");
         LALSPrintTimeSeries(&segmentPad2, "segmentPad2.dat");
        }
       }
      /* resample */
     
      /* resample first data stream */
      if (resampleRate1 != sampleRate1)
       {
        if (verbose_flag)
         fprintf(stdout, "Resampling first data stream to %d Hz...\n", resampleRate1);
       LAL_CALL(LALResampleREAL4TimeSeries(&status,&segmentPad1,&resampleParams1),&status);	
       }

      /* resample second data stream */
      if (resampleRate2 != sampleRate2)
       {
        if (verbose_flag)
         fprintf(stdout, "Resampling second data stream to %d Hz...\n",resampleRate2);
        LAL_CALL(LALResampleREAL4TimeSeries(&status,&segmentPad2,&resampleParams2),&status);
        }
    
      /* high pass filtering */
      if (high_pass_flag)
       {               
        if (verbose_flag)
         fprintf(stdout, "High pass filter data streams...\n");
        LAL_CALL(LALButterworthREAL4TimeSeries(&status,&segmentPad1,&highpassParam),&status);
        LAL_CALL(LALButterworthREAL4TimeSeries(&status,&segmentPad2,&highpassParam),&status);
       } 
            
      /* throw away pad data */   
      segment1.epoch = segment2.epoch = gpsStartTime;
      for (i=0;i<segmentLength1;i++)
       segment1.data->data[i]=calibFactor*segmentPad1.data->data[i+padData*resampleRate1];
      for (i=0;i<segmentLength2;i++)
       segment2.data->data[i]=calibFactor*segmentPad2.data->data[i+padData*resampleRate2];
  
      /*clean up */
      LAL_CALL(LALDestroyVector(&status,&(segmentPad1.data)),&status);
      LAL_CALL(LALDestroyVector(&status,&(segmentPad2.data)),&status);   

      /* print */
      if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
       {
        LALSPrintTimeSeries(&segment1, "segment1.dat");
	LALSPrintTimeSeries(&segment2, "segment2.dat");
       }

      /* store in memory */
      for (i=0;i<segmentLength1;i ++)
       seg1[numSegments-1]->data[i] = segment1.data->data[i];
      for (i=0;i<segmentLength2;i ++)
       seg2[numSegments-1]->data[i] = segment2.data->data[i]; 

   
    /* PSD estimation */

    /* shift segments */
      if (verbose_flag)
       fprintf(stdout, "Shift PSDs\n interval %d out of %d\n",lInter, numIntervals);
    
      for (lSeg=0;lSeg<numSegments-1;lSeg++)
       {
        for (i=0;i<filterLength;i++)
	 {
	  PSD1[lSeg]->data[i] = PSD1[lSeg+1]->data[i];
	  PSD2[lSeg]->data[i] = PSD2[lSeg+1]->data[i];
         }
       }

    if (verbose_flag)
     fprintf(stdout, "Estimating extra PSDs...\n");

    /* compute PSDs */
    PSDTemp1.epoch = PSDTemp2.epoch = gpsStartTime;
    LAL_CALL(LALREAL4AverageSpectrum(&status,&PSDTemp1,&segment1, 
             &specparPSD1),&status);
    LAL_CALL(LALREAL4AverageSpectrum(&status,&PSDTemp2,&segment2, 
             &specparPSD2),&status);
   
    /* print */
    if ((test_flag)&&(lInter==testInter)&&(lSeg==testSeg))
     {
      LALSPrintFrequencySeries(&PSDTemp1, "PSDTemp1.dat");
      LALSPrintFrequencySeries(&PSDTemp2, "PSDTemp2.dat");
     }

    /* reduce to the optimal filter frequency range and store*/
    if (verbose_flag)
         fprintf(stdout, "Reduce to optimal filter frequency range...\n");
    for (i=0;i<filterLength;i++)
     {
       PSD1[numSegments-1]->data[i] = PSDTemp1.data->data[i + numPointInf];
       PSD2[numSegments-1]->data[i] = PSDTemp2.data->data[i + numPointInf];
     }
    

    /* increment epoch */
    gpsStartPadTime.gpsSeconds=gpsStartPadTime.gpsSeconds+segmentDuration;
    gpsStartTime.gpsSeconds=gpsStartTime.gpsSeconds+segmentDuration; 

    /* compute average PSD, excluding middle segment*/
    if (verbose_flag)
     fprintf(stdout, "Sum PSDs over segments...\n");

    for (lSeg=0;lSeg<numSegments;lSeg++)
     {
      if(lSeg==midSegment)
	 { 
          if (verbose_flag)
	    fprintf(stdout, "skip middle segment...\n");}
         else
	  {
           for (i=0;i<filterLength;i++)
            {
             avPSD1.data->data[i]=avPSD1.data->data[i]+PSD1[lSeg]->data[i];
             avPSD2.data->data[i]=avPSD2.data->data[i]+PSD2[lSeg]->data[i];
	    }
	  }
     }

    if (verbose_flag)
     fprintf(stdout, "Compute average PSD...\n");
    avPSD1.epoch = avPSD2.epoch = gpsMidSegTime;
    for (i=0;i<filterLength;i++)
     {
      avPSD1.data->data[i]=avPSD1.data->data[i]/(REAL4)(numSegments-1);
      avPSD2.data->data[i]=avPSD2.data->data[i]/(REAL4)(numSegments-1);
     }

    /* print */
    if ((test_flag)&&(lInter==testInter))
     {
      LALSPrintFrequencySeries(&avPSD1, "avPSD1.dat");
      LALSPrintFrequencySeries(&avPSD2, "avPSD2.dat");
     }

     }
   
  /* analyse middle segment */             
  /*gpsStartTime.gpsSeconds = startTime+(lInter+midSegment)*segmentDuration;*/
  segment1.epoch = segment2.epoch = gpsMidSegTime;          
  if (verbose_flag) 
   fprintf(stdout, "analysing segment at GPS %d\n", gpsMidSegTime.gpsSeconds);
  for (i=0;i<segmentLength1;i++)
   segment1.data->data[i] = seg1[midSegment]->data[i];
  for (i=0;i<segmentLength2;i++)
   segment2.data->data[i] = seg2[midSegment]->data[i];
	     
  /* zero pad and fft */
  if (verbose_flag)  
   fprintf(stdout, "zeropadding and fft...\n");
  hBarTildeTemp1.epoch = hBarTildeTemp2.epoch = gpsMidSegTime;
  LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeTemp1, &segment1, 
            &zeroPadParams1), &status );
  LAL_CALL( LALSZeroPadAndFFT(&status, &hBarTildeTemp2, &segment2, 
            &zeroPadParams2), &status );   

  /* print */
  if ((test_flag)&&(lInter==testInter))
   {
    LALCPrintFrequencySeries(&hBarTildeTemp1, "hBarTildeTemp1.dat");
    LALCPrintFrequencySeries(&hBarTildeTemp2, "hBarTildeTemp2.dat");
   }
       
  /* get hBarTild1 and hbarTild2 of the same length */
  hBarTilde1.epoch = hBarTilde2.epoch = gpsMidSegTime;
  if (verbose_flag)
   fprintf(stdout, "get hBarTild1 and hbarTild2 of the same length...\n");
  for (i=0;i<(INT4)hBarLength;i++)
   {
    hBarTilde1.data->data[i] =  hBarTildeTemp1.data->data[i];
    hBarTilde2.data->data[i] =  hBarTildeTemp2.data->data[i];
   }  

  /* print */
  if ((test_flag)&&(lInter==testInter))
   {
    LALCPrintFrequencySeries(&hBarTilde1, "hBarTilde1.dat");
    LALCPrintFrequencySeries(&hBarTilde2, "hBarTilde2.dat");
   }
 
   /* compute the cross correlation product */
   if (verbose_flag)
   fprintf(stdout, "compute the cross correlation product...\n");
   LAL_CALL( LALCCVectorMultiplyConjugate(&status, h1StarH2Temp.data,
             hBarTilde2.data, hBarTilde1.data), &status );

   /* coarse grain to optimal filter frequency resolution and range */
   if (verbose_flag)
   fprintf(stdout, "coarse grain to optimal filter frequency resolution and range...\n");
   LAL_CALL( LALCCoarseGrainFrequencySeries(&status, &h1StarH2,
              &h1StarH2Temp, &freqParams), &status );

  /* print */
  if ((test_flag)&&(lInter==testInter))
   {
    LALCPrintFrequencySeries(&h1StarH2Temp, "h1StarH2Temp.dat");
    LALCPrintFrequencySeries(&h1StarH2, "h1StarH2.dat");
   }

   /* write into frames */
   if (verbose_flag)
    fprintf(stdout, "write average PSDs into frames\n"); 
   LALFrWriteREAL4FrequencySeries( &status, &avPSD1, &frPSD1,0);
   LALFrWriteREAL4FrequencySeries( &status, &avPSD2, &frPSD2,0);
   
    /* write into frames */
    if (verbose_flag)
     fprintf(stdout, "write cross correlation product into frame\n"); 
    LALFrWriteCOMPLEX8FrequencySeries( &status, &h1StarH2, &frCC12,0);

   /* increment epoch */
  gpsMidSegTime.gpsSeconds=gpsMidSegTime.gpsSeconds+segmentDuration;  

   }
            
	    
  /* close frame caches */
  LAL_CALL( LALFrClose( &status, &frStream1), &status );
  LAL_CALL( LALFrClose( &status, &frStream2), &status );
  
  /* cleanup */
  LAL_CALL( LALDestroyRealFFTPlan(&status,&(specparPSD1.plan)),&status);
  LAL_CALL( LALDestroyRealFFTPlan(&status,&(specparPSD2.plan)),&status);
  XLALDestroyREAL4Window(dataWindow1);
  XLALDestroyREAL4Window(dataWindow2);
  LAL_CALL(LALDestroyRealFFTPlan(&status,&fftDataPlan1),&status);
  LAL_CALL(LALDestroyRealFFTPlan(&status,&fftDataPlan2),&status);
  LAL_CALL( LALDestroyVector(&status, &(segment1.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(segment2.data)), &status );

  XLALDestroyREAL4Window(hannWindow1);
  XLALDestroyREAL4Window(hannWindow2);
  LAL_CALL( LALCDestroyVector(&status, &(hBarTildeTemp1.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(hBarTildeTemp2.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(hBarTilde1.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(hBarTilde2.data)), &status );
  LAL_CALL( LALCDestroyVector(&status, &(h1StarH2.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(PSDTemp1.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(PSDTemp2.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(avPSD1.data)), &status );
  LAL_CALL( LALDestroyVector(&status, &(avPSD2.data)), &status );

  return 0;
 }
	    
	
 /* parse command line options */
void parseOptions(INT4 argc, CHAR *argv[])
 {
  int c = -1;

  while(1)
   {
    static struct option long_options[] =
     {
      /* options that set a flag */
      {"double", no_argument, &double_flag,1},
      {"double-high-pass-filter", no_argument, &double_high_pass_flag, 1},
      {"recenter", no_argument, &recenter_flag, 1},
      {"apply-mask", no_argument, &apply_mask_flag, 1},
      {"high-pass-filter", no_argument, &high_pass_flag, 1},
      {"overlap-hann", no_argument, &overlap_hann_flag, 1},
      {"verbose", no_argument, &verbose_flag, 1},
      {"test", no_argument, &test_flag, 1},

      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"gps-start-time", required_argument, 0, 't'},
      {"gps-end-time", required_argument, 0, 'T'},
      {"ifo1", required_argument, 0, 'i'},
      {"ifo2", required_argument, 0, 'I'},
      {"channel1", required_argument, 0, 'c'},
      {"channel2", required_argument, 0, 'C'},
      {"frame-cache1", required_argument, 0, 'd'},
      {"frame-cache2", required_argument, 0, 'D'},
      {"segment-duration", required_argument, 0, 'l'},
      {"sample-rate1", required_argument, 0, 's'},
      {"sample-rate2", required_argument, 0, 'S'},
      {"resample-rate1", required_argument, 0, 'r'},
      {"resample-rate2", required_argument, 0, 'R'},
      {"f-min", required_argument, 0, 'f'},
      {"f-max", required_argument, 0, 'F'},
      {"hann-duration", required_argument, 0, 'w'},
      {"hpf-frequency", required_argument, 0, 'k'},
      {"hpf-attenuation", required_argument, 0, 'p'},
      {"hpf-order", required_argument, 0, 'P'},
      {"mask-bin", required_argument, 0, 'b'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
     };

    /* getopt_long stores the option here */
    int option_index = 0;

    c = getopt_long(argc, argv, 
                  "ht:T:i:I:c:C:d:D:l:s:S:r:R:f:F:w:k:p:P:b:v",
 		   long_options, &option_index);

    if (c == -1)
     {
      /* end of options, break loop */
      break;
     }

    switch(c)
     {
      case 0:
             /* If this option set a flag, do nothing else now. */
             if (long_options[option_index].flag != 0)
              break;
             printf ("option %s", long_options[option_index].name);
             if (optarg)
              printf (" with arg %s", optarg);
             printf ("\n");
             break;

      case 'h':
               /* HELP!!! */
               displayUsage(0);
               break;

      case 't':
               /* start time */
	       startTime = atoi(optarg);
	       break;

      case 'T':
	       /* stop time */
	       stopTime = atoi(optarg);
	       break;

      case 'i':
	       /* ifo for first stream */
	       strncpy(ifo1, optarg, LALNameLength);
               break;

      case 'I':
	       /* ifo for first stream */
	       strncpy(ifo2, optarg, LALNameLength);
               break;
     
      case 'c':
	       /* ifo for first stream */
	       strncpy(channel1, optarg, LALNameLength);
               break;

      case 'C':
	       /* ifo for first stream */
	       strncpy(channel2, optarg, LALNameLength);
               break;


      case 'd':
               /* data cache one */
               strncpy(frameCache1, optarg, 200);
               break;

      case 'D':
               /* data cache two */
               strncpy(frameCache2, optarg, 200);
               break;

      case 'l':
	       /* duration */
	       segmentDuration = atoi(optarg);
	       break;
      case 's':
               /* sample rate */
               sampleRate1 = atoi(optarg);
               break;

      case 'S':
              /* sample rate */
              sampleRate2 = atoi(optarg);
              break;

      case 'r':
	      /* resampling */
	      resampleRate1 = atoi(optarg);
	      break;
     
      case 'R':
	      /* resampling */
	      resampleRate2 = atoi(optarg);
	      break;

      case 'f':
	       /* minimal frequency */
	       fMin = atoi(optarg);
	       break;

      case 'F':
	       /* maximal frequency */
	       fMax = atoi(optarg);
	       break;

      case 'w':
	       /* hann window duration */
	       hannDuration = atoi(optarg);
	       break;

      case 'k':
	       /* high pass knee filter frequency  */
	       highPassFreq= atof(optarg);
	       break;
                          
      case 'p':
	       /* high pass filter attenuation  */
	       highPassAt = atof(optarg);
	       break;

      case 'P':
	       /* high pass filter order  */
	       highPassOrder = atoi(optarg);
   

      case 'b':
              /* bin for frequency mask */
              maskBin = atoi(optarg);
              break;
        
     case 'v':
	     /* display version info and exit */
	     fprintf(stdout, "Standalone SGWB Search Engine\n" CVS_ID "\n");
	     exit(0);
	     break;

     default:
     displayUsage(1);
     }
    }

   if (optind < argc)
    {
     displayUsage(1);
    }

  return;
}

/* display program usage */
void displayUsage(INT4 exitcode)
 {
  fprintf(stderr, "Usage: pipeline [options]\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, " -h                    print this message\n");
  fprintf(stderr, " -V                    display version\n");
  fprintf(stderr, " --verbose             verbose mode\n");
  fprintf(stderr, " --test                print intermediate results\n");
  fprintf(stderr, " --double              read double precision data in frames\n"); 
   fprintf(stderr, " --double-high-pass-filter              high pass data before casting to single precision\n"); 
  fprintf(stderr, " --recenter            recenter jobs\n");
  fprintf(stderr, " -t                    GPS start time\n");
  fprintf(stderr, " -T                    GPS stop time\n");
  fprintf(stderr, " -i                    ifo for first stream\n");
  fprintf(stderr, " -I                    ifo for second stream\n");
  fprintf(stderr, " -c                    channel for first stream\n");
  fprintf(stderr, " -C                    channel for second stream\n");
  fprintf(stderr, " -d                    cache file for first stream\n");
  fprintf(stderr, " -D                    cache file for second stream\n");
  fprintf(stderr, " -l                    segment duration\n");
  fprintf(stderr, " -s                    sample rate for first stream\n");
  fprintf(stderr, " -S                    sample rate for second stream\n");
  fprintf(stderr, " -r                    resample rate for first stream\n");
  fprintf(stderr, " -R                    resample rate for second stream\n");
  fprintf(stderr, " -f                    minimal frequency\n");
  fprintf(stderr, " -F                    maximal frequency\n");
  fprintf(stderr, " -- high-pass-filter   apply high pass filter\n");
  fprintf(stderr, " -k                    high pass filter knee frequency\n");
  fprintf(stderr, " -p                    high pass filter attenuation\n");
  fprintf(stderr, " -P                    high pass filter order\n");        
  fprintf(stderr, " --overlap-hann        use overlap window\n");             
  fprintf(stderr, " -w                    hann duration\n");
  fprintf(stderr, " --apply-mask          apply frequency masking\n");
  fprintf(stderr, " -b                    number of bin for frequency mask\n");
  fprintf(stderr, " -O                    directory for output files\n");
  fprintf(stderr, " -z                    debugging level\n");     
  exit(exitcode);
}
