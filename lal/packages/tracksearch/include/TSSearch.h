/*
*  Copyright (C) 2007 Cristina Valeria Torres, Kipp Cannon
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

/********************************** <lalVerbatim file="TrackSearchHV">
Author: Torres, Cristina
$Id$
**************************************************** </lalVerbatim> */

#ifndef _TSSEARCH_H
#define _TSSEARCH_H

#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Interpolate.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdlib.h>
#include <lal/PrintFTSeries.h>
#include <lal/PrintVector.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/TSSearch.h>
#include <lal/TimeFreq.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/TrackSearch.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <stdlib.h>
#include <string.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

  NRCSID (TSSEARCHH, "$Id$");

  /******** <lalErrTable file="TSSearchHErrTab"> ********/
#define TSSEARCHH_ENULLP       1
#define TSSEARCHH_EPOSARG      2
#define TSSEARCHH_EPOW2        4
#define TSSEARCHH_EMALLOC      8
#define TSSEARCHH_EINCOMP      16
#define TSSEARCHH_EORDER       32
#define TSSEARCHH_ENONNULL     64
#define TSSEARCHH_ETILES       65
#define TSSEARCHH_EDELF        128


#define TSSEARCHH_MSGENULLP    "Null pointer"
#define TSSEARCHH_MSGEPOSARG   "Arguments must be non-negative"
#define TSSEARCHH_MSGEPOW2     "Length of supplied data must be a power of 2"
#define TSSEARCHH_MSGEMALLOC   "Malloc failure"
#define TSSEARCHH_MSGEINCOMP   "Incompatible arguments"
#define TSSEARCHH_MSGEORDER    "Routines called in illegal order"
#define TSSEARCHH_MSGENONNULL  "Null pointer expected"
#define TSSEARCHH_MSGETILES    "Malloc failed while assigning memory for a tile"
#define TSSEARCHH_MSGEDELF     "Inconsistent deltaF in spectrum and data"
  /******** </lalErrTable> ********/

  /*
   * Enumeration type to hold diagnostic flag information
   */
  typedef enum tagTSDiagnosticType
    {
      quiet, verbose, printFiles, all
    }TSDiagnosticType;


  /*
   * Structure to hold a collection of data segments which may be
   * overlapped by n points
   */
  typedef struct
  tagTSSegmentVector
  {
    UINT4               length;  /* Number of segments long */
    UINT4               SegBufferPoints; /*
					  * Number of points before
					  * and after which are buffer
					  * to remove TFR artifacts.
					  */
    REAL4TimeSeries   **dataSeg; /*
				  * Structure for individual data
				  * segments
				  */
  }TSSegmentVector;

  /*
   * Structure that holds all possible parameters for the tracksearch
   * library functions collectively not just LALTracksearch
   */
  typedef struct
  tagTSSearchParams
  {
    BOOLEAN           tSeriesAnalysis; /*Yes we prep for tseries processing*/
    BOOLEAN           searchMaster; /*DO NOT USE*/
    BOOLEAN           haveData;/*DO NOT USE */
    UINT4            *numSlaves;/* DO NOT USE*/
    LIGOTimeGPS       GPSstart; /* GPS start time of entire stretch */
    UINT4             TimeLengthPoints; /* Product of NumSeg&SegLenthPoints*/
    UINT4             discardTLP;/* Points that need to be
				  * discarded given
				  * input map overlap
				  */
    UINT4             SegLengthPoints;/*Data Seg Length*/
    UINT4             colsToClip;/*
				  * Should be colsToClip before and
				  * after region of interest.
				  */
    UINT4             SegBufferPoints;/*
				       *Extra points that are later
				       *clipped from TFR making correct
				       *time boundaries to avoid TFR
				       *edge effects, should be at
				       *least 4 cols worth of
				       *points. Two before and two
				       *after part of TFR we want to
				       * keep
				       */
    UINT4             NumSeg;/* Number of segments length TLP */
    REAL8             SamplingRate; /*
				     * Sampling rate of input data
				     * should be resampled at
				     */
    REAL8             SamplingRateOriginal; /*
					     * Original sampling rate
					     * of input time series
					     * data
					     */
    LIGOTimeGPS       Tlength;/*Data set time length*/
    TimeFreqRepType   TransformType;/*Type of TF rep to make */
    INT4              LineWidth;/* Sigma-convolution kernel width*/
    REAL4             StartThresh;/*Lh-2ndDerv Start threshold*/
    REAL4             LinePThresh;/*Ll-2ndDerv Member threshold*/
    REAL4             autoThresh;/* Z score auto Lh threshold value*/
    REAL4             relaThresh;/* Relative LinePThresh given found Lh*/
    BOOLEAN           autoLambda;/*Bool flag to track auto Lh state */
    INT4              MinLength;/*Minimum length of a curve*/
    REAL4             MinPower;/*Minimum Power in a curve*/
    REAL4             MinSNR;/*minimum estimate SNR to accept*/
    UINT4             overlapFlag;/*
				   * Number of points to overlap the
				   * individual data segments by. This
				   * is apart from implicit Segment
				   * buffer for TFR clipping option
				   */
    UINT4             whiten;/*Flags typ of whitening to do*/
    AvgSpecMethod     avgSpecMethod;/*Type of PSD averaging to do*/
    WindowType        avgSpecWindow;/*Type of PSD averaging window*/
    UINT4             multiResolution;/*FlagMultiResolutionRun*/
    UINT4             FreqBins; /*Number of bins to use*/
    UINT4             TimeBins; /*Number of bins to use*/
    INT4              windowsize; /*Number of points in window*/
    WindowType        window; /*Window to use in TF map creation*/
    UINT4             numEvents; /*Does map have features*/
    CHAR             *channelName; /*Data Channel Name */
    CHAR             *dataDirPath; /*Path to data frames */
    CHAR             *singleDataCache; /*Explicit name to 1 data cache*/
    CHAR             *detectorPSDCache;/*Explicit cache for PSD*/
    CHAR             *channelNamePSD;/*DO NOT USE*/
    CHAR             *calFrameCache;/*Cache file for cal frames*/
    BOOLEAN           calibrate;/*Calibration flag Y/N */
    CHAR              calibrateIFO[3];/*3LetterName of IFO*/
    CHAR             *calCatalog;/*Holds calib coeffs*/
    TSSegmentVector  *dataSegVec; /*Vector of NumSeg of data */
    UINT4             currentSeg; /*Denotes current chosen seg */
    INT4              makenoise; /*SeedFlag to fake lalapps data*/
    CHAR             *auxlabel; /*For labeling out etc testing*/
    BOOLEAN           joinCurves; /*Flag joins 1 sigma gap curves */
    TSDiagnosticType  verbosity; /* 1 verbose,2 aux files,0 quiet */
    BOOLEAN           printPGM; /* Create output PGMs (BWdefault) */
    CHAR             *pgmColorMapFile; /* User specifiable colormap AsciiPGM*/
    CHAR             *injectMapCache;/* filename to text set of maps*/
    CHAR             *injectSingleMap;/* explicit map file read */
    UINT4             smoothAvgPSD;/*(>0) True apply running median to AvgPSD*/
    REAL4             highPass;/*Real f value to high pass filter with*/
    REAL4             lowPass;/*Real f value to low pass filter with*/
    REAL4             listLinesToRemove[512];/*static list of lines 2 remove*/
    UINT4             numLinesToRemove;/*Tally of harms to remove*/
    UINT4             maxHarmonics;/*Max num of harms to remove*/
  }TSSearchParams;


  /*
   * This is a structure which gives detailed information about a
   * signal candidate.  CURRENTLY NOT USED.
   */
  typedef struct
  tagTrackSearchEvent
  {
    LIGOTimeGPS                   mapStartTime;/*GPS map start*/
    LIGOTimeGPS                   mapStopTime;/*GPS map stop*/
    REAL4                         samplingRate;/*Input data sample rate*/
    REAL4Vector                  *fvalues;/*Pointer to freq indexs*/
    REAL4Vector                  *pvalues;/*Pixel power values */
    REAL4Vector                  *tvalues;/*Pointer to time indexs*/
    REAL8                         peakPixelPower;/*Peak pixel value*/
    REAL8                         power;/*Integrated curve power*/
    UINT4                         FreqBins;/*Number of freq bins in map*/
    UINT4                         FstartPixel;/*Pixel Fstart location IMAGE*/
    UINT4                         FstopPixel;/*Pixel FStop location IMAGE*/
    UINT4                         TimeBins;/*Number of time bins in map*/
    UINT4                         TstartPixel;/*Pixel Tstart location IMAGE*/
    UINT4                         TstopPixel;/*Pixel TStop location IMAGE*/
    UINT4                         dateString;/*GPS candidate start*/
    UINT4                         durationPoints;/*Time points used in map*/
    UINT4                         fftLength;/* Num points in transform */
    UINT4                         junction; /* 1 yes 0 no*/
    UINT4                         overlap; /* Points for overlaped fft */
    UINT4                         peakPixelF;/*Peak pixel coord F*/
    UINT4                         peakPixelT;/*Peak pixel coord T*/
    UINT4                         whiten;  /* 0 No 1 yes 2 Overwhiten */
    UINT4                         windowsize;/* FFT window length */
    WindowType                    window;/* Type of window to make map*/
    struct tagTrackSearchEvent   *nextEvent; /* Pointer for linked listing */
  } TrackSearchEvent;

#ifdef  __cplusplus

#endif  /* C++ protection. */

#endif
