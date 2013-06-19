/*
*  Copyright (C) 2007 Cristina Valeria Torres
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

/**
      \author Torres. C
      \file
*/

#ifndef TRACKSEARCH_H
#define TRACKSEARCH_H

#include "config.h"

#include <math.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <getopt.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FrameCalibration.h>
#include <lal/LALFrStream.h>
#include <lal/GenerateBurst.h>
#include <lal/IIRFilter.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALError.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/SeqFactories.h>
#include <lal/TSData.h>
#include <lal/TSSearch.h>
#include <lal/TimeFreq.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/CLR.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/LALRunningMedian.h>
#include <lal/RngMedBias.h>

#include <lalapps.h>
#include <processtable.h>

#include "tracksearchToolbox.h"
#include "FrameDataConvert.h"
#include "tracksearchAverager.h"

#define maxFilenameLength 2048

/*
 * LALapps Structures
 */
typedef struct
tagTSappsInjectParams
{
  REAL8    startTimeOffset; /*Seconds to wait before first inject*/
  UINT4    numOfInjects; /*Number of total injects to try and do*/
  REAL8    injectSpace; /*Number of seconds between injects*/
  REAL8    scaleFactor; /*Multiplicative scale factor sets ~SNR*/
  REAL8    sampleRate; /*Sample rate to make data at if diffs from file*/
  CHARVector *injectTxtFile; /*Path&Name of file to inject*/
}TSappsInjectParams;

/*
 * lalapps functions
 */
/*
 * SemiPrivate functions
 */
void
LALappsTrackSearchPrepareData(
			      REAL4TimeSeries*,
			      REAL4TimeSeries*,
			      TSSegmentVector*,
			      TSSearchParams);

int
LALappsQuickHeterodyneTimeSeries(REAL4TimeSeries*,
				 REAL8);

int
LALappsSmoothWithRunningMedian(REAL4Vector*,
			       UINT4);

void
LALappsTrackSearchCalibrate( REAL4TimeSeries*,
			     TSSearchParams);

void
LALappsTrackSearchBandPassing(
			       REAL4TimeSeries*,
			       TSSearchParams);

void
LALappsTracksearchRemoveHarmonics(
				   REAL4TimeSeries*,
				   TSSearchParams);

void LALappsTracksearchRemoveHarmonicsFromSegments(
						   REAL4TimeSeries*,
						   TSSegmentVector*,
						   TSSearchParams);
void
LALappsTrackSearchPerformInjection(
				    REAL4TimeSeries*,
				    REAL4TimeSeries*,
				    TSSearchParams);

void
LALappsTrackSearchWhitenSegments( REAL4TimeSeries*,
				  TSSegmentVector*,
				  TSSearchParams);

void
LALappsTrackSearchInitialize(
			     int argc,
			     char* argv[],
			     TSSearchParams*,
			     TSappsInjectParams*,
			     CHARVector**,
			     CHARVector**);
void
LALappsGetFrameData(TSSearchParams*,
		    REAL4TimeSeries*,
		    CHARVector*,
		    CHAR*);

void
LALappsGetAsciiData(
		    TSSearchParams*,
		    REAL4TimeSeries*,
		    CHARVector*);

void
LALappsDoTrackSearch(TimeFreqRep*,
		     TrackSearchParams,
		     TrackSearchMapMarkingParams ,
		     TSSearchParams);

void
LALappsDoTSeriesSearch(REAL4TimeSeries*,
		       TSSearchParams,
		       INT4);

void
LALappsDoTimeSeriesAnalysis(TSSearchParams,
			    TSappsInjectParams,
			    CHAR*,
			    CHARVector*);

void
LALappsDoTSAMapSearch(TSAMap*,
		      TSSearchParams*,
		      INT4);


void
LALappsDoTSAMapAnalysis(LALStatus*,
			TSSearchParams);


void
LALappsWriteCurveList(LALStatus*,
		      CHAR*,
		      TrackSearchOut,
		      TSSearchParams*);

void
LALappsWriteSearchConfig(
			 const CHAR*,
			 TSSearchParams);

void
LALappsWriteSearchResults(
			  const CHAR*,
			  TrackSearchOut);

void
LALappsWriteBreveResults(
			 const CHAR*,
			 TrackSearchOut);

void
LALappsCreateInjectableData(REAL4TimeSeries**,
			    TSappsInjectParams);
/*
 * Following routine need to be fixed
 * they limp along and are prone to error
 * if use is not careful
 * CONSIDER PLACING IN TSDATA.C FROM LALLIBS
 */
int
LALappsCreateCurveDataSection(Curve**);

void
LALappsDestroyCurveDataSection(Curve**,
			       INT4);
/*
 * Private functions
 */
void Dump_Search_Data(TSSearchParams,TrackSearchOut,CHAR*);
void QuickDump_Data(TSSearchParams,TrackSearchOut,CHAR*);
void fakeDataGeneration(REAL4TimeSeries*,INT4,INT4);

#endif

