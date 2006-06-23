/**** <lalVerbatim file="TSDatgenHV"> *********
Author: Torres. C
$ID: tracksearch.h,v 1.0 2004/04/14 02:00:00 charlie Exp $
***** </lalVerbatim> **********************************/

#ifndef TRACKSEARCH_H
#define TRACKSEARCH_H

#include <config.h>
#include <fcntl.h>
#include <getopt.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/GenerateBurst.h>
#include <lal/IIRFilter.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALError.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/PrintFTSeries.h>
#include <lal/Random.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/SeqFactories.h>
#include <lal/StdBurstSearch.h>
#include <lal/TSData.h>
#include <lal/TSSearch.h>
#include <lal/TimeFreq.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lalapps.h>
#include <math.h>
#include <processtable.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <tracksearchToolbox.h>
#include <tracksearchAverager.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/LALRunningMedian.h>
#include <lal/RngMedBias.h>

#define maxFilenameLength 2048

/*
 * lalapps functions
 */
/*
 * SemiPrivate functions
 */
void
LALappsTrackSearchPrepareData(
			      LALStatus*,      
			      REAL4TimeSeries*,
			      TSSegmentVector*,
			      TSSearchParams);

void 
LALappsTrackSearchInitialize(LALStatus*,
			     int argc,
			     char* argv[],
			     TSSearchParams*,
			     CHARVector**,
			     CHARVector**);
void
LALappsGetFrameData(LALStatus*,
		    TSSearchParams*,
		    REAL4TimeSeries*,
		    CHARVector*,
		    CHAR*);

void
LALappsGetAsciiData(LALStatus*,
		    TSSearchParams*,
		    REAL4TimeSeries*,
		    CHARVector*);

void
LALappsDoTrackSearch(LALStatus*,
		     TimeFreqRep*,
		     TrackSearchParams,
		     TrackSearchMapMarkingParams ,
		     TSSearchParams);

void
LALappsDoTSeriesSearch(LALStatus*,
		       REAL4TimeSeries*,
		       TSSearchParams,
		       INT4);

void
LALappsDoTimeSeriesAnalysis(LALStatus*,
			    TSSearchParams,
			    CHAR*,
			    CHARVector*);

void
LALappsDoTSAMapSearch(LALStatus*,
		      TSAMap*,
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

/*
 * Following routine need to be fixed
 * they limp along and are prone to error
 * if use is not careful
 * CONSIDER PLACING IN TSDATA.C FROM LALLIBS
 */
void
LALappsCreateCurveDataSection(LALStatus*,
			      Curve**);

void
LALappsDestroyCurveDataSection(LALStatus*,
				Curve**,
				INT4);
/*
 * Private functions
 */
void Dump_Search_Data(TSSearchParams,TrackSearchOut,CHAR*);
void QuickDump_Data(TSSearchParams,TrackSearchOut,CHAR*);
void fakeDataGeneration(LALStatus*,REAL4TimeSeries*,INT4,INT4);

#endif
