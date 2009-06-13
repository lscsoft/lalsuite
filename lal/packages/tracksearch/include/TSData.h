/*
*  Copyright (C) 2007 Cristina Valeria Torres, Jolien Creighton, Kipp Cannon
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
 * File TSData.h
 * Author: Torres, C
 */
#if 0
<lalVerbatim file="TSDataHV">
Author C. Torres
</lalVerbatim>

<lalLaTeX file="TSDataH">
\section{Header \texttt{TSData.h}}
\label{s:TSData.h}

Provides an intermediate level of functions and structures for testing
and production use of the tracksearch libraries.

\begin{verbatim}
#include <lal/TSData.h>
\end{verbatim}

\subsection*{Error Conditions}
\input{TSDataHErrTab}

</lalLaTeX>
#endif

#ifndef _TSDATA_H
#define _TSDATA_H

#include <lal/Date.h>
#include <lal/FrameCache.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameData.h>
#include <lal/FrameStream.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/LALMoment.h>
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/Matrix.h>
#include <lal/RealFFT.h>
#include <lal/TSSearch.h>
#include <lal/TimeSeries.h>
#include <string.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID (TSDATAH, "$Id$");

/******** <lalErrTable file="TSDataHErrTab"> ********/
#define TSDATA_ENULL    1
#define TSDATA_ENNUL    2
#define TSDATA_EALOC    4
#define TSDATA_ESEGZ    8
#define TSDATA_ENUMZ    16
#define TSDATA_ESUBR    32
#define TSDATA_EWHIT    64
#define TSDATA_ERESP    128
#define TSDATA_EINTP    256
#define TSDATA_EINVA    512

#define TSDATA_MSGENULL "Null pointer"
#define TSDATA_MSGENNUL "Non-null pointer"
#define TSDATA_MSGEALOC "Memory allocation error"
#define TSDATA_MSGESEGZ "Invalid number of segments"
#define TSDATA_MSGENUMZ "Invalid number of points in segment"
#define TSDATA_MSGESUBR "Condition Data Internal Subroutine Failure"
#define TSDATA_MSGEWHIT "Response function does not have enough frequency information to whiten data"
#define TSDATA_MSGERESP "Response function start frequency not equal to 0"
#define TSDATA_MSGEINTP "Not enough points for interpolate function."
#define TSDATA_MSGEINVA "Inconsistent Argument(s)"
/******** </lalErrTable> ********/

/*
 * All fields are depricated minus dataSegmentPoints and numberDataSegments
 */
typedef struct
tagTSCreateParams
{
  UINT4      dataSegmentPoints;     /* Fixed Length Varies Values  */
  UINT4      responseSegmentPoints; /* Fixed Copy for each segment */
  UINT4      spectraSegmentPoints;  /* Fixed Copy for each segment */
  UINT4      numberDataSegments;    /* Number of data segments to analyze */
  UINT4      SegBufferPoints; /*Number of data points to buffer seg with*/
}TSCreateParams;


  /*
   * Struture to allow TSDatagen to make fake signals files/frames
   */
typedef struct
tagTSDatagen
{
  REAL8      Amp_I;
  REAL8      Amp_F;
  REAL8      F_start;
  REAL8      F_end;
  REAL8      Fsample;
  REAL8      Noise_Amp;
}TSDatagen;


/*
 * Struture used by the line connection subroutine
 * we want to take all event candidates and match them via param SIGMA
 */
typedef struct
tagTSConnectParams
{
  TrackSearchParams     *generalParams;
  TimeFreqRep           *map;
}TSConnectParams;


/*
 * Struture to whiten the data including renormalization information
 * we plan on following EPSearch definition
 * This structure might be dropped
 */
typedef struct
tagTSWhitenParams
{
  BOOLEAN         renormalize;
  REAL4           mean;
  REAL4           variance;
  UINT4           whitenLevel;
}TSWhitenParams;

/*
 * Routine to determine the best Lh and set Ll given Lrelative.
 */
void
LALTracksearchFindLambdaMean(
			     LALStatus                *status,
			     TimeFreqRep               map,
			     TSSearchParams           *searchParams
			     );
void
LALTracksearchFindLambdaMedian(
			       LALStatus                *status,
			       TimeFreqRep               map,
			       TSSearchParams           *searchParams
			       );
/*
 * Routine to break up time series input and make a collection of
 * segments which overlap by the overlap(points) parameter
 */
void
LALCreateTSDataSegmentVector (
			      LALStatus                  *status,
			      TSSegmentVector           **vector,
			      TSCreateParams             *params
			      );

/*
 * Routine to deallocate this collection of segments
 */
void
LALDestroyTSDataSegmentVector (
			       LALStatus                  *status,
			       TSSegmentVector            *vector
			       );

/*
 * This routine will handle a great deal of data conditioning
 * Whitening and calibration will occur in this function
 */
void TrackSearchPrep(
		     LALStatus                      *status,
		     REAL4TimeSeries                *TSSearchData,
		     COMPLEX8Vector                 *TSSearchResponse,
		     REAL8FrequencySeries           *TSSearchSpectrum,
		     TSSegmentVector                *PreparedData,
		     TSSearchParams                  params
		     );
/*
 * This is a less functional version of TrackSearchPrep which is
 * greatly simplified in light of analysis pipeline design
 */
void LALTrackSearchDataSegmenter(
				 LALStatus           *status,
				 REAL4TimeSeries     *TSSearchData,
				 TSSegmentVector     *PreparedData,
				 TSSearchParams       params);
/*
 * This routine applies our thresholds on length and power
 * It expects to allocate a thresholded candidate list
 */
void
LALTrackSearchApplyThreshold(
			     LALStatus         *status,
			     TrackSearchOut    *curveinfo,
			     TrackSearchOut    *dataProduct,
			     TSSearchParams     params
			     );

/*
 * Routine connects candidates who begin or end within SIGMA of
 * one another
 */
void
LALTrackSearchConnectSigma(
			   LALStatus                   *status,
			   TrackSearchOut              *curveinfo,
			   TimeFreqRep                  map,
			   TrackSearchParams            params
			   );
/*
 * This routine can be run alone but is called via TrackSearchPrep to
 * whiten the time series signal as done in Power
 */
void
LALTrackSearchWhitenREAL4TimeSeries(
				    LALStatus              *status,
				    REAL4TimeSeries        *signal,
				    REAL4FrequencySeries   *signalPSD,
				    TSWhitenParams          params
				    );

/*
 * This function does simple manipulation to avoid
 * wasting CPU time on FFTs
 */

void
LALTrackSearchWhitenCOMPLEX8FrequencySeries(
					    LALStatus                *status,
					    COMPLEX8FrequencySeries  *fSeries,
					    REAL4FrequencySeries     *PSD,
					    UINT4                     level
					    );


/*
 * This routine can be run alone but is called via TrackSearchPrep to
 * use a response curve for the segment epoch and calibrate that
 * time series segment DEAD FUNCTION
 */
void
LALTrackSearchCalibrateREAL4TimeSeries(LALStatus               *status,
				       REAL4TimeSeries         *signal,
				       COMPLEX8FrequencySeries *response);

/*
 * This routine will calibrate the fourier data given a corresponding
 * transfer function working just in fourier domain
 */
void
LALTrackSearchCalibrateCOMPLEX8FrequencySeries(
					       LALStatus                 *status,
					       COMPLEX8FrequencySeries   *fSeries,
					       COMPLEX8FrequencySeries   *response
					       );

/*
 * This routine is unfinished it was an attempt to mimic
 * Matlab's interp1 routine
 */
void
LALSVectorPolynomialInterpolation(
				  LALStatus        *status,
				  REAL4Sequence    *newDomain,
				  REAL4Sequence    *newRange,
				  REAL4Sequence    *Domain,
				  REAL4Sequence    *Range
				  );
/*
 * Noncompliant code
 * Local function not meant for general use
 */
void connect2Segments(
		      TimeFreqRep    map,
		      Curve          *curveA,
		      Curve          *curveB
		      );

void WriteMap(
	      TimeFreqRep       map,
	      REAL4Vector       signal
	      );

void DumpTFImage(
		 REAL4         **image,
		 const CHAR           *filename,
		 INT4           height,
		 INT4           width,
		 BOOLEAN        killNeg
		 );

void DumpTFImageCHAR(
		 CHAR         **image,
		 const CHAR           *filename,
		 INT4           height,
		 INT4           width,
		 BOOLEAN        killNeg
		 );

void DumpREAL8KernelMask(
			 REAL8      *kernel,
			 const CHAR       *filename,
			 INT4        ksize
			 );

void cleanLinkedList(
		     TrackSearchOut      *inList,
		     TrackSearchOut      *outList
		     );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif


