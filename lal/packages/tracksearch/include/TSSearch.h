/********************************** <lalVerbatim file="TrackSearchHV">
Author: Torres, C
$Id$
**************************************************** </lalVerbatim> */

#ifndef _TSSEARCH_H
#define _TSSEARCH_H

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/TimeFreq.h>
#include <lal/TrackSearch.h>
#include <lal/LALRCSID.h>
#include <lal/Window.h>

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


typedef struct
tagTSDataSegment
{
  REAL4TimeSeries               *TSSearchData;
  COMPLEX8TimeSeries            *TSSearchResponse;
  REAL4TimeSeries               *TSSearchSpectrum;
  UINT4                          SegNum; /*Segment Number*/
}
TSDataSegment;

typedef struct 
tagTSSegmentVector
{
  UINT4             length; /* Number of segments long */
  TSDataSegment    *dataSeg; /* Structure for individual data segments */
}TSSegmentVector;

typedef struct
tagTSSearchParams
{
  BOOLEAN                       searchMaster;
  BOOLEAN                       haveData;
  UINT4                        *numSlaves;          
  LIGOTimeGPS                   GPSstart;
  UINT4                         TimeLengthPoints;/*Each Data Seg Lengths*/
  UINT4                         NumSeg;/* Number of segments length TLP */
  LIGOTimeGPS                   Tlength;/*Data set time length*/
  TimeFreqRepType               TransformType;
  INT4                          LineWidth;
  REAL4                         StartThresh;
  REAL4                         LinePThresh;
  INT4                          MinLength;
  REAL4                         MinPower;
  UINT4                         FreqBins; /*Number of bins to use*/
  UINT4                         TimeBins; /*Number of bins to use*/
  INT4                          windowsize; /*Number of points in window*/
  WindowType                    window; /*Window to use*/
  UINT4                         numEvents; /*Does map have features*/
  CHAR                         *channelName; /*Data Channel Name */
  TSSegmentVector               dataSegVec; /*Vector of NumSeg of data */
  UINT4                         currentSeg; /*Denotes current chosen seg */
}TSSearchParams;

#ifdef  __cplusplus

#endif  /* C++ protection. */

#endif
