/**** <lalVerbatim file="PowerInitSearchHV"> ***********
Author: Brady, P and Brown, D
$Id$ 
***** </lalVerbatim> ***********************************/

#ifndef _INITSEARCH_H
#define _INITSEARCH_H

/* #include <lal/LALStdlib.h> */
#include "Power.h"

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (INITSEARCHH, "$Id$");

/******** <lalErrTable file="PowerInitSearchHErrTab"> ********/
#define INITSEARCHH_ENULL     1
#define INITSEARCHH_ENNUL     2
#define INITSEARCHH_EALOC     3
#define INITSEARCHH_EARGS     4
#define INITSEARCHH_ENUMZ     5
#define INITSEARCHH_ESEGZ     6
#define INITSEARCHH_EOVLP     7
#define INITSEARCHH_EOVPF     8
#define INITSEARCHH_EMFBZ     9
#define INITSEARCHH_EMTBZ     10
#define INITSEARCHH_EFLOW     11
#define INITSEARCHH_EDELF     12
#define INITSEARCHH_ELTFZ     13
#define INITSEARCHH_ESIGM     14
#define INITSEARCHH_EALPH     15
#define INITSEARCHH_EFREE     16
#define INITSEARCHH_ENLAL     17
#define INITSEARCHH_ENDAT     18
#define INITSEARCHH_EDUTY     19
#define INITSEARCHH_EAMAX     20 
#define INITSEARCHH_EE2MS     21 
#define INITSEARCHH_ECHNL     22 
#define INITSEARCHH_ESIM      23


#define INITSEARCHH_MSGENULL     "Null pointer"
#define INITSEARCHH_MSGENNUL     "Non-null pointer"
#define INITSEARCHH_MSGEALOC     "Memory allocation error"
#define INITSEARCHH_MSGEARGS     "Wrong number of arguments"
#define INITSEARCHH_MSGENUMZ     "Data segment length is zero or negative"
#define INITSEARCHH_MSGESEGZ     "Number of data segments is zero or negative"
#define INITSEARCHH_MSGEOVLP     "Overlap of data segments is negative"
#define INITSEARCHH_MSGEOVPF     "Overlap of TF tiles is negative"
#define INITSEARCHH_MSGEMFBZ     "Smallest extent in freq of TF tile is zero" 
#define INITSEARCHH_MSGEMTBZ     "Smallest extent in time of TF tile is zero"
#define INITSEARCHH_MSGEFLOW     "Lowest frequency to be searched is <= 0"
#define INITSEARCHH_MSGEDELF     "Freq resolution of 1st time-freq plane is <=0"
#define INITSEARCHH_MSGELTFZ     "Length (Nf) of 1st TF plane (with Nt=1) is <=0"
#define INITSEARCHH_MSGESIGM     "Threshold number of sigma is <=1"
#define INITSEARCHH_MSGEALPH     "Default alpha value is out of range"
#define INITSEARCHH_MSGEFREE     "Memory free error"
#define INITSEARCHH_MSGENLAL     "Tried to allocate to non-null pointer"
#define INITSEARCHH_MSGENDAT     "No data read"
#define INITSEARCHH_MSGEDUTY     "Number of segments sent to slave is zero"
#define INITSEARCHH_MSGEAMAX     "The threshold value of alpha is negative"
#define INITSEARCHH_MSGEE2MS     "Number of events out of range[1..99]"
#define INITSEARCHH_MSGECHNL     "Channel name not valid"
#define INITSEARCHH_MSGESIM      "Invalid simulation type: 0, 1, or 2"
/******** </lalErrTable> ********/

/*
 *
 * libpower.so should be invoked with the following arguments to wrapperAPI
 *
 * argv[0]  = "-filterparams"
 * argv[1]  = numPoints         Number of data points in a segment (2 * 2^p + 2)
 * argv[2]  = numSegments       Number of overlapping segments
 * argv[3]  = ovrlap            Number of points overlap between segments
 * argv[4]  = overlapFactor     Amount of overlap between neighboring TF tiles
 * argv[5]  = minFreqBins       Smallest extent in freq of TF tiles to search
 * argv[6]  = minTimeBins       Smallest extent in time of TF tiles to search
 * argv[7]  = flow              Lowest frequency in Hz to be searched
 * argv[8]  = deltaF            Frequency resolution of first TF plane
 * argv[9]  = length            Length (N_F) of first TF plane (with N_T = 1)
 * argv[10] = numSigmaMin       threshold number of sigma
 * argv[11] = alphaDefault      default alpha value for tiles with sigma < numSigmaMin
 * argv[12] = segDutyCycle      Number of segments sent to slave for analysis
 * argv[13] = alphaThreshold    Indentify events with alpha less than this
 * argv[14] = events2Master     Number of events out or range[1..99]
 * argv[15] = channelName       Channel name
 */

#ifdef  __cplusplus
}
#endif

#endif /* _INITSEARCH_H  */


